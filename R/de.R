# Live differential expression. presto::wilcoxauc needs an in-memory
# genes x cells matrix, so a DE run reconstructs the contrast's expression
# matrix from the Parquet store (the one runtime operation that materialises a
# chunk of the matrix — hence on-demand, behind a Compute button).

# Above this many tested cells (with no max_cells cap) scroll_de emits a heads-up
# that the in-RAM matrix may be large — the one place runtime memory scales.
.SCROLL_DE_WARN_CELLS <- 100000L

#' Live differential expression for a contrast (presto / Wilcoxon)
#'
#' Compares `ident1` against `ident2` or against the rest of the cells when
#' `ident2` is `NULL`/`"rest"`. Either side may name **several** levels of
#' `group_col` (they are pooled into one group). Reconstructs the contrast's
#' expression matrix from the store and runs [presto::wilcoxauc()].
#'
#' @param data A scroll data handle (`con`, `cells`, `manifest`).
#' @param assay Assay to test.
#' @param group_col Metadata column(s) defining the groups. Naming several columns
#'   compares their **interaction levels** (e.g. `c("genotype", "timepoint")` gives
#'   groups like `"KO | d7"`); `ident1`/`ident2` then select which combined levels
#'   form each side.
#' @param ident1 One or more levels of `group_col` forming the group of interest.
#' @param ident2 One or more comparison levels, or `NULL`/`"rest"` for
#'   one-vs-rest (every cell not in `ident1`).
#' @param min_pct Keep genes expressed in at least this fraction of either side.
#' @param cells Cells to test over (default `data$cells`); pass a subset to run
#'   DE within an active cell-subset filter.
#' @param max_cells Optional cap on cells **per group**. When set (a positive
#'   number), each side is randomly down-sampled to at most `max_cells` cells
#'   before testing, which bounds DE time and peak memory on large contrasts
#'   (a standard marker-detection shortcut; cf. Seurat's `max.cells.per.ident`).
#'   `NULL` (default) uses every cell. The subsample is deterministic (fixed
#'   seed) so repeated runs match, and the session RNG is left untouched.
#' @param progress Optional `function(fraction, detail)` called at the read and
#'   test stages, for a UI progress bar. `NULL` (default) is a no-op.
#' @return A data.frame of results, ranked by adjusted p-value, with columns
#'   `gene`, `logFC`, `avg_log2FC`, `auc`, `pct.1`, `pct.2`, `p_val`, `p_val_adj`.
#'   Two fold-change columns are reported because presto and Seurat define it
#'   differently: **`logFC`** is presto's value -- a natural-log "mean of log"
#'   difference (`mean(log-data in group1) - mean(...group2)`); **`avg_log2FC`**
#'   matches `Seurat::FindMarkers` -- the log2 of the mean of the *un-logged*
#'   normalized counts, `log2((sum(expm1(x1))+1)/n1) - log2((sum(expm1(x2))+1)/n2)`
#'   (Seurat v5 `FoldChange`, pseudocount 1). For zero-inflated data `avg_log2FC` is
#'   typically larger in magnitude than `logFC`.
#' @details With the default quantized build (`quantize = TRUE` in
#'   [scroll_build()]), expression values below ~`max/510` round to zero. Such a
#'   floored value is still stored (as an explicit zero), so `pct.1`/`pct.2` are
#'   unaffected; instead those cells become tied at zero, which makes the Wilcoxon
#'   p-values approximate, and `avg_log2FC` close to but not bit-identical to
#'   Seurat's (max abs diff ~0.04, median ~0.005 in practice; see
#'   `dev/validate_de_vs_findmarkers.R`). Build with `quantize = FALSE` for exact
#'   statistics -- `avg_log2FC`, `pct`, and the p-value ranking then match
#'   `Seurat::FindMarkers` to scroll's output rounding (Seurat's default `wilcox`
#'   test itself dispatches to presto when it is installed, so the two run the
#'   same Wilcoxon).
#' @export
scroll_de <- function(data, assay, group_col, ident1, ident2 = NULL, min_pct = 0.1,
                      cells = NULL, max_cells = NULL, progress = NULL) {
  if (!requireNamespace("presto", quietly = TRUE))
    stop("Live DE needs the 'presto' package.", call. = FALSE)
  pr <- function(f, d) if (is.function(progress)) progress(f, d)
  cells <- cells %||% data$cells
  ident1 <- as.character(ident1); ident1 <- ident1[nzchar(ident1)]
  if (!length(ident1)) stop("Pick at least one level for group 1.", call. = FALSE)
  # per-cell label ("group1" / "group2" / "rest", NA = not in the contrast); the
  # single source of truth shared with the panel's live contrast preview.
  lab <- .scroll_contrast_labels(cells, group_col, ident1, ident2)
  # `.gidx` rides on the cells handle; `kept` indexes the tested cells.
  gidx <- if (!is.null(cells$.gidx)) cells$.gidx else seq_len(nrow(cells))
  kept <- which(!is.na(lab))
  labels <- lab[kept]
  # optional per-group cap: down-sample each side to <= max_cells so a large
  # contrast stays bounded in time + peak RAM (Wilcoxon on a subsample is a
  # standard marker shortcut). Deterministic + RNG-neutral (see .scroll_cap_groups).
  if (!is.null(max_cells) && is.finite(max_cells)) {
    sel <- .scroll_cap_groups(labels, max_cells)
    kept <- kept[sel]; labels <- labels[sel]
  }
  if (sum(labels == "group1") < 3 || sum(labels != "group1") < 3)
    stop("Each side of the contrast needs at least 3 cells.", call. = FALSE)
  # Heads-up: DE materialises an all-genes x tested-cells matrix in RAM (the one
  # runtime operation whose memory scales with the matrix). On a very large
  # uncapped contrast that can be sizable — point the caller at max_cells.
  if (is.null(max_cells) && length(labels) > .SCROLL_DE_WARN_CELLS)
    message(sprintf(
      "scroll_de: testing %s cells with no max_cells cap; the in-memory matrix ",
      format(length(labels), big.mark = ",")),
      "may be large. Set max_cells to bound DE time and peak RAM.")
  # `bc` = barcodes (matrix dimnames); `key` = the store's cell key used to query +
  # join (v2 int global index / v1 barcode).
  bc  <- cells$cell[kept]
  key <- if (isTRUE(data$manifest$cell_index)) gidx[kept] else bc

  # dict = TRUE returns `feature` dictionary-encoded (a factor): its levels are the
  # store's OWN feature names and as.integer() gives the matrix row index directly.
  # This avoids sort(unique()) + match() over the tens of millions of repeated
  # gene-name strings (the dominant cost of a large DE) — arrow builds the
  # dictionary during the scan. Indexing by the store's names (not the manifest
  # list) also stays correct when yaml has mangled a non-ASCII feature name, and
  # zero-expression genes are absent here so restricting to expressed features is
  # correct too.
  pr(0.15, "Reading expression...")
  long <- scroll_query_cells(data$con, assay, key, dict = TRUE)
  val   <- scroll_dequantize(long$value, data$manifest, assay)
  feats <- levels(long$feature)
  i     <- as.integer(long$feature)
  # Map each stored `cell` (a global index / barcode) to its matrix column via a
  # reverse-index gather instead of match() over tens of millions of ids, and
  # drop `long` BEFORE allocating the matrix so the long table and the sparse
  # matrix never coexist — together this roughly halves peak RAM and speeds the
  # matrix build on a large contrast.
  revmap <- if (is.numeric(key)) {              # v2: keys are 1-based global indices
              r <- integer(max(key)); r[key] <- seq_along(key); r[long$cell]
            } else match(long$cell, key)         # v1 barcode store: fall back to match
  rm(long)
  X <- Matrix::sparseMatrix(i = i, j = revmap, x = val,
                            dims = c(length(feats), length(bc)), dimnames = list(feats, bc))
  rm(i, revmap, val)

  pr(0.7, "Testing genes (Wilcoxon)...")
  res <- presto::wilcoxauc(X, labels)
  # Seurat-compatible avg_log2FC on the SAME matrix presto tested. presto's own
  # `logFC` is a natural-log "mean of log" difference (mean(x1) - mean(x2)); Seurat's
  # FindMarkers reports a "log of mean" instead: log2 of the mean of the un-logged
  # normalized counts, with pseudocount 1 added to the group SUM (Seurat v5
  # FoldChange). We surface BOTH: `logFC` (presto, unchanged) and `avg_log2FC`.
  # expm1() recovers the normalized counts and is applied in place (X is not used
  # afterwards); the per-group sums are a sparse mat-vec, so this stays O(nnz)
  # without densifying the matrix.
  n1 <- sum(labels == "group1"); n2 <- length(labels) - n1
  X@x <- expm1(X@x)
  s1 <- as.vector(X %*% as.numeric(labels == "group1"))
  s2 <- as.vector(X %*% as.numeric(labels != "group1"))
  log2fc <- stats::setNames(log2((s1 + 1) / n1) - log2((s2 + 1) / n2), rownames(X))
  res <- res[res$group == "group1", , drop = FALSE]
  res <- res[pmax(res$pct_in, res$pct_out) >= min_pct * 100, , drop = FALSE]
  res <- res[order(res$padj, -abs(res$logFC)), , drop = FALSE]
  data.frame(
    gene = res$feature,
    logFC = round(res$logFC, 3),
    avg_log2FC = round(unname(log2fc[res$feature]), 3),
    auc = round(res$auc, 3),
    pct.1 = round(res$pct_in / 100, 3),
    pct.2 = round(res$pct_out / 100, 3),
    p_val = signif(res$pval, 3),
    p_val_adj = signif(res$padj, 3),
    row.names = NULL, stringsAsFactors = FALSE
  )
}

# Per-cell contrast label aligned to the rows of `cells`: "group1" (a cell whose
# grouping value is in `ident1`), "group2" (in `ident2`), "rest" (one-vs-rest: any
# other annotated cell), or NA for a cell not in the contrast (un-annotated, or
# outside both idents in a two-group contrast). `group_cols` may name several
# metadata columns, which are combined into interaction levels ("a | b"), so the
# idents are combined levels. This is what `scroll_de()` tests over AND what the
# DE panel's live preview colours, so the picture always matches the computation.
# Empty `ident1` => everything NA (nothing selected).
.scroll_contrast_labels <- function(cells, group_cols, ident1, ident2 = NULL) {
  ident1 <- as.character(ident1); ident1 <- ident1[nzchar(ident1)]
  if (!length(ident1)) return(rep(NA_character_, nrow(cells)))
  ident2 <- as.character(ident2); ident2 <- ident2[nzchar(ident2)]
  g <- .scroll_combo_levels(cells, group_cols)   # 1 col => that col; many => "a | b"
  in1 <- !is.na(g) & g %in% ident1
  if (!length(ident2) || "rest" %in% ident2)
    ifelse(in1, "group1", ifelse(!is.na(g), "rest", NA_character_))
  else
    ifelse(in1, "group1", ifelse(!is.na(g) & g %in% ident2, "group2", NA_character_))
}

# Evaluate `expr` with the RNG seeded to `seed` (for reproducible sampling),
# restoring the caller's `.Random.seed` afterward so the draw never perturbs the
# session's RNG stream. `seed = NULL` runs `expr` against the current RNG state
# untouched. Shared by the DE cap and the pseudobulk pseudo-replicate draws.
.scroll_with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  old <- if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
           get(".Random.seed", .GlobalEnv) else NULL
  set.seed(seed)
  on.exit(if (!is.null(old)) assign(".Random.seed", old, .GlobalEnv)
          else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
            rm(".Random.seed", envir = .GlobalEnv))
  force(expr)
}

# Positions (into `labels`) keeping at most `max_cells` cells per group. The
# subsample is deterministic (fixed seed) so repeated DE runs match, and the
# caller's RNG stream is preserved (save + restore `.Random.seed`).
.scroll_cap_groups <- function(labels, max_cells) {
  .scroll_with_seed(1L, {
    sel <- lapply(split(seq_along(labels), labels), function(ii)
      if (length(ii) > max_cells) sample(ii, max_cells) else ii)
    sort(unlist(sel, use.names = FALSE))
  })
}
