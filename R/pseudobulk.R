# Pseudobulk differential expression. Cells are aggregated into sample-level
# count profiles (summed raw counts, via arrow) and compared with the
# edgeR / limma-voom pipeline, so the statistics are replicate-based rather than
# per-cell. Needs a counts store (scroll_build(counts = TRUE)).

# Combined-interaction level for each cell (e.g. "sgIcos | 0"); NA if any of the
# aggregate columns is NA for that cell.
.scroll_combo_levels <- function(cells, cols) {
  parts <- lapply(cols, function(c) as.character(cells[[c]]))
  combo <- do.call(paste, c(parts, sep = " | "))
  combo[Reduce(`|`, lapply(parts, is.na))] <- NA
  combo
}

# Group label ("group1" / "group2" / NA) for each cell given its combined level.
# NA = not in the contrast (its combo is NA, or it is outside both idents in a
# two-group contrast). Shared by scroll_pseudobulk_de() and the panel's live
# preview so the preview matches what the compute aggregates. Empty ident1 => NA.
.scroll_combo_group <- function(combo, ident1, ident2 = NULL) {
  ident1 <- as.character(ident1); ident1 <- ident1[nzchar(ident1)]
  if (!length(ident1)) return(rep(NA_character_, length(combo)))
  ident2 <- as.character(ident2); ident2 <- ident2[nzchar(ident2)]
  in1 <- !is.na(combo) & combo %in% ident1
  if (!length(ident2) || "rest" %in% ident2)
    ifelse(in1, "group1", ifelse(!is.na(combo), "group2", NA_character_))
  else
    ifelse(in1, "group1", ifelse(!is.na(combo) & combo %in% ident2, "group2", NA_character_))
}

.rbind_or_empty <- function(out)
  if (length(out)) do.call(rbind, out) else
    data.frame(cell = character(), psample = character(), group = character())

# Draw n_pseudo pseudo-replicates of cells_per_pseudo cells from a cell vector.
.scroll_pseudo_reps <- function(cellv, tag, group, min_cells, n_pseudo, cells_per_pseudo) {
  if (length(cellv) < min_cells) return(list())
  per <- min(cells_per_pseudo, length(cellv))
  lapply(seq_len(n_pseudo), function(k)
    data.frame(cell = sample(cellv, per), psample = paste0(tag, " :: pseudo", k),
               group = group, stringsAsFactors = FALSE))
}

# Build the cell -> pseudobulk-sample map, tagged group1/group2.
# - replicate_col NULL / "no_replicate": pseudo-replicate the two GROUPS directly
#   (pool each group's cells, sample cells_per_pseudo cells n_pseudo times).
# - a real replicate column: one sample per replicate level per combined level
#   (>= min_cells each, >= 2 of them), with a per-level pseudo fallback otherwise.
# Returns list(mapping = data.frame(cell, psample, group), pseudo = logical).
.scroll_pseudobulk_map <- function(df, replicate_col, min_cells, n_pseudo, cells_per_pseudo) {
  no_rep <- is.null(replicate_col) || !nzchar(replicate_col) ||
    identical(replicate_col, "no_replicate")

  if (no_rep) {                              # pseudo-replicate the two groups
    out <- c(.scroll_pseudo_reps(df$cell[df$grp == "group1"], "group1", "group1",
                                 min_cells, n_pseudo, cells_per_pseudo),
             .scroll_pseudo_reps(df$cell[df$grp == "group2"], "group2", "group2",
                                 min_cells, n_pseudo, cells_per_pseudo))
    return(list(mapping = .rbind_or_empty(out), pseudo = TRUE))
  }

  out <- list(); pseudo <- FALSE             # real replicates, per combined level
  for (L in unique(df$combo)) {
    sub <- df[df$combo == L, , drop = FALSE]; g <- sub$grp[[1]]
    reps <- split(sub$cell, sub$rep)
    reps <- reps[vapply(reps, length, integer(1)) >= min_cells]
    if (length(reps) >= 2) {
      out <- c(out, lapply(names(reps), function(r)
        data.frame(cell = reps[[r]], psample = paste(L, r, sep = " :: "),
                   group = g, stringsAsFactors = FALSE)))
    } else {                                 # this level lacks replicates
      made <- .scroll_pseudo_reps(sub$cell, L, g, min_cells, n_pseudo, cells_per_pseudo)
      if (length(made)) pseudo <- TRUE
      out <- c(out, made)
    }
  }
  list(mapping = .rbind_or_empty(out), pseudo = pseudo)
}

#' Pseudobulk differential expression (edgeR / limma-voom)
#'
#' Aggregates cells into pseudobulk samples by summing raw counts, then compares
#' two groups with the limma-voom pipeline. Groups are built from
#' **combined-interaction levels**: the observed combinations of `aggregate_cols`
#' (e.g. `"sgIcos | 0"`); `ident1`/`ident2` select which combinations form each
#' group (either may name several). Pseudobulk samples come from `replicate_col`
#' when real replicates exist, otherwise from random pseudo-replicates.
#'
#' Requires a counts store — build with [scroll_build()] `counts = TRUE`.
#'
#' @param data A scroll data handle (`con`, `cells`, `manifest`).
#' @param assay Assay to test.
#' @param aggregate_cols Metadata columns whose combinations define the groups.
#' @param ident1 Combined levels forming the group of interest.
#' @param ident2 Combined levels for the comparison group, or `NULL`/`"rest"` for
#'   one-vs-rest (all cells not in `ident1`, restricted to observed combinations).
#' @param replicate_col Metadata column defining biological replicates. `NULL` or
#'   `"no_replicate"` pseudo-replicates the two groups directly (pooling each
#'   group's cells) instead of using a real replicate column.
#' @param min_cells Minimum cells for a pseudobulk sample to be kept.
#' @param n_pseudo,cells_per_pseudo Pseudo-replicate fallback: number of
#'   pseudo-replicates and cells sampled per pseudo-replicate.
#' @param cells Cells to test over (default `data$cells`).
#' @param seed Integer seed for the pseudo-replicate draw, so a compute is
#'   reproducible and the caller's RNG stream is left untouched. `NULL` uses the
#'   current RNG state (used internally to vary draws across stability runs).
#' @param progress Optional `function(fraction, detail)` called at the aggregate
#'   and model-fit stages, for a UI progress bar. `NULL` (default) is a no-op.
#' @return A data.frame (`gene`, `logFC`, `avg_expr`, `p_val`, `p_val_adj`) ranked
#'   by adjusted p-value; positive `logFC` is up in `ident1`. The attribute
#'   `"pseudo"` is `TRUE` when pseudo-replicates were used.
#' @export
scroll_pseudobulk_de <- function(data, assay, aggregate_cols, ident1, ident2 = NULL,
                                 replicate_col = NULL, min_cells = 10,
                                 n_pseudo = 3, cells_per_pseudo = 50, cells = NULL,
                                 seed = 1L, progress = NULL) {
  pr <- function(f, d) if (is.function(progress)) progress(f, d)
  if (!requireNamespace("edgeR", quietly = TRUE) ||
      !requireNamespace("limma", quietly = TRUE))
    stop("Pseudobulk DE needs the 'edgeR' and 'limma' packages.", call. = FALSE)
  if (!isTRUE(data$manifest$has_counts))
    stop("This project has no counts store; rebuild with scroll_build(counts = TRUE).",
         call. = FALSE)
  cells <- cells %||% data$cells
  aggregate_cols <- as.character(aggregate_cols)
  if (!length(aggregate_cols)) stop("Pick at least one aggregate-by column.", call. = FALSE)
  ident1 <- as.character(ident1); ident1 <- ident1[nzchar(ident1)]
  ident2 <- as.character(ident2); ident2 <- ident2[nzchar(ident2)]
  if (!length(ident1)) stop("Pick at least one Group 1 level.", call. = FALSE)

  combo <- .scroll_combo_levels(cells, aggregate_cols)
  grp <- .scroll_combo_group(combo, ident1, ident2)   # "group1"/"group2"/NA per cell

  rep_vals <- if (!is.null(replicate_col) && replicate_col %in% names(cells))
                as.character(cells[[replicate_col]]) else NA_character_
  # `cell` here is the store's cell key that scroll_aggregate_counts joins on:
  # the int32 global index for a v2 store, the barcode string for v1.
  cell_key <- if (isTRUE(data$manifest$cell_index))
                (if (!is.null(cells$.gidx)) cells$.gidx else seq_len(nrow(cells))) else cells$cell
  df <- data.frame(cell = cell_key, combo = combo, grp = grp, rep = rep_vals,
                   stringsAsFactors = FALSE)
  df <- df[!is.na(df$combo) & !is.na(df$grp), , drop = FALSE]
  if (!nrow(df)) stop("No cells match the selected groups.", call. = FALSE)

  # seed the pseudo-replicate draw (reproducible) without perturbing the session
  # RNG; stability passes seed = NULL / a per-run seed to vary the draws.
  pm <- .scroll_with_seed(seed,
    .scroll_pseudobulk_map(df, replicate_col, min_cells, n_pseudo, cells_per_pseudo))
  mapping <- pm$mapping
  samp <- unique(mapping[, c("psample", "group")])
  if (sum(samp$group == "group1") < 2 || sum(samp$group == "group2") < 2)
    stop("Each group needs >= 2 pseudobulk samples after the min-cell filter. ",
         "Lower Min cells, choose a replicate column, or add pseudo-replicates.",
         call. = FALSE)

  pr(0.2, "Aggregating pseudobulk counts...")
  agg <- scroll_aggregate_counts(data$con, assay, mapping[, c("cell", "psample")])
  # store's own feature names (encoding-robust; the manifest yaml can mangle a
  # non-ASCII feature name into an escaped form that would not match() the store).
  feats <- sort(unique(agg$feature))
  ps <- samp$psample
  M <- matrix(0L, nrow = length(feats), ncol = length(ps), dimnames = list(feats, ps))
  M[cbind(match(agg$feature, feats), match(agg$psample, ps))] <- as.integer(agg$count)

  pr(0.7, "Fitting model (edgeR / limma-voom)...")
  # group1 as the tested coefficient: positive logFC = up in ident1
  group <- factor(samp$group[match(ps, samp$psample)], levels = c("group2", "group1"))
  dge <- edgeR::DGEList(counts = M, group = group)
  keep <- edgeR::filterByExpr(dge, group = group)
  dge <- edgeR::calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
  design <- stats::model.matrix(~ group)
  fit <- limma::eBayes(limma::lmFit(limma::voom(dge, design), design))
  tt <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "P")

  res <- data.frame(
    gene = rownames(tt),
    logFC = round(tt$logFC, 3),
    avg_expr = round(tt$AveExpr, 3),
    p_val = signif(tt$P.Value, 3),
    p_val_adj = signif(tt$adj.P.Val, 3),
    row.names = NULL, stringsAsFactors = FALSE)
  attr(res, "pseudo") <- pm$pseudo
  attr(res, "n_samples") <- c(group1 = sum(samp$group == "group1"),
                              group2 = sum(samp$group == "group2"))
  res
}

# Run scroll_pseudobulk_de `runs` times (fresh pseudo-replicates each) and return
# the successful per-run result data.frames. A distinct per-run seed keeps the
# draws varied AND the whole stability computation reproducible, and (like a
# single compute) leaves the caller's RNG stream untouched.
.scroll_pseudobulk_runs <- function(..., runs) {
  out <- lapply(seq_len(runs),
                function(i) tryCatch(scroll_pseudobulk_de(..., seed = i),
                                     error = function(e) NULL))
  out <- Filter(Negate(is.null), out)
  if (length(out) < 2)
    stop("Stability needs >= 2 successful runs (", length(out), " succeeded); ",
         "check the contrast and min-cell settings.", call. = FALSE)
  out
}

# Summarise a list of per-run DE results into per-gene stability statistics: how
# often each gene clears the lfc + padj cutoffs (selection frequency), its median
# effect, and how consistently the effect points the same way.
.scroll_stability_aggregate <- function(runs, lfc = 1, padj = 0.05) {
  k <- length(runs)
  all <- do.call(rbind, runs)
  sig <- abs(all$logFC) >= lfc & all$p_val_adj <= padj
  ix <- split(seq_len(nrow(all)), all$gene)
  lf <- all$logFC; pj <- all$p_val_adj
  res <- data.frame(
    gene = names(ix),
    sel_freq = round(vapply(ix, function(j) sum(sig[j]) / k, numeric(1)), 3),
    median_logFC = round(vapply(ix, function(j) stats::median(lf[j]), numeric(1)), 3),
    sign_agree = round(vapply(ix, function(j) {
      s <- sign(lf[j]); max(mean(s > 0), mean(s < 0)) }, numeric(1)), 2),
    median_padj = signif(vapply(ix, function(j) stats::median(pj[j]), numeric(1)), 3),
    n_tested = vapply(ix, length, integer(1)),
    row.names = NULL, stringsAsFactors = FALSE)
  res <- res[order(-res$sel_freq, res$median_padj), , drop = FALSE]
  attr(res, "runs") <- k
  # aliases so view_volcano / the shared table still work
  res$logFC <- res$median_logFC; res$p_val_adj <- res$median_padj
  res
}

#' Stability of pseudobulk DE across random pseudo-replicate draws
#'
#' Runs [scroll_pseudobulk_de()] `runs` times with fresh random pseudo-replicates
#' and reports, per gene, how often it clears the `lfc` + `padj` cutoffs
#' (`sel_freq`), its `median_logFC`, sign agreement, and median adjusted p. Use it
#' to see which genes come up **consistently** rather than by a lucky draw.
#'
#' Note: this checks robustness to the random sampling only; it does **not** fix
#' the anti-conservative bias of pseudo-replication (fabricated replicates
#' understate biological variance). Treat `sel_freq` as a ranking heuristic, not a
#' p-value, and prefer real replicates when available.
#'
#' @inheritParams scroll_pseudobulk_de
#' @param runs Number of random re-runs.
#' @param lfc,padj logFC and adjusted-p cutoffs defining a "hit" in each run.
#' @return A data.frame with `gene`, `sel_freq`, `median_logFC`, `sign_agree`,
#'   `median_padj`, `n_tested`, ranked by selection frequency.
#' @export
scroll_pseudobulk_stability <- function(data, assay, aggregate_cols, ident1, ident2 = NULL,
                                        replicate_col = "no_replicate", min_cells = 10,
                                        n_pseudo = 3, cells_per_pseudo = 50, cells = NULL,
                                        runs = 25, lfc = 1, padj = 0.05) {
  runs_list <- .scroll_pseudobulk_runs(
    data, assay, aggregate_cols, ident1, ident2, replicate_col = replicate_col,
    min_cells = min_cells, n_pseudo = n_pseudo, cells_per_pseudo = cells_per_pseudo,
    cells = cells, runs = runs)
  .scroll_stability_aggregate(runs_list, lfc, padj)
}
