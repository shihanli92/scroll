# Live differential expression. presto::wilcoxauc needs an in-memory
# genes x cells matrix, so a DE run reconstructs the contrast's expression
# matrix from the Parquet store (the one runtime operation that materialises a
# chunk of the matrix — hence on-demand, behind a Compute button).

#' Live differential expression for a contrast (presto / Wilcoxon)
#'
#' Compares `ident1` against `ident2` or against the rest of the cells when
#' `ident2` is `NULL`/`"rest"`. Either side may name **several** levels of
#' `group_col` (they are pooled into one group). Reconstructs the contrast's
#' expression matrix from the store and runs [presto::wilcoxauc()].
#'
#' @param data A scroll data handle (`con`, `cells`, `manifest`).
#' @param assay Assay to test.
#' @param group_col Metadata column defining the groups.
#' @param ident1 One or more levels of `group_col` forming the group of interest.
#' @param ident2 One or more comparison levels, or `NULL`/`"rest"` for
#'   one-vs-rest (every cell not in `ident1`).
#' @param min_pct Keep genes expressed in at least this fraction of either side.
#' @param cells Cells to test over (default `data$cells`); pass a subset to run
#'   DE within an active cell-subset filter.
#' @return A data.frame of results, ranked by adjusted p-value.
#' @details With the default quantized build (`quantize = TRUE` in
#'   [scroll_build()]), expression values below ~`max/510` round to zero, so the
#'   fraction-expressing columns (`pct.1`/`pct.2`) slightly under-count cells with
#'   very low expression and p-values are approximate. Build with
#'   `quantize = FALSE` for exact statistics.
#' @export
scroll_de <- function(data, assay, group_col, ident1, ident2 = NULL, min_pct = 0.1,
                      cells = NULL) {
  if (!requireNamespace("presto", quietly = TRUE))
    stop("Live DE needs the 'presto' package.", call. = FALSE)
  cells <- cells %||% data$cells
  g <- as.character(cells[[group_col]])
  ident1 <- as.character(ident1); ident1 <- ident1[nzchar(ident1)]
  ident2 <- as.character(ident2); ident2 <- ident2[nzchar(ident2)]
  if (!length(ident1)) stop("Pick at least one level for group 1.", call. = FALSE)

  in1 <- !is.na(g) & g %in% ident1
  one_vs_rest <- !length(ident2) || "rest" %in% ident2
  if (one_vs_rest) {
    keep <- !is.na(g)                       # drop un-annotated cells (else NA labels)
    labels <- ifelse(in1, "group1", "rest")
  } else {
    in2 <- !is.na(g) & g %in% ident2
    keep <- in1 | in2                       # in1 wins if a level is in both
    labels <- ifelse(in1, "group1", "group2")
  }
  # `bc` = barcodes (matrix dimnames); `key` = the store's cell key used to query +
  # join (v2 int global index / v1 barcode). `.gidx` rides on the cells handle.
  gidx <- if (!is.null(cells$.gidx)) cells$.gidx else seq_len(nrow(cells))
  bc  <- cells$cell[keep]
  key <- if (isTRUE(data$manifest$cell_index)) gidx[keep] else bc
  labels <- labels[keep]
  if (sum(labels == "group1") < 3 || sum(labels != "group1") < 3)
    stop("Each side of the contrast needs at least 3 cells.", call. = FALSE)

  # dict = TRUE returns `feature` dictionary-encoded (a factor): its levels are the
  # store's OWN feature names and as.integer() gives the matrix row index directly.
  # This avoids sort(unique()) + match() over the tens of millions of repeated
  # gene-name strings (the dominant cost of a large DE) — arrow builds the
  # dictionary during the scan. Indexing by the store's names (not the manifest
  # list) also stays correct when yaml has mangled a non-ASCII feature name, and
  # zero-expression genes are absent here so restricting to expressed features is
  # correct too.
  long <- scroll_query_cells(data$con, assay, key, dict = TRUE)
  long$value <- scroll_dequantize(long$value, data$manifest, assay)
  feats <- levels(long$feature)
  X <- Matrix::sparseMatrix(
    i = as.integer(long$feature), j = match(long$cell, key), x = long$value,
    dims = c(length(feats), length(bc)), dimnames = list(feats, bc))

  res <- presto::wilcoxauc(X, labels)
  res <- res[res$group == "group1", , drop = FALSE]
  res <- res[pmax(res$pct_in, res$pct_out) >= min_pct * 100, , drop = FALSE]
  res <- res[order(res$padj, -abs(res$logFC)), , drop = FALSE]
  data.frame(
    gene = res$feature,
    logFC = round(res$logFC, 3),
    auc = round(res$auc, 3),
    pct.1 = round(res$pct_in / 100, 3),
    pct.2 = round(res$pct_out / 100, 3),
    p_val = signif(res$pval, 3),
    p_val_adj = signif(res$padj, 3),
    row.names = NULL, stringsAsFactors = FALSE
  )
}
