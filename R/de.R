# Live differential expression. presto::wilcoxauc needs an in-memory
# genes x cells matrix, so a DE run reconstructs the contrast's expression
# matrix from the duckdb store (the one runtime operation that materialises a
# chunk of the matrix — hence on-demand, behind a Compute button).

#' Live differential expression for a contrast (presto / Wilcoxon)
#'
#' Compares `ident1` against `ident2` (another level) or against the rest of the
#' cells when `ident2` is `NULL`/`"rest"`. Reconstructs the contrast's expression
#' matrix from the store and runs [presto::wilcoxauc()].
#'
#' @param data A scroll data handle (`con`, `cells`, `manifest`).
#' @param assay Assay to test.
#' @param group_col Metadata column defining the groups.
#' @param ident1 The group of interest.
#' @param ident2 The comparison group, or `NULL`/`"rest"` for one-vs-rest.
#' @param min_pct Keep genes expressed in at least this fraction of either side.
#' @return A data.frame of results, ranked by adjusted p-value.
#' @details With the default quantized build (`quantize = TRUE` in
#'   [scroll_build()]), expression values below ~`max/510` round to zero, so the
#'   fraction-expressing columns (`pct.1`/`pct.2`) slightly under-count cells with
#'   very low expression and p-values are approximate. Build with
#'   `quantize = FALSE` for exact statistics.
#' @export
scroll_de <- function(data, assay, group_col, ident1, ident2 = NULL, min_pct = 0.1) {
  if (!requireNamespace("presto", quietly = TRUE))
    stop("Live DE needs the 'presto' package.", call. = FALSE)
  cells <- data$cells
  g <- as.character(cells[[group_col]])

  one_vs_rest <- is.null(ident2) || !nzchar(ident2) || identical(ident2, "rest")
  if (one_vs_rest) {
    keep <- !is.na(g)                       # drop un-annotated cells (else NA labels)
    labels <- ifelse(!is.na(g) & g == ident1, ident1, "rest")
  } else {
    keep <- !is.na(g) & g %in% c(ident1, ident2)
    labels <- g
  }
  ccells <- cells$cell[keep]
  labels <- labels[keep]
  if (sum(labels == ident1) < 3 || sum(labels != ident1) < 3)
    stop("Each side of the contrast needs at least 3 cells.", call. = FALSE)

  long <- scroll_query_cells(data$con, assay, ccells)
  long$value <- scroll_dequantize(long$value, data$manifest, assay)
  feats <- .scroll_features_of(data$manifest, assay)
  X <- Matrix::sparseMatrix(
    i = match(long$feature, feats), j = match(long$cell, ccells), x = long$value,
    dims = c(length(feats), length(ccells)), dimnames = list(feats, ccells))

  res <- presto::wilcoxauc(X, labels)
  res <- res[res$group == ident1, , drop = FALSE]
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
