# --- Gene-signature scoring ---------------------------------------------------
# Runtime module scoring from a gene list, bounded to the signature's genes (plus
# AddModuleScore's control genes). Methods:
#   mean / scaled    mean of the genes' log-norm expression / of per-gene z-scores
#   addmodulescore   Seurat's AddModuleScore: signature mean minus the mean of an
#                    expression-matched control set, drawn exactly as Seurat draws it
#                    from the build-time gene means (stats/<assay>/genes.parquet)
#   ucell / aucell   rank-based scores from the within-cell ranks baked with
#                    scroll_build(ranks = TRUE) (the expr `rank2` column): UCell's
#                    Mann-Whitney U, and AUCell's area under the recovery curve (ties
#                    averaged rather than random)
# Projects built before the gene-mean bake fall back: AddModuleScore scans the store
# for gene means once (close to, but not exactly, Seurat). Without ranks, the
# rank-based scores are not offered (and stop with a rebuild hint if called).

# Per-cell score aligned to `cells` for the bounded methods. `long` = data$queryN(assay,
# genes) (dequantized sparse long); absent cell => 0. "scaled" z-scores each gene across
# the shown cells before averaging.
.scroll_signature_score <- function(cells, genes, long, method = c("mean", "scaled")) {
  method <- match.arg(method)
  n <- nrow(cells)
  if (!length(genes) || n == 0) return(rep(0, n))
  mat <- vapply(genes, function(g) {
    v <- if (!is.null(long) && nrow(long))
           long[long$feature == g, c("cell", "value"), drop = FALSE] else NULL
    .scroll_expr_vector(cells, v)                       # dense, zero-filled, aligned to `cells`
  }, numeric(n))
  if (!is.matrix(mat)) mat <- matrix(mat, nrow = n)     # n == 1 guard
  if (identical(method, "scaled")) {                    # z per gene over the shown cells
    mat <- apply(mat, 2, function(x) {
      s <- stats::sd(x); if (is.na(s) || s == 0) x * 0 else (x - mean(x)) / s
    })
    if (!is.matrix(mat)) mat <- matrix(mat, nrow = n)
  }
  rowMeans(mat)
}

# Per-gene mean expression across ALL cells (dequantized), for an assay. One grouped-sum
# pass over the whole expr store, memoized on the connection handle so it is paid once.
# Dequantization is linear (and 0 -> 0), so the sum of stored values dequantizes to the
# sum of normalized values; dividing by n_cells (zeros implicit) gives the mean.
.scroll_gene_means <- function(con, manifest, assay) {
  if (is.null(con$gene_means)) con$gene_means <- new.env(parent = emptyenv())
  if (!is.null(con$gene_means[[assay]])) return(con$gene_means[[assay]])
  ds  <- .scroll_dataset(con, assay)
  agg <- dplyr::collect(dplyr::summarise(dplyr::group_by(ds, .data$feature),
                                         s = sum(.data$value, na.rm = TRUE)))
  s   <- scroll_dequantize(as.numeric(agg$s), manifest, assay) / manifest$n_cells
  feats <- .scroll_features_of(manifest, assay)
  full  <- stats::setNames(numeric(length(feats)), feats)        # absent (all-zero) genes -> 0
  full[agg$feature] <- s
  con$gene_means[[assay]] <- full
  full
}

# ---- build-time statistics --------------------------------------------------

# An assay's stats/<assay>/<file> table, read once per connection (NULL when the
# project was built without it).
.scroll_stats_table <- function(con, manifest, assay, file) {
  if (!isTRUE(manifest$assays[[assay]]$stats)) return(NULL)
  if (is.null(con$stats)) con$stats <- new.env(parent = emptyenv())
  key <- paste(assay, file)
  if (is.null(con$stats[[key]])) {
    .scroll_assay_dir(con$dir, assay)                   # validates the assay name
    path <- file.path(con$dir, "stats", assay, file)
    if (!file.exists(path)) return(NULL)
    con$stats[[key]] <- as.data.frame(arrow::read_parquet(path, mmap = FALSE))
  }
  con$stats[[key]]
}

# Per-gene means over all cells, in the assay's row order (named). From the build-time
# bake when present (exact: Matrix::rowMeans of the original data); else the one-off
# store scan above.
.scroll_signature_gene_means <- function(con, manifest, assay) {
  g <- .scroll_stats_table(con, manifest, assay, "genes.parquet")
  if (!is.null(g)) return(structure(stats::setNames(g$mean, g$feature), exact = TRUE))
  structure(.scroll_gene_means(con, manifest, assay), exact = FALSE)
}

# ---- per-cell aggregates, summed inside arrow ----------------------------------

# Per-cell mean of `features` (dequantized; absent = 0), aligned to `cells`. The sum
# runs in arrow's engine, so a large control set never materializes per-gene rows.
.scroll_cell_mean <- function(cells, features, data, assay) {
  if (!length(features) || !nrow(cells)) return(rep(0, nrow(cells)))
  ds <- .scroll_dataset(data$con, assay)
  agg <- dplyr::collect(dplyr::summarise(
    dplyr::group_by(dplyr::filter(ds, .data$feature %in% !!features), .data$cell),
    value = sum(.data$value, na.rm = TRUE)))
  sums <- data.frame(cell = agg$cell,
                     value = scroll_dequantize(as.numeric(agg$value), data$manifest, assay))
  .scroll_expr_vector(cells, sums) / length(features)
}

# ---- AddModuleScore ------------------------------------------------------------

# The control genes Seurat's AddModuleScore draws for `genes` (its Assay method, step
# for step): sort the gene means, add a tiny seeded jitter, cut into `nbin`
# equal-count bins, then sample `ctrl` genes from each signature gene's bin, in the
# signature's order. With identical means, this is the same set Seurat uses. Only
# where Seurat itself would error (a bin smaller than `ctrl`) does it take the whole
# bin. The app's global RNG stream is left untouched.
.scroll_seurat_controls <- function(avg, genes, nbin = 24L, ctrl = 100L, seed = 1L) {
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    old <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
  } else on.exit(rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  if (!is.null(seed)) set.seed(seed)
  avg <- avg[order(avg)]
  cut <- tryCatch(
    ggplot2::cut_number(avg + stats::rnorm(length(avg)) / 1e30, n = nbin,
                        labels = FALSE, right = FALSE),
    error = function(e) stop(sprintf(
      "Too few distinct gene means to form %d expression bins; lower the bin count.", nbin),
      call. = FALSE))
  names(cut) <- names(avg)
  ctrl.use <- character(0)
  for (g in genes) {
    pool <- cut[which(cut == cut[[g]])]
    ctrl.use <- c(ctrl.use,
                  if (length(pool) <= 1L) names(pool)
                  else names(sample(pool, size = min(ctrl, length(pool)), replace = FALSE)))
  }
  unique(ctrl.use)
}

# Seurat's AddModuleScore: mean(signature) - mean(controls), per cell (so a cell's
# score does not depend on the active subset, as in Seurat). Exact when the project
# carries build-time gene means; the result's "exact" attribute says which.
.scroll_module_score <- function(cells, genes, data, assay = NULL,
                                 nbin = 24L, ctrl = 100L, seed = 1L, progress = NULL) {
  pr <- function(f, d) if (is.function(progress)) progress(f, d)
  assay <- assay %||% data$manifest$default_assay
  genes <- intersect(genes, .scroll_features_of(data$manifest, assay))
  if (!length(genes) || !nrow(cells)) return(rep(0, nrow(cells)))
  exact <- isTRUE(data$manifest$assays[[assay]]$stats)
  # set the message BEFORE each slow step so the bar reflects what is running.
  pr(0.10, if (exact) "Reading gene means..." else "Scanning gene background (one-off)...")
  avg <- .scroll_signature_gene_means(data$con, data$manifest, assay)
  pr(0.40, "Sampling control genes...")
  ctrl.use <- .scroll_seurat_controls(avg, genes, nbin, ctrl, seed)
  pr(0.55, "Scoring the signature...")
  sig <- .scroll_cell_mean(cells, genes, data, assay)
  pr(0.80, "Scoring the controls...")
  ctl <- .scroll_cell_mean(cells, ctrl.use, data, assay)
  pr(1, "")
  structure(sig - ctl, exact = isTRUE(attr(avg, "exact")))
}

# ---- rank-based scores (UCell / AUCell) ------------------------------------------

.scroll_has_ranks <- function(manifest, assay) isTRUE(manifest$assays[[assay]]$ranks)

# The signature genes' within-cell ranks as a cells x genes matrix (ranks, not x2):
# stored ranks from the expr `rank2` column; a gene a cell doesn't express takes that
# cell's zero-tie rank, from the per-cell positive / negative counts.
.scroll_rank_matrix <- function(cells, genes, data, assay) {
  m <- data$manifest
  if (!.scroll_has_ranks(m, assay))
    stop("Rank-based scores need within-cell ranks: rebuild the project with ",
         "scroll_build(..., ranks = TRUE) (scroll 0.3.7 or later).", call. = FALSE)
  G <- m$assays[[assay]]$n_features
  ds <- .scroll_dataset(data$con, assay)
  long <- dplyr::collect(dplyr::select(dplyr::filter(ds, .data$feature %in% !!genes),
                                       "feature", "cell", "rank2"))
  cs <- .scroll_stats_table(data$con, m, assay, "cells.parquet")
  key <- .scroll_cell_key(cells, data.frame(cell = integer(0)))
  i <- match(key, cs$cell)
  zero <- (2 * cs$n_pos[i] + (G - cs$n_pos[i] - cs$n_neg[i]) + 1) / 2
  zero[is.na(zero)] <- (G + 1) / 2                      # a cell with nothing stored
  vapply(genes, function(g) {
    sub <- long[long$feature == g, , drop = FALSE]
    r <- sub$rank2[match(key, sub$cell)] / 2
    ifelse(is.na(r), zero, r)
  }, numeric(nrow(cells)))
}

# UCell (Andreatta & Carmona 2021): 1 - U / (n * max_rank), with ranks at or beyond
# `max_rank` counted as `max_rank` (as UCell::ScoreSignatures_UCell does).
.scroll_ucell_score <- function(cells, genes, data, assay = NULL, max_rank = 1500) {
  assay <- assay %||% data$manifest$default_assay
  genes <- intersect(genes, .scroll_features_of(data$manifest, assay))
  if (!length(genes) || !nrow(cells)) return(rep(0, nrow(cells)))
  max_rank <- min(max_rank, data$manifest$assays[[assay]]$n_features)
  if (length(genes) > max_rank) stop("The signature has more genes than the max rank.", call. = FALSE)
  r <- pmin(.scroll_rank_matrix(cells, genes, data, assay), max_rank)
  if (!is.matrix(r)) r <- matrix(r, nrow = nrow(cells))
  n <- length(genes); rmin <- n * (n + 1) / 2
  1 - (rowSums(r) - rmin) / (n * max_rank - rmin)
}

# AUCell-style (Aibar et al. 2017): the area under each cell's recovery curve over its
# top `auc_max_rank` genes, normalized by the best possible area. The curve's area
# telescopes to sum(auc_max_rank - rank) over signature genes ranked below the
# threshold, so no per-cell sort is needed. AUCell breaks rank ties at random; ties
# here are averaged, so values can differ from AUCell's slightly.
.scroll_aucell_score <- function(cells, genes, data, assay = NULL, auc_max_rank = NULL) {
  assay <- assay %||% data$manifest$default_assay
  genes <- intersect(genes, .scroll_features_of(data$manifest, assay))
  if (!length(genes) || !nrow(cells)) return(rep(0, nrow(cells)))
  thr <- auc_max_rank %||% ceiling(0.05 * data$manifest$assays[[assay]]$n_features)
  r <- .scroll_rank_matrix(cells, genes, data, assay)
  if (!is.matrix(r)) r <- matrix(r, nrow = nrow(cells))
  best <- seq_len(length(genes)); best <- best[best < thr]
  max_auc <- sum(thr - best)
  if (max_auc <= 0) return(rep(0, nrow(cells)))
  rowSums(pmax(thr - r, 0)) / max_auc
}

# ---- up / down signatures ---------------------------------------------------------

# Combine a signature's up- and down-gene scores (each from the same method): UCell's
# own rule -- up minus w_neg x down, floored at 0 -- for UCell (so it matches
# UCell::ScoreSignatures_UCell with "gene-" entries), and plain up minus w_neg x down
# for the others (AddModuleScore, AUCell, Mean and Scaled have no down genes of their
# own; scoring each list separately and subtracting is the usual way to add them).
.scroll_signed_score <- function(up, down = NULL, method = "mean", w_neg = 1) {
  if (is.null(down)) return(up)
  s <- up - w_neg * down
  if (identical(method, "ucell")) s <- pmax(s, 0)
  s
}

