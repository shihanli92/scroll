# --- Gene-signature scoring ---------------------------------------------------
# Runtime module scoring from a gene list. `mean` / `scaled` are bounded to the
# signature's genes (a `queryN` read). `addmodulescore` reproduces Seurat's
# AddModuleScore -- the mean of the signature minus the mean of an expression-matched
# control set -- which needs each gene's average expression across all cells. That
# average is computed once by a full-store grouped-sum scan and cached on the connection
# handle (the same class of on-demand full scan as live DE, paid once, not per interaction).

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

# Equal-frequency binning: a base-R stand-in for ggplot2::cut_number() (which
# AddModuleScore uses) -- n bins by quantiles of x, one integer bin per element, names
# preserved. Collapsed (tied) breaks reduce the effective bin count, as expected when
# many genes share a mean (e.g. all-zero genes).
.scroll_cut_number <- function(x, n) {
  br <- unique(stats::quantile(x, probs = seq(0, 1, length.out = n + 1L),
                               na.rm = TRUE, names = FALSE))
  b <- if (length(br) < 2L) rep(1L, length(x))
       else as.integer(cut(x, breaks = br, include.lowest = TRUE, labels = FALSE))
  stats::setNames(b, names(x))
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

# Seurat-style AddModuleScore at runtime: mean(signature) - mean(control), where the
# control set is `ctrl` genes drawn per signature gene from the same expression bin
# (seed-stable). Binning + control genes come from the full-data gene means, and the two
# per-cell means are per cell, so a cell's score is independent of the active subset
# (filter-invariant), matching Seurat.
.scroll_module_score <- function(cells, genes, data, assay = NULL,
                                 nbin = 24L, ctrl = 100L, seed = 1L, progress = NULL) {
  pr <- function(f, d) if (is.function(progress)) progress(f, d)
  assay <- assay %||% data$manifest$default_assay
  genes <- intersect(genes, .scroll_features_of(data$manifest, assay))
  if (!length(genes) || !nrow(cells)) return(rep(0, nrow(cells)))
  # set the message BEFORE each slow step so the bar reflects what is running.
  pr(0.10, "Scanning gene background (one-off)...")       # the full-store mean scan (cached after)
  avg  <- .scroll_gene_means(data$con, data$manifest, assay)
  pr(0.55, "Sampling control genes...")
  bins <- .scroll_cut_number(avg, min(nbin, length(unique(avg))))
  # seed the control sampling without disturbing the app's global RNG stream
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    old <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
  }
  set.seed(seed)
  ctrl.use <- character(0)
  for (g in genes) {
    pool <- names(bins)[bins == bins[[g]]]
    ctrl.use <- c(ctrl.use, sample(pool, min(ctrl, length(pool))))
  }
  ctrl.use <- unique(ctrl.use)
  pr(0.70, "Querying signature expression...")
  sig <- .scroll_signature_score(cells, genes,    data$queryN(assay, genes),    "mean")
  pr(0.85, "Querying control expression...")
  ctl <- .scroll_signature_score(cells, ctrl.use, data$queryN(assay, ctrl.use), "mean")
  pr(1, "")
  sig - ctl
}
