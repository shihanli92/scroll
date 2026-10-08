# The Signature panel: live per-cell module scoring from a gene list, shown on the
# embedding or as a violin (R/signature.R, R/panel-signature.R).

test_that(".scroll_signature_score computes mean and scaled correctly", {
  cells <- data.frame(cell = c("c1", "c2", "c3", "c4"), .gidx = 1:4,
                      stringsAsFactors = FALSE)
  # v2-style store: `long$cell` is the integer global row index (== cells$.gidx).
  long <- data.frame(
    feature = c("A", "A", "B"),
    cell    = c(1L, 2L, 1L),          # gene A nonzero in cells 1,2; gene B in cell 1 only
    value   = c(2,  4,  6), stringsAsFactors = FALSE)

  # Mean: per-gene dense vectors A=[2,4,0,0], B=[6,0,0,0]; rowMeans -> [4,2,0,0].
  # Absent cells (3,4) are zero (sparse convention).
  expect_equal(scroll:::.scroll_signature_score(cells, c("A", "B"), long, "mean"),
               c(4, 2, 0, 0))

  # Scaled: each gene z-scored across the shown cells, then averaged. Per-gene column
  # means are ~0, so the score is centred; check a couple of structural properties.
  sc <- scroll:::.scroll_signature_score(cells, c("A", "B"), long, "scaled")
  expect_length(sc, 4)
  expect_equal(mean(sc), 0, tolerance = 1e-8)          # average of two centred z columns
  expect_true(which.max(sc) == 1)                       # cell 1 (high in both) scores highest

  # An all-zero (or single-value) gene has sd 0 -> contributes 0, never NaN.
  expect_equal(scroll:::.scroll_signature_score(cells, "C", long, "scaled"), c(0, 0, 0, 0))
  # a barcode-keyed (v1) long resolves the same way
  long_v1 <- data.frame(feature = "A", cell = c("c1", "c2"), value = c(2, 4),
                        stringsAsFactors = FALSE)
  expect_equal(scroll:::.scroll_signature_score(cells, "A", long_v1, "mean"), c(2, 4, 0, 0))
})

test_that("the signature panel is a built-in, placed after FeaturePlot", {
  ids <- vapply(scroll:::.scroll_builtin_panels(), function(p) p$id, character(1))
  expect_true("signature" %in% ids)
  expect_equal(which(ids == "signature"), which(ids == "featureplot") + 1L)
})

test_that("signature_server renders the UMAP and violin views from a gene list", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  feats <- scroll:::.scroll_features_of(data$manifest, data$manifest$default_assay)
  red   <- scroll:::.scroll_view_embeddings(data$manifest, NULL)[[1]]
  cat1  <- scroll:::.scroll_cat_cols(data$manifest)[[1]]

  panel <- Filter(function(p) identical(p$id, "signature"), scroll:::.scroll_builtin_panels())[[1]]
  shiny::testServer(panel$server, args = list(data = data), {
    session$setInputs(sig = feats[1:3], method = "mean", view = "umap", reduction = red,
                      palette = "grey-purple", size = 0.7, clip = c(0, 100),
                      order = TRUE, legend = TRUE, raster = FALSE, aspect = 1)
    expect_error(plot_r())                                 # nothing until Calculate is clicked
    session$setInputs(compute = 1)
    expect_s3_class(plot_r(), "ggplot")

    # CSV export is one score per cell
    csv <- csv_r()
    expect_setequal(names(csv), c("cell", "signature_score"))
    expect_equal(nrow(csv), nrow(data$cells))

    # switching to the violin view reuses the computed score (no recompute needed)
    session$setInputs(view = "violin", group = cat1)
    expect_s3_class(plot_r(), "ggplot")

    # scaled scoring after a re-Calculate
    session$setInputs(view = "umap", method = "scaled", compute = 2)
    expect_s3_class(plot_r(), "ggplot")
  })
})

test_that(".scroll_gene_means scans the store and matches the object's per-gene means", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  means <- scroll:::.scroll_gene_means(data$con, data$manifest, "RNA")
  feats <- scroll:::.scroll_features_of(data$manifest, "RNA")
  expect_setequal(names(means), feats)
  # self-consistency: the grouped-sum scan equals per-gene query aggregation
  for (g in c("CD3D", "CD8A", "NKG7")) {
    q <- data$query1("RNA", g)                            # dequantized (cell, value)
    expect_equal(unname(means[[g]]), sum(q$value) / data$manifest$n_cells, tolerance = 1e-6)
  }
  # and it tracks the source object's true log-norm means (within quantization)
  obj <- make_test_object()
  truth <- Matrix::rowMeans(SeuratObject::GetAssayData(obj, layer = "data"))
  expect_equal(means[names(truth)], truth[names(truth)], tolerance = 0.02, ignore_attr = TRUE)
})

test_that(".scroll_module_score is seed-stable and filter-invariant per cell", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  sig <- c("CD3D", "CD8A", "NKG7")
  s1 <- scroll:::.scroll_module_score(data$cells, sig, data, "RNA", nbin = 4, ctrl = 3, seed = 1)
  s2 <- scroll:::.scroll_module_score(data$cells, sig, data, "RNA", nbin = 4, ctrl = 3, seed = 1)
  expect_equal(s1, s2)                                    # deterministic given a seed
  expect_length(s1, nrow(data$cells))
  # a cell's score does not depend on which other cells are shown
  sub <- data$cells[1:60, , drop = FALSE]
  s_sub <- scroll:::.scroll_module_score(sub, sig, data, "RNA", nbin = 4, ctrl = 3, seed = 1)
  expect_equal(as.numeric(s_sub), as.numeric(s1[1:60]))
})

# A small object with a real, cell-varying signature (`sig`) for the parity tests.
sig_object <- function(ng = 300L, nc = 150L, lambda = 0.8, seed = 42) {
  set.seed(seed)
  genes <- paste0("Gene", seq_len(ng))
  cm <- matrix(rpois(ng * nc, lambda), nrow = ng,
               dimnames = list(genes, paste0("c", seq_len(nc))))
  sig <- genes[1:8]
  hot <- seq_len(nc / 2)                                  # give the signature real signal
  cm[sig, hot] <- cm[sig, hot] + rpois(length(sig) * length(hot), 5)
  obj <- SeuratObject::CreateSeuratObject(counts = Matrix::Matrix(cm, sparse = TRUE),
                                          min.cells = 0, min.features = 0)
  obj <- SeuratObject::SetAssayData(obj, layer = "data",
                                    new.data = Matrix::Matrix(.lognorm(cm, 1e4), sparse = TRUE))
  emb <- matrix(rnorm(nc * 2), ncol = 2, dimnames = list(colnames(obj), c("UMAP_1", "UMAP_2")))
  obj[["umap"]] <- SeuratObject::CreateDimReducObject(embeddings = emb, key = "UMAP_", assay = "RNA")
  list(obj = obj, sig = sig)
}
sig_build <- function(obj, quantize = FALSE, ranks = TRUE, ...) {
  dir <- tempfile("sig")
  suppressMessages(scroll_build(obj, dir, assays = "RNA", embeddings = "umap",
                                quantize = quantize, ranks = ranks, verbose = FALSE, ...))
  scroll:::.scroll_load(dir)
}

test_that("AddModuleScore matches Seurat::AddModuleScore exactly (same control genes)", {
  skip_if_not_installed("Seurat")
  for (cfg in list(list(ng = 300L, nbin = 10L, ctrl = 20L),
                   list(ng = 3000L, nbin = 24L, ctrl = 100L))) {   # Seurat's defaults
    x <- sig_object(cfg$ng)
    ref <- suppressWarnings(Seurat::AddModuleScore(x$obj, features = list(x$sig), name = "s",
                                                   nbin = cfg$nbin, ctrl = cfg$ctrl, seed = 1))$s1
    data <- sig_build(x$obj)
    scr <- scroll:::.scroll_module_score(data$cells, x$sig, data, "RNA",
                                         nbin = cfg$nbin, ctrl = cfg$ctrl, seed = 1)
    expect_true(attr(scr, "exact"))
    # identical control set; the residue is only the store's float32 values
    expect_equal(as.numeric(scr), unname(ref[match(data$cells$cell, colnames(x$obj))]),
                 tolerance = 1e-5)
    scroll_disconnect(data$con)
  }
})

test_that("the control draw leaves the session's RNG untouched", {
  set.seed(99); a <- stats::runif(1)
  set.seed(99)
  scroll:::.scroll_seurat_controls(stats::setNames(seq(0, 1, length.out = 50), paste0("g", 1:50)),
                                   c("g3", "g40"), nbin = 5, ctrl = 4, seed = 1)
  expect_identical(stats::runif(1), a)
})

test_that("without build-time gene means, AddModuleScore falls back to a store scan", {
  x <- sig_object()
  data <- sig_build(x$obj)
  on.exit(scroll_disconnect(data$con), add = TRUE)
  data$manifest$assays$RNA$stats <- NULL                  # as an older build
  scr <- scroll:::.scroll_module_score(data$cells, x$sig, data, "RNA", nbin = 10, ctrl = 20)
  expect_false(attr(scr, "exact"))
  expect_gt(stats::cor(scr, as.numeric(scroll:::.scroll_module_score(
    data$cells, x$sig, sig_build(x$obj), "RNA", nbin = 10, ctrl = 20))), 0.9)
})

test_that("UCell matches UCell::ScoreSignatures_UCell exactly, quantized or not", {
  skip_if_not_installed("UCell")
  x <- sig_object(ng = 2000L, nc = 120L, lambda = 0.3, seed = 7)
  dat <- SeuratObject::GetAssayData(x$obj, layer = "data")
  for (q in c(FALSE, TRUE)) {
    data <- sig_build(x$obj, quantize = q)
    for (mr in c(1500, 100)) {
      u <- suppressWarnings(UCell::ScoreSignatures_UCell(dat, features = list(s = x$sig),
                                                         maxRank = mr, ncores = 1))
      scr <- scroll:::.scroll_ucell_score(data$cells, x$sig, data, "RNA", max_rank = mr)
      expect_equal(scr, unname(u[match(data$cells$cell, rownames(u)), 1]), tolerance = 1e-12)
    }
    scroll_disconnect(data$con)
  }
})

test_that("UCell and AUCell-style match a dense reference built from base::rank", {
  x <- sig_object(ng = 400L, nc = 60L, lambda = 0.4, seed = 3)
  dat <- as.matrix(SeuratObject::GetAssayData(x$obj, layer = "data"))
  R <- apply(dat, 2, function(v) rank(-v, ties.method = "average"))   # genes x cells
  rs <- t(R[x$sig, ]); n <- length(x$sig)
  data <- sig_build(x$obj, quantize = TRUE)
  on.exit(scroll_disconnect(data$con), add = TRUE)
  key <- match(data$cells$cell, rownames(rs))
  # UCell, max rank 50
  u <- 1 - (rowSums(pmin(rs, 50)) - n * (n + 1) / 2) / (n * 50 - n * (n + 1) / 2)
  expect_equal(scroll:::.scroll_ucell_score(data$cells, x$sig, data, "RNA", max_rank = 50),
               unname(u[key]))
  # AUCell-style: recovery-curve area over the top ceiling(5%) ranks, normalized
  thr <- ceiling(0.05 * nrow(dat)); best <- seq_len(n); best <- best[best < thr]
  a <- vapply(seq_len(nrow(rs)), function(i) {
    v <- sort(rs[i, ][rs[i, ] < thr]); sum(diff(c(v, thr)) * seq_along(v))
  }, numeric(1)) / sum(diff(c(best, thr)) * seq_along(best))
  expect_equal(scroll:::.scroll_aucell_score(data$cells, x$sig, data, "RNA"), a[key])
})

test_that("rank-based scores need a project built with ranks", {
  x <- sig_object()
  data <- sig_build(x$obj, ranks = FALSE)
  on.exit(scroll_disconnect(data$con), add = TRUE)
  expect_null(data$manifest$assays$RNA$ranks)
  expect_error(scroll:::.scroll_ucell_score(data$cells, x$sig, data, "RNA"), "ranks = TRUE")
  expect_true(isTRUE(data$manifest$assays$RNA$stats))    # gene means are always baked
})

test_that("signature_server requires at least one gene", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  panel <- Filter(function(p) identical(p$id, "signature"), scroll:::.scroll_builtin_panels())[[1]]
  shiny::testServer(panel$server, args = list(data = data), {
    session$setInputs(sig = character(0), method = "mean", view = "umap", compute = 1)
    expect_error(plot_r())                               # validate() -> "Add one or more genes..."
  })
})

test_that("up / down signatures: up minus weight x down, floored at 0 for UCell", {
  ss <- scroll:::.scroll_signed_score
  expect_equal(ss(c(0.5, 0.2), NULL, "ucell"), c(0.5, 0.2))           # no down genes: unchanged
  expect_equal(ss(c(0.5, 0.2), c(0.1, 0.4), "ucell"), c(0.4, 0))      # floored at 0
  expect_equal(ss(c(0.5, 0.2), c(0.1, 0.4), "mean"), c(0.4, -0.2))    # others may go negative
  expect_equal(ss(c(0.5, 0.2), c(0.1, 0.4), "aucell", w_neg = 0.5), c(0.45, 0))
})

test_that("UCell with down genes matches UCell::ScoreSignatures_UCell's 'gene-' entries", {
  skip_if_not_installed("UCell")
  x <- sig_object(ng = 2000L, nc = 120L, lambda = 0.3, seed = 7)
  dat <- SeuratObject::GetAssayData(x$obj, layer = "data")
  up <- x$sig[1:3]; down <- setdiff(rownames(dat), x$sig)[1:4]
  data <- sig_build(x$obj)
  on.exit(scroll_disconnect(data$con), add = TRUE)
  for (w in c(1, 0.5)) {
    u <- suppressWarnings(UCell::ScoreSignatures_UCell(dat, ncores = 1, w_neg = w,
           features = list(s = c(paste0(up, "+"), paste0(down, "-")))))
    scr <- scroll:::.scroll_signed_score(
      scroll:::.scroll_ucell_score(data$cells, up, data, "RNA"),
      scroll:::.scroll_ucell_score(data$cells, down, data, "RNA"), "ucell", w_neg = w)
    expect_equal(scr, unname(u[match(data$cells$cell, rownames(u)), 1]), tolerance = 1e-12)
  }
})

test_that("signature_server scores up and down genes and rejects a gene in both", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  feats <- scroll:::.scroll_features_of(data$manifest, data$manifest$default_assay)
  red <- scroll:::.scroll_view_embeddings(data$manifest, NULL)[[1]]
  panel <- Filter(function(p) identical(p$id, "signature"), scroll:::.scroll_builtin_panels())[[1]]
  shiny::testServer(panel$server, args = list(data = data), {
    session$setInputs(sig = feats[1:3], sig_down = feats[4:5], method = "mean", w_neg = 1,
                      view = "umap", reduction = red, palette = "grey-purple", clip = c(0, 100),
                      order = TRUE, legend = TRUE, raster = FALSE, aspect = 1, compute = 1)
    d <- score_r()
    up <- scroll:::.scroll_cell_mean(d$cells, feats[1:3], data, data$manifest$default_assay)
    dn <- scroll:::.scroll_cell_mean(d$cells, feats[4:5], data, data$manifest$default_assay)
    expect_equal(d$values$value, as.numeric(up - dn))
    expect_match(scroll:::.scroll_sig_label(d$genes, d$down), "3 up, 2 down")
    expect_s3_class(plot_r(), "ggplot")
    session$setInputs(sig_down = feats[c(1, 4)], compute = 2)            # feats[1] in both
    expect_error(score_r(), "in both")
  })
})
