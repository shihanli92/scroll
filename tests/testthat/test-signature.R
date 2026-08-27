# The Signature panel: live per-cell module scoring from a gene list, shown on the
# embedding or as a violin. Runtime-only (Mean / Scaled), no build artifact.

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
    expect_s3_class(plot_r(), "ggplot")

    # CSV export is one score per cell
    csv <- csv_r()
    expect_setequal(names(csv), c("cell", "signature_score"))
    expect_equal(nrow(csv), nrow(data$cells))

    # scaled scoring still renders
    session$setInputs(method = "scaled")
    expect_s3_class(plot_r(), "ggplot")

    # violin view, grouped by a categorical column
    session$setInputs(view = "violin", group = cat1)
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
  expect_equal(s_sub, s1[1:60])
})

test_that("runtime AddModuleScore agrees with Seurat::AddModuleScore", {
  skip_if_not_installed("Seurat")
  set.seed(42)
  ng <- 300L; nc <- 150L
  genes <- paste0("Gene", seq_len(ng))
  cm <- matrix(rpois(ng * nc, 0.8), nrow = ng,
               dimnames = list(genes, paste0("c", seq_len(nc))))
  sig <- genes[1:8]
  hot <- seq_len(nc / 2)                                  # give the signature real, cell-varying signal
  cm[sig, hot] <- cm[sig, hot] + rpois(length(sig) * length(hot), 5)
  obj <- SeuratObject::CreateSeuratObject(counts = Matrix::Matrix(cm, sparse = TRUE),
                                          min.cells = 0, min.features = 0)
  obj <- SeuratObject::SetAssayData(obj, layer = "data",
                                    new.data = Matrix::Matrix(.lognorm(cm, 1e4), sparse = TRUE))
  emb <- matrix(rnorm(nc * 2), ncol = 2, dimnames = list(colnames(obj), c("UMAP_1", "UMAP_2")))
  obj[["umap"]] <- SeuratObject::CreateDimReducObject(embeddings = emb, key = "UMAP_", assay = "RNA")

  obj <- Seurat::AddModuleScore(obj, features = list(sig), name = "sigscore",
                                nbin = 10, ctrl = 20, seed = 1)
  seurat <- obj$sigscore1

  dir <- file.path(tempdir(), "scroll-sig-parity")
  suppressMessages(scroll_build(obj, dir, assays = "RNA", embeddings = "umap",
                                quantize = FALSE, overwrite = TRUE))
  data <- scroll:::.scroll_load(dir)
  on.exit(scroll_disconnect(data$con), add = TRUE)
  scr <- scroll:::.scroll_module_score(data$cells, sig, data, "RNA", nbin = 10, ctrl = 20, seed = 1)

  seurat <- seurat[match(data$cells$cell, colnames(obj))]  # align by barcode
  expect_gt(stats::cor(scr, seurat), 0.9)                  # near-identical; control RNG differs slightly
})

test_that("signature_server requires at least one gene", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  panel <- Filter(function(p) identical(p$id, "signature"), scroll:::.scroll_builtin_panels())[[1]]
  shiny::testServer(panel$server, args = list(data = data), {
    session$setInputs(sig = character(0), method = "mean", view = "umap")
    expect_error(plot_r())                               # validate() -> "Add one or more genes..."
  })
})
