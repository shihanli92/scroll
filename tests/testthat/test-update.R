test_that("scroll_update adds a reduction + metadata without touching expr", {
  dir <- file.path(tempdir(), "scroll-update-proj")
  obj <- make_test_object()
  suppressMessages(scroll_build(obj, dir, assays = "RNA", overwrite = TRUE))

  exprfiles <- list.files(file.path(dir, "expr"), recursive = TRUE, full.names = TRUE)
  before <- file.info(exprfiles)[, c("size", "mtime")]

  # a new 3-D reduction + a numeric metadata column
  set.seed(2)
  pca <- matrix(rnorm(ncol(obj) * 3), ncol = 3,
                dimnames = list(colnames(obj), paste0("PC_", 1:3)))
  obj[["pca"]] <- SeuratObject::CreateDimReducObject(embeddings = pca, key = "PC_", assay = "RNA")
  obj$score <- runif(ncol(obj))

  suppressMessages(scroll_update(dir, obj, embeddings = "pca", meta_cols = "score"))

  man <- scroll_manifest(dir)
  expect_true("pca" %in% names(man$embeddings))
  expect_equal(man$embeddings$pca$dims, 3)
  expect_equal(man$embeddings$pca$n_covered, ncol(obj))
  expect_equal(man$meta$score$type, "numeric")

  cells <- as.data.frame(arrow::read_parquet(file.path(dir, "cells.parquet")))
  expect_true(all(c("pca_1", "pca_2", "pca_3", "score") %in% names(cells)))
  expect_equal(nrow(cells), ncol(obj))

  # expr store byte-identical (no re-export)
  after <- file.info(exprfiles)[, c("size", "mtime")]
  expect_identical(before$size, after$size)
  expect_identical(before$mtime, after$mtime)

  # the query layer still works against the untouched store
  con <- scroll_connect(dir); on.exit(scroll_disconnect(con))
  expect_s3_class(scroll_query_feature(con, "RNA", "CD3D"), "data.frame")
})

test_that("scroll_update registers a subset view (partial coverage + scoped meta)", {
  dir <- file.path(tempdir(), "scroll-update-sub")
  obj <- make_test_object()
  suppressMessages(scroll_build(obj, dir, assays = "RNA", overwrite = TRUE))

  tcells <- colnames(obj)[obj$celltype == "T"]
  sub <- subset(obj, cells = tcells)
  emb <- matrix(rnorm(length(tcells) * 2), ncol = 2,
                dimnames = list(tcells, c("UMAP_1", "UMAP_2")))
  sub[["umap_t"]] <- SeuratObject::CreateDimReducObject(embeddings = emb, key = "umapt_", assay = "RNA")
  sub$tsub <- factor(sample(c("a", "b"), length(tcells), replace = TRUE))
  obj2 <- scroll_add_subset(obj, sub, "tcell", embeddings = c(t_umap = "umap_t"),
                            meta = "tsub", label = "T cells")

  suppressMessages(scroll_update(dir, obj2, embeddings = "t_umap", meta_cols = "tsub"))

  man <- scroll_manifest(dir)
  expect_true("tcell" %in% names(man$subsets))
  expect_equal(man$subsets$tcell$primary_embedding, "t_umap")
  expect_lt(man$embeddings$t_umap$n_covered, man$n_cells)   # partial coverage
  expect_equal(man$meta$tsub$scope, "tcell")
})

test_that("scroll_update preserves non-ASCII feature names in the manifest", {
  # Regression: scroll_update rewrites the manifest, so it must write it with
  # unicode = TRUE (as the build does); otherwise a non-ASCII feature name is
  # re-escaped to a `<U+XXXX>` form that no longer matches the store. yaml only
  # keeps UTF-8 under a UTF-8 locale, so pin one for this test (skip if none).
  old_loc <- Sys.getlocale("LC_CTYPE")
  on.exit(suppressWarnings(Sys.setlocale("LC_CTYPE", old_loc)), add = TRUE)
  set_ok <- ""
  for (loc in c("en_US.UTF-8", "C.UTF-8", "en_US.utf8"))
    if (nzchar(set_ok <- suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) break
  skip_if_not(nzchar(set_ok), "no UTF-8 locale available")

  dir <- file.path(tempdir(), "scroll-update-unicode")
  set.seed(3)
  genes <- c("CD3D", "Fc\u03b5RIa", "Na\u00efve1")   # non-ASCII feature names
  n <- 40
  cm <- matrix(rpois(length(genes) * n, 0.8), nrow = length(genes),
               dimnames = list(genes, paste0("cell", seq_len(n))))
  norm <- log1p(sweep(cm, 2, pmax(colSums(cm), 1), "/") * 1e4)
  obj <- SeuratObject::CreateSeuratObject(counts = Matrix::Matrix(cm, sparse = TRUE),
                                          min.cells = 0, min.features = 0)
  obj <- SeuratObject::SetAssayData(obj, layer = "data",
                                    new.data = Matrix::Matrix(norm, sparse = TRUE))
  emb <- matrix(rnorm(n * 2), ncol = 2,
                dimnames = list(colnames(obj), c("UMAP_1", "UMAP_2")))
  obj[["umap"]] <- SeuratObject::CreateDimReducObject(embeddings = emb, key = "UMAP_", assay = "RNA")

  suppressMessages(scroll_build(obj, dir, assays = "RNA", overwrite = TRUE))
  feats0 <- unlist(scroll_manifest(dir)$assays$RNA$features)
  expect_true(all(genes %in% feats0))              # build round-trips UTF-8

  obj$score <- runif(n)
  suppressMessages(scroll_update(dir, obj, meta_cols = "score"))
  feats1 <- unlist(scroll_manifest(dir)$assays$RNA$features)
  expect_true(all(genes %in% feats1))              # update must NOT re-escape them
  expect_false(any(grepl("<U\\+", feats1)))
})

test_that("scroll_update errors on unknown refs and no-ops when empty", {
  dir <- file.path(tempdir(), "scroll-update-proj")   # built above
  obj <- make_test_object()
  expect_error(scroll_update(dir, obj, embeddings = "nope"), "no reduction")
  expect_error(scroll_update(dir, obj, meta_cols = "nope"), "no metadata")
  expect_equal(suppressMessages(scroll_update(dir, obj)), dir)   # nothing to do
})
