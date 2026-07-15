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

test_that("scroll_update errors on unknown refs and no-ops when empty", {
  dir <- file.path(tempdir(), "scroll-update-proj")   # built above
  obj <- make_test_object()
  expect_error(scroll_update(dir, obj, embeddings = "nope"), "no reduction")
  expect_error(scroll_update(dir, obj, meta_cols = "nope"), "no metadata")
  expect_equal(suppressMessages(scroll_update(dir, obj)), dir)   # nothing to do
})
