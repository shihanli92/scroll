test_that("scroll_build emits the expected artifacts", {
  dir <- test_project()
  expect_true(file.exists(file.path(dir, "cells.parquet")))
  expect_true(file.exists(file.path(dir, "manifest.yaml")))
  expect_true(file.exists(file.path(dir, "config.yaml")))
  expect_true(file.exists(file.path(dir, "app.R")))
  expect_true(dir.exists(file.path(dir, "expr", "RNA")))

  cells <- as.data.frame(arrow::read_parquet(file.path(dir, "cells.parquet")))
  expect_true(all(c("cell", "condition", "celltype", "umap_1", "umap_2") %in% names(cells)))
  expect_equal(nrow(cells), 120)
})

test_that("manifest records assays, embeddings, and defaults", {
  man <- scroll_manifest(test_project())
  expect_equal(man$default_assay, "RNA")
  expect_equal(man$default_embedding, "umap")
  expect_true(man$quantize)
  expect_true("CD3D" %in% unlist(man$assays$RNA$features))
  expect_equal(man$embeddings$umap$dims, 2)
})

test_that("feature partitions exist for exported features", {
  dir <- test_project()
  parts <- list.files(file.path(dir, "expr", "RNA"))
  expect_true("feature=MS4A1" %in% parts)
})
