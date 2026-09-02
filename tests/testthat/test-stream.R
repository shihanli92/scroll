test_that("scroll_build_stream assembles two sources into one store", {
  dir <- file.path(tempdir(), "scroll-stream-ok")
  a <- make_test_object(n = 40, seed = 1)
  b <- make_test_object(n = 30, seed = 2)
  colnames(a) <- paste0("a_", colnames(a))   # disjoint barcodes across sources
  colnames(b) <- paste0("b_", colnames(b))
  reader <- function(s) if (s == "A") a else b

  suppressMessages(scroll_build_stream(
    dir, sources = c("A", "B"), reader = reader, assays = "RNA",
    embeddings = "umap", meta_cols = c("condition", "celltype"), overwrite = TRUE))

  man <- scroll_manifest(dir)
  expect_equal(man$n_cells, ncol(a) + ncol(b))

  con <- scroll_connect(dir); on.exit(scroll_disconnect(con))
  q <- scroll_query_feature(con, "RNA", "CD3D")     # spans both sources
  expect_s3_class(q, "data.frame")
  expect_gt(nrow(q), 0)
  # the two sources' parts are compacted into one feature-sorted file per assay
  expect_identical(list.files(file.path(dir, "expr", "RNA")), "part-0.parquet")
})

test_that("scroll_build_stream errors, naming the source, on a missing declared column", {
  dir <- file.path(tempdir(), "scroll-stream-bad")
  a <- make_test_object(n = 40, seed = 1)
  b <- make_test_object(n = 30, seed = 2)
  colnames(a) <- paste0("a_", colnames(a))
  colnames(b) <- paste0("b_", colnames(b))
  b$condition <- NULL                                # source B lacks a declared col
  reader <- function(s) if (s == "A") a else b

  expect_error(
    suppressMessages(scroll_build_stream(
      dir, sources = c("A", "B"), reader = reader, assays = "RNA",
      embeddings = "umap", meta_cols = c("condition", "celltype"), overwrite = TRUE)),
    "'B'.*condition")
})
