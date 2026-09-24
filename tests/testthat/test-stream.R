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

test_that("scroll_build_stream(counts = TRUE) matches a single build and runs pseudobulk", {
  obj <- make_test_object()                          # 120 cells, RNA with a counts layer
  d1 <- file.path(tempdir(), "scroll-stream-counts-single")
  suppressMessages(scroll_build(obj, d1, assays = "RNA", embeddings = "umap",
    meta_cols = c("condition", "celltype"), quantize = FALSE, counts = TRUE,
    overwrite = TRUE))
  d2 <- file.path(tempdir(), "scroll-stream-counts")
  suppressMessages(scroll_build_stream(d2, list(1:60, 61:120), function(idx) obj[, idx],
    assays = "RNA", embeddings = "umap", meta_cols = c("condition", "celltype"),
    counts = TRUE, id_of = function(x) paste0("h", x[[1]]), overwrite = TRUE,
    verbose = FALSE))

  expect_true(isTRUE(scroll_manifest(d2)$has_counts))
  expect_length(list.files(file.path(d2, "counts", "RNA")), 2L)   # one part per source
  norm <- function(x) { x <- as.data.frame(x); x[order(x$feature, x$cell), ] |> `rownames<-`(NULL) }
  single <- norm(arrow::read_parquet(file.path(d1, "counts", "RNA.parquet")))
  stream <- norm(dplyr::collect(arrow::open_dataset(file.path(d2, "counts", "RNA"))))
  expect_equal(stream, single)                       # same global cell index, same counts

  skip_if_not_installed("edgeR"); skip_if_not_installed("limma")
  run <- function(d) {
    data <- scroll:::.scroll_load(d); on.exit(scroll_disconnect(data$con))
    scroll_pseudobulk_de(data, "RNA", aggregate_cols = "celltype", ident1 = "T",
                         ident2 = "B", replicate_col = "condition", min_cells = 5)
  }
  r1 <- run(d1); r2 <- run(d2)
  expect_equal(r2$gene, r1$gene)
  expect_equal(r2$logFC, r1$logFC, tolerance = 1e-8)
})

test_that("an append run must keep the project's counts setting", {
  obj <- make_test_object()
  dir <- file.path(tempdir(), "scroll-stream-counts-mismatch")
  reader <- function(idx) obj[, idx]
  args <- list(outdir = dir, reader = reader, assays = "RNA", embeddings = "umap",
               meta_cols = c("condition", "celltype"),
               id_of = function(x) paste0("h", x[[1]]), verbose = FALSE)
  suppressMessages(do.call(scroll_build_stream,
    c(args, list(sources = list(1:60), counts = TRUE, overwrite = TRUE))))
  expect_error(
    suppressMessages(do.call(scroll_build_stream,
      c(args, list(sources = list(1:60, 61:120), counts = FALSE)))),
    "built with counts = TRUE")
})
