# v2 store: int32 cell-index + float32 values, and the store-agnostic join.

test_that("v2 build: manifest flags, int32 cell, one file per assay", {
  skip_if_not_installed("SeuratObject")
  dir <- file.path(tempdir(), "scroll-v2-flags")
  suppressMessages(scroll_build(make_test_object(), dir, quantize = FALSE, overwrite = TRUE))
  m <- scroll_manifest(dir)
  expect_equal(m$store_version, 2L)
  expect_true(isTRUE(m$cell_index))
  adir <- file.path(dir, "expr", "RNA")
  t <- arrow::open_dataset(adir)
  expect_equal(t$schema$GetFieldByName("cell")$type$ToString(), "int32")
  expect_equal(t$schema$GetFieldByName("value")$type$ToString(), "float")   # float32
  # `feature` must be large_utf8 (int64 offsets): on a big assay the per-nonzero
  # gene-name column tops 2 GB, which overflows plain utf8's int32 offsets.
  expect_equal(t$schema$GetFieldByName("feature")$type$ToString(), "large_string")
  expect_identical(list.files(adir, pattern = "\\.parquet$"), "part-0.parquet")  # single file
})

test_that("v2 query/join parity: reproduces the true expression, whole + subset", {
  skip_if_not_installed("SeuratObject")
  obj <- make_test_object()
  truth <- as.numeric(SeuratObject::GetAssayData(obj, assay = "RNA", layer = "data")["MS4A1", ])

  # float32 (unquantized) store: near-exact parity
  dirf <- file.path(tempdir(), "scroll-v2-f32")
  suppressMessages(scroll_build(obj, dirf, quantize = FALSE, overwrite = TRUE))
  data <- scroll:::.scroll_load(dirf); on.exit(scroll_disconnect(data$con))
  expect_true(".gidx" %in% names(data$cells))
  v <- data$query1("RNA", "MS4A1")
  expect_true(is.integer(v$cell))                                   # int32 -> R integer
  got <- scroll:::.scroll_expr_vector(data$cells, v)
  expect_equal(got, truth, tolerance = 1e-5)                        # float32 ~lossless

  # subset: the vector still aligns to the filtered rows via the global index
  set.seed(1); sub <- data$cells[sample(nrow(data$cells), 40), , drop = FALSE]
  expect_equal(scroll:::.scroll_expr_vector(sub, v), truth[sub$.gidx], tolerance = 1e-5)
})

test_that("scroll_build_stream matches a single build of the concatenation", {
  skip_if_not_installed("SeuratObject")
  obj <- make_test_object()                          # 120 cells, RNA + umap
  d1 <- file.path(tempdir(), "scroll-stream-single")
  suppressMessages(scroll_build(obj, d1, assays = "RNA", embeddings = "umap",
    meta_cols = c("condition", "celltype"), quantize = FALSE, overwrite = TRUE))

  # stream two halves (same cell order => same global indices as the single build)
  reader <- function(idx) obj[, idx]
  d2 <- file.path(tempdir(), "scroll-stream-multi"); unlink(d2, recursive = TRUE)
  suppressMessages(scroll_build_stream(d2, list(1:60, 61:120), reader,
    assays = "RNA", embeddings = "umap", meta_cols = c("condition", "celltype"),
    id_of = function(x) paste0("h", x[[1]]), verbose = FALSE))

  m1 <- scroll_manifest(d1); m2 <- scroll_manifest(d2)
  expect_equal(m2$n_cells, m1$n_cells)
  expect_true(isTRUE(m2$cell_index))
  expect_setequal(m2$assays$RNA$features, m1$assays$RNA$features)
  expect_setequal(m2$meta$celltype$levels, m1$meta$celltype$levels)

  c1 <- scroll_connect(d1); c2 <- scroll_connect(d2)
  on.exit({ scroll_disconnect(c1); scroll_disconnect(c2) })
  g <- m1$assays$RNA$features[[1]]
  q1 <- scroll_query_feature(c1, "RNA", g); q1 <- q1[order(q1$cell), ]
  q2 <- scroll_query_feature(c2, "RNA", g); q2 <- q2[order(q2$cell), ]
  expect_equal(q2$cell, q1$cell)                     # identical global indices
  expect_equal(q2$value, q1$value, tolerance = 1e-6) # identical values
})

test_that(".scroll_expr_vector is store-agnostic (v2 int key + v1 barcode key)", {
  cells <- data.frame(cell = c("A", "B", "C", "D"), .gidx = 1:4, stringsAsFactors = FALSE)
  expect_equal(scroll:::.scroll_expr_vector(cells, data.frame(cell = c(3L, 1L), value = c(9, 5))),
               c(5, 0, 9, 0))                                       # v2: integer key
  expect_equal(scroll:::.scroll_expr_vector(
                 cells, data.frame(cell = c("C", "A"), value = c(9, 5), stringsAsFactors = FALSE)),
               c(5, 0, 9, 0))                                       # v1: barcode key (back-compat)
  expect_equal(scroll:::.scroll_expr_vector(cells, data.frame(cell = integer(), value = numeric())),
               c(0, 0, 0, 0))                                       # empty -> all zero
})
