test_that("scroll_query_feature reads a single feature and dequantizes", {
  dir <- test_project()
  man <- scroll_manifest(dir)
  con <- scroll_connect(dir)
  on.exit(scroll_disconnect(con))

  hit <- scroll_query_feature(con, "RNA", "MS4A1")
  expect_true(all(c("cell", "value") %in% names(hit)))
  expect_gt(nrow(hit), 0)
  # stored as quantized uint8
  expect_true(all(hit$value >= 0 & hit$value <= 255))

  deq <- scroll_dequantize(hit$value, man, "RNA")
  expect_true(all(deq >= 0 & deq <= man$assays$RNA$max + 1e-9))
})

test_that("scroll_query_features reads several partitions at once", {
  dir <- test_project()
  con <- scroll_connect(dir)
  on.exit(scroll_disconnect(con))
  hit <- scroll_query_features(con, "RNA", c("CD3D", "MS4A1", "NOT_A_GENE"))
  expect_true(all(c("feature", "cell", "value") %in% names(hit)))
  expect_setequal(unique(hit$feature), c("CD3D", "MS4A1"))
  expect_equal(nrow(scroll_query_features(con, "RNA", character(0))), 0)
})

test_that("querying a missing feature returns zero rows, not an error", {
  dir <- test_project()
  con <- scroll_connect(dir)
  on.exit(scroll_disconnect(con))
  expect_equal(nrow(scroll_query_feature(con, "RNA", "NOT_A_GENE")), 0)
})

test_that("scroll_query_feature matches a direct file read (and repeats identically)", {
  dir <- test_project()
  con <- scroll_connect(dir); on.exit(scroll_disconnect(con))
  part <- file.path(dir, "expr", "RNA", "part-0.parquet")     # single per-assay file
  direct <- as.data.frame(arrow::read_parquet(part))
  direct <- direct[direct$feature == "MS4A1", , drop = FALSE] # ...filtered to the feature
  ord <- function(d) d[order(d$cell), c("cell", "value")]
  hit <- scroll_query_feature(con, "RNA", "MS4A1")
  expect_equal(ord(hit), ord(direct), ignore_attr = TRUE)
  # the cached Dataset in the handle returns identical rows on a repeat lookup
  expect_equal(ord(scroll_query_feature(con, "RNA", "MS4A1")), ord(hit), ignore_attr = TRUE)
})

test_that("scroll_aggregate_counts sums a many-to-many mapping correctly", {
  dir <- test_project()                                       # built with counts = TRUE
  con <- scroll_connect(dir); on.exit(scroll_disconnect(con))
  mapping <- data.frame(cell = c(1L, 2L, 1L), psample = c("A", "A", "B"))  # v2 int index; c1 in A and B
  agg <- scroll_aggregate_counts(con, "RNA", mapping)
  expect_true(all(c("feature", "psample", "count") %in% names(agg)))

  # independent expected sums from a direct read of the counts store
  raw <- as.data.frame(arrow::read_parquet(file.path(dir, "counts", "RNA.parquet")))
  expected <- do.call(rbind, lapply(unique(mapping$psample), function(ps) {
    sub <- raw[raw$cell %in% mapping$cell[mapping$psample == ps], ]
    a <- stats::aggregate(value ~ feature, sub, sum)
    data.frame(feature = a$feature, psample = ps, count = a$value)
  }))
  key <- function(d) d[order(d$feature, d$psample), c("feature", "psample", "count")]
  expect_equal(key(agg), key(expected), ignore_attr = TRUE)
})

test_that("dequantized values match the source matrix within the quant bound", {
  skip_if_not_installed("SeuratObject")
  obj <- make_test_object()
  dir <- file.path(tempdir(), "scroll-qc-proj")
  suppressMessages(scroll_build(obj, dir, overwrite = TRUE))
  man <- scroll_manifest(dir)
  con <- scroll_connect(dir); on.exit(scroll_disconnect(con))

  truth <- SeuratObject::GetAssayData(obj, assay = "RNA", layer = "data")["MS4A1", ]
  truth <- truth[truth > 0]
  hit <- scroll_query_feature(con, "RNA", "MS4A1")
  hit$value <- scroll_dequantize(hit$value, man, "RNA")

  bc <- colnames(obj)[hit$cell]                # v2 int index -> barcode (matrix col order)
  expect_setequal(bc, names(truth))
  m <- merge(data.frame(cell = bc, value = hit$value),
             data.frame(cell = names(truth), truth = as.numeric(truth)), by = "cell")
  expect_lt(max(abs(m$value - m$truth)), man$assays$RNA$max / 255 + 1e-9)
})
