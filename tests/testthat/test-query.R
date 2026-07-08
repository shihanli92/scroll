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

  expect_setequal(hit$cell, names(truth))
  m <- merge(hit, data.frame(cell = names(truth), truth = as.numeric(truth)), by = "cell")
  expect_lt(max(abs(m$value - m$truth)), man$assays$RNA$max / 255 + 1e-9)
})
