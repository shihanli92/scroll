test_that("scroll_query_cells returns all features for the given cells", {
  dir <- test_project()
  con <- scroll_connect(dir); on.exit(scroll_disconnect(con))
  cells <- as.data.frame(arrow::read_parquet(file.path(dir, "cells.parquet")))$cell[1:20]
  hit <- scroll_query_cells(con, "RNA", cells)
  expect_true(all(c("feature", "cell", "value") %in% names(hit)))
  expect_true(all(hit$cell %in% cells))
})

test_that("scroll_de runs presto and returns ranked marker results", {
  skip_if_not_installed("presto")
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  ident1 <- sort(unique(as.character(data$cells$celltype)))[[1]]
  res <- scroll_de(data, "RNA", "celltype", ident1 = ident1)     # one-vs-rest
  expect_true(all(c("gene", "logFC", "auc", "pct.1", "pct.2", "p_val", "p_val_adj") %in% names(res)))
  expect_gt(nrow(res), 0)
  expect_false(is.unsorted(res$p_val_adj))                       # ranked by adjusted p
})

test_that("de_server computes a table on click", {
  skip_if_not_installed("presto")
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  ct <- sort(unique(as.character(data$cells$celltype)))[[1]]
  shiny::testServer(scroll:::de_server, args = list(data = data), {
    session$setInputs(group = "celltype", ident1 = ct, ident2 = "rest",
                      minpct = 10, topn = 20, compute = 1)
    expect_false(is.null(output$table))
  })
})
