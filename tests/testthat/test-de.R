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

test_that("scroll_de one-vs-rest drops NA-group cells (no NA labels to presto)", {
  skip_if_not_installed("presto")
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  data$cells$celltype[1:10] <- NA                       # un-annotated cells
  ident1 <- sort(unique(stats::na.omit(as.character(data$cells$celltype))))[[1]]
  res <- scroll_de(data, "RNA", "celltype", ident1 = ident1)   # would error on NA labels
  expect_gt(nrow(res), 0)
})

test_that("scroll_de pools multiple levels per group", {
  skip_if_not_installed("presto")
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  cts <- sort(unique(as.character(data$cells$celltype)))   # e.g. B, NK, T
  # two pooled levels as group 1 vs the third
  res <- scroll_de(data, "RNA", "celltype", ident1 = cts[1:2], ident2 = cts[3])
  expect_true(all(c("gene", "logFC", "pct.1", "pct.2") %in% names(res)))
  expect_gt(nrow(res), 0)
  # pooled group 1 vs rest
  expect_gt(nrow(scroll_de(data, "RNA", "celltype", ident1 = cts[1:2])), 0)
})

test_that("view_volcano renders from a DE result (and empty input)", {
  set.seed(1)
  de <- data.frame(gene = paste0("G", 1:50), logFC = stats::rnorm(50, 0, 1.5),
                   auc = runif(50), pct.1 = runif(50), pct.2 = runif(50),
                   p_val = 10^(-runif(50, 0, 20)), p_val_adj = 10^(-runif(50, 0, 20)))
  expect_s3_class(view_volcano(de, list(lfc = 1, padj = 0.05, label_n = 10)), "ggplot")
  expect_s3_class(view_volcano(data.frame()), "ggplot")   # empty -> placeholder, no error
})

test_that("de_server computes table + volcano from one contrast", {
  skip_if_not_installed("presto")
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  ct <- sort(unique(as.character(data$cells$celltype)))[[1]]
  shiny::testServer(scroll:::de_server, args = list(data = data), {
    session$setInputs(group = "celltype", ident1 = ct, ident2 = "rest",
                      minpct = 10, topn = 20, lfc = 1, padj = 0.05, labeln = 10,
                      aspect = 1, compute = 1)
    # a single Compute feeds both the DT table and the volcano plot
    expect_false(is.null(output$table))
    expect_false(is.null(output$plot))
  })
})

test_that("de_server can export all genes (not just the top-N)", {
  skip_if_not_installed("presto")
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  ct <- sort(unique(as.character(data$cells$celltype)))[[1]]
  shiny::testServer(scroll:::de_server, args = list(data = data), {
    session$setInputs(group = "celltype", ident1 = ct, ident2 = "rest",
                      minpct = 0, topn = 5, lfc = 1, padj = 0.05, labeln = 10,
                      aspect = 1, export_all = FALSE, compute = 1)
    full <- nrow(de_df())
    expect_gt(full, 5)                         # full result has more than the cap
    expect_lte(nrow(table_rows()), 5)          # default CSV is the displayed top-N
    session$setInputs(export_all = TRUE)       # "Export all genes" -> full de_df()
    exported <- if (isTRUE(input$export_all)) de_df() else table_rows()
    expect_equal(nrow(exported), full)
  })
})

test_that("de_server captures a too-few-cells contrast as a friendly message", {
  skip_if_not_installed("presto")
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  cts <- sort(unique(as.character(data$cells$celltype)))
  # a 2-cell subset (all one ident) trips presto's >=3-per-side floor
  tiny <- utils::head(which(as.character(data$cells$celltype) == cts[[1]]), 2)
  small <- reactive(data$cells[tiny, , drop = FALSE])
  shiny::testServer(scroll:::de_server, args = list(data = data, cells_r = small), {
    session$setInputs(group = "celltype", ident1 = cts[[1]], ident2 = "rest",
                      minpct = 10, topn = 20, lfc = 1, padj = 0.05, labeln = 10,
                      aspect = 1, compute = 1)
    expect_false(is.null(result()$err))          # error captured, not raised
    expect_match(result()$err, "at least 3 cells")
  })
})
