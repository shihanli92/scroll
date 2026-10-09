# Wide CSVs for GraphPad Prism (R/csv-wide.R).

test_that("per-cell values become one column per group, padded with blanks", {
  w <- scroll:::.scroll_wide_columns(list(x = c(1, 2, 3, 4, NA)), c("b", "a", "b", "10", "a"))
  expect_equal(names(w), c("10", "a", "b"))                 # natural order
  expect_equal(w$b, c(1, 3))
  expect_equal(w$a, c(2, NA))                               # NA value dropped, then padded
  w2 <- scroll:::.scroll_wide_columns(list(g1 = 1:2, g2 = 3:4), c("a", "b"), order = "b")
  expect_equal(names(w2), c("g1: b", "g1: a", "g2: b", "g2: a"))   # manual order, per gene
})

test_that("the wide CSV writer handles a table and titled blocks", {
  f <- tempfile(fileext = ".csv")
  scroll:::.scroll_write_wide(data.frame(a = c(1, NA), b = 2:3, check.names = FALSE), f)
  expect_equal(readLines(f), c('"a","b"', '1,2', ',3'))
  scroll:::.scroll_write_wide(list("Mean" = data.frame(gene = "X", A = 1),
                                   "% expressing" = data.frame(gene = "X", A = 50)), f)
  expect_equal(readLines(f), c('"Mean"', '"gene","A"', '"X",1', "", '"% expressing"', '"gene","A"', '"X",50'))
})

test_that("panels build their wide tables", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  lv <- scroll:::.scroll_meta_levels(data, "celltype")
  # proportions: rows = groups, columns = categories, percents sum to 100
  w <- scroll:::.scroll_proportions_wide(data$cells, "celltype", "condition")
  expect_equal(w$celltype, lv)
  expect_equal(unname(rowSums(w[, -1])), rep(100, length(lv)))
  # dotplot: genes x groups, two blocks
  el <- data$queryN("RNA", c("CD3D", "NKG7"))
  w <- scroll:::.scroll_dotplot_wide(data$cells, c("CD3D", "NKG7"), "celltype", el)
  expect_named(w, c("Mean expression (log-normalised)", "% expressing"))
  expect_equal(w[[1]]$gene, c("CD3D", "NKG7"))
  expect_equal(names(w[[2]])[-1], lv)
  expect_true(all(w[[2]][, -1] >= 0 & w[[2]][, -1] <= 100))
  # violin server: wide CSV with one column per group
  shiny::testServer(scroll:::violin_server, args = list(data = data), {
    session$setInputs(group = "celltype", feature = "CD3D", split = "", order = "", onscreen = TRUE)
    w <- wide_r()
    expect_equal(names(w), lv)
    expect_equal(sum(!is.na(as.matrix(w))), nrow(data$cells))
  })
})
