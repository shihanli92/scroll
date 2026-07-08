# The scroll x reactive coupling, tested as pure logic (no Shiny/Quarto):
# given the active closeread section + feature-search value + toggle state, what
# does the single sticky resolve to?

test_that("no feature: section drives the view", {
  data <- scroll:::.scroll_data(test_project())
  res <- scroll:::.scroll_resolve("overview", "", NULL, data)
  expect_equal(res$view, "umap_colorby")
  expect_null(res$feature_values)
  expect_equal(res$kind, "plot")
})

test_that("active feature overrides colouring on a scatter section", {
  data <- scroll:::.scroll_data(test_project())
  res <- scroll:::.scroll_resolve("overview", "MS4A1", NULL, data)
  expect_false(is.null(res$feature_values))
  expect_equal(res$feature_name, "MS4A1")
})

test_that("dotplot section resolves an expr_long payload for its panel", {
  data <- scroll:::.scroll_data(test_project())
  res <- scroll:::.scroll_resolve("markers", "", NULL, data)
  expect_equal(res$view, "dotplot")
  expect_true(is.data.frame(res$expr_long))
  expect_true(all(c("feature", "cell", "value") %in% names(res$expr_long)))
  # a searched gene is appended to the panel
  res2 <- scroll:::.scroll_resolve("markers", "NKG7", NULL, data)
  expect_true("NKG7" %in% unlist(res2$params$features))
})

test_that("subset toggle restricts the cells all views see", {
  data <- scroll:::.scroll_data(test_project())
  lvl <- as.character(data$cells$celltype[[1]])
  st <- list(subset_col = "celltype", subset_val = lvl)
  res <- scroll:::.scroll_resolve("overview", "", st, data)
  expect_true(all(as.character(res$cells$celltype) == lvl))
  expect_lt(nrow(res$cells), nrow(data$cells))
})

test_that("null / unknown active section falls back to the first section", {
  data <- scroll:::.scroll_data(test_project())
  first <- data$config$sections[[1]]$view
  expect_equal(scroll:::.scroll_resolve(NULL, "", NULL, data)$view, first)
  expect_equal(scroll:::.scroll_resolve("nope", "", NULL, data)$view, first)
})

test_that("de_table resolves a table read from de/", {
  dir <- test_project()
  # drop a contrast table into de/
  de_dir <- file.path(dir, "de")
  dir.create(de_dir, showWarnings = FALSE)
  arrow::write_parquet(data.frame(gene = c("A", "B"), avg_log2FC = c(1, -1)),
                       file.path(de_dir, "T_vs_B.parquet"))
  data <- scroll:::.scroll_data(dir)
  # inject a de_table section
  data$config$sections <- c(data$config$sections,
    list(list(id = "de", view = "de_table", params = list(contrast = "T_vs_B"))))
  res <- scroll:::.scroll_resolve("de", "", NULL, data)
  expect_equal(res$kind, "table")
  expect_true(is.data.frame(res$de_data))
  expect_equal(nrow(res$de_data), 2)
})
