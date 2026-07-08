# The scroll x reactive coupling, tested as pure logic (no Shiny/Quarto):
# given the active closeread section + the feature-search value, what does the
# single sticky resolve to?

test_that("no feature: section drives the view", {
  data <- scroll:::.scroll_data(test_project())
  res <- scroll:::.scroll_resolve("overview", "", data)
  expect_equal(res$view, "umap_colorby")
  expect_null(res$feature_values)
})

test_that("active feature overrides colouring on any section", {
  data <- scroll:::.scroll_data(test_project())
  res <- scroll:::.scroll_resolve("overview", "MS4A1", data)
  expect_false(is.null(res$feature_values))
  expect_equal(res$feature_name, "MS4A1")
  expect_true(all(c("cell", "value") %in% names(res$feature_values)))
})

test_that("null / unknown active section falls back to the first section", {
  data <- scroll:::.scroll_data(test_project())
  first <- data$config$sections[[1]]$view
  expect_equal(scroll:::.scroll_resolve(NULL, "", data)$view, first)
  expect_equal(scroll:::.scroll_resolve("does-not-exist", "", data)$view, first)
})

test_that("a missing feature does not override the section", {
  data <- scroll:::.scroll_data(test_project())
  res <- scroll:::.scroll_resolve("overview", "NOT_A_GENE", data)
  expect_null(res$feature_values)
})
