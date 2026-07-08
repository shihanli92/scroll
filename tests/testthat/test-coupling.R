# Per-section coupling, tested as pure logic (no Shiny/Quarto): given a section
# + the global feature search + toggle state, what does that section's sticky
# resolve to?

sec_by_id <- function(dir, id) scroll:::.scroll_sections_by_id(scroll_config(dir))[[id]]

test_that("no feature: the section's own view drives the sticky", {
  data <- scroll:::.scroll_data(test_project())
  ctx <- scroll:::.scroll_section_ctx(sec_by_id(test_project(), "overview"), "", NULL, data)
  expect_equal(ctx$view, "umap_colorby")
  expect_null(ctx$feature_values)
  expect_equal(ctx$kind, "plot")
})

test_that("an active feature recolours a scatter section", {
  data <- scroll:::.scroll_data(test_project())
  ctx <- scroll:::.scroll_section_ctx(sec_by_id(test_project(), "overview"), "MS4A1", NULL, data)
  expect_false(is.null(ctx$feature_values))
  expect_equal(ctx$feature_name, "MS4A1")
})

test_that("dotplot section resolves an expr_long payload; search appends the gene", {
  data <- scroll:::.scroll_data(test_project())
  ctx <- scroll:::.scroll_section_ctx(sec_by_id(test_project(), "markers"), "", NULL, data)
  expect_equal(ctx$view, "dotplot")
  expect_true(is.data.frame(ctx$expr_long))
  ctx2 <- scroll:::.scroll_section_ctx(sec_by_id(test_project(), "markers"), "NKG7", NULL, data)
  expect_true("NKG7" %in% unlist(ctx2$params$features))
})

test_that("subset toggle restricts the cells the section sees", {
  data <- scroll:::.scroll_data(test_project())
  lvl <- as.character(data$cells$celltype[[1]])
  st <- list(subset_col = "celltype", subset_val = lvl)
  ctx <- scroll:::.scroll_section_ctx(sec_by_id(test_project(), "overview"), "", st, data)
  expect_true(all(as.character(ctx$cells$celltype) == lvl))
  expect_lt(nrow(ctx$cells), nrow(data$cells))
})

test_that("de_table section resolves a table read from de/", {
  dir <- test_project()
  de_dir <- file.path(dir, "de"); dir.create(de_dir, showWarnings = FALSE)
  arrow::write_parquet(data.frame(gene = c("A", "B"), avg_log2FC = c(1, -1)),
                       file.path(de_dir, "T_vs_B.parquet"))
  data <- scroll:::.scroll_data(dir)
  section <- list(view = "de_table", params = list(contrast = "T_vs_B"))
  ctx <- scroll:::.scroll_section_ctx(section, "", NULL, data)
  expect_equal(ctx$kind, "table")
  expect_true(is.data.frame(ctx$de_data))
  expect_equal(nrow(ctx$de_data), 2)
})
