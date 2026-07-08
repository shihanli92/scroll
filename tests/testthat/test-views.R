read_cells <- function(dir) as.data.frame(arrow::read_parquet(file.path(dir, "cells.parquet")))

test_that("scatter views render ggplots", {
  cells <- read_cells(test_project())
  expect_s3_class(view_umap_colorby(cells, list(embedding = "umap", color_by = "celltype")), "ggplot")
  vals <- data.frame(cell = cells$cell[1:10], value = runif(10))
  expect_s3_class(view_feature_plot(cells, list(embedding = "umap", feature = "CD3D"), vals), "ggplot")
})

test_that("dotplot aggregates fraction and mean per group", {
  cells <- read_cells(test_project())
  feats <- c("CD3D", "MS4A1")
  el <- data.frame(
    feature = c("CD3D", "CD3D", "MS4A1"),
    cell = cells$cell[1:3],
    value = c(1, 2, 3))
  p <- view_dotplot(cells, list(group_by = "celltype", features = feats), el)
  expect_s3_class(p, "ggplot")
  # grid is complete: every feature x group combination present
  expect_equal(nrow(p$data), length(feats) * length(unique(cells$celltype)))
})

test_that("violin and proportions render ggplots", {
  cells <- read_cells(test_project())
  vals <- data.frame(cell = cells$cell[1:20], value = runif(20))
  expect_s3_class(view_violin(cells, list(group_by = "celltype", feature = "CD3D"), vals), "ggplot")
  expect_s3_class(view_proportions(cells, list(group_by = "condition", fill_by = "celltype")), "ggplot")
})

test_that("de_table returns a data.frame (and a message when absent)", {
  df <- data.frame(gene = c("A", "B"), avg_log2FC = c(1.2, -0.5))
  expect_s3_class(view_de_table(df, list(top = 1)), "data.frame")
  expect_equal(nrow(view_de_table(df, list(top = 1))), 1)
  expect_true("message" %in% names(view_de_table(NULL, list(contrast = "x"))))
})

test_that("render_view dispatches by ctx and honours feature override", {
  cells <- read_cells(test_project())
  ctx <- list(view = "umap_colorby", params = list(embedding = "umap", color_by = "condition"),
              cells = cells)
  expect_s3_class(render_view("umap_colorby", ctx), "ggplot")

  vals <- data.frame(cell = cells$cell[1:5], value = runif(5))
  ctx2 <- c(ctx, list(feature_values = vals, feature_name = "CD3D"))
  # feature override routes a colorby section through the feature plot
  expect_true(".expr" %in% names(render_view("umap_colorby", ctx2)$data))

  expect_equal(scroll_view_kind("de_table"), "table")
  expect_equal(scroll_view_kind("violin"), "plot")
  expect_error(render_view("nope", ctx), "Unknown view type")
})

test_that("view cores honour palette / size / order / scale controls", {
  cells <- read_cells(test_project())
  expect_s3_class(
    view_umap_colorby(cells, list(embedding = "umap", color_by = "celltype"),
                      state = list(palette = "Okabe-Ito", point_size = 1.2, alpha = 0.5)),
    "ggplot")
  vals <- data.frame(cell = cells$cell[1:10], value = runif(10))
  expect_s3_class(
    view_feature_plot(cells, list(embedding = "umap", feature = "CD3D"), vals,
                      state = list(palette = "magma", order = FALSE, point_size = 1)),
    "ggplot")
  el <- data.frame(feature = c("CD3D", "CD3D", "MS4A1"), cell = cells$cell[1:3],
                   value = c(1, 2, 3))
  expect_s3_class(
    view_dotplot(cells, list(group_by = "celltype", features = c("CD3D", "MS4A1")), el,
                 state = list(scale = TRUE, dot_size = c(2, 8), palette = "viridis")),
    "ggplot")
})

test_that("discrete colors are deterministic by level name", {
  a <- scroll:::.scroll_discrete_colors(c("B", "T", "NK"))
  b <- scroll:::.scroll_discrete_colors(c("NK", "B"))   # subset, different order
  expect_equal(a[["B"]], b[["B"]])                       # B keeps its color
  expect_equal(a[["NK"]], b[["NK"]])
})

test_that("toggles: split_by facets and overrides group_by", {
  cells <- read_cells(test_project())
  p <- view_umap_colorby(cells, list(embedding = "umap", color_by = "celltype"),
                         state = list(split_by = "condition"))
  expect_s3_class(p, "ggplot")
  # violin group column follows the split_by toggle
  v <- view_violin(cells, list(feature = "CD3D"), NULL, state = list(split_by = "condition"))
  expect_true(all(v$data$group %in% as.character(cells$condition)))
})
