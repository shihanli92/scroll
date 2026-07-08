test_that("views render ggplots", {
  dir <- test_project()
  cells <- as.data.frame(arrow::read_parquet(file.path(dir, "cells.parquet")))

  p1 <- view_umap_colorby(cells, list(embedding = "umap", color_by = "celltype"))
  expect_s3_class(p1, "ggplot")

  vals <- data.frame(cell = cells$cell[1:10], value = runif(10))
  p2 <- view_feature_plot(cells, list(embedding = "umap", feature = "CD3D"), vals)
  expect_s3_class(p2, "ggplot")
})

test_that("render_view dispatches and feature override colours by expression", {
  dir <- test_project()
  cells <- as.data.frame(arrow::read_parquet(file.path(dir, "cells.parquet")))

  expect_s3_class(
    render_view("umap_colorby", cells, list(embedding = "umap", color_by = "condition")),
    "ggplot")

  vals <- data.frame(cell = cells$cell[1:5], value = runif(5))
  p <- render_view("umap_colorby", cells, list(embedding = "umap", color_by = "condition"),
                   feature_values = vals, feature_name = "CD3D")
  # feature override routes through the feature plot regardless of the section's view
  expect_true(".expr" %in% names(p$data))
})

test_that("unknown view type errors", {
  dir <- test_project()
  cells <- as.data.frame(arrow::read_parquet(file.path(dir, "cells.parquet")))
  expect_error(render_view("nope", cells, list()), "Unknown view type")
})
