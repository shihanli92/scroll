# Spatial support: the baked tissue-image asset, the manifest images block +
# embedding kind, panel gating via the `when` predicate, and rendering.

test_that("spatial_spec captures defaults", {
  s <- spatial_spec()
  expect_s3_class(s, "scroll_spatial_spec")
  expect_identical(s$name, "spatial")
  expect_true(s$image_asset)
  expect_null(s$image)
})

test_that("baked image asset round-trips (write -> read)", {
  dir <- spatial_test_project()
  man <- scroll_manifest(dir)
  expect_false(is.null(man$images))
  expect_true(file.exists(file.path(dir, man$images$spatial$file)))
  expect_identical(man$embeddings$spatial$kind, "spatial")

  data <- scroll:::.scroll_load(dir)
  on.exit(scroll_disconnect(data$con), add = TRUE)
  img <- scroll:::.scroll_spatial_image(data, "spatial")
  expect_false(is.null(img))
  expect_identical(img$width, 60L)
  expect_identical(img$height, 50L)
  expect_s3_class(img$raster, "raster")
  expect_identical(scroll:::.scroll_spatial_embeddings(data$manifest), "spatial")
})

test_that("the Spatial panel surfaces only for a project with an images block", {
  sp_man <- scroll_manifest(spatial_test_project())
  keep <- vapply(scroll:::.scroll_assemble_panels(sp_man), `[[`, "", "id")
  expect_true("spatial" %in% keep)

  rna_man <- scroll_manifest(test_project())          # no images block
  drop <- vapply(scroll:::.scroll_assemble_panels(rna_man), `[[`, "", "id")
  expect_false("spatial" %in% drop)
})

test_that("the Spatial panel renders by metadata and by gene, with/without image", {
  data <- scroll:::.scroll_load(spatial_test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  p <- Filter(function(x) identical(x$id, "spatial"),
              scroll:::.scroll_assemble_panels(data$manifest))[[1]]
  gene <- scroll:::.scroll_features_of(data$manifest, data$manifest$default_assay)[[1]]
  shiny::testServer(p$server, args = list(data = data), {
    session$setInputs(mode = "Metadata", meta = "celltype", size = 1.4, image = TRUE)
    expect_error(force(output$plot), NA)
    session$setInputs(mode = "Gene", gene = gene, image = FALSE)
    expect_error(force(output$plot), NA)
  })
})

test_that("brush zooms the Spatial panel and double-click resets", {
  data <- scroll:::.scroll_load(spatial_test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  p <- Filter(function(x) identical(x$id, "spatial"),
              scroll:::.scroll_assemble_panels(data$manifest))[[1]]
  shiny::testServer(p$server, args = list(data = data), {
    session$setInputs(mode = "Metadata", meta = "celltype", size = 1.4, image = TRUE)
    force(output$plot)
    session$setInputs(brush = list(xmin = 10, xmax = 30, ymin = 10, ymax = 25))
    b <- ggplot2::ggplot_build(plot_r())          # coord limits follow the brush
    expect_equal(b$layout$panel_params[[1]]$x.range, c(10, 30), tolerance = 1e-6)
    session$setInputs(dblclick = list(x = 1, y = 1))
    expect_false(isTRUE(all.equal(
      ggplot2::ggplot_build(plot_r())$layout$panel_params[[1]]$x.range, c(10, 30))))
  })
})

test_that("a spatial project without a tissue image still surfaces the panel", {
  # imaging platforms (Xenium/CosMx) carry cell centroids but no H&E raster: the
  # spatial embedding exists (kind: spatial) with no `images` block.
  m <- list(scroll_version = "0.0",
            embeddings = list(fov = list(dims = 2, kind = "spatial")),
            assays = list(RNA = list(features = c("g1"))),
            meta = list(ct = list(type = "categorical", levels = c("a", "b"))))
  ids <- vapply(scroll:::.scroll_assemble_panels(m), `[[`, "", "id")
  expect_true("spatial" %in% ids)
})

test_that("imaging-based (Xenium-style FOV) build extracts centroids, no image", {
  skip_if_not(exists("CreateFOV", where = asNamespace("SeuratObject")),
              "SeuratObject too old for CreateFOV")
  dir <- file.path(tempdir(), "scroll-xenium-proj")
  suppressWarnings(suppressMessages(scroll_build(
    make_xenium_object(), dir, assays = "Xenium", meta_cols = "celltype",
    spatial = spatial_spec(), overwrite = TRUE)))
  man <- scroll_manifest(dir)
  expect_identical(man$embeddings$spatial$kind, "spatial")
  expect_identical(as.integer(man$embeddings$spatial$n_covered), 200L)
  expect_null(man$images)                       # no H&E raster for Xenium
  expect_false(dir.exists(file.path(dir, "spatial")))

  data <- scroll:::.scroll_load(dir)
  on.exit(scroll_disconnect(data$con), add = TRUE)
  ids <- vapply(scroll:::.scroll_assemble_panels(data$manifest), `[[`, "", "id")
  expect_true("spatial" %in% ids)               # panel surfaces without an image
  p <- Filter(function(x) identical(x$id, "spatial"),
              scroll:::.scroll_assemble_panels(data$manifest))[[1]]
  shiny::testServer(p$server, args = list(data = data), {
    session$setInputs(mode = "Metadata", meta = "celltype", size = 1)
    expect_error(force(output$plot), NA)
  })
})

test_that("real Visium build path extracts coords + bakes the image", {
  skip_if_not_installed("stxBrain.SeuratData")
  skip_if_not_installed("Seurat")
  suppressMessages(library(Seurat))
  obj <- SeuratData::LoadData("stxBrain", type = "anterior1")
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  dir <- file.path(tempdir(), "scroll-stx-test")
  suppressMessages(scroll_build(obj, dir, assays = "Spatial",
                                meta_cols = "region", spatial = spatial_spec(),
                                overwrite = TRUE))
  man <- scroll_manifest(dir)
  expect_identical(man$embeddings$spatial$kind, "spatial")
  expect_identical(as.integer(man$embeddings$spatial$n_covered), as.integer(ncol(obj)))
  expect_true(file.exists(file.path(dir, man$images$spatial$file)))
  cells <- as.data.frame(arrow::read_parquet(file.path(dir, "cells.parquet")))
  # coordinates fall inside the baked image extent (points overlay the tissue)
  expect_true(all(cells$spatial_1 >= 0 & cells$spatial_1 <= man$images$spatial$width))
  expect_true(all(cells$spatial_2 >= 0 & cells$spatial_2 <= man$images$spatial$height))
})
