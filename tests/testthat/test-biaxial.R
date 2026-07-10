# Biaxial panel: pairwise scatters of numeric metadata columns, coloured by a
# categorical selection.

test_that("view_biaxial builds one facet per column pair", {
  cells <- as.data.frame(arrow::read_parquet(
    file.path(test_project(), "cells.parquet")))
  p <- view_biaxial(cells, list(features = c("nCount_RNA", "nFeature_RNA"),
                                color_by = "celltype"))
  expect_s3_class(p, "ggplot")
  # one pair -> one facet level; every plotted row keeps both coordinates
  expect_setequal(unique(p$data$pair), "nCount_RNA vs nFeature_RNA")
  expect_false(any(is.na(p$data$.x)) || any(is.na(p$data$.y)))
})

test_that("view_biaxial makes all unique pairs from >2 columns", {
  cells <- as.data.frame(arrow::read_parquet(
    file.path(test_project(), "cells.parquet")))
  # duplicate a numeric column so we have three to pair (3 choose 2 = 3 facets)
  cells$q3 <- cells$nCount_RNA / 2
  p <- view_biaxial(cells, list(features = c("nCount_RNA", "nFeature_RNA", "q3"),
                                color_by = "condition"))
  expect_length(unique(p$data$pair), 3)
})

test_that("view_biaxial errors on <2 features or a bad colour column", {
  cells <- as.data.frame(arrow::read_parquet(
    file.path(test_project(), "cells.parquet")))
  expect_error(view_biaxial(cells, list(features = "nCount_RNA", color_by = "celltype")),
               "at least two")
  expect_error(view_biaxial(cells, list(features = c("nCount_RNA", "nFeature_RNA"),
                                        color_by = "nope")), "color_by")
})

test_that(".scroll_default_biaxial prefers hashtag/antibody columns", {
  expect_equal(scroll:::.scroll_default_biaxial(
    c("nCount_RNA", "hto_Hashtag-1", "hto_Hashtag-2", "percent.mt")),
    c("hto_Hashtag-1", "hto_Hashtag-2"))
  # no HTO/ADT columns -> first few numerics
  expect_equal(scroll:::.scroll_default_biaxial(c("a", "b", "c", "d")), c("a", "b", "c"))
})

test_that("biaxial is a registered built-in panel", {
  ids <- vapply(scroll:::.scroll_builtin_panels(), function(p) p$id, "")
  expect_true("biaxial" %in% ids)
})

test_that("biaxial_server renders a plot", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  shiny::testServer(scroll:::biaxial_server, args = list(data = data), {
    session$setInputs(features = c("nCount_RNA", "nFeature_RNA"), colorby = "celltype",
                      palette = "Tableau 10", size = 0.5, alpha = 0.6, legend = TRUE, aspect = 1)
    expect_false(is.null(output$plot))
  })
})
