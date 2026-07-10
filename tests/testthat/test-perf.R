# Performance refactor: rasterized point layer (vector by default / export),
# query LRU memoization, and the data/cosmetic reactive split (cosmetic changes
# must not re-query duckdb).

test_that(".scroll_point_layer is a vector geom_point by default", {
  lyr <- scroll:::.scroll_point_layer(size = 0.6)
  expect_true(inherits(lyr$geom, "GeomPoint"))
})

test_that(".scroll_point_layer rasterizes only when asked and scattermore present", {
  # absent scattermore, the raster request still falls back to a vector geom_point
  if (!requireNamespace("scattermore", quietly = TRUE)) {
    expect_true(inherits(scroll:::.scroll_point_layer(raster = TRUE)$geom, "GeomPoint"))
  } else {
    lyr <- scroll:::.scroll_point_layer(raster = TRUE)
    expect_true(inherits(lyr$geom, "GeomScattermore"))
  }
})

test_that("scatter views default to vector points (no accidental raster default)", {
  cells <- as.data.frame(arrow::read_parquet(file.path(test_project(), "cells.parquet")))
  p <- view_umap_colorby(cells, list(embedding = "umap", color_by = "celltype"))
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], "")
  expect_true(all(geoms == "GeomPoint"))
})

test_that(".scroll_use_raster honours the toggle and the threshold", {
  thr <- scroll:::.scroll_raster_threshold
  expect_true(scroll:::.scroll_use_raster(TRUE, thr + 1))    # on + large -> raster
  expect_false(scroll:::.scroll_use_raster(TRUE, thr - 1))   # on + small -> vector
  expect_false(scroll:::.scroll_use_raster(FALSE, thr + 1))  # off -> vector even if large
  expect_true(scroll:::.scroll_use_raster(NULL, thr + 1))    # unset defaults to on
})

test_that(".scroll_lru get / set / evict / MRU", {
  lru <- scroll:::.scroll_lru(max = 2L)
  expect_null(lru$get("a"))
  lru$set("a", 1); lru$set("b", 2)
  expect_equal(lru$get("a"), 1)          # touch 'a' -> most-recently-used
  lru$set("c", 3)                        # evicts the LRU key ('b'), not 'a'
  expect_equal(lru$get("a"), 1)
  expect_null(lru$get("b"))
  expect_equal(lru$get("c"), 3)
})

test_that("query results are memoized (repeat lookup does not re-hit duckdb)", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  calls <- 0L
  testthat::with_mocked_bindings(
    scroll_query_feature = function(con, assay, feature) {
      calls <<- calls + 1L
      data.frame(cell = character(), value = numeric())
    },
    {
      data$query1("RNA", "CD3D"); data$query1("RNA", "CD3D")
      expect_equal(calls, 1L)            # second call served from the LRU
      data$query1("RNA", "CD8A")
      expect_equal(calls, 2L)            # a different gene does hit duckdb
    },
    .package = "scroll")
})

test_that("cosmetic changes do not re-query; data changes do (featureplot)", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  calls <- 0L
  orig <- data$query1
  data$query1 <- function(assay, feature) { calls <<- calls + 1L; orig(assay, feature) }
  shiny::testServer(scroll:::featureplot_server, args = list(data = data), {
    session$setInputs(reduction = "umap", feature = "CD3D", assay = "RNA",
                      palette = "grey-purple", size = 0.7, order = TRUE, legend = TRUE,
                      clip = c(0, 100), split = "", aspect = 1)
    session$flushReact(); force(output$plot)
    n1 <- calls
    expect_gte(n1, 1L)
    session$setInputs(size = 1.5, palette = "viridis")   # cosmetic only
    session$elapse(200); session$flushReact(); force(output$plot)
    expect_equal(calls, n1)                               # no extra query
    session$setInputs(feature = "CD8A")                   # data change
    session$flushReact(); force(output$plot)
    expect_gt(calls, n1)                                  # re-queried
    # the export reactive is always vector points
    expect_true(inherits(export_r()$layers[[1]]$geom, "GeomPoint"))
  })
})

test_that("featureplot colours a numeric metadata column without querying", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  calls <- 0L
  orig <- data$query1
  data$query1 <- function(assay, feature) { calls <<- calls + 1L; orig(assay, feature) }
  shiny::testServer(scroll:::featureplot_server, args = list(data = data), {
    session$setInputs(reduction = "umap", feature = "", metacol = "nCount_RNA", assay = "RNA",
                      palette = "grey-purple", size = 0.7, order = TRUE, legend = TRUE,
                      clip = c(0, 100), split = "", aspect = 1)
    session$flushReact(); force(output$plot)
    expect_equal(calls, 0L)                          # metadata path never queries duckdb
    expect_true(inherits(export_r()$layers[[1]]$geom, "GeomPoint"))
    # the plotted colour values are the metadata column, not zero-filled expression
    expect_equal(sort(plot_r()$data$.expr), sort(cells_r()$nCount_RNA))
  })
})

test_that("view_dotplot accepts a precomputed assembly", {
  cells <- as.data.frame(arrow::read_parquet(file.path(test_project(), "cells.parquet")))
  feats <- c("CD3D", "CD8A")
  a <- scroll:::.scroll_dotplot_assemble(cells, feats, "celltype", NULL)
  p <- view_dotplot(NULL, list(features = feats, group_by = "celltype"),
                    NULL, list(), assembly = a)
  expect_s3_class(p, c("ggplot", "patchwork", "gg"))
})

test_that("view_biaxial accepts a precomputed pair data.frame", {
  cells <- as.data.frame(arrow::read_parquet(file.path(test_project(), "cells.parquet")))
  params <- list(features = c("nCount_RNA", "nFeature_RNA"), color_by = "celltype")
  df <- scroll:::.scroll_biaxial_df(cells, params)
  p <- view_biaxial(NULL, params, list(), df = df)
  expect_s3_class(p, "ggplot")
  expect_setequal(unique(p$data$pair), "nCount_RNA vs nFeature_RNA")
})
