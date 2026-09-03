# Rendered-image caching (bindCache) on the scatter panels. The panels take an
# opt-in `cache`; passing an explicit cachem store exercises the cached path
# (cache = NULL, the default everywhere else, is the unchanged plain renderPlot).

test_that("featureplot: returning to a rendered state serves from cache (no re-query)", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("cachem")
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  calls <- 0L; orig <- data$query1
  data$query1 <- function(assay, feature) { calls <<- calls + 1L; orig(assay, feature) }

  store <- cachem::cache_mem()
  shiny::testServer(scroll:::featureplot_server, args = list(data = data, cache = store), {
    set <- function(f) session$setInputs(
      reduction = "umap", feature = f, assay = "RNA", onscreen = TRUE,
      palette = "grey-purple", size = 0.7, order = TRUE, legend = TRUE,
      clip = c(0, 100), split = "", aspect = 1, raster = TRUE)
    set("CD3D"); session$flushReact(); force(output$plot); n1 <- calls
    expect_gte(n1, 1L)
    set("CD8A"); session$flushReact(); force(output$plot)          # new state -> renders
    expect_gt(calls, n1); n2 <- calls
    set("CD3D"); session$flushReact(); force(output$plot)          # previously rendered state
    expect_equal(calls, n2)                                        # cache hit -> no re-query
  })
})

test_that("dimplot cache key tracks the controls (and is stable otherwise)", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("cachem")
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  shiny::testServer(scroll:::dimplot_server,
                    args = list(data = data, cache = cachem::cache_mem()), {
    session$setInputs(reduction = "umap", palette = "Tableau 10", size = 0.6, alpha = 0.85,
                      labels = TRUE, legend = TRUE, split = "", raster = TRUE, aspect = 1,
                      onscreen = TRUE)
    session$setInputs(colorby = "celltype"); session$elapse(200); k1 <- key_r()
    session$setInputs(colorby = "condition"); session$elapse(200); k2 <- key_r()
    expect_false(identical(k1, k2))                       # colour column change -> new key
    session$setInputs(colorby = "celltype"); session$elapse(200)
    expect_equal(key_r(), k1)                             # back to the same state -> same key
    force(output$plot)                                    # renders through bindCache, no error
    expect_false(is.null(output$plot))
  })
})
