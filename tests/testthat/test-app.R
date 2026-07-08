# App-layer logic: the data handle and the DimPlot module server (no browser).

test_that(".scroll_load exposes cells, manifest, config and query helpers", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  expect_true(all(c("cells", "manifest", "config", "con", "query1", "queryN") %in% names(data)))
  hit <- data$query1("RNA", "MS4A1")
  expect_true(all(c("cell", "value") %in% names(hit)))
})

test_that("dimplot_server renders a plot from its controls", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  shiny::testServer(scroll:::dimplot_server, args = list(data = data), {
    session$setInputs(reduction = "umap", colorby = "celltype", palette = "Tableau 10",
                      size = 0.6, alpha = 0.85, labels = TRUE, split = "")
    expect_false(is.null(output$plot))
    # switching to a numeric color-by swaps the palette family (no error)
    session$setInputs(colorby = "nCount_RNA")
    expect_false(is.null(output$plot))
  })
})
