test_that("scroll_preview_panel mounts a single registered panel", {
  on.exit(scroll_reset_panels(), add = TRUE)
  dir <- test_project()

  register_plot_panel(
    "prev_demo", label = "Preview demo", title = "Preview demo",
    controls = list(scroll_input_column("grp", "Group", type = "categorical")),
    plot = function(cells, input, data) {
      ggplot2::ggplot(cells, ggplot2::aes(.data[[input$grp]])) + ggplot2::geom_bar()
    })

  app <- scroll_preview_panel("prev_demo", dir)
  expect_s3_class(app, "shiny.appobj")
})

test_that("scroll_preview_panel errors clearly on an unknown id", {
  dir <- test_project()
  expect_error(scroll_preview_panel("does_not_exist", dir), "No panel with id")
  expect_error(scroll_preview_panel(c("a", "b"), dir), "single")
})

test_that("scroll_preview_panel can preview a built-in panel", {
  dir <- test_project()
  expect_s3_class(scroll_preview_panel("dotplot", dir), "shiny.appobj")
})
