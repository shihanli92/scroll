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

test_that("scroll_preview_panel never strands the warm-up overlay", {
  # A project configured for the startup warm-up (prewarm_views + subset views) runs
  # NO warm-up when previewed as a single panel, so its overlay must be absent -- else
  # it would hang until the 60s JS failsafe. Copy the fixture so we don't mutate it.
  dir <- file.path(tempdir(), "scroll-preview-warm")
  if (dir.exists(dir)) unlink(dir, recursive = TRUE)
  file.copy(subset_test_project(), tempdir(), recursive = TRUE)
  file.rename(file.path(tempdir(), basename(subset_test_project())), dir)
  cfg <- yaml::read_yaml(file.path(dir, "config.yaml"))
  cfg$prewarm_views <- 3L                       # would gate the overlay on in scroll_app
  yaml::write_yaml(cfg, file.path(dir, "config.yaml"))

  app <- scroll_preview_panel("dimplot", dir)
  expect_s3_class(app, "shiny.appobj")
  ui   <- environment(app$httpHandler)$ui
  html <- as.character(if (is.function(ui)) ui(list()) else ui)
  expect_false(grepl("scroll-warm-overlay", html))
})
