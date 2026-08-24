# Global theme controls: .scroll_ggtheme() builds a trailing theme() override, wired
# into every plot (RNA views + modality panels) and a no-op at the Default selection.

test_that(".scroll_ggtheme is a no-op by default and overrides only what changed", {
  expect_null(scroll:::.scroll_ggtheme(NULL))
  expect_null(scroll:::.scroll_ggtheme(list(font = "", legend = "", grid = "")))

  th <- scroll:::.scroll_ggtheme(list(legend = "bottom", grid = "off", font = "16"))
  expect_s3_class(th, "theme")
  expect_identical(th$legend.position, "bottom")
  expect_s3_class(th$panel.grid.major, "element_blank")

  # border / axes / background overrides
  th2 <- scroll:::.scroll_ggtheme(list(border = "on", axes = "hide", bg = "grey"))
  expect_s3_class(th2$panel.border, "element_rect")
  expect_s3_class(th2$axis.text, "element_blank")
  expect_identical(th2$panel.background$fill, "grey95")
  expect_s3_class(scroll:::.scroll_ggtheme(list(border = "off"))$panel.border, "element_blank")

  # axis-title toggle + x-label rotation
  th3 <- scroll:::.scroll_ggtheme(list(axis_titles = "hide", angle = "45"))
  expect_s3_class(th3$axis.title, "element_blank")
  expect_equal(th3$axis.text.x$angle, 45)

  p <- ggplot2::ggplot(mtcars, ggplot2::aes(.data$mpg, .data$wt)) + ggplot2::geom_point()
  expect_no_error(ggplot2::ggplot_build(p + scroll:::.scroll_ggtheme(NULL)))   # +NULL safe
  b <- ggplot2::ggplot_build(p + th)
  expect_identical(b$plot$theme$legend.position, "bottom")
})

test_that("the active-theme reactive maps the sidebar inputs (blank = default)", {
  shiny::reactiveConsole(TRUE); on.exit(shiny::reactiveConsole(FALSE))
  input <- list(scroll_theme_font = "", scroll_theme_legend = "bottom", scroll_theme_grid = "")
  gs <- scroll:::.scroll_active_theme(input)()
  expect_null(gs$font); expect_identical(gs$legend, "bottom"); expect_null(gs$grid)
})

test_that("panels accept a theme_r and apply it (RNA scatter + modality)", {
  data <- scroll:::.scroll_load(vdj_test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  theme_r <- shiny::reactive(list(legend = "bottom", grid = "off", font = "14"))

  shiny::testServer(scroll:::dimplot_server, args = list(data = data, theme_r = theme_r), {
    session$setInputs(reduction = "umap", colorby = "celltype", palette = "Tableau 10",
                      size = 0.6, alpha = 0.85, labels = TRUE, legend = TRUE, split = "",
                      aspect = 1, highlight = character(0), raster = FALSE)
    b <- ggplot2::ggplot_build(export_r())
    expect_identical(b$plot$theme$legend.position, "bottom")     # global override wins
  })

  cp <- Filter(function(x) identical(x$id, "clone_overview"),
               scroll:::.scroll_assemble_panels(data$manifest))[[1]]
  shiny::testServer(cp$server, args = list(data = data, theme_r = theme_r), {
    session$setInputs(view = "Rank-abundance", clone_col = "clone_id", group = "group",
                      group_levels = character(0), palette = "Tableau 10")
    b <- ggplot2::ggplot_build(plot_r())
    expect_identical(b$plot$theme$legend.position, "bottom")     # themed via the builder
  })
})
