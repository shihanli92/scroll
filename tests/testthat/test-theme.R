# Global theme controls: .scroll_ggtheme() builds a trailing theme() override, wired
# into every plot (RNA views + modality panels) and a no-op at the Default selection.

test_that(".scroll_ggtheme is a no-op by default and overrides only what changed", {
  expect_null(scroll:::.scroll_ggtheme(NULL))
  expect_null(scroll:::.scroll_ggtheme(list(font = "", legend = "", grid = "")))

  th <- scroll:::.scroll_ggtheme(list(legend = "bottom", grid_major = "off", font = "16"))
  expect_s3_class(th, "theme")
  expect_identical(th$legend.position, "bottom")
  expect_s3_class(th$panel.grid.major, "element_blank")
  expect_equal(th$text$size, 16)

  # border / axes / background overrides (bg now a raw picker colour)
  th2 <- scroll:::.scroll_ggtheme(list(border = "on", axes = "hide", bg = "#f0f0f0"))
  expect_s3_class(th2$panel.border, "element_rect")
  expect_s3_class(th2$axis.text, "element_blank")
  expect_identical(th2$panel.background$fill, "#f0f0f0")
  expect_s3_class(scroll:::.scroll_ggtheme(list(border = "off"))$panel.border, "element_blank")

  # axis-title toggle + x/y label rotation + axis line/ticks
  th3 <- scroll:::.scroll_ggtheme(list(axis_titles = "hide", angle = "45", yangle = "90",
                                       axis_line = "show", axis_ticks = "hide"))
  expect_s3_class(th3$axis.title, "element_blank")
  expect_equal(th3$axis.text.x$angle, 45)
  expect_equal(th3$axis.text.y$angle, 90)
  expect_s3_class(th3$axis.line, "element_line")
  expect_s3_class(th3$axis.ticks, "element_blank")

  # the rest of the expansive surface; colours are now raw picker values (hex/transparent)
  th4 <- scroll:::.scroll_ggtheme(list(
    font_family = "serif", text_colour = "#333333", title_style = "bold", title_size = "20",
    legend_dir = "horizontal", legend_title = "hide", legend_key = "none",
    grid_minor = "off", grid_colour = "#cccccc", border_colour = "#999999",
    plot_bg = "transparent", strip_bg = "#eeeeee", margin = "roomy", strip_text_size = "12"))
  expect_identical(th4$text$family, "serif")
  expect_identical(th4$text$colour, "#333333")
  expect_identical(th4$plot.title$face, "bold"); expect_equal(th4$plot.title$size, 20)
  expect_identical(th4$legend.direction, "horizontal")
  expect_s3_class(th4$legend.title, "element_blank")
  expect_s3_class(th4$panel.grid.minor, "element_blank")
  expect_identical(th4$panel.grid.major$colour, "#cccccc")    # gridline colour recolours majors
  expect_identical(th4$panel.border$colour, "#999999")
  expect_identical(th4$plot.background$fill, "transparent")
  expect_identical(th4$strip.background$fill, "#eeeeee")
  expect_false(is.null(th4$plot.margin)); expect_true(inherits(th4$plot.margin, "unit"))
  expect_equal(th4$strip.text$size, 12)

  # base line element: thickness + colour (inherited by grid / axis lines / ticks)
  th5 <- scroll:::.scroll_ggtheme(list(line_size = "thick", line_colour = "#555555"))
  expect_s3_class(th5$line, "element_line")
  expect_equal(as.numeric(th5$line$linewidth), 0.9)
  expect_identical(th5$line$colour, "#555555")

  # axis colour: one picker shared by the axis line + ticks
  th6 <- scroll:::.scroll_ggtheme(list(axis_colour = "#ff0000"))
  expect_identical(th6$axis.line$colour, "#ff0000")
  expect_identical(th6$axis.ticks$colour, "#ff0000")

  p <- ggplot2::ggplot(mtcars, ggplot2::aes(.data$mpg, .data$wt)) + ggplot2::geom_point()
  expect_no_error(ggplot2::ggplot_build(p + scroll:::.scroll_ggtheme(NULL)))   # +NULL safe
  b <- ggplot2::ggplot_build(p + th)
  expect_identical(b$plot$theme$legend.position, "bottom")
})

test_that("theme snapshot maps the sidebar inputs (blank = no override)", {
  input <- list(scroll_theme_font = "", scroll_theme_legend = "bottom",
                scroll_theme_grid_major = "off")
  gs <- scroll:::.scroll_theme_values(input)
  expect_null(gs$font); expect_identical(gs$legend, "bottom")
  expect_identical(gs$grid_major, "off")
})

test_that("theme + filters are deferred: Apply commits, Reset clears", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  ui <- as.character(scroll:::.scroll_controls_ui(data))
  expect_true(any(grepl("scroll_theme_apply", ui)))       # theme Apply button
  expect_true(any(grepl("scroll_filter_apply", ui)))      # filter Apply button

  # theme_rv only changes on Apply, and Reset clears it back to the no-op default
  shiny::testServer(function(input, output, session) {
    theme_rv <- shiny::reactiveVal(list())
    scroll:::.scroll_bind_theme(input, session, data, theme_rv)
  }, {
    session$setInputs(scroll_theme_legend = "bottom")
    expect_length(theme_rv(), 0)                              # not applied yet
    session$setInputs(scroll_theme_apply = 1)
    expect_identical(theme_rv()$legend, "bottom")             # applied on click
    session$setInputs(scroll_theme_reset = 1)
    expect_length(theme_rv(), 0)                              # cleared on reset
  })
})

test_that("panels accept a theme_r and apply it (RNA scatter + modality)", {
  data <- scroll:::.scroll_load(vdj_test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  theme_r <- shiny::reactive(list(legend = "bottom", grid_major = "off", font = "14"))

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
