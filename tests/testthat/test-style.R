# Per-plot Style sheet (R/style.R): the toolbar button, the sheet shell (rendered
# outside the cards), the per-panel style value, and apply-to-all / reset.

test_that("the plot toolbar puts the Style button before the export slider", {
  html <- as.character(scroll:::.scroll_plot_area(shiny::NS("dimplot")))
  expect_match(html, "scroll-style-btn")
  expect_match(html, 'data-sheet="dimplot-sheet"', fixed = TRUE)
  expect_match(html, 'data-open="dimplot-style_open"', fixed = TRUE)
  expect_lt(regexpr("scroll-style-btn", html), regexpr("scroll-size-wrap", html))
})

test_that("sheets render outside the layout grid and the cards", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  panels <- scroll:::.scroll_assemble_panels(data$manifest)
  body <- scroll:::.scroll_body(data, "t", panels)
  # the sheets container is the body's own child, a sibling of .scroll-layout -- so
  # no sheet sits inside .scroll-layout / .scroll-content / a card
  expect_identical(body[[3]]$attribs$class, "scroll-sheets")
  expect_match(body[[2]]$attribs$class, "scroll-layout")
  expect_false(grepl("scroll-sheet", as.character(body[[2]])))
  html <- as.character(body)
  # one sheet per panel, addressed through the panel namespace
  n <- lengths(regmatches(html, gregexpr('class="scroll-sheet"', html, fixed = TRUE)))
  expect_equal(n, length(panels))
  expect_match(html, 'id="dimplot-sheet"', fixed = TRUE)
  expect_match(html, 'id="dimplot-style_body"', fixed = TRUE)
  # the rail no longer carries theme inputs
  expect_false(grepl("scroll_theme_", html))
})

test_that("theme_controls: false drops the theme apply-all from the sheet", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  sec <- scroll:::.scroll_builtin_panels()[[1]]
  on_html <- as.character(scroll:::.scroll_style_sheet(sec, data))
  expect_match(on_html, "style_all")
  data$config$theme_controls <- FALSE
  expect_false(grepl("style_all", as.character(scroll:::.scroll_style_sheet(sec, data))))
})

test_that("sheet inputs feed the panel's style; reset clears it", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  rv <- shiny::reactiveVal(list())
  shiny::testServer(scroll:::.scroll_style_server,
                    args = list(data = data, rv = rv, all = list(rv)), {
    session$setInputs(st_theme_legend = "bottom")        # not inserted yet: ignored
    session$elapse(300)
    expect_length(rv(), 0)
    session$setInputs(style_open = 1)                     # first open inserts the inputs
    vals <- stats::setNames(as.list(rep("", length(scroll:::.scroll_theme_keys()))),
                            paste0("st_theme_", scroll:::.scroll_theme_keys()))
    vals$st_theme_legend <- "bottom"; vals$st_theme_grid_major <- "off"
    do.call(session$setInputs, vals)
    session$elapse(300)
    expect_identical(rv()$theme, list(legend = "bottom", grid_major = "off"))
    session$setInputs(style_reset = 1)
    expect_length(rv(), 0)
  })
})

test_that("a half-bound sheet never blanks a style set elsewhere", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  rv <- shiny::reactiveVal(list(theme = list(legend = "top")))
  shiny::testServer(scroll:::.scroll_style_server,
                    args = list(data = data, rv = rv, all = list(rv)), {
    session$setInputs(style_open = 1, st_theme_font = "")   # only one control bound so far
    session$elapse(300)
    expect_identical(rv()$theme, list(legend = "top"))
  })
})

test_that("apply-to-all copies this panel's theme to every panel", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  a <- shiny::reactiveVal(list(theme = list(legend = "bottom", font = "16")))
  b <- shiny::reactiveVal(list())
  c <- shiny::reactiveVal(list(theme = list(border = "on"), other = 1))
  shiny::testServer(scroll:::.scroll_style_server,
                    args = list(data = data, rv = a, all = list(a = a, b = b, c = c)), {
    session$setInputs(style_all = 1)
    expect_identical(b()$theme, list(font = "16", legend = "bottom"))
    expect_identical(c()$theme, list(font = "16", legend = "bottom"))   # theme replaced
    expect_identical(c()$other, 1)                                      # other families kept
  })
})

test_that("mounting threads style_r and theme_r per panel", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  got <- new.env()
  panels <- list(
    list(id = "p1", label = "P1", server = function(id, data, cells_r, style_r)
      got$s1 <- style_r),
    list(id = "p2", label = "P2", server = function(id, data, cells_r, theme_r)
      got$t2 <- theme_r))
  rvs <- list(p1 = shiny::reactiveVal(list(theme = list(legend = "top"))),
              p2 = shiny::reactiveVal(list(theme = list(legend = "left"))))
  shiny::testServer(function(input, output, session) {
    scroll:::.scroll_mount_panels(data, panels, shiny::reactive(data$cells),
                                  shiny::reactive(NULL), rvs)
  }, {
    expect_identical(got$s1()$theme$legend, "top")
    expect_identical(got$t2()$legend, "left")                # each panel its own theme
  })
})

test_that(".scroll_apply_style sets labels on a ggplot and keeps attributes", {
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3, g = c("a", "b", "c")),
                       ggplot2::aes(x, y, colour = g)) + ggplot2::geom_point()
  attr(p, "scroll_source") <- data.frame(a = 1)
  st <- list(labels = list(title = "T", subtitle = "S", x = "X", legend = "Group"))
  q <- scroll:::.scroll_apply_style(p, st)
  expect_identical(q$labels$title, "T")
  expect_identical(q$labels$subtitle, "S")
  expect_identical(q$labels$x, "X")
  expect_identical(q$labels$colour, "Group")
  expect_identical(attr(q, "scroll_source"), data.frame(a = 1))
  # nothing set -> unchanged object
  expect_identical(scroll:::.scroll_apply_style(p, list()), p)
  expect_identical(scroll:::.scroll_apply_style(p, list(labels = list(title = ""))), p)
})

test_that(".scroll_apply_style annotates a patchwork grid and skips aplot", {
  skip_if_not_installed("patchwork")
  mk <- function() ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
    ggplot2::geom_point()
  pw <- patchwork::wrap_plots(mk(), mk())
  q <- scroll:::.scroll_apply_style(pw, list(labels = list(title = "Grid", x = "X")))
  expect_s3_class(q, "patchwork")
  expect_identical(q$patches$annotation$title, "Grid")
  fake_aplot <- structure(list(), class = "aplot")
  expect_identical(scroll:::.scroll_apply_style(fake_aplot, list(labels = list(title = "t"))),
                   fake_aplot)
})

test_that("labels typed in the sheet reach the panel style", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  rv <- shiny::reactiveVal(list())
  shiny::testServer(scroll:::.scroll_style_server,
                    args = list(data = data, rv = rv, all = list(rv)), {
    session$setInputs(style_open = 1)
    vals <- stats::setNames(as.list(rep("", 6)), paste0("st_labels_", scroll:::.scroll_label_keys()))
    vals$st_labels_title <- "My plot"
    do.call(session$setInputs, vals)
    session$elapse(700)
    expect_identical(rv()$labels, list(title = "My plot"))
  })
})

test_that("dimplot applies its style: title on the export, key changes", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  style_r <- shiny::reactiveVal(list())
  shiny::testServer(scroll:::dimplot_server, args = list(data = data, style_r = style_r), {
    session$setInputs(reduction = "umap", colorby = "celltype", palette = "Tableau 10",
                      size = 0.6, alpha = 0.85, labels = TRUE, legend = TRUE, split = "",
                      aspect = 1, highlight = character(0), raster = FALSE)
    k0 <- key_r()
    style_r(list(labels = list(title = "Cells")))
    expect_identical(export_r()$labels$title, "Cells")
    expect_false(identical(key_r(), k0))
  })
})

test_that("a Compute-gated builder panel restyles without re-computing", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  calls <- 0L
  ps <- scroll:::.scroll_plot_panel_uiserver(
    plot = function(cells, input, data) { calls <<- calls + 1L
      ggplot2::ggplot(cells, ggplot2::aes(umap_1, umap_2)) + ggplot2::geom_point() },
    controls = list(), compute = TRUE)
  style_r <- shiny::reactiveVal(list())
  shiny::testServer(ps$server, args = list(data = data, style_r = style_r), {
    session$setInputs(scroll_compute = 1)
    expect_null(output$plot$error)
    n <- calls
    style_r(list(labels = list(title = "Restyled")))
    session$flushReact()
    expect_identical(calls, n)                               # no recompute
  })
})

test_that("theme colour pickers start at the plot's colours; the default = no override", {
  expect_true(scroll:::.scroll_is_no_colour("#FFFFFF00"))            # legacy marker
  x <- scroll:::.scroll_theme_compact(list(text_colour = "#000000", grid_colour = "#ffffff",
                                           border_colour = "#FF0000", strip_bg = "#FFFFFF00",
                                           legend = "top"))
  expect_identical(x, list(legend = "top", border_colour = "#FF0000"))
  expect_identical(scroll:::.scroll_theme_compact(list(text_colour = "#C7444400")),
                   list(text_colour = "#C74444"))                     # opaque, not ignored
  expect_null(scroll:::.scroll_ggtheme(list(text_colour = "#FFFFFF00")))   # never white text
  skip_if_not_installed("colourpicker")
  html <- as.character(scroll:::.scroll_theme_inputs(shiny::NS("p")))
  expect_match(html, 'id="p-st_theme_text_colour"[^>]*data-init-value="#000000"')
  expect_match(html, 'id="p-st_theme_grid_colour"[^>]*data-init-value="#FFFFFF"')
  expect_false(grepl('data-allow-alpha="true"', html))                # opaque pickers
})

test_that("a field still being typed is not reset when another family syncs", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  rv <- shiny::reactiveVal(list())
  shiny::testServer(scroll:::.scroll_style_server,
                    args = list(data = data, rv = rv, all = list(rv)), {
    session$setInputs(style_open = 1)
    th <- stats::setNames(as.list(rep("", length(scroll:::.scroll_theme_keys()))),
                          paste0("st_theme_", scroll:::.scroll_theme_keys()))
    lb <- stats::setNames(as.list(rep("", 6)), paste0("st_labels_", scroll:::.scroll_label_keys()))
    do.call(session$setInputs, c(th, lb))
    session$elapse(700)
    # user types a title AND changes the legend: the theme snapshot (150 ms) lands
    # before the labels snapshot (500 ms); the title must survive to be stored
    session$setInputs(st_labels_title = "Kept", st_theme_legend = "bottom")
    session$elapse(200)
    expect_identical(rv()$theme$legend, "bottom")
    expect_identical(input$st_labels_title, "Kept")
    session$elapse(500)
    expect_identical(rv()$labels$title, "Kept")
  })
})

test_that("the legend title only targets the aesthetics a plot maps", {
  p <- ggplot2::ggplot(data.frame(x = 1, y = 1, g = "a"), ggplot2::aes(x, y, colour = g)) +
    ggplot2::geom_point()
  expect_identical(scroll:::.scroll_legend_aes(p), "colour")
  expect_no_message(q <- scroll:::.scroll_apply_style(p, list(labels = list(legend = "G"))))
  expect_null(q$labels$fill)
})

test_that("look controls live in each built-in's Style sheet, not its column", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  moved <- list(dimplot = c("palette", "size", "alpha", "labels", "legend", "raster", "aspect"),
                featureplot = c("palette", "size", "clip", "order", "legend", "raster", "aspect"),
                dotplot = c("palette", "dotrange", "aspect"),
                heatmap = c("palette", "legend", "aspect"),
                violin = c("palette", "vwidth", "legend", "jitter", "aspect"),
                proportions = c("palette", "barwidth", "outline", "labels", "horizontal", "aspect"))
  kept <- list(dimplot = c("reduction", "colorby", "highlight", "split"),
               featureplot = c("feature", "blend", "reduction", "split"),
               dotplot = c("markers", "group", "display", "scale", "cluster"),
               heatmap = c("markers", "group", "cellcap", "scale", "cluster"),
               violin = c("feature", "group", "split", "stack"),
               proportions = c("group", "fill", "position"))
  specs <- scroll:::.scroll_builtin_panels()
  for (p in names(moved)) {
    sec <- Filter(function(s) identical(s$id, p), specs)[[1]]
    col <- as.character(sec$ui(p, data))
    sheet <- as.character(sec$style_ui(p, data))
    for (i in moved[[p]]) {
      id <- sprintf('id="%s-%s"', p, i)
      expect_false(grepl(id, col, fixed = TRUE), info = paste(p, i, "in column"))
      expect_true(grepl(id, sheet, fixed = TRUE), info = paste(p, i, "in sheet"))
    }
    for (i in kept[[p]])
      expect_true(grepl(sprintf('id="%s-%s"', p, i), col, fixed = TRUE), info = paste(p, i, "kept"))
  }
})

test_that("scroll_style_input routes a builder control to the sheet", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  us <- scroll:::.scroll_plot_panel_uiserver(
    plot = function(cells, input, data) ggplot2::ggplot(),
    controls = list(scroll_input_column("grp", "Group", "categorical"),
                    scroll_style_input(scroll_input_slider("sz", "Size", 0.1, 3, 1, 0.1))),
    compute = FALSE)
  col <- as.character(us$ui("cp", data)); sheet <- as.character(us$style_ui("cp", data))
  expect_match(col, 'id="cp-grp"', fixed = TRUE)
  expect_false(grepl('id="cp-sz"', col, fixed = TRUE))
  expect_match(sheet, 'id="cp-sz"', fixed = TRUE)
  expect_error(scroll_style_input(list()), "scroll_input")
  # a builder panel with no sheet controls has no style_ui
  expect_null(scroll:::.scroll_plot_panel_uiserver(function(...) NULL, list(), FALSE)$style_ui)
})

test_that("scales: limits keep the coord class, transforms keep labels/expansion", {
  p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point() +
    ggplot2::scale_y_continuous(labels = scales::percent,
                                expand = ggplot2::expansion(mult = c(0, 0.05)))
  q <- scroll:::.scroll_apply_style(p, list(scales = list(ytrans = "sqrt", ybreaks = "6", xmin = "2")))
  expect_identical(q$coordinates$limits$x, c(2, NA))
  sc <- q$scales$get_scales("y")
  expect_match(sc$trans$name, "sqrt")
  expect_identical(sc$n.breaks, 6L)
  expect_equal(sc$expand, ggplot2::expansion(mult = c(0, 0.05)))
  expect_match(ggplot2::ggplot_build(q)$layout$panel_params[[1]]$y$get_labels()[2], "%")
  # original untouched (scales cloned, not mutated in place)
  expect_null(p$scales$get_scales("y")$n.breaks)
  # an equal-aspect embedding stays equal-aspect; a flipped plot stays flipped
  e <- scroll:::.scroll_apply_style(p + ggplot2::coord_fixed(2), list(scales = list(ymax = "30")))
  expect_identical(e$coordinates$ratio, 2)
  f <- scroll:::.scroll_apply_style(p + ggplot2::coord_flip(), list(scales = list(xmin = "1")))
  expect_s3_class(f$coordinates, "CoordFlip")
  # a transform on a discrete axis is skipped, not an error at draw time
  d <- scroll:::.scroll_apply_style(
    ggplot2::ggplot(mtcars, ggplot2::aes(factor(cyl), mpg)) + ggplot2::geom_boxplot(),
    list(scales = list(xtrans = "log10")))
  expect_no_error(ggplot2::ggplot_build(d))
  # facet scales become free; junk limits are ignored
  w <- scroll:::.scroll_apply_style(p + ggplot2::facet_wrap(~cyl),
                                    list(scales = list(facet = "free_y", ymin = "abc")))
  expect_identical(w$facet$params$free, list(x = FALSE, y = TRUE))
  expect_no_error(ggplot2::ggplot_build(w))
})

test_that("the sheet only offers the scale options a panel declares", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  specs <- scroll:::.scroll_builtin_panels()
  caps <- function(id) scroll:::.scroll_style_caps(Filter(function(s) identical(s$id, id), specs)[[1]])
  expect_identical(scroll:::.scroll_scale_keys(caps("dimplot")),
                   c("xmin", "xmax", "xtrans", "ymin", "ymax", "ytrans"))
  expect_length(scroll:::.scroll_scale_keys(caps("dotplot")), 0)
  expect_true(all(c("flip", "ytrans", "ybreaks") %in% scroll:::.scroll_scale_keys(caps("violin"))))
  expect_true("facet" %in% scroll:::.scroll_scale_keys(caps("biaxial")))
  html <- as.character(scroll:::.scroll_scale_inputs(shiny::NS("d"), list(), caps("dimplot")))
  expect_match(html, "Reversed"); expect_false(grepl("Log10", html))
  expect_null(scroll:::.scroll_scale_inputs(shiny::NS("d"), list(), list()))
  # no scales family at all for a panel without caps
  expect_null(scroll:::.scroll_sheet_families(data, list())$scales)
})

test_that("scale inputs typed in the sheet reach the style", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  rv <- shiny::reactiveVal(list())
  caps <- list(limits = "y", trans = list(y = c("sqrt", "reverse")))
  shiny::testServer(scroll:::.scroll_style_server,
                    args = list(data = data, rv = rv, all = list(rv), caps = caps), {
    session$setInputs(style_open = 1)
    session$setInputs(st_scales_ymin = "0", st_scales_ymax = "", st_scales_ytrans = "sqrt")
    session$elapse(700)
    expect_identical(rv()$scales, list(ymin = "0", ytrans = "sqrt"))
  })
})

test_that("a colour picked from the transparent default becomes opaque, not 'no change'", {
  expect_null(scroll:::.scroll_colour_value("#FFFFFF00"))
  expect_identical(scroll:::.scroll_colour_value("#C7444400"), "#C74444")
  expect_identical(scroll:::.scroll_colour_value("#C74444"), "#C74444")
  expect_identical(scroll:::.scroll_colour_value("#C7444480"), "#C7444480")   # real alpha kept
})

test_that("the text colour reaches axis, legend and strip text (not just the title)", {
  p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg, colour = factor(cyl))) + ggplot2::geom_point() +
    scroll:::.scroll_base_theme() + scroll:::.scroll_ggtheme(list(text_colour = "#D62839"))
  th <- p$theme
  for (e in c("axis.text", "axis.title", "legend.text", "legend.title", "strip.text"))
    expect_identical(ggplot2::calc_element(e, ggplot2::theme_grey() + th)$colour, "#D62839", info = e)
})
