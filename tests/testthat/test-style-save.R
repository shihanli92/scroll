# Saved plot styles (R/style-save.R): <project>/style.yaml, read at session start and
# written by the Style sheet's Save button.

style_data <- function() {
  data <- scroll:::.scroll_load(test_project())
  dir <- tempfile("style"); dir.create(dir)
  data$dir <- dir
  data
}

test_that("style.yaml round-trips and is read defensively", {
  dir <- tempfile("style"); dir.create(dir)
  st <- list(dimplot = list(theme = list(legend = "bottom"), labels = list(title = "UMAP"),
                            layers = list(Cells = list(size = "1.2")),
                            plot = list(palette = "Set2", aspect = 1.5, labels = FALSE),
                            manual = list(col = c(B = "#1F77B4"))))
  scroll:::.scroll_write_styles(dir, st)
  got <- scroll:::.scroll_read_styles(dir)
  expect_identical(got$dimplot$theme, list(legend = "bottom"))
  expect_identical(got$dimplot$labels$title, "UMAP")
  expect_identical(got$dimplot$layers, list(Cells = list(size = "1.2")))
  expect_identical(got$dimplot$plot$palette, "Set2")
  expect_identical(got$dimplot$plot$aspect, 1.5)
  expect_false(got$dimplot$plot$labels)
  expect_identical(got$dimplot$manual$col, c(B = "#1F77B4"))
  expect_match(readLines(file.path(dir, "style.yaml"))[1], "^# scroll plot styles")

  # junk is dropped; YAML expressions are never evaluated
  writeLines(c("scroll_style: 1", "panels:", "  dimplot:",
               "    theme: {legend: bottom, bogus: 1, font: [1, 2]}",
               "    labels: {title: !expr stop('evaluated')}",
               "    plot: {palette: Set2, 'bad id': 1, deep: {a: 1}}",
               "    manual: {col: {B: notacolour, C: '#000000'}}",
               "  other: 5"), file.path(dir, "style.yaml"))
  got <- scroll:::.scroll_read_styles(dir)
  expect_identical(names(got), "dimplot")
  expect_identical(got$dimplot$theme, list(legend = "bottom"))
  expect_false(identical(got$dimplot$labels$title, "evaluated"))
  expect_identical(names(got$dimplot$plot), "palette")
  expect_identical(got$dimplot$manual$col, c(C = "#000000"))

  # a malformed file warns and reads as "no styles"
  writeLines("panels: [unclosed", file.path(dir, "style.yaml"))
  expect_warning(expect_length(scroll:::.scroll_read_styles(dir), 0), "ignoring")
  expect_length(scroll:::.scroll_read_styles(tempfile()), 0)          # absent file
})

test_that("a panel's look-control ids are found from its style_ui", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  sec <- Filter(function(s) identical(s$id, "dimplot"), scroll:::.scroll_builtin_panels())[[1]]
  ids <- scroll:::.scroll_style_plot_ids(sec, data)
  expect_true(all(c("palette", "labels", "legend", "raster", "aspect") %in% ids))
  expect_length(scroll:::.scroll_style_plot_ids(list(id = "x"), data), 0)
})

test_that("a saved style seeds the panel; Save writes changed plots and merges the file", {
  data <- style_data()
  on.exit(scroll_disconnect(data$con))
  # another session saved a different panel; this one starts from a saved dimplot style
  scroll:::.scroll_write_styles(data$dir, list(
    dimplot = list(theme = list(legend = "top"), scales = list(xmin = "-5", flip = "flip"),
                   plot = list(aspect = 2), manual = list(col = c(B = "#123456"))),
    violin = list(labels = list(title = "Kept"))))
  ctx <- scroll:::.scroll_style_ctx(data)
  rv <- shiny::reactiveVal(list())
  shiny::testServer(scroll:::.scroll_style_server,
                    args = list(id = "dimplot", data = data, rv = rv, all = list(rv),
                                caps = list(limits = c("x", "y")), ctx = ctx,
                                plot_ids = c("aspect", "palette")), {
    expect_identical(rv()$theme, list(legend = "top"))
    expect_identical(rv()$scales, list(xmin = "-5"))        # flip isn't in this panel's caps
    expect_identical(session$userData$scroll_manual_saved[[session$ns("col")]], c(B = "#123456"))
    session$setInputs(style_save = 1)                         # nothing changed yet
    expect_identical(ctx$saved$dimplot$theme, list(legend = "top"))

    # a look control is saved only once changed in the open sheet
    session$setInputs(palette = "Set1", aspect = 2, style_save = 2)
    expect_null(scroll:::.scroll_read_styles(data$dir)$dimplot$plot$palette)
    expect_identical(scroll:::.scroll_read_styles(data$dir)$dimplot$plot$aspect, 2)  # kept as saved
    session$setInputs(style_open = 1)
    rv(list(theme = list(legend = "bottom"), labels = list(title = "Mine")))
    session$setInputs(palette = "Set2", style_save = 3)
    disk <- scroll:::.scroll_read_styles(data$dir)
    expect_identical(disk$dimplot$theme, list(legend = "bottom"))
    expect_identical(disk$dimplot$labels$title, "Mine")
    expect_identical(disk$dimplot$plot$palette, "Set2")
    expect_identical(disk$violin$labels$title, "Kept")       # other panels untouched

    session$setInputs(style_reset = 1)                        # back to defaults ...
    expect_length(rv(), 0)
    expect_length(session$userData$scroll_manual_saved, 0)   # saved Manual colours too
    expect_identical(scroll:::.scroll_read_styles(data$dir)$dimplot$theme$legend, "bottom")
    session$setInputs(palette = "Set1", aspect = NULL, style_save = 4)   # ... kept once saved
    expect_null(scroll:::.scroll_read_styles(data$dir)$dimplot)
    expect_false(is.null(scroll:::.scroll_read_styles(data$dir)$violin))
  })
})

test_that("Save reports an unwritable project instead of failing", {
  skip_on_os("windows")
  skip_if(identical(Sys.info()[["user"]], "root"))
  data <- style_data()
  on.exit({ Sys.chmod(data$dir, "755"); scroll_disconnect(data$con) })
  Sys.chmod(data$dir, "555")
  ctx <- scroll:::.scroll_style_ctx(data)
  rv <- shiny::reactiveVal(list())
  shiny::testServer(scroll:::.scroll_style_server,
                    args = list(id = "dimplot", data = data, rv = rv, all = list(rv), ctx = ctx), {
    rv(list(theme = list(legend = "bottom")))
    expect_no_error(session$setInputs(style_save = 1))
    expect_false(file.exists(file.path(data$dir, "style.yaml")))
  })
  expect_error(scroll:::.scroll_save_styles(data, ctx), "can't write")
})

test_that("config style_save: false hides Save; the sheet shows it otherwise", {
  data <- style_data()
  on.exit(scroll_disconnect(data$con))
  sec <- Filter(function(s) identical(s$id, "dimplot"), scroll:::.scroll_builtin_panels())[[1]]
  expect_match(as.character(scroll:::.scroll_style_sheet(sec, data)), 'id="dimplot-style_save"',
               fixed = TRUE)
  data$config$style_save <- FALSE
  expect_false(grepl("style_save", as.character(scroll:::.scroll_style_sheet(sec, data))))
})

test_that("Manual colours are saved by level name and seed the pickers", {
  ud <- new.env()
  session <- list(userData = ud, ns = function(x) paste0("dimplot-", x))
  ud$scroll_manual_levels <- list(`dimplot-col` = list(levels = c("A", "B"),
                                  defaults = c(A = "#111111", B = "#222222")))
  input <- list(col_1 = "#111111", col_2 = "#FF0000")          # B changed, A left alone
  expect_identical(scroll:::.scroll_manual_snapshot(input, session, list(col = c(Z = "#00FF00"))),
                   list(col = c(B = "#FF0000", Z = "#00FF00")))
  # a picker render seeds from saved colours (only for levels it shows)
  shiny::testServer(function(input, output, session) {
    session$userData$scroll_manual_saved <- list(`proxy-col` = c(B = "#ABCDEF", Q = "#000000"))
    output$m <- shiny::renderUI(scroll:::.scroll_manual_ui(shiny::NS("proxy"), c("A", "B")))
  }, {
    html <- output$m$html
    expect_match(as.character(html), "#ABCDEF", fixed = TRUE)
    expect_identical(session$userData$scroll_manual_levels$`proxy-col`$levels, c("A", "B"))
  })
})

test_that("an overwrite rebuild keeps the project's style.yaml", {
  dir <- tempfile("proj")
  obj <- make_test_object()
  build <- function() suppressMessages(utils::capture.output(
    scroll_build(obj, dir, overwrite = TRUE, verbose = FALSE)))
  build()
  scroll:::.scroll_write_styles(dir, list(dimplot = list(theme = list(legend = "top"))))
  build()
  expect_identical(scroll:::.scroll_read_styles(dir)$dimplot$theme$legend, "top")
})
