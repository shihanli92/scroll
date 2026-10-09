# The global filter rail (right sidebar): auto-generated per-column filters that
# narrow the cells every panel sees, composed into the active-cells reactive.

test_that("filter specs auto-generate from the manifest (categorical + numeric)", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  specs <- scroll:::.scroll_filter_specs(data)
  types <- vapply(specs, `[[`, "", "type")
  expect_true("categorical" %in% types && "numeric" %in% types)
  expect_true(all(vapply(specs, function(s) s$col %in% names(data$cells), logical(1))))
  # ids are unique, and stable per column (they survive the list changing)
  expect_equal(length(unique(vapply(specs, `[[`, "", "id"))), length(specs))
  expect_identical(specs[[1]]$id, scroll:::.scroll_filter_id(specs[[1]]$col))
})

test_that("categorical and numeric filters narrow the cell table (AND, untouched = no-op)", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  specs <- scroll:::.scroll_filter_specs(data)

  # no filters applied -> the very same table (no copy)
  expect_identical(scroll:::.scroll_filter_cells(data$cells, list()), data$cells)

  cs <- Filter(function(s) s$type == "categorical", specs)[[1]]
  lvl <- as.character(stats::na.omit(data$cells[[cs$col]])[1])
  out <- scroll:::.scroll_filter_cells(data$cells, list(c(cs, list(value = lvl))))
  expect_true(all(as.character(out[[cs$col]]) == lvl))
  expect_lt(nrow(out), nrow(data$cells))

  ns <- Filter(function(s) s$type == "numeric", specs)[[1]]
  mid <- (ns$min + ns$max) / 2
  out2 <- scroll:::.scroll_filter_cells(data$cells, list(c(ns, list(value = c(ns$min, mid)))))
  v <- out2[[ns$col]]
  expect_true(all(is.na(v) | v <= mid))
  expect_lte(nrow(out2), nrow(data$cells))
})

test_that("config filters = FALSE disables the rail; a character vector curates it", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  data$config$filters <- FALSE
  expect_length(scroll:::.scroll_filter_specs(data), 0)
  expect_null(scroll:::.scroll_filters_ui(data))

  data$config$filters <- "celltype"
  specs <- scroll:::.scroll_filter_specs(data)
  expect_length(specs, 1)
  expect_identical(specs[[1]]$col, "celltype")
})

test_that("filters are deferred behind Apply and cleared on Reset", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  cs <- Filter(function(s) s$type == "categorical", scroll:::.scroll_filter_specs(data))[[1]]
  lvl <- as.character(stats::na.omit(data$cells[[cs$col]])[1])
  shiny::testServer(function(input, output, session) {
    filt_rv <- shiny::reactiveVal(list())
    scroll:::.scroll_bind_filters(input, output, session, data, filt_rv)
  }, {
    do.call(session$setInputs, stats::setNames(list(lvl), cs$id))  # set the control
    expect_length(filt_rv(), 0)                                    # ...but not applied yet
    session$setInputs(scroll_filter_apply = 1)
    expect_identical(filt_rv()[[cs$id]]$value, lvl)                # committed on Apply
    session$setInputs(scroll_filter_reset = 1)
    expect_length(filt_rv(), 0)                                    # cleared on Reset
  })
})

test_that("the control rail renders a collapsible Filters section (no theme)", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  aside <- as.character(scroll:::.scroll_controls_ui(data))
  expect_true(any(grepl("scroll-filters", aside)))       # the sidebar wrapper
  expect_true(any(grepl("<details", aside)))             # collapsible sections
  expect_true(any(grepl("Filters", aside)))
  expect_false(any(grepl("theme", aside)))               # styling moved to the per-plot sheet
  expect_true(any(grepl("Reset all", aside)))
})

test_that("the filters follow the View: a subset's own columns appear in its view", {
  data <- scroll:::.scroll_load(subset_test_project())
  on.exit(scroll_disconnect(data$con))
  cols <- function(v) vapply(scroll:::.scroll_filter_specs(data, v), `[[`, "", "col")
  expect_false("tsub" %in% cols(NULL))                 # whole dataset: no subset columns
  expect_true(all(c("tsub", "tscore") %in% cols("tcell")))
  shiny::testServer(function(input, output, session) {
    view <- shiny::reactiveVal(NULL)
    filt_rv <- shiny::reactiveVal(list())
    scroll:::.scroll_bind_filters(input, output, session, data, filt_rv, view)
  }, {
    view("tcell"); session$flushReact()
    expect_match(as.character(output$scroll_filter_controls$html), "tsub")
    id <- scroll:::.scroll_filter_id("tsub")
    lv <- scroll:::.scroll_meta_levels(data, "tsub")[1]
    do.call(session$setInputs, stats::setNames(list(lv), id))
    session$setInputs(scroll_filter_apply = 1)
    expect_named(filt_rv(), id)
    view(NULL); session$flushReact()                   # leaving the view drops its filter
    expect_length(filt_rv(), 0)
  })
})

test_that("session columns (New column, Add as column) can be filtered on", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  panels <- Filter(function(p) p$id == "dimplot", scroll:::.scroll_builtin_panels())
  shiny::testServer(function(input, output, session) {
    reg <- scroll:::.scroll_wire(input, output, session, data, panels, prewarm = FALSE)
  }, {
    lv <- scroll:::.scroll_meta_levels(data, "celltype")
    session$setInputs(scroll_newcol_from = "celltype", scroll_newcol_name = "grp",
                      scroll_newcol_map = sprintf("[%s]: a, *: b", lv[1]), scroll_newcol_add = 1)
    expect_match(as.character(output$scroll_filter_controls$html), "grp")
    id <- scroll:::.scroll_filter_id("grp")
    do.call(session$setInputs, stats::setNames(list("a"), id))
    session$setInputs(scroll_filter_apply = 1)
    expect_true(all(as.character(reg$cells()$celltype) == lv[1]))
  })
})
