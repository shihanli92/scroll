# The global filter rail (right sidebar): auto-generated per-column filters that
# narrow the cells every panel sees, composed into the active-cells reactive.

test_that("filter specs auto-generate from the manifest (categorical + numeric)", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  specs <- scroll:::.scroll_filter_specs(data)
  types <- vapply(specs, `[[`, "", "type")
  expect_true("categorical" %in% types && "numeric" %in% types)
  expect_true(all(vapply(specs, function(s) s$col %in% names(data$cells), logical(1))))
  # ids are unique and line up 1:1 with the spec order
  expect_equal(length(unique(vapply(specs, `[[`, "", "id"))), length(specs))
})

test_that("categorical and numeric filters narrow the cell table (AND, untouched = no-op)", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  specs <- scroll:::.scroll_filter_specs(data)

  # empty inputs -> unchanged
  expect_equal(nrow(scroll:::.scroll_filter_cells(data$cells, data, list())), nrow(data$cells))

  cs <- Filter(function(s) s$type == "categorical", specs)[[1]]
  lvl <- as.character(stats::na.omit(data$cells[[cs$col]])[1])
  out <- scroll:::.scroll_filter_cells(data$cells, data, stats::setNames(list(lvl), cs$id))
  expect_true(all(as.character(out[[cs$col]]) == lvl))
  expect_lt(nrow(out), nrow(data$cells))

  ns <- Filter(function(s) s$type == "numeric", specs)[[1]]
  mid <- (ns$min + ns$max) / 2
  out2 <- scroll:::.scroll_filter_cells(data$cells, data, stats::setNames(list(c(ns$min, mid)), ns$id))
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
    scroll:::.scroll_bind_filters(input, session, data, filt_rv)
  }, {
    do.call(session$setInputs, stats::setNames(list(lvl), cs$id))  # set the control
    expect_length(filt_rv(), 0)                                    # ...but not applied yet
    session$setInputs(scroll_filter_apply = 1)
    expect_identical(filt_rv()[[cs$id]], lvl)                      # committed on Apply
    session$setInputs(scroll_filter_reset = 1)
    expect_length(filt_rv(), 0)                                    # cleared on Reset
  })
})

test_that("the control rail renders collapsible Theme + Filters sections", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  aside <- as.character(scroll:::.scroll_controls_ui(data))
  expect_true(any(grepl("scroll-filters", aside)))       # the sidebar wrapper
  expect_true(any(grepl("<details", aside)))             # collapsible sections
  expect_true(any(grepl("Theme", aside)) && any(grepl("Filters", aside)))
  expect_true(any(grepl("Reset all", aside)))
})
