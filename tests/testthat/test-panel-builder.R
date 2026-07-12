# The declarative panel builder (register_plot_panel + scroll_input_* + the
# Layer-2 helpers). Each test resets the session-global registry on exit.

test_that("register_plot_panel assembles into the registry with position/number", {
  on.exit(scroll_reset_panels()); scroll_reset_panels()
  register_plot_panel("cc", plot = function(cells, input, data) ggplot2::ggplot(),
                      controls = list(scroll_input_column("g", "G", "categorical")),
                      label = "CC", after = "de")
  p <- scroll:::.scroll_assemble_panels()
  ids <- vapply(p, `[[`, "", "id")
  expect_true("cc" %in% ids)
  expect_equal(ids[which(ids == "de") + 1L], "cc")          # inserted after 'de'
  cc <- p[[which(ids == "cc")]]
  expect_equal(cc$label, "CC")
  expect_match(cc$num, "^[0-9]{2}$")                         # numbered by position
})

test_that("built panel: placeholder, required message, inline error, render, downloads", {
  on.exit(scroll_reset_panels()); scroll_reset_panels()
  register_plot_panel("cc",
    controls = list(
      scroll_input_column("grp", "Group", "categorical"),
      scroll_input_levels("lv", "Levels", from = "grp", required = TRUE)),
    plot = function(cells, input, data) {
      if ("STOP" %in% input$lv) stop("boom")
      ggplot2::ggplot(cells, ggplot2::aes(.data[[input$grp]])) + ggplot2::geom_bar()
    })
  p <- Filter(function(x) identical(x$id, "cc"), scroll:::.scroll_assemble_panels())[[1]]
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)

  shiny::testServer(p$server, args = list(data = data), {
    msg <- function() tryCatch({ force(output$plot); "" },
                               error = function(e) conditionMessage(e))
    session$setInputs(grp = "celltype")
    expect_match(msg(), "click Compute")                    # pre-compute placeholder
    session$setInputs(lv = character(0), scroll_compute = 1)
    expect_match(msg(), "Select Levels")                    # required field, not a blank
    session$setInputs(lv = "STOP", scroll_compute = 2)
    expect_match(msg(), "boom")                             # plot error surfaced inline
    session$setInputs(lv = c("T", "B"), scroll_compute = 3)
    expect_equal(msg(), "")                                 # valid -> renders, no message
    expect_false(is.null(output$png)); expect_false(is.null(output$pdf))
  })
})

test_that("compute = FALSE makes the plot live (no Compute gate)", {
  on.exit(scroll_reset_panels()); scroll_reset_panels()
  register_plot_panel("live1", compute = FALSE,
    controls = list(scroll_input_column("grp", "Group", "categorical")),
    plot = function(cells, input, data)
      ggplot2::ggplot(cells, ggplot2::aes(.data[[input$grp]])) + ggplot2::geom_bar())
  p <- Filter(function(x) identical(x$id, "live1"), scroll:::.scroll_assemble_panels())[[1]]
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  shiny::testServer(p$server, args = list(data = data), {
    session$setInputs(grp = "celltype")
    expect_error(force(output$plot), NA)                    # renders immediately
  })
})

test_that("a panel needing an absent column type renders an empty-state message", {
  on.exit(scroll_reset_panels()); scroll_reset_panels()
  register_plot_panel("needsnum", plot = function(cells, input, data) ggplot2::ggplot(),
                      controls = list(scroll_input_column("x", "X", "numeric")))
  p <- Filter(function(x) identical(x$id, "needsnum"), scroll:::.scroll_assemble_panels())[[1]]
  # a fake handle whose only metadata column is categorical -> no numeric columns
  fake <- list(manifest = list(meta = list(ct = list(type = "categorical", levels = c("A", "B"))),
                               assays = list(RNA = list(features = c("g1"))), default_assay = "RNA"))
  h <- as.character(p$ui("needsnum", fake))
  expect_match(h, "needs a numeric metadata column")
})

test_that("scroll_columns reads categorical/numeric columns from the manifest", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  expect_true(all(c("condition", "celltype") %in% scroll_columns(data, "categorical")))
  expect_false("condition" %in% scroll_columns(data, "numeric"))
  expect_true(length(scroll_columns(data, "any")) >= length(scroll_columns(data, "categorical")))
})

test_that("scroll_render_plot surfaces errors inline and wires downloads (standalone)", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  srv <- function(id, data) shiny::moduleServer(id, function(input, output, session) {
    scroll_render_plot(output, id, function() {
      if (isTRUE(input$bad)) stop("nope")
      ggplot2::ggplot(data$cells, ggplot2::aes(condition)) + ggplot2::geom_bar()
    })                                                       # event = NULL -> live
  })
  shiny::testServer(srv, args = list(data = data), {
    session$setInputs(bad = FALSE)
    expect_error(force(output$plot), NA)                    # renders
    expect_false(is.null(output$png))
    session$setInputs(bad = TRUE)
    expect_error(force(output$plot), "nope")                # error surfaced, not silent
  })
})
