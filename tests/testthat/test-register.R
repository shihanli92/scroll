# The public panel-registration extension point. Each test that registers a
# panel resets the (session-global) registry on exit so state never leaks.

count_ui <- function(id, data) {
  ns <- shiny::NS(id)
  cols <- names(Filter(function(x) identical(x$type, "categorical"), data$manifest$meta))
  shiny::tagList(shiny::selectInput(ns("grp"), "Group", cols),
                 shiny::plotOutput(ns("plot")))
}
count_server <- function(id, data, cells_r = shiny::reactive(data$cells)) {
  shiny::moduleServer(id, function(input, output, session) {
    output$plot <- shiny::renderPlot({
      shiny::req(input$grp); barplot(table(cells_r()[[input$grp]]))
    })
  })
}

test_that("default assembly is the six built-ins, numbered by position", {
  scroll_reset_panels()
  p <- scroll:::.scroll_assemble_panels()
  expect_equal(vapply(p, `[[`, "", "id"),
               c("dimplot", "featureplot", "dotplot", "violin", "proportions", "de"))
  expect_equal(vapply(p, `[[`, "", "num"), sprintf("%02d", 1:6))
})

test_that("register_panel appends a custom panel, and it renders in the page", {
  on.exit(scroll_reset_panels())
  scroll_reset_panels()
  register_panel("counts", count_ui, count_server, label = "Counts",
                 title = "Cells per group")
  p <- scroll:::.scroll_assemble_panels()
  expect_length(p, 7)
  expect_equal(p[[7]]$id, "counts")
  expect_equal(p[[7]]$num, "07")            # numbered after the six built-ins

  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  h <- as.character(scroll:::.scroll_page(data, "T", p))
  expect_match(h, "counts")                 # section + rail anchor present
  expect_match(h, "Cells per group")
  shiny::testServer(count_server, args = list(data = data), {
    session$setInputs(grp = "celltype")
    expect_false(is.null(output$plot))      # custom server wires up
  })
})

test_that("register_panel positions with `before` (top of the app)", {
  on.exit(scroll_reset_panels())
  scroll_reset_panels()
  register_panel("qc", count_ui, count_server, before = "dimplot")
  ids <- vapply(scroll:::.scroll_assemble_panels(), `[[`, "", "id")
  expect_equal(ids[1:2], c("qc", "dimplot"))    # inserted before the first built-in
  expect_equal(ids[[1]], "qc")                  # ...i.e. at the very top
})

test_that("register_panel positions with `after` and replaces by id in place", {
  on.exit(scroll_reset_panels())
  scroll_reset_panels()
  register_panel("mid", count_ui, count_server, after = "dimplot")
  ids <- vapply(scroll:::.scroll_assemble_panels(), `[[`, "", "id")
  expect_equal(ids[1:2], c("dimplot", "mid"))   # inserted right after dimplot

  # re-registering an id overrides in place rather than duplicating
  register_panel("dimplot", count_ui, count_server, label = "MyDim")
  p <- scroll:::.scroll_assemble_panels()
  expect_length(p, 7)                            # 6 built-ins + mid; dimplot swapped
  expect_equal(p[[1]]$id, "dimplot")
  expect_equal(p[[1]]$label, "MyDim")
})

test_that("re-registering the same custom id replaces (no duplicate)", {
  on.exit(scroll_reset_panels())
  scroll_reset_panels()
  register_panel("counts", count_ui, count_server, label = "A")
  register_panel("counts", count_ui, count_server, label = "B")
  p <- scroll:::.scroll_assemble_panels()
  expect_length(p, 7)
  expect_equal(p[[7]]$label, "B")
})

test_that("scroll_reset_panels clears custom registrations", {
  register_panel("tmp", count_ui, count_server)
  expect_gt(length(scroll:::.scroll_assemble_panels()), 6)
  scroll_reset_panels()
  expect_length(scroll:::.scroll_assemble_panels(), 6)
})

test_that("register_panel validates the id and the functions", {
  on.exit(scroll_reset_panels())
  expect_error(register_panel("1bad", count_ui, count_server), "single name")
  expect_error(register_panel("has space", count_ui, count_server), "single name")
  expect_error(register_panel(c("a", "b"), count_ui, count_server), "single name")
  expect_error(register_panel("ok", "notfun", count_server), "must be functions")
})
