# The public panel-registration extension point. Each test that registers a
# panel resets the (session-global) registry on exit so state never leaks.

nbuiltin <- length(scroll:::.scroll_builtin_panels())   # number of built-in panels

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

test_that("default assembly is the built-ins, numbered by position", {
  scroll_reset_panels()
  p <- scroll:::.scroll_assemble_panels()
  ids <- vapply(p, `[[`, "", "id")
  # the 8 always-on RNA panels lead, in order; the modality panels (VDJ) follow
  # and are only *gated out* when a manifest is supplied (see gating test below).
  expect_equal(ids[1:8],
               c("dimplot", "featureplot", "biaxial", "dotplot", "violin",
                 "proportions", "de", "pseudobulk"))
  expect_true(all(c("clone_overview", "gene_usage", "cdr3_length", "diversity") %in% ids))
  expect_equal(vapply(p, `[[`, "", "num"), sprintf("%02d", seq_along(p)))
})

test_that("modality panels are gated out for a manifest without their block", {
  scroll_reset_panels()
  rna_man <- scroll_manifest(test_project())          # no vdj block
  ids <- vapply(scroll:::.scroll_assemble_panels(rna_man), `[[`, "", "id")
  expect_equal(ids,
               c("dimplot", "featureplot", "biaxial", "dotplot", "violin",
                 "proportions", "de", "pseudobulk"))
})

test_that("built-in panels drop out when the data can't support them", {
  scroll_reset_panels()
  ids_of <- function(m) vapply(scroll:::.scroll_assemble_panels(m), `[[`, "", "id")

  # a project whose only categorical column has a single level (e.g. a one-region
  # spatial slide): no contrast -> DE and Pseudobulk vanish; one categorical ->
  # Proportions vanishes; DimPlot/FeaturePlot always remain.
  one_region <- list(
    scroll_version = "0.0",             # marks this as a single manifest, not a list
    embeddings = list(spatial = list(dims = 2, kind = "spatial")),
    assays = list(RNA = list(features = c("g1", "g2"))),
    meta = list(region = list(type = "categorical", levels = "anterior"),
                nCount = list(type = "numeric"), nFeature = list(type = "numeric")),
    has_counts = FALSE,
    images = list(spatial = list(embedding = "spatial")))
  ids <- ids_of(one_region)
  expect_true(all(c("dimplot", "featureplot", "spatial") %in% ids))
  expect_false(any(c("de", "pseudobulk", "proportions") %in% ids))
  expect_true("biaxial" %in% ids)        # two numeric columns present

  # add a second, multi-level categorical -> DE returns (a contrast now exists),
  # Proportions returns (two categoricals); still no counts -> Pseudobulk stays out.
  two_region <- one_region
  two_region$meta$region$levels <- c("cortex", "hippocampus")
  two_region$meta$layer <- list(type = "categorical", levels = c("L1", "L2"))
  ids2 <- ids_of(two_region)
  expect_true(all(c("de", "proportions") %in% ids2))
  expect_false("pseudobulk" %in% ids2)   # gated on has_counts
})

test_that("register_panel appends a custom panel, and it renders in the page", {
  on.exit(scroll_reset_panels())
  scroll_reset_panels()
  register_panel("counts", count_ui, count_server, label = "Counts",
                 title = "Cells per group")
  p <- scroll:::.scroll_assemble_panels()
  expect_length(p, nbuiltin + 1L)
  expect_equal(p[[length(p)]]$id, "counts")     # appended after the built-ins
  expect_equal(p[[length(p)]]$num, sprintf("%02d", nbuiltin + 1L))

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
  expect_length(p, nbuiltin + 1L)               # built-ins + mid; dimplot swapped
  expect_equal(p[[1]]$id, "dimplot")
  expect_equal(p[[1]]$label, "MyDim")
})

test_that("re-registering the same custom id replaces (no duplicate)", {
  on.exit(scroll_reset_panels())
  scroll_reset_panels()
  register_panel("counts", count_ui, count_server, label = "A")
  register_panel("counts", count_ui, count_server, label = "B")
  p <- scroll:::.scroll_assemble_panels()
  expect_length(p, nbuiltin + 1L)
  expect_equal(p[[length(p)]]$label, "B")
})

test_that("scroll_reset_panels clears custom registrations", {
  register_panel("tmp", count_ui, count_server)
  expect_gt(length(scroll:::.scroll_assemble_panels()), nbuiltin)
  scroll_reset_panels()
  expect_length(scroll:::.scroll_assemble_panels(), nbuiltin)
})

test_that("register_panel validates the id and the functions", {
  on.exit(scroll_reset_panels())
  expect_error(register_panel("1bad", count_ui, count_server), "single name")
  expect_error(register_panel("has space", count_ui, count_server), "single name")
  expect_error(register_panel(c("a", "b"), count_ui, count_server), "single name")
  expect_error(register_panel("ok", "notfun", count_server), "must be functions")
})
