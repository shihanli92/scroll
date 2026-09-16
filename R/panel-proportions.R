# --- Proportions panel --------------------------------------------------------

proportions_ui <- function(id, data) {
  ns <- NS(id)
  cats <- .scroll_cat_cols(data$manifest)
  if (!length(cats)) return(.scroll_empty_panel("Needs categorical metadata columns (none found)."))
  x_default <- if (length(cats) >= 2) cats[[2]] else cats[[1]]
  bslib::layout_columns(
    col_widths = c(3, 9), class = "scroll-panel",
    div(
      class = "scroll-controls",
      .scroll_group("Composition",
        selectInput(ns("group"), "Group by (x)", stats::setNames(cats, cats), selected = x_default),
        # multi-select: pick >1 column to fill by their "a | b" interaction (like DimPlot colour-by)
        selectizeInput(ns("fill"), "Fill by", stats::setNames(cats, cats), selected = cats[[1]],
                       multiple = TRUE, options = list(plugins = list("remove_button")))),
      .scroll_group("Appearance",
        selectInput(ns("palette"), "Palette", .scroll_cat_palettes()),
        uiOutput(ns("manual")),
        bslib::input_switch(ns("normalize"), "Normalize to 100%", TRUE),
        bslib::input_switch(ns("legend"), "Legend", TRUE)),
      .scroll_group("Layout", .scroll_aspect_input(ns))
    ),
    .scroll_plot_area(ns, "460px", csv = TRUE)
  )
}

proportions_server <- function(id, data, cells_r = reactive(data$cells),
                               view_r = reactive(NULL), theme_r = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    m <- data$manifest
    .scroll_bind_view_cats(input, session, view_r, m, "group")   # single-select x
    # Fill by is multi-select (interaction), so bind it like DimPlot's colour-by:
    # keep the in-view categorical columns, defaulting to the first if none survive.
    observeEvent(view_r(), {
      cats_v <- .scroll_cat_cols(m, view_r())
      sel <- intersect(input$fill, cats_v); if (!length(sel) && length(cats_v)) sel <- cats_v[[1]]
      updateSelectizeInput(session, "fill",
                           choices = stats::setNames(cats_v, cats_v), selected = sel)
    }, ignoreNULL = FALSE)
    # per-fill color pickers when the palette is "Manual": levels of the (possibly
    # composite) Fill-by, computed on the active cells for >1 column.
    lvl_r <- reactive({ req(input$fill)
      if (length(input$fill) > 1L)
        sort(unique(stats::na.omit(.scroll_combo_levels(cells_r(), input$fill))))
      else .scroll_meta_levels(data, input$fill) })
    output$manual <- renderUI(
      if (identical(input$palette, "Manual")) .scroll_manual_ui(session$ns, lvl_r()))
    manual_colors <- reactive(
      if (identical(input$palette, "Manual")) .scroll_manual_colors(input, lvl_r()))
    data_r <- reactive({
      req(input$group, input$fill)
      list(cells = cells_r(), group_by = input$group, fill_by = input$fill)
    })
    cosmetic_r <- .scroll_cosmetic(reactive(
      list(theme = theme_r(), palette = input$palette, normalize = isTRUE(input$normalize),
           legend = isTRUE(input$legend), aspect = input$aspect,
           manual_colors = manual_colors())))
    plot_r <- .scroll_lazy_plot(input, function() {
      d <- data_r()
      view_proportions(d$cells, list(group_by = d$group_by, fill_by = d$fill_by), cosmetic_r())
    })
    output$plot <- renderPlot(plot_r())
    csv_r <- reactive({ d <- data_r(); .scroll_proportions_source(d$cells, d$group_by, d$fill_by) })
    .scroll_plot_downloads(output, plot_r, id, csv_r = csv_r)
  })
}

