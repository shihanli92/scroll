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
                       multiple = TRUE, options = list(plugins = list("remove_button"))),
        # with >1 fill column: combine into one interaction plot, or facet one per
        # column -- as a grid or stacked rows (like the stacked violin)
        conditionalPanel("input['fill'] && input['fill'].length > 1", ns = ns,
          selectInput(ns("filllayout"), "Fill columns",
                      c("Combine levels" = "combine", "Facet (grid)" = "grid",
                        "Facet (stacked rows)" = "stacked"), selected = "combine"))),
      .scroll_group("Appearance",
        selectInput(ns("position"), "Bars",
                    c("Fill (100%)" = "fill", "Stack (counts)" = "stack",
                      "Grouped (dodge)" = "dodge"), selected = "fill"),
        sliderInput(ns("barwidth"), "Bar width", 0.3, 1, 0.8, 0.05),
        sliderInput(ns("outline"), "Bar outline", 0, 1, 0.2, 0.1),
        selectInput(ns("palette"), "Palette", .scroll_cat_palettes()),
        uiOutput(ns("manual")),
        bslib::input_switch(ns("legend"), "Legend", TRUE)),
      .scroll_group("Order & labels",
        selectInput(ns("xorder"), "Order groups",
                    c("Alphabetical" = "alpha", "Total count" = "total",
                      "By fill level" = "level", "Reverse" = "reverse")),
        conditionalPanel("input['xorder'] == 'level'", ns = ns,
          selectInput(ns("orderlevel"), "Order by level", choices = NULL)),
        selectInput(ns("fillorder"), "Order fill",
                    c("Alphabetical" = "alpha", "Abundance" = "abundance", "Reverse" = "reverse")),
        selectInput(ns("labels"), "Segment labels",
                    c("None" = "none", "Count" = "count", "Percent" = "percent")),
        conditionalPanel("input['labels'] != 'none'", ns = ns,
          sliderInput(ns("labelmin"), "Hide labels below (%)", 0, 50, 0, 1),
          sliderInput(ns("labelsize"), "Label size", 1.5, 6, 2.8, 0.1)),
        bslib::input_switch(ns("totals"), "Show group totals", FALSE),
        bslib::input_switch(ns("horizontal"), "Horizontal bars", FALSE)),
      .scroll_group("Layout", .scroll_aspect_input(ns))
    ),
    .scroll_plot_area(ns, csv = TRUE)
  )
}

proportions_server <- function(id, data, cells_r = reactive(data$cells),
                               view_r = reactive(NULL), theme_r = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    m <- data$manifest
    .scroll_bind_view_cats(input, session, view_r, m, "group")                  # single-select x
    .scroll_bind_view_cats(input, session, view_r, m, "fill", multiple = TRUE)  # multi-select interaction
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
    # "order groups by fill level" picker tracks the current fill levels
    observeEvent(lvl_r(), {
      lv <- lvl_r(); cur <- isolate(input$orderlevel)
      updateSelectInput(session, "orderlevel", choices = lv,
                        selected = if (!is.null(cur) && cur %in% lv) cur else lv[[1]])
    })
    data_r <- reactive({
      req(input$group, input$fill)
      list(cells = cells_r(), group_by = input$group, fill_by = input$fill)
    })
    cosmetic_r <- .scroll_cosmetic(reactive(
      list(theme = theme_r(), palette = input$palette, position = input$position,
           fill_layout = input$filllayout %||% "combine", bar_width = input$barwidth %||% 0.8,
           outline = input$outline %||% 0.2, x_order = input$xorder %||% "alpha",
           x_order_level = input$orderlevel, fill_order = input$fillorder %||% "alpha",
           labels = input$labels %||% "none", label_min = input$labelmin %||% 0,
           label_size = input$labelsize %||% 2.8, totals = isTRUE(input$totals),
           horizontal = isTRUE(input$horizontal),
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

