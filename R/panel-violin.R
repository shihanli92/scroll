# --- Violin panel -------------------------------------------------------------

violin_ui <- function(id, data) {
  ns <- NS(id)
  m <- data$manifest
  assays <- .scroll_assays_of(m); cats <- .scroll_cat_cols(m); nums <- .scroll_num_cols(m)
  if (!length(cats)) return(.scroll_empty_panel("Needs a categorical metadata column (none found)."))
  # bare control id: conditionalPanel(ns = ns) prepends the module prefix itself
  when <- function(val) sprintf("input['%s'] == '%s'", "source", val)
  bslib::layout_columns(
    col_widths = c(3, 9), class = "scroll-panel",
    div(
      class = "scroll-controls",
      .scroll_group("Value",
        if (length(nums))
          radioButtons(ns("source"), NULL, c("Gene", "Metadata"), inline = TRUE)
        else NULL,
        conditionalPanel(
          if (length(nums)) when("Gene") else "true", ns = ns,
          selectizeInput(ns("feature"), "Gene", choices = NULL, multiple = FALSE,
                         options = list(placeholder = "Search a gene...", maxOptions = 50)),
          if (length(assays) > 1)
            selectInput(ns("assay"), "Assay", assays, selected = m$default_assay)),
        if (length(nums))
          conditionalPanel(when("Metadata"), ns = ns,
            selectInput(ns("metacol"), "Numeric column",
                        stats::setNames(nums, nums), selected = nums[[1]]))),
      .scroll_group("Grouping",
        selectInput(ns("group"), "Group by", stats::setNames(cats, cats), selected = cats[[1]])),
      .scroll_group("Appearance",
        selectInput(ns("palette"), "Palette", .scroll_cat_palettes()),
        uiOutput(ns("manual")),
        bslib::input_switch(ns("jitter"), "Show points", FALSE),
        bslib::input_switch(ns("legend"), "Legend", FALSE)),
      .scroll_group("Layout", .scroll_aspect_input(ns))
    ),
    .scroll_plot_area(ns, "460px")
  )
}

violin_server <- function(id, data, cells_r = reactive(data$cells),
                          view_r = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    m <- data$manifest
    .scroll_bind_view_cats(input, session, view_r, m, "group")
    # numeric-column picker (Metadata mode) is view-aware too, so scoped numeric
    # columns (e.g. module scores) surface only in the view where they're defined
    if (length(.scroll_num_cols(m))) .scroll_bind_view_nums(input, session, view_r, m, "metacol")
    assay <- reactive(input$assay %||% m$default_assay)
    # repopulate the gene list for the active assay; drop a selection that does
    # not exist in the newly chosen assay (else it silently queries empty)
    observeEvent(assay(), {
      feats <- .scroll_features_of(m, assay())
      cur <- isolate(input$feature)
      keep <- if (!is.null(cur) && cur %in% feats) cur else character(0)
      updateSelectizeInput(session, "feature", choices = feats, server = TRUE, selected = keep)
    })
    # per-group color pickers when the palette is "Manual" (levels of Group by)
    lvl_r <- reactive({ req(input$group); .scroll_meta_levels(data, input$group) })
    output$manual <- renderUI(
      if (identical(input$palette, "Manual")) .scroll_manual_ui(session$ns, lvl_r()))
    manual_colors <- reactive(
      if (identical(input$palette, "Manual")) .scroll_manual_colors(input, lvl_r()))
    # DATA reactive: cells + query (cosmetic changes no longer re-hit the store).
    # A numeric metadata column (source = "Metadata") plots directly; otherwise a
    # queried feature.
    data_r <- reactive({
      req(input$group)
      cells <- cells_r()
      if (identical(input$source, "Metadata")) {
        col <- .scroll_nz(input$metacol)
        validate(need(!is.null(col), "Pick a numeric metadata column."))
        list(cells = cells, feature = col, group_by = input$group,
             values = NULL, value_col = col)
      } else {
        feat <- .scroll_nz(input$feature)
        validate(need(!is.null(feat), "Search for a gene to plot its distribution."))
        list(cells = cells, feature = feat, group_by = input$group,
             values = data$query1(assay(), feat), value_col = NULL)
      }
    })
    cosmetic_r <- .scroll_cosmetic(reactive(
      list(palette = input$palette, jitter = isTRUE(input$jitter),
           legend = isTRUE(input$legend), aspect = input$aspect,
           manual_colors = manual_colors())))
    plot_r <- reactive({
      d <- data_r()
      view_violin(d$cells, list(feature = d$feature, group_by = d$group_by,
                                value_col = d$value_col),
                  d$values, cosmetic_r())
    })
    output$plot <- renderPlot(plot_r())
    .scroll_plot_downloads(output, plot_r, id)
  })
}

