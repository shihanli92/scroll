# --- Ridge plot panel -------------------------------------------------------------

# Seurat's RidgePlot: one density ridge per group for a gene (or a numeric column).
# Same controls as the Violin panel (several genes -> a grid), minus split / stacked.

ridge_ui <- function(id, data) {
  ns <- NS(id)
  m <- data$manifest
  assays <- .scroll_assays_of(m); cats <- .scroll_cat_cols(m); nums <- .scroll_num_cols(m)
  if (!length(cats)) return(.scroll_empty_panel("Needs a categorical metadata column (none found)."))
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
          selectizeInput(ns("feature"), "Genes", choices = NULL, multiple = TRUE,
                         options = list(placeholder = "Add genes, or paste a list...",
                                        maxOptions = 50, plugins = list("remove_button"))),
          .scroll_paste_handler(ns("feature"), ns("feature_paste")),
          if (length(assays) > 1)
            selectInput(ns("assay"), "Assay", assays, selected = m$default_assay),
          # zeros pile up as one tall spike at 0 for most genes
          bslib::input_switch(ns("nonzero"), "Expressing cells only", FALSE)),
        if (length(nums))
          conditionalPanel(when("Metadata"), ns = ns,
            selectInput(ns("metacol"), "Numeric column",
                        stats::setNames(nums, nums), selected = nums[[1]]))),
      .scroll_group("Grouping",
        selectizeInput(ns("group"), "Group by", stats::setNames(cats, cats), selected = cats[[1]],
                       multiple = TRUE, options = list(plugins = list("remove_button"))))
    ),
    .scroll_plot_area(ns, csv = TRUE)
  )
}

# The look controls, shown in the panel's Style sheet (same input ids).
ridge_style_ui <- function(id, data) {
  ns <- NS(id)
  tagList(
    selectInput(ns("palette"), "Palette", .scroll_cat_palettes()),
    uiOutput(ns("manual")),
    bslib::input_switch(ns("legend"), "Legend", FALSE),
    selectInput(ns("ridgemode"), "Layout",
                c("Ridges (overlapping)" = "ridges", "Separate rows" = "separate",
                  "Overlaid (one axis)" = "overlay")),
    conditionalPanel("input['ridgemode'] == 'ridges'", ns = ns,
      sliderInput(ns("ridgescale"), "Ridge height (overlap)", 0.5, 4, 1.4, 0.1)),
    bslib::input_switch(ns("ltgroup"), "Line type by group", FALSE),
    uiOutput(ns("linetypes")),
    .scroll_order_input(ns("order")),
    .scroll_aspect_input(ns))
}

ridge_server <- function(id, data, cells_r = reactive(data$cells),
                         view_r = reactive(NULL), theme_r = reactive(NULL), style_r = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    m <- data$manifest
    .scroll_bind_view_cats(input, session, view_r, m, "group", multiple = TRUE)
    if (length(.scroll_num_cols(m))) .scroll_bind_view_nums(input, session, view_r, m, "metacol")
    assay <- reactive(input$assay %||% m$default_assay)
    .scroll_bind_gene_box(input, session, m, assay, "feature", paste_id = "feature_paste")
    lvl_r <- reactive({ req(input$group)
      .scroll_natural_sort(.scroll_group_levels(cells_r(), input$group)) })
    .scroll_bind_order(session, "order", lvl_r)
    output$manual <- renderUI(
      if (identical(input$palette, "Manual")) .scroll_manual_ui(session$ns, lvl_r()))
    manual_colors <- reactive(
      if (identical(input$palette, "Manual")) .scroll_manual_colors(input, lvl_r()))
    # a line type per group (few groups only: there are six line types)
    lt_n <- length(.SCROLL_LINETYPES)
    output$linetypes <- renderUI({
      if (!isTRUE(input$ltgroup)) return(NULL)
      lv <- lvl_r()
      if (length(lv) > lt_n)
        return(helpText(sprintf("Line types by group need %d or fewer groups (this grouping has %d).",
                                lt_n, length(lv))))
      tagList(lapply(seq_along(lv), function(k)
        div(class = "scroll-filter",
            selectInput(session$ns(paste0("lt_", k)), lv[[k]], .SCROLL_LINETYPES,
                        selected = isolate(input[[paste0("lt_", k)]]) %||% .SCROLL_LINETYPES[[k]]))))
    })
    linetypes_r <- reactive({
      lv <- lvl_r()
      if (!isTRUE(input$ltgroup) || length(lv) > lt_n) return(NULL)
      stats::setNames(vapply(seq_along(lv), function(k)
        input[[paste0("lt_", k)]] %||% .SCROLL_LINETYPES[[k]], ""), lv)
    })
    data_r <- reactive({
      req(input$group)
      cells <- cells_r()
      if (identical(input$source, "Metadata")) {
        col <- .scroll_nz(input$metacol)
        validate(need(!is.null(col), "Pick a numeric metadata column."))
        return(list(cells = cells, mode = "single", feature = col, group_by = input$group,
                    values = NULL, value_col = col, nonzero = FALSE))
      }
      feats <- input$feature; feats <- feats[!is.na(feats) & nzchar(feats)]
      validate(need(length(feats) >= 1, "Search for a gene to plot its distribution."),
               need(length(feats) <= 12, "Select at most 12 genes for the grid."))
      nz <- isTRUE(input$nonzero)
      if (length(feats) > 1)
        return(list(cells = cells, mode = "multi", features = feats, group_by = input$group,
                    values = data$queryN(assay(), feats), value_col = NULL, nonzero = nz))
      list(cells = cells, mode = "single", feature = feats[1], group_by = input$group,
           values = data$query1(assay(), feats[1]), value_col = NULL, nonzero = nz)
    })
    cosmetic_r <- .scroll_cosmetic(reactive(
      list(theme = theme_r(), palette = input$palette, legend = isTRUE(input$legend),
           ridge_scale = input$ridgescale %||% 1.4, ridge_mode = input$ridgemode %||% "ridges",
           linetypes = linetypes_r(), group_order = .scroll_order_value(input$order),
           aspect = input$aspect, manual_colors = manual_colors())))
    plot_r <- .scroll_lazy_plot(input, function() {
      d <- data_r()
      if (identical(d$mode, "multi"))
        return(view_ridge_multi(d$cells, list(group_by = d$group_by, features = d$features,
                                              nonzero = d$nonzero), d$values, cosmetic_r()))
      view_ridge(d$cells, list(feature = d$feature, group_by = d$group_by,
                               value_col = d$value_col, nonzero = d$nonzero),
                 d$values, cosmetic_r())
    }, style_r)
    output$plot <- renderPlot(plot_r())
    csv_r <- reactive({
      d <- data_r()
      if (identical(d$mode, "multi")) {
        out <- data.frame(cell = d$cells$cell, stringsAsFactors = FALSE)
        out[[paste(d$group_by, collapse = " | ")]] <- .scroll_combo_levels(d$cells, d$group_by)
        for (g in d$features)
          out[[g]] <- .scroll_expr_vector(d$cells, d$values[d$values$feature == g, c("cell", "value"), drop = FALSE])
        return(out)
      }
      .scroll_violin_source(d$cells, d$group_by, d$feature, d$values, d$value_col)
    })
    .scroll_plot_downloads(output, plot_r, id, csv_r = csv_r)
  })
}
