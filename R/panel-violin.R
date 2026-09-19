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
          selectizeInput(ns("feature"), "Genes", choices = NULL, multiple = TRUE,
                         options = list(placeholder = "Add genes, or paste a list...",
                                        maxOptions = 50, plugins = list("remove_button"))),
          .scroll_paste_handler(ns("feature"), ns("feature_paste")),
          if (length(assays) > 1)
            selectInput(ns("assay"), "Assay", assays, selected = m$default_assay)),
        if (length(nums))
          conditionalPanel(when("Metadata"), ns = ns,
            selectInput(ns("metacol"), "Numeric column",
                        stats::setNames(nums, nums), selected = nums[[1]]))),
      .scroll_group("Grouping",
        # multi-select: group by the "a | b" interaction of several columns (like DimPlot colour-by)
        selectizeInput(ns("group"), "Group by", stats::setNames(cats, cats), selected = cats[[1]],
                       multiple = TRUE, options = list(plugins = list("remove_button"))),
        # split -> side-by-side coloured violins within each group (Seurat split.by)
        selectInput(ns("split"), "Split by", c("None" = "", stats::setNames(cats, cats)))),
      .scroll_group("Violins",
        sliderInput(ns("vwidth"), "Violin width", 0.3, 1.2, 0.9, 0.05),
        # stacked = one compact row per gene (Seurat stacked violin); points off then
        bslib::input_switch(ns("stack"), "Stacked (one row per gene)", FALSE)),
      .scroll_group("Appearance",
        selectInput(ns("palette"), "Palette", .scroll_cat_palettes()),
        uiOutput(ns("manual")),
        bslib::input_switch(ns("legend"), "Legend", FALSE)),
      .scroll_group("Points",
        bslib::input_switch(ns("jitter"), "Show points", FALSE),
        # point controls only apply when points are shown and not stacked
        conditionalPanel("input['jitter'] == true && input['stack'] != true", ns = ns,
          sliderInput(ns("psize"), "Point size", 0.1, 2, 0.3, 0.1),
          sliderInput(ns("palpha"), "Point opacity", 0.05, 1, 0.3, 0.05),
          sliderInput(ns("pfrac"), "Subsample points (%)", 1, 100, 100, 1))),
      .scroll_group("Layout", .scroll_aspect_input(ns))
    ),
    .scroll_plot_area(ns, csv = TRUE)
  )
}

violin_server <- function(id, data, cells_r = reactive(data$cells),
                          view_r = reactive(NULL), theme_r = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    m <- data$manifest
    # Group by is multi-select (interaction) -> bind it like DimPlot colour-by;
    # Split by is single-select and view-aware via the shared binder.
    observeEvent(view_r(), {
      cats_v <- .scroll_cat_cols(m, view_r())
      sel <- intersect(input$group, cats_v); if (!length(sel) && length(cats_v)) sel <- cats_v[[1]]
      updateSelectizeInput(session, "group",
                           choices = stats::setNames(cats_v, cats_v), selected = sel)
    }, ignoreNULL = FALSE)
    .scroll_bind_view_cats(input, session, view_r, m, "split", prepend = c("None" = ""))
    # numeric-column picker (Metadata mode) is view-aware too, so scoped numeric
    # columns (e.g. module scores) surface only in the view where they're defined
    if (length(.scroll_num_cols(m))) .scroll_bind_view_nums(input, session, view_r, m, "metacol")
    assay <- reactive(input$assay %||% m$default_assay)
    # repopulate the gene list for the active assay; drop a selection that does
    # not exist in the newly chosen assay (else it silently queries empty)
    observeEvent(assay(), {
      feats <- .scroll_features_of(m, assay())
      cur <- isolate(input$feature)
      updateSelectizeInput(session, "feature", choices = feats, server = TRUE,
                           selected = intersect(cur, feats))
    })
    # a delimited gene list pasted into the gene box (see .scroll_paste_handler)
    observeEvent(input$feature_paste, {
      feats <- .scroll_features_of(m, assay())
      parsed <- .scroll_parse_gene_list(input$feature_paste, feats)
      req(length(parsed$ok) > 0 || length(parsed$missing) > 0)
      sel <- unique(c(isolate(input$feature), parsed$ok))
      updateSelectizeInput(session, "feature", choices = feats, server = TRUE, selected = sel)
      if (length(parsed$missing))
        showNotification(paste("Not in this assay:", paste(parsed$missing, collapse = ", ")),
                         type = "warning", duration = 6)
    })
    # per-group color pickers when the palette is "Manual": levels of the fill
    # variable -- the split column when splitting, else the (possibly composite) group.
    lvl_r <- reactive({ req(input$group)
      col <- if (nzchar(.scroll_nz(input$split) %||% "")) input$split else input$group
      if (length(col) > 1L)
        sort(unique(stats::na.omit(.scroll_combo_levels(cells_r(), col))))
      else .scroll_meta_levels(data, col) })
    output$manual <- renderUI(
      if (identical(input$palette, "Manual")) .scroll_manual_ui(session$ns, lvl_r()))
    manual_colors <- reactive(
      if (identical(input$palette, "Manual")) .scroll_manual_colors(input, lvl_r()))
    # DATA reactive: cells + query (cosmetic changes no longer re-hit the store).
    # A numeric metadata column (source = "Metadata") plots directly; otherwise a
    # queried feature.
    data_r <- reactive({
      req(input$group)
      cells <- cells_r(); split_by <- .scroll_nz(input$split); stacked <- isTRUE(input$stack)
      if (identical(input$source, "Metadata")) {
        col <- .scroll_nz(input$metacol)
        validate(need(!is.null(col), "Pick a numeric metadata column."))
        return(list(cells = cells, mode = "single", feature = col, group_by = input$group,
                    split_by = split_by, values = NULL, value_col = col))
      }
      feats <- input$feature; feats <- feats[!is.na(feats) & nzchar(feats)]
      validate(need(length(feats) >= 1, "Search for a gene to plot its distribution."),
               need(length(feats) <= 12, "Select at most 12 genes for the grid."))
      # stacked -> one compact row per gene; multi -> a grid; else a single violin
      if (stacked || length(feats) > 1)
        return(list(cells = cells, mode = if (stacked) "stacked" else "multi",
                    features = feats, group_by = input$group, split_by = split_by,
                    values = data$queryN(assay(), feats), value_col = NULL))
      list(cells = cells, mode = "single", feature = feats[1], group_by = input$group,
           split_by = split_by, values = data$query1(assay(), feats[1]), value_col = NULL)
    })
    cosmetic_r <- .scroll_cosmetic(reactive(
      list(theme = theme_r(), palette = input$palette, jitter = isTRUE(input$jitter),
           legend = isTRUE(input$legend), aspect = input$aspect,
           violin_width = input$vwidth %||% 0.9, point_size = input$psize %||% 0.3,
           point_alpha = input$palpha %||% 0.3, point_frac = (input$pfrac %||% 100) / 100,
           manual_colors = manual_colors())))
    plot_r <- .scroll_lazy_plot(input, function() {
      d <- data_r()
      if (identical(d$mode, "stacked"))
        return(view_violin_stacked(d$cells, list(group_by = d$group_by, features = d$features,
                                                 split_by = d$split_by), d$values, cosmetic_r()))
      if (identical(d$mode, "multi"))
        return(view_violin_multi(d$cells, list(group_by = d$group_by, features = d$features,
                                               split_by = d$split_by), d$values, cosmetic_r()))
      view_violin(d$cells, list(feature = d$feature, group_by = d$group_by,
                                value_col = d$value_col, split_by = d$split_by),
                  d$values, cosmetic_r())
    })
    output$plot <- renderPlot(plot_r())
    csv_r <- reactive({
      d <- data_r()
      if (identical(d$mode, "multi") || identical(d$mode, "stacked")) {
        out <- data.frame(cell = d$cells$cell, stringsAsFactors = FALSE)
        out[[paste(d$group_by, collapse = " | ")]] <- .scroll_combo_levels(d$cells, d$group_by)
        if (length(d$split_by))
          out[[paste(d$split_by, collapse = " | ")]] <- .scroll_combo_levels(d$cells, d$split_by)
        for (g in d$features)
          out[[g]] <- .scroll_expr_vector(d$cells, d$values[d$values$feature == g, c("cell", "value"), drop = FALSE])
        return(out)
      }
      .scroll_violin_source(d$cells, d$group_by, d$feature, d$values, d$value_col, d$split_by)
    })
    .scroll_plot_downloads(output, plot_r, id, csv_r = csv_r)
  })
}

