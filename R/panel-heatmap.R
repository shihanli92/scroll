# --- Heatmap panel ------------------------------------------------------------
# Genes x groups tile heatmap (aggregated; the cell dimension collapses at
# aggregation, so it is safe at any cell count) with an optional genes x CELLS
# mode that randomly subsamples cells to a safe cap. Mirrors the DotPlot panel.

heatmap_ui <- function(id, data) {
  ns <- NS(id)
  m <- data$manifest
  assays <- .scroll_assays_of(m); cats <- .scroll_cat_cols(m)
  if (!length(cats)) return(.scroll_empty_panel("Needs a categorical metadata column (none found)."))
  bslib::layout_columns(
    col_widths = c(3, 9), class = "scroll-panel",
    div(
      class = "scroll-controls",
      .scroll_group("Genes",
        # maxOptions must exceed the default gene count, else a large server-side
        # pre-selection (e.g. 100 top-variable genes) fails to bind to the input.
        selectizeInput(ns("markers"), "Genes", choices = NULL, multiple = TRUE,
                       options = list(placeholder = "Add genes, or paste a list...",
                                      maxOptions = 500, plugins = list("remove_button"))),
        .scroll_paste_handler(ns("markers"), ns("marker_paste"))),
      .scroll_group("Grouping",
        # multi-select: group by the "a | b" interaction of several columns
        selectizeInput(ns("group"), "Group by", stats::setNames(cats, cats), selected = cats[[1]],
                       multiple = TRUE, options = list(plugins = list("remove_button"))),
        if (length(assays) > 1)
          selectInput(ns("assay"), "Assay", assays, selected = m$default_assay)),
      .scroll_group("Mode",
        selectInput(ns("mode"), "Columns",
                    c("Groups (aggregated)" = "groups", "Cells (subsampled)" = "cells")),
        conditionalPanel("input['mode'] == 'cells'", ns = ns,
          sliderInput(ns("cellcap"), "Max cells shown", 500, 20000, 5000, 500),
          # order cells within each group so similar cells sit together, no tree
          selectInput(ns("cellorder"), "Cell order",
                      c("Grouped" = "group", "By PC1 (similar together)" = "pc1"))),
        conditionalPanel("input['mode'] == 'groups'", ns = ns,
          selectInput(ns("stat"), "Value",
                      c("Mean expression" = "mean", "% expressing" = "frac")))),
      .scroll_group("Values",
        bslib::input_switch(ns("scale"), "Scale per gene (z-score)", TRUE),
        conditionalPanel("input['scale']", ns = ns,
          sliderInput(ns("clip"), "Clip z at ±", 0.5, 5, 2.5, 0.5))),
      .scroll_group("Appearance",
        selectInput(ns("palette"), "Palette", .scroll_continuous_palettes, selected = "RdBu"),
        bslib::input_switch(ns("legend"), "Legend", TRUE)),
      .scroll_group("Layout",
        selectInput(ns("cluster"), "Cluster (hclust)",
                    c("Off" = "off", "Rows" = "rows", "Columns" = "columns", "Both" = "both")),
        .scroll_aspect_input(ns))
    ),
    .scroll_plot_area(ns, csv = TRUE)
  )
}

heatmap_server <- function(id, data, cells_r = reactive(data$cells),
                           view_r = reactive(NULL), theme_r = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    m <- data$manifest
    # multi-select group-by, view-aware (like the Violin / Composition panels)
    observeEvent(view_r(), {
      cats_v <- .scroll_cat_cols(m, view_r())
      sel <- intersect(input$group, cats_v); if (!length(sel) && length(cats_v)) sel <- cats_v[[1]]
      updateSelectizeInput(session, "group",
                           choices = stats::setNames(cats_v, cats_v), selected = sel)
    }, ignoreNULL = FALSE)
    assay <- reactive(input$assay %||% m$default_assay)
    # default gene set: a project's `heatmap_markers:` (e.g. top variable genes)
    # when set, else the curated `markers:` shared with DotPlot.
    defaults <- intersect(unlist(data$config$heatmap_markers %||% data$config$markers),
                          .scroll_features_of(m, m$default_assay))
    observeEvent(assay(), {
      feats <- .scroll_features_of(m, assay())
      cur <- isolate(input$markers) %||% defaults
      updateSelectizeInput(session, "markers", choices = feats, server = TRUE,
                           selected = intersect(cur, feats))
    })
    observeEvent(input$marker_paste, {                 # pasted delimited gene list
      feats <- .scroll_features_of(m, assay())
      parsed <- .scroll_parse_gene_list(input$marker_paste, feats)
      req(length(parsed$ok) > 0 || length(parsed$missing) > 0)
      sel <- unique(c(isolate(input$markers), parsed$ok))
      updateSelectizeInput(session, "markers", choices = feats, server = TRUE, selected = sel)
      if (length(parsed$missing))
        showNotification(paste("Not in this assay:", paste(parsed$missing, collapse = ", ")),
                         type = "warning", duration = 6)
    })
    observeEvent(input$markers, {
      if (length(input$markers) > 150L)                  # heatmaps read many rows fine
        showNotification(sprintf(
          "%d genes selected - the heatmap may be slow to compute.",
          length(input$markers)), type = "warning", duration = 5)
    }, ignoreInit = TRUE)
    # DATA reactive: query + assemble. mode/stat/scale/cluster/cellcap change the
    # matrix or ordering, so they are DATA inputs; palette/clip/legend are cosmetic.
    data_r <- reactive({
      req(input$group)
      feats <- input$markers
      validate(need(length(feats) > 0, "Add one or more genes to build the heatmap."))
      cells <- cells_r(); el <- data$queryN(assay(), feats); mode <- input$mode %||% "groups"
      assembly <- if (identical(mode, "cells"))
        .scroll_heatmap_cells_assemble(cells, feats, input$group, el,
          scale = if (isTRUE(input$scale)) "zscore" else "none",
          cluster = input$cluster %||% "off", cap = input$cellcap %||% 5000,
          cell_order = input$cellorder %||% "group")
      else .scroll_heatmap_assemble(cells, feats, input$group, el,
          stat = input$stat %||% "mean",
          scale = if (isTRUE(input$scale)) "zscore" else "none",
          cluster = input$cluster %||% "off")
      list(cells = cells, group_by = input$group, features = feats,
           expr_long = el, assembly = assembly)
    })
    cosmetic_r <- .scroll_cosmetic(reactive(
      list(theme = theme_r(), palette = input$palette, clip = input$clip %||% 2.5,
           legend = isTRUE(input$legend), aspect = input$aspect)))
    plot_r <- .scroll_lazy_plot(input, function() {
      d <- data_r()
      view_heatmap(d$cells, list(group_by = d$group_by, features = d$features),
                   NULL, cosmetic_r(), assembly = d$assembly)
    })
    output$plot <- renderPlot(plot_r())
    csv_r <- reactive({ d <- data_r()
      .scroll_heatmap_source(d$cells, d$features, d$group_by, d$expr_long) })
    .scroll_plot_downloads(output, plot_r, id, csv_r = csv_r)
  })
}
