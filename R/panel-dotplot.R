# --- DotPlot panel ------------------------------------------------------------

dotplot_ui <- function(id, data) {
  ns <- NS(id)
  m <- data$manifest
  assays <- .scroll_assays_of(m); cats <- .scroll_cat_cols(m)
  if (!length(cats)) return(.scroll_empty_panel("Needs a categorical metadata column (none found)."))
  bslib::layout_columns(
    col_widths = c(3, 9), class = "scroll-panel",
    div(
      class = "scroll-controls",
      .scroll_group("Genes",
        selectizeInput(ns("markers"), "Marker genes", choices = NULL, multiple = TRUE,
                       options = list(placeholder = "Add genes...", maxOptions = 50,
                                      plugins = list("remove_button")))),
      .scroll_group("Grouping",
        selectInput(ns("group"), "Group by", stats::setNames(cats, cats), selected = cats[[1]]),
        if (length(assays) > 1)
          selectInput(ns("assay"), "Assay", assays, selected = m$default_assay)),
      .scroll_group("Appearance",
        bslib::input_switch(ns("scale"), "Scale expression (z-score)", TRUE),
        selectInput(ns("palette"), "Palette", .scroll_continuous_palettes, selected = "magma"),
        sliderInput(ns("dotrange"), "Dot size", 0, 10, c(1, 6), 0.5)),
      .scroll_group("Layout",
        selectInput(ns("cluster"), "Cluster (hclust)",
                    c("Off" = "off", "Rows" = "rows", "Columns" = "columns", "Both" = "both")),
        .scroll_aspect_input(ns))
    ),
    .scroll_plot_area(ns, "520px")
  )
}

dotplot_server <- function(id, data, cells_r = reactive(data$cells),
                           view_r = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    m <- data$manifest
    .scroll_bind_view_cats(input, session, view_r, m, "group")
    assay <- reactive(input$assay %||% m$default_assay)
    defaults <- intersect(unlist(data$config$markers), .scroll_features_of(m, m$default_assay))
    # repopulate markers for the active assay, keeping only those present in it
    observeEvent(assay(), {
      feats <- .scroll_features_of(m, assay())
      cur <- isolate(input$markers) %||% defaults
      updateSelectizeInput(session, "markers", choices = feats, server = TRUE,
                           selected = intersect(cur, feats))
    })
    # DATA reactive: query + aggregation + hclust. `scale` and `cluster` change
    # the aggregation/clustering, so they are DATA inputs (not cosmetic); palette
    # and dot size are cosmetic. No rasterization (dots = features x groups).
    data_r <- reactive({
      req(input$group)
      feats <- input$markers
      validate(need(length(feats) > 0, "Add one or more marker genes to build the panel."))
      cells <- cells_r()
      list(cells = cells, group_by = input$group, features = feats,
           assembly = .scroll_dotplot_assemble(
             cells, feats, input$group, data$queryN(assay(), feats),
             scale = isTRUE(input$scale), cluster = input$cluster))
    })
    cosmetic_r <- .scroll_cosmetic(reactive(
      list(palette = input$palette, dot_size = input$dotrange, aspect = input$aspect)))
    plot_r <- reactive({
      d <- data_r()
      view_dotplot(d$cells, list(group_by = d$group_by, features = d$features),
                   NULL, cosmetic_r(), assembly = d$assembly)
    })
    output$plot <- renderPlot(plot_r())
    .scroll_plot_downloads(output, plot_r, id)
  })
}

