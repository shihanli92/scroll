# --- DotPlot panel ------------------------------------------------------------

# Parse a pasted gene list (space/comma/tab/newline separated) and match it to a
# feature set: exact first, then case-insensitively. Returns the matched features
# in pasted order (`ok`) plus the unmatched tokens (`missing`).
.scroll_parse_gene_list <- function(text, feats) {
  toks <- strsplit(text %||% "", "[[:space:],;]+", perl = TRUE)[[1]]
  toks <- unique(toks[nzchar(toks)])
  if (!length(toks)) return(list(ok = character(0), missing = character(0)))
  lut <- stats::setNames(feats, toupper(feats))
  matched <- ifelse(toks %in% feats, toks, unname(lut[toupper(toks)]))
  list(ok = unique(matched[!is.na(matched)]), missing = toks[is.na(matched)])
}

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
                                      plugins = list("remove_button"))),
        textAreaInput(ns("marker_paste"), NULL, rows = 2,
                      placeholder = "...or paste a list (space, comma, tab or newline separated)"),
        actionButton(ns("marker_set"), "Set from list",
                     class = "btn-outline-secondary btn-sm")),
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
    .scroll_plot_area(ns, "520px", csv = TRUE)
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
    # paste a whitespace/comma/tab/newline-separated gene list -> the selection,
    # in the pasted order. Unknown symbols are matched case-insensitively, then
    # reported. This replaces the current selection (a curated marker panel is
    # usually pasted whole).
    observeEvent(input$marker_set, {
      feats <- .scroll_features_of(m, assay())
      parsed <- .scroll_parse_gene_list(input$marker_paste, feats)
      req(length(parsed$ok) > 0 || length(parsed$missing) > 0)
      updateSelectizeInput(session, "markers", choices = feats, server = TRUE,
                           selected = parsed$ok)
      if (length(parsed$missing))
        showNotification(paste("Not in this assay:", paste(parsed$missing, collapse = ", ")),
                         type = "warning", duration = 6)
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
    csv_r <- reactive(.scroll_dotplot_source(data_r()$assembly))
    .scroll_plot_downloads(output, plot_r, id, csv_r = csv_r)
  })
}

