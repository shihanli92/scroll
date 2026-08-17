# --- DimPlot panel ------------------------------------------------------------

dimplot_ui <- function(id, data) {
  ns <- NS(id)
  m <- data$manifest
  reductions <- .scroll_view_embeddings(m, NULL)   # whole-dataset (full) embeddings
  cats <- .scroll_cat_cols(m)
  first_cat <- if (length(cats)) cats[[1]] else .scroll_num_cols(m)[[1]]
  bslib::layout_columns(
    col_widths = c(3, 9), class = "scroll-panel",
    div(
      class = "scroll-controls",
      .scroll_group("Embedding",
        selectInput(ns("reduction"), "Reduction", reductions,
                    selected = .scroll_default(data, "default_embedding", reductions[[1]])),
        selectInput(ns("colorby"), "Color by", .scroll_colorby_choices(m), selected = first_cat)),
      .scroll_group("Groups",
        selectizeInput(ns("highlight"), "Highlight", choices = NULL, multiple = TRUE,
                       options = list(placeholder = "All groups"))),
      .scroll_group("Appearance",
        selectInput(ns("palette"), "Palette", names(.scroll_discrete_palettes)),
        uiOutput(ns("manual")),
        sliderInput(ns("size"), "Point size", 0.1, 5, 0.6, 0.1),
        sliderInput(ns("alpha"), "Opacity", 0.1, 1, 0.85, 0.05),
        bslib::input_switch(ns("labels"), "Cluster labels", TRUE),
        bslib::input_switch(ns("legend"), "Legend", TRUE),
        bslib::input_switch(ns("raster"), "Rasterize (fast)", TRUE)),
      .scroll_group("Layout",
        selectInput(ns("split"), "Split by",
                    c("None" = "", stats::setNames(cats, cats))),
        .scroll_aspect_input(ns))
    ),
    .scroll_plot_area(ns, "460px", csv = TRUE)
  )
}

dimplot_server <- function(id, data, cells_r = reactive(data$cells),
                           view_r = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    m <- data$manifest
    is_cat <- reactive(identical(m$meta[[input$colorby]]$type, "categorical"))
    levels_of <- reactive(if (is_cat()) .scroll_meta_levels(data, input$colorby) else character(0))

    # When the active view changes: restrict reductions to that view's embeddings
    # (subset views default to the primary sub-embedding) and add the subset's
    # scoped columns to Color-by. Split-by follows the same scoped-column rule.
    observeEvent(view_r(), {
      reds <- .scroll_view_embeddings(m, view_r())
      sel_red <- if (is.null(view_r()))
        .scroll_default(data, "default_embedding", reds[[1]]) else reds[[1]]
      if (!sel_red %in% reds) sel_red <- reds[[1]]
      updateSelectInput(session, "reduction", choices = reds, selected = sel_red)
      cb <- .scroll_colorby_choices(m, view_r())
      flat <- unlist(cb, use.names = FALSE)
      cur <- input$colorby
      updateSelectInput(session, "colorby", choices = cb,
                        selected = if (!is.null(cur) && cur %in% flat) cur else flat[[1]])
    }, ignoreNULL = FALSE)
    .scroll_bind_view_cats(input, session, view_r, m, "split", prepend = c("None" = ""))

    observeEvent(input$colorby, {
      if (is_cat()) {
        updateSelectInput(session, "palette",
                          choices = c(names(.scroll_discrete_palettes), "Manual"), selected = "Tableau 10")
        updateSelectizeInput(session, "highlight", choices = levels_of(), selected = character(0))
      } else {
        updateSelectInput(session, "palette",
                          choices = .scroll_continuous_palettes, selected = "viridis")
        updateSelectizeInput(session, "highlight", choices = character(0), selected = character(0))
      }
    })

    # per-group color pickers, shown only when the palette is "Manual"
    output$manual <- renderUI({
      if (!is_cat() || !identical(input$palette, "Manual")) return(NULL)
      .scroll_manual_ui(session$ns, levels_of())
    })
    manual_colors <- reactive({
      if (!is_cat() || !identical(input$palette, "Manual")) return(NULL)
      .scroll_manual_colors(input, levels_of())
    })

    # DATA reactive (cells + params) invalidates only on data-input changes.
    data_r <- reactive({
      req(input$reduction, input$colorby)
      cells <- cells_r()
      list(cells = cells, embedding = input$reduction, color_by = input$colorby,
           n = nrow(cells))
    })
    # COSMETIC reactive, debounced; restyle-only inputs (incl. highlight/manual).
    cosmetic_r <- .scroll_cosmetic(reactive(
      list(palette = input$palette, point_size = input$size, alpha = input$alpha,
           show_labels = isTRUE(input$labels), legend = isTRUE(input$legend),
           split_by = .scroll_nz(input$split), aspect = input$aspect,
           highlight = input$highlight, manual_colors = manual_colors())))
    build <- function(raster) {
      d <- data_r(); st <- cosmetic_r(); st$raster <- raster
      view_umap_colorby(d$cells, list(embedding = d$embedding, color_by = d$color_by), st)
    }
    plot_r   <- reactive(build(.scroll_use_raster(input$raster, data_r()$n)))  # rasterized on screen
    export_r <- reactive(build(FALSE))                                         # vector for downloads
    output$plot <- renderPlot(plot_r())
    csv_r <- reactive({ d <- data_r(); .scroll_dimplot_source(d$cells, d$embedding, d$color_by) })
    .scroll_plot_downloads(output, export_r, id, csv_r = csv_r)
  })
}

