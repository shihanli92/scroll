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
        # multi-select: pick >1 categorical column to colour by their "a | b" interaction
        selectizeInput(ns("colorby"), "Color by", .scroll_colorby_choices(m),
                       selected = first_cat, multiple = TRUE,
                       options = list(placeholder = "Pick column(s)"))),
      .scroll_group("Groups",
        selectizeInput(ns("highlight"), "Highlight", choices = NULL, multiple = TRUE,
                       options = list(placeholder = "All groups")),
        # background colour for non-highlighted cells (only relevant while highlighting)
        conditionalPanel(
          condition = sprintf("input['%s'] && input['%s'].length > 0",
                              ns("highlight"), ns("highlight")),
          colourpicker::colourInput(ns("bg_color"), "Background colour", value = "grey85"))),
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
                           view_r = reactive(NULL), theme_r = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    m <- data$manifest
    cats0 <- .scroll_cat_cols(m)
    first_cat <- if (length(cats0)) cats0[[1]] else .scroll_num_cols(m)[[1]]

    # When the active view changes: restrict reductions to that view's embeddings
    # (subset views default to the primary sub-embedding) and add the subset's
    # scoped columns to Color-by. Split-by follows the same scoped-column rule.
    # Effective reduction / colour column as deduped reactiveVals. The plot reads
    # these, not input$reduction/input$colorby, so the selectInput updates a View
    # switch triggers (which round-trip through the client) do NOT each re-render:
    # the values are resolved server-side and set once. The view observer runs at
    # high priority so they are updated before the plot renders (avoiding a stale
    # first draw), and manual dropdown changes feed the same reactiveVals (a set to
    # the current value is a no-op, so the echoed round-trip does not re-render).
    red_rv <- reactiveVal(.scroll_default(data, "default_embedding",
                                          .scroll_view_embeddings(m, NULL)[[1]]))
    cb_rv  <- reactiveVal(first_cat)                    # holds 1+ column names (composite if >1)
    observeEvent(view_r(), {
      reds <- .scroll_view_embeddings(m, view_r())
      sel_red <- if (is.null(view_r()))
        .scroll_default(data, "default_embedding", reds[[1]]) else reds[[1]]
      if (!sel_red %in% reds) sel_red <- reds[[1]]
      cb <- .scroll_colorby_choices(m, view_r()); flat <- unlist(cb, use.names = FALSE)
      cur <- cb_rv(); sel_cb <- if (length(cur) && all(cur %in% flat)) cur else flat[[1]]
      red_rv(sel_red); cb_rv(sel_cb)
      updateSelectInput(session, "reduction", choices = reds, selected = sel_red)
      updateSelectizeInput(session, "colorby", choices = cb, selected = sel_cb)
    }, ignoreNULL = FALSE, priority = 100)
    observeEvent(input$reduction, red_rv(input$reduction), ignoreInit = TRUE)
    # keep the last valid selection if the box is cleared, so color_by is never empty
    observeEvent(input$colorby, if (length(input$colorby)) cb_rv(input$colorby),
                 ignoreInit = TRUE, ignoreNULL = FALSE)
    # composite (>1 column) is categorical; a single column follows its manifest type
    is_cat <- reactive({ cb <- cb_rv()
      length(cb) > 1L || identical(m$meta[[cb]]$type, "categorical") })
    levels_of <- reactive({ cb <- cb_rv()
      if (!is_cat()) character(0)
      else if (length(cb) > 1L)                          # composite levels: computed on active cells
        sort(unique(stats::na.omit(.scroll_combo_levels(cells_r(), cb))))
      else .scroll_meta_levels(data, cb) })
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

    # Levels that get a manual colour picker: when highlighting, only the
    # highlighted groups (you can't recolour cells that aren't shown) — this also
    # keeps a big/composite column's Manual pickers to the handful in focus; when
    # not highlighting, every level.
    manual_levels <- reactive({
      hl <- intersect(input$highlight %||% character(0), levels_of())
      if (length(hl)) hl else levels_of()
    })
    # per-group color pickers, shown only when the palette is "Manual"
    output$manual <- renderUI({
      if (!is_cat() || !identical(input$palette, "Manual")) return(NULL)
      .scroll_manual_ui(session$ns, manual_levels())
    })
    manual_colors <- reactive({
      if (!is_cat() || !identical(input$palette, "Manual")) return(NULL)
      .scroll_manual_colors(input, manual_levels())
    })

    # DATA reactive (cells + params) invalidates only on data-input changes.
    data_r <- reactive({
      req(red_rv(), cb_rv())
      cells <- cells_r()
      list(cells = cells, embedding = red_rv(), color_by = cb_rv(),
           n = nrow(cells))
    })
    # COSMETIC reactive, debounced; restyle-only inputs (incl. highlight/manual).
    cosmetic_r <- .scroll_cosmetic(reactive(
      list(theme = theme_r(), palette = input$palette, point_size = input$size, alpha = input$alpha,
           show_labels = isTRUE(input$labels), legend = isTRUE(input$legend),
           split_by = .scroll_nz(input$split), aspect = input$aspect,
           highlight = input$highlight, bg_color = input$bg_color %||% "grey85",
           manual_colors = manual_colors())))
    build <- function(raster) {
      d <- data_r(); st <- cosmetic_r(); st$raster <- raster
      view_umap_colorby(d$cells, list(embedding = d$embedding, color_by = d$color_by), st)
    }
    plot_r   <- .scroll_lazy_plot(input, function() build(.scroll_use_raster(input$raster, data_r()$n)))  # rasterized on screen; recomputed only on-screen
    export_r <- reactive(build(FALSE))                                         # vector for downloads
    output$plot <- renderPlot(plot_r())
    csv_r <- reactive({ d <- data_r(); .scroll_dimplot_source(d$cells, d$embedding, d$color_by) })
    .scroll_plot_downloads(output, export_r, id, csv_r = csv_r)
  })
}

