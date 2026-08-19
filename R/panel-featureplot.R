# --- FeaturePlot panel --------------------------------------------------------

featureplot_ui <- function(id, data) {
  ns <- NS(id)
  m <- data$manifest
  assays <- .scroll_assays_of(m); cats <- .scroll_cat_cols(m)
  bslib::layout_columns(
    col_widths = c(3, 9), class = "scroll-panel",
    div(
      class = "scroll-controls",
      .scroll_group("Feature",
        selectizeInput(ns("feature"), "Gene", choices = NULL, multiple = FALSE,
                       options = list(placeholder = "Search a gene...", maxOptions = 50)),
        if (length(.scroll_num_cols(m)))
          selectInput(ns("metacol"), "...or numeric column",
                      c("(use gene)" = "", stats::setNames(.scroll_num_cols(m), .scroll_num_cols(m)))),
        if (length(assays) > 1)
          selectInput(ns("assay"), "Assay", assays, selected = m$default_assay)),
      .scroll_group("Co-expression",
        bslib::input_switch(ns("blend"), "Blend two genes", FALSE),
        conditionalPanel(
          condition = sprintf("input['%s']", ns("blend")),
          selectizeInput(ns("feature2"), "Second gene", choices = NULL, multiple = FALSE,
                         options = list(placeholder = "Search a gene...", maxOptions = 50)),
          sliderInput(ns("blend_threshold"), "Blend threshold", 0, 1, 0.5, 0.05),
          colourpicker::colourInput(ns("blend_c1"), "First gene colour", "#FF0000"),
          colourpicker::colourInput(ns("blend_c2"), "Second gene colour", "#00FF00"))),
      .scroll_group("Embedding",
        selectInput(ns("reduction"), "Reduction", .scroll_view_embeddings(m, NULL),
                    selected = .scroll_default(data, "default_embedding", .scroll_view_embeddings(m, NULL)[[1]]))),
      .scroll_group("Appearance",
        selectInput(ns("palette"), "Palette", .scroll_continuous_palettes, selected = "grey-purple"),
        sliderInput(ns("size"), "Point size", 0.1, 5, 0.7, 0.1),
        sliderInput(ns("clip"), "Color quantiles (%)", 0, 100, c(0, 100), 1),
        bslib::input_switch(ns("order"), "Expressing cells on top", TRUE),
        bslib::input_switch(ns("legend"), "Legend", TRUE),
        bslib::input_switch(ns("raster"), "Rasterize (fast)", TRUE)),
      .scroll_group("Layout",
        selectInput(ns("split"), "Split by", c("None" = "", stats::setNames(cats, cats))),
        .scroll_aspect_input(ns))
    ),
    .scroll_plot_area(ns, "560px", csv = TRUE)
  )
}

featureplot_server <- function(id, data, cells_r = reactive(data$cells),
                               view_r = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    m <- data$manifest
    assay <- reactive(input$assay %||% m$default_assay)
    observeEvent(view_r(), {
      reds <- .scroll_view_embeddings(m, view_r())
      sel_red <- if (is.null(view_r()))
        .scroll_default(data, "default_embedding", reds[[1]]) else reds[[1]]
      if (!sel_red %in% reds) sel_red <- reds[[1]]
      updateSelectInput(session, "reduction", choices = reds, selected = sel_red)
      # numeric-column choices follow the view (scoped numerics appear in-view)
      nums <- .scroll_num_cols(m, view_r()); cur <- input$metacol
      updateSelectInput(session, "metacol",
                        choices = c("(use gene)" = "", stats::setNames(nums, nums)),
                        selected = if (!is.null(cur) && cur %in% nums) cur else "")
    }, ignoreNULL = FALSE)
    .scroll_bind_view_cats(input, session, view_r, m, "split", prepend = c("None" = ""))
    # repopulate the gene list for the active assay; drop a selection that does
    # not exist in the newly chosen assay (else it silently queries empty)
    observeEvent(assay(), {
      feats <- .scroll_features_of(m, assay())
      keepsel <- function(cur) if (!is.null(cur) && cur %in% feats) cur else character(0)
      updateSelectizeInput(session, "feature", choices = feats, server = TRUE,
                           selected = keepsel(isolate(input$feature)))
      updateSelectizeInput(session, "feature2", choices = feats, server = TRUE,
                           selected = keepsel(isolate(input$feature2)))
    })
    # DATA reactive: cells + the (cached) expression query. Cosmetic drags do not
    # invalidate it, so they never re-hit the store.
    data_r <- reactive({
      req(input$reduction)
      cells <- cells_r()
      # blend mode: two genes, co-expression. Takes precedence over the numeric
      # column (blend is gene-vs-gene) but needs both genes chosen.
      if (isTRUE(input$blend)) {
        f1 <- .scroll_nz(input$feature); f2 <- .scroll_nz(input$feature2)
        validate(need(!is.null(f1) && !is.null(f2),
                      "Pick two genes to blend their co-expression."))
        return(list(cells = cells, embedding = input$reduction, blend = TRUE,
                    feature = f1, feature2 = f2,
                    values = data$query1(assay(), f1),
                    values2 = data$query1(assay(), f2), n = nrow(cells)))
      }
      metacol <- .scroll_nz(input$metacol)
      # a numeric metadata column is coloured like expression (same continuous
      # controls: quantile clip, order-on-top, palette) but read from `cells`,
      # not queried; picking one takes precedence over the gene selector.
      if (!is.null(metacol) && metacol %in% names(cells)) {
        vals <- data.frame(cell = cells$cell,
                           value = suppressWarnings(as.numeric(cells[[metacol]])),
                           stringsAsFactors = FALSE)
        return(list(cells = cells, embedding = input$reduction, feature = metacol,
                    values = vals, n = nrow(cells)))
      }
      feat <- .scroll_nz(input$feature)
      validate(need(!is.null(feat),
                    "Search for a gene, or pick a numeric column, to colour the embedding."))
      list(cells = cells, embedding = input$reduction, feature = feat,
           values = data$query1(assay(), feat), n = nrow(cells))
    })
    cosmetic_r <- .scroll_cosmetic(reactive(
      list(palette = input$palette, point_size = input$size,
           order = isTRUE(input$order), legend = isTRUE(input$legend),
           clip = input$clip / 100, split_by = .scroll_nz(input$split),
           aspect = input$aspect)))
    build <- function(raster) {
      d <- data_r(); st <- cosmetic_r(); st$raster <- raster
      if (isTRUE(d$blend))
        return(view_feature_blend(
          d$cells,
          list(embedding = d$embedding, feature1 = d$feature, feature2 = d$feature2,
               blend_threshold = input$blend_threshold %||% 0.5,
               colors = c(input$blend_c1 %||% "#FF0000", input$blend_c2 %||% "#00FF00")),
          d$values, d$values2, st))
      view_feature_plot(d$cells, list(embedding = d$embedding, feature = d$feature),
                        d$values, st)
    }
    plot_r   <- reactive(build(.scroll_use_raster(input$raster, data_r()$n)))
    export_r <- reactive(build(FALSE))
    output$plot <- renderPlot(plot_r())
    csv_r <- reactive({
      d <- data_r()
      out <- .scroll_featureplot_source(d$cells, d$embedding, d$feature, d$values)
      if (isTRUE(d$blend))
        out[[d$feature2]] <- .scroll_expr_vector(.scroll_embedding_xy(d$cells, d$embedding), d$values2)
      out
    })
    .scroll_plot_downloads(output, export_r, id, csv_r = csv_r)
  })
}

