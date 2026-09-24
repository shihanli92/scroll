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
        selectizeInput(ns("feature"), "Genes", choices = NULL, multiple = TRUE,
                       options = list(placeholder = "Search genes (pick 2+ for a grid)...",
                                      maxOptions = 50)),
        if (length(.scroll_num_cols(m)))
          selectInput(ns("metacol"), "...and/or numeric columns",
                      stats::setNames(.scroll_num_cols(m), .scroll_num_cols(m)),
                      multiple = TRUE),
        if (length(assays) > 1)
          selectInput(ns("assay"), "Assay", assays, selected = m$default_assay)),
      .scroll_group("Co-expression",
        bslib::input_switch(ns("blend"), "Blend (exactly 2 genes)", FALSE)),
      .scroll_group("Embedding",
        selectInput(ns("reduction"), "Reduction", .scroll_view_embeddings(m, NULL),
                    selected = .scroll_default(data, "default_embedding", .scroll_view_embeddings(m, NULL)[[1]]))),
      .scroll_group("Layout",
        selectInput(ns("split"), "Split by", c("None" = "", stats::setNames(cats, cats))))
    ),
    .scroll_plot_area(ns, csv = TRUE)
  )
}

# The look controls, shown in the panel's Style sheet (same input ids).
featureplot_style_ui <- function(id, data) {
  ns <- NS(id)
  tagList(
    selectInput(ns("palette"), "Palette", .scroll_continuous_palettes, selected = "grey-purple"),
    sliderInput(ns("clip"), "Color quantiles (%)", 0, 100, c(0, 100), 1),
    bslib::input_switch(ns("order"), "Expressing cells on top", TRUE),
    bslib::input_switch(ns("legend"), "Legend", TRUE),
    bslib::input_switch(ns("raster"), "Rasterize (fast)", TRUE),
    .scroll_aspect_input(ns),
    # blend colours (the Blend switch itself stays with the data controls)
    conditionalPanel(
      condition = sprintf("input['%s']", ns("blend")),
      sliderInput(ns("blend_threshold"), "Blend threshold", 0, 1, 0.5, 0.05),
      colourpicker::colourInput(ns("blend_c1"), "First gene colour", "#FF0000"),
      colourpicker::colourInput(ns("blend_c2"), "Second gene colour", "#00FF00")))
}

featureplot_server <- function(id, data, cells_r = reactive(data$cells),
                               view_r = reactive(NULL), theme_r = reactive(NULL), style_r = reactive(NULL),
                               cache_key_r = reactive(NULL), cache = NULL) {
  moduleServer(id, function(input, output, session) {
    m <- data$manifest
    assay <- reactive(input$assay %||% m$default_assay)
    # Effective reduction as a deduped reactiveVal so the view-driven selectInput
    # update does not add a second render (see panel-dimplot.R for the rationale).
    red_rv <- reactiveVal(.scroll_default(data, "default_embedding",
                                          .scroll_view_embeddings(m, NULL)[[1]]))
    observeEvent(view_r(), {
      reds <- .scroll_view_embeddings(m, view_r())
      sel_red <- if (is.null(view_r()))
        .scroll_default(data, "default_embedding", reds[[1]]) else reds[[1]]
      if (!sel_red %in% reds) sel_red <- reds[[1]]
      red_rv(sel_red)
      updateSelectInput(session, "reduction", choices = reds, selected = sel_red)
      # numeric-column choices follow the view (scoped numerics appear in-view)
      nums <- .scroll_num_cols(m, view_r()); cur <- input$metacol
      updateSelectInput(session, "metacol",
                        choices = stats::setNames(nums, nums),
                        selected = intersect(cur, nums))
    }, ignoreNULL = FALSE, priority = 100)
    observeEvent(input$reduction, red_rv(input$reduction), ignoreInit = TRUE)
    .scroll_bind_view_cats(input, session, view_r, m, "split", prepend = c("None" = ""))
    .scroll_bind_gene_box(input, session, m, assay, "feature")   # assay-aware gene list
    # DATA reactive: cells + the (cached) expression query. Cosmetic drags do not
    # invalidate it, so they never re-hit the store.
    data_r <- reactive({
      req(red_rv())
      cells <- cells_r()
      clean <- function(x) { x <- x[!is.na(x) & nzchar(x)]; x }
      feats    <- clean(input$feature)                       # genes (queried)
      metacols <- clean(input$metacol); metacols <- metacols[metacols %in% names(cells)]  # numeric columns (from cells)
      # blend mode: co-expression of exactly two selected genes (numeric columns
      # are ignored here — blend is gene-vs-gene).
      if (isTRUE(input$blend)) {
        validate(need(length(feats) == 2,
                      "Blend needs exactly two genes selected."))
        return(list(cells = cells, embedding = red_rv(), blend = TRUE,
                    feature = feats[1], feature2 = feats[2],
                    values = data$query1(assay(), feats[1]),
                    values2 = data$query1(assay(), feats[2]), n = nrow(cells)))
      }
      # genes and numeric columns display side by side, one panel each. Build one
      # long (feature, cell, value) table over both: genes come from the store,
      # numeric columns from `cells` (keyed on the same global row index).
      items <- c(feats, metacols)
      validate(need(length(items) >= 1,
                    "Search for a gene, or pick a numeric column, to colour the embedding."),
               need(length(items) <= 12, "Select at most 12 features for the grid."))
      vl <- if (length(feats) == 1) {              # single gene keeps the cached query1 hot path
        q <- data$query1(assay(), feats)
        data.frame(feature = feats, cell = q$cell, value = q$value, stringsAsFactors = FALSE)
      } else if (length(feats) > 1) data$queryN(assay(), feats) else NULL
      if (length(metacols)) {
        key <- if (!is.null(cells$.gidx)) cells$.gidx else seq_len(nrow(cells))
        mvl <- do.call(rbind, lapply(metacols, function(mc)
          data.frame(feature = mc, cell = key,
                     value = suppressWarnings(as.numeric(cells[[mc]])),
                     stringsAsFactors = FALSE)))
        vl <- if (is.null(vl)) mvl else rbind(vl, mvl)
      }
      list(cells = cells, embedding = red_rv(), multi = TRUE,
           features = items, values = vl, n = nrow(cells))
    })
    cosmetic_r <- .scroll_cosmetic(reactive(
      list(theme = theme_r(), palette = input$palette, point_size = input$size,
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
      if (isTRUE(d$multi))
        return(view_feature_multi(
          d$cells, list(embedding = d$embedding, features = d$features), d$values, st))
      view_feature_plot(d$cells, list(embedding = d$embedding, feature = d$feature),
                        d$values, st)
    }
    plot_r   <- .scroll_lazy_plot(input, function() build(.scroll_use_raster(input$raster, data_r()$n)),
                                 style_r)
    export_r <- reactive(.scroll_apply_style(build(FALSE), style_r()))
    # cache key: active cells + reduction/assay/features/blend + raster + cosmetics
    # (+ onscreen, so an off-screen lazy render is never cached under a shown key)
    key_r <- reactive(list(cache_key_r(), isTRUE(input$onscreen %||% TRUE),
                           red_rv(), assay(), input$feature, input$metacol,
                           isTRUE(input$blend), input$blend_threshold,
                           input$blend_c1, input$blend_c2,
                           .scroll_use_raster(input$raster, nrow(cells_r())),
                           cosmetic_r(), style_r()))
    .scroll_render_cached(output, plot_r, key_r, cache)
    csv_r <- reactive({
      d <- data_r()
      # multi-gene grid: embedding coords + one expression column per gene
      if (isTRUE(d$multi)) {
        df <- .scroll_embedding_xy(d$cells, d$embedding)
        out <- df[, intersect(c("cell", paste0(d$embedding, c("_1", "_2"))), names(df)), drop = FALSE]
        for (g in d$features)
          out[[g]] <- .scroll_expr_vector(df, d$values[d$values$feature == g, c("cell", "value"), drop = FALSE])
        return(out)
      }
      out <- .scroll_featureplot_source(d$cells, d$embedding, d$feature, d$values)
      if (isTRUE(d$blend))
        out[[d$feature2]] <- .scroll_expr_vector(.scroll_embedding_xy(d$cells, d$embedding), d$values2)
      out
    })
    .scroll_plot_downloads(output, export_r, id, csv_r = csv_r)
  })
}

