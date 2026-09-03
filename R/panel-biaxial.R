# --- Biaxial panel ------------------------------------------------------------

# Prefer hashtag / antibody CLR columns for the default axes, else the first few
# numerics — so an HTO or CITE-seq dataset opens on its biaxial signal plots.
.scroll_default_biaxial <- function(nums) {
  hit <- grep("^hto_|hashtag|adt_|_adt$|^ab_", nums, value = TRUE, ignore.case = TRUE)
  if (length(hit) >= 2) hit else utils::head(nums, 3)
}

biaxial_ui <- function(id, data) {
  ns <- NS(id)
  m <- data$manifest
  nums <- .scroll_num_cols(m); cats <- .scroll_cat_cols(m)
  if (length(nums) < 2)
    return(.scroll_empty_panel("Needs at least two numeric metadata columns (none found)."))
  if (!length(cats))
    return(.scroll_empty_panel("Needs a categorical column to colour by (none found)."))
  color_default <- if ("hto" %in% cats) "hto" else cats[[1]]
  assays <- .scroll_assays_of(m)
  # bare control id: conditionalPanel(ns = ns) prepends the module prefix itself
  when <- function(v) sprintf("input['%s'] == '%s'", "source", v)
  bslib::layout_columns(
    col_widths = c(3, 9), class = "scroll-panel",
    div(
      class = "scroll-controls",
      .scroll_group("Axes",
        radioButtons(ns("source"), NULL, c("Metadata", "Genes"), inline = TRUE),
        conditionalPanel(when("Metadata"), ns = ns,
          selectizeInput(ns("features"), "Numeric columns", choices = nums,
                         selected = .scroll_default_biaxial(nums), multiple = TRUE,
                         options = list(placeholder = "Pick 2+ numeric columns"))),
        conditionalPanel(when("Genes"), ns = ns,
          if (length(assays) > 1)
            selectInput(ns("gassay"), "Assay", assays, selected = m$default_assay),
          selectizeInput(ns("genes"), "Genes", choices = NULL, multiple = TRUE,
                         options = list(placeholder = "Pick 2+ genes", maxOptions = 50)))),
      .scroll_group("Colour",
        selectInput(ns("colorby"), "Colour by", stats::setNames(cats, cats),
                    selected = color_default)),
      .scroll_group("Appearance",
        selectInput(ns("palette"), "Palette", .scroll_cat_palettes()),
        uiOutput(ns("manual")),
        sliderInput(ns("size"), "Point size", 0.1, 3, 0.5, 0.1),
        sliderInput(ns("alpha"), "Opacity", 0.1, 1, 0.6, 0.05),
        bslib::input_switch(ns("legend"), "Legend", TRUE),
        bslib::input_switch(ns("raster"), "Rasterize (fast)", TRUE)),
      .scroll_group("Layout", .scroll_aspect_input(ns),
        numericInput(ns("ncol"), "Facet columns (blank = auto)", value = NA, min = 1, step = 1),
        numericInput(ns("nrow"), "Facet rows (blank = auto)", value = NA, min = 1, step = 1))
    ),
    .scroll_plot_area(ns, "460px", csv = TRUE)
  )
}

biaxial_server <- function(id, data, cells_r = reactive(data$cells),
                           view_r = reactive(NULL), theme_r = reactive(NULL),
                           cache_key_r = reactive(NULL), cache = NULL) {
  moduleServer(id, function(input, output, session) {
    m <- data$manifest
    .scroll_bind_view_cats(input, session, view_r, m, "colorby")
    # numeric axis choices are view-aware too (scoped numeric columns appear only in
    # their view). Multi-select: keep the current pair when still valid, else fall
    # back to the default axes for the new column set.
    observeEvent(view_r(), {
      nums <- .scroll_num_cols(m, view_r())
      cur <- isolate(input$features); keep <- cur[cur %in% nums]
      sel <- if (length(keep) >= 2) keep else .scroll_default_biaxial(nums)
      updateSelectizeInput(session, "features", choices = nums, selected = sel)
    }, ignoreNULL = FALSE)
    # per-level color pickers when the palette is "Manual" (levels of Colour by)
    lvl_r <- reactive({ req(input$colorby); .scroll_meta_levels(data, input$colorby) })
    output$manual <- renderUI(
      if (identical(input$palette, "Manual")) .scroll_manual_ui(session$ns, lvl_r()))
    manual_colors <- reactive(
      if (identical(input$palette, "Manual")) .scroll_manual_colors(input, lvl_r()))
    gassay <- reactive(input$gassay %||% m$default_assay)
    # repopulate the gene list for the active assay (Genes mode)
    observeEvent(gassay(), {
      feats <- .scroll_features_of(m, gassay())
      cur <- isolate(input$genes)
      updateSelectizeInput(session, "genes", choices = feats, server = TRUE,
                           selected = cur[cur %in% feats])
    })
    # DATA reactive: build the (potentially large) pairwise long df once; cosmetic
    # drags no longer re-expand it. Axes are either numeric metadata columns or
    # queried genes (attached to `cells` as columns, then paired like any numeric).
    data_r <- reactive({
      req(input$colorby)
      cells <- cells_r()
      if (identical(input$source, "Genes")) {
        genes <- input$genes
        validate(need(length(genes) >= 2, "Pick at least two genes."))
        for (g in genes) cells[[g]] <- .scroll_expr_vector(cells, data$query1(gassay(), g))
        feats <- genes
      } else {
        feats <- input$features
        validate(need(length(feats) >= 2, "Pick at least two numeric columns."))
      }
      params <- list(features = feats, color_by = input$colorby)
      # `cells` here carries the gene-augmented columns (Genes mode) so the CSV
      # source builder can read the plotted feature values per cell.
      list(params = params, df = .scroll_biaxial_df(cells, params), cells = cells)
    })
    cosmetic_r <- .scroll_cosmetic(reactive(
      list(theme = theme_r(), palette = input$palette, point_size = input$size,
           alpha = input$alpha, legend = isTRUE(input$legend), aspect = input$aspect,
           ncol = input$ncol, nrow = input$nrow, manual_colors = manual_colors())))
    build <- function(raster) {
      d <- data_r(); st <- cosmetic_r(); st$raster <- raster
      view_biaxial(NULL, d$params, st, df = d$df)
    }
    plot_r   <- .scroll_lazy_plot(input, function() build(.scroll_use_raster(input$raster, nrow(data_r()$df))))
    export_r <- reactive(build(FALSE))
    # cache key: active cells + axes/source/colour + raster + cosmetics (+ onscreen)
    key_r <- reactive(list(cache_key_r(), isTRUE(input$onscreen %||% TRUE),
                           input$source, input$features, input$genes, gassay(),
                           input$colorby,
                           .scroll_use_raster(input$raster, nrow(cells_r())),
                           cosmetic_r()))
    .scroll_render_cached(output, plot_r, key_r, cache)
    csv_r <- reactive({ d <- data_r(); .scroll_biaxial_source(d$cells, d$params$features, d$params$color_by) })
    .scroll_plot_downloads(output, export_r, id, csv_r = csv_r)
  })
}

