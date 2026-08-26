# The Spatial (tissue-map) panel + its runtime accessors. Build-time coordinate/
# image extraction lives in R/spatial.R; these read the baked image asset and render
# the panel (brush-zoom, image underlay). A dedicated ui/server pair, not the builder.

# ---- runtime: read the baked image ------------------------------------------

# Load the baked raster + extent for a spatial embedding, or NULL.
.scroll_spatial_image <- function(data, embedding) {
  blk <- data$manifest$images[[embedding]]
  if (is.null(blk) || is.null(blk$file)) return(NULL)
  .scroll_read_asset(file.path(data$dir, blk$file), "rds")
}

# The spatial embeddings declared in a manifest (marked kind: spatial).
.scroll_spatial_embeddings <- function(m) {
  if (is.null(m$embeddings)) return(character())
  names(Filter(function(e) identical(e$kind, "spatial"), m$embeddings))
}

# ---- the Spatial panel ------------------------------------------------------

# Build the tissue-map ggplot: points at the spatial coordinates coloured by a
# gene/metadata vector, optionally over the baked tissue raster, at a fixed aspect
# ratio. `zoom` (a list(x=, y=) of axis limits, or NULL) restricts the view — the
# brush/zoom hook feeds it. Pure, so it renders identically on screen and on export.
.scroll_spatial_plot <- function(df, img, values, is_num, size, lab, zoom = NULL,
                                 cpalette = "viridis", dpalette = "Tableau 10") {
  df$.col <- values
  p <- ggplot2::ggplot(df, ggplot2::aes(.data$.x, .data$.y))
  if (!is.null(img))
    p <- p + ggplot2::annotation_raster(img$raster, xmin = 0, xmax = img$width,
                                        ymin = 0, ymax = img$height, interpolate = TRUE)
  p <- p + .scroll_point_layer(ggplot2::aes(color = .data$.col), size = size,
                               raster = isTRUE(nrow(df) > .scroll_raster_threshold))
  p <- p + if (is_num) .scroll_continuous_scale(cpalette %||% "viridis", name = lab)
           else ggplot2::scale_color_manual(
             values = .scroll_discrete_colors(df$.col, dpalette %||% "Tableau 10"), name = lab)
  xl <- if (!is.null(zoom)) zoom$x else if (!is.null(img)) c(0, img$width) else range(df$.x)
  yl <- if (!is.null(zoom)) zoom$y else if (!is.null(img)) c(0, img$height) else range(df$.y)
  p + ggplot2::coord_fixed(xlim = xl, ylim = yl, expand = FALSE) +
    ggplot2::labs(title = lab) + ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                   legend.position = "right")
}

# The spatial coordinates + colour vector for the current controls, or a stop().
# Honours the selected embedding (a project may bake several FOVs / slides).
.scroll_spatial_frame <- function(cells, input, data) {
  emb_all <- .scroll_spatial_embeddings(data$manifest)
  if (!length(emb_all)) stop("No spatial embedding in this project.")
  emb <- if (!is.null(input$embedding) && input$embedding %in% emb_all)
    input$embedding else emb_all[[1]]
  df <- .scroll_embedding_xy(cells, emb)
  if (!nrow(df)) stop("No cells with spatial coordinates.")
  if (identical(input$mode %||% "Gene", "Metadata")) {
    col <- input$meta
    if (is.null(col) || !col %in% names(df)) stop("Pick a metadata column.")
    list(df = df, values = df[[col]], is_num = is.numeric(df[[col]]), lab = col, emb = emb)
  } else {
    feat <- input$gene
    if (is.null(feat) || !nzchar(feat)) stop("Type a gene to colour by.")
    v <- .scroll_expr_vector(df, data$query1(data$manifest$default_assay, feat))
    list(df = df, values = v, is_num = TRUE, lab = feat, emb = emb)
  }
}

# A dedicated (non-declarative) panel so it can own brush-to-zoom + dblclick-reset
# on the plot. Gated on a spatial embedding existing — so it also surfaces for
# imaging platforms (Xenium/CosMx) that carry cell centroids but no H&E raster.
spatial_ui <- function(id, data) {
  ns <- shiny::NS(id)
  m <- data$manifest
  emb <- .scroll_spatial_embeddings(m)
  if (!length(emb)) return(.scroll_empty_panel("No spatial embedding in this project."))
  has_img <- length(m$images) > 0                          # any FOV/slide carries an image
  meta_cols <- c(.scroll_cat_cols(m), .scroll_num_cols(m))
  ctl <- list(
    # embedding selector only when several tissue maps were baked (multi-FOV/slide)
    if (length(emb) > 1) selectInput(ns("embedding"), "Tissue map", stats::setNames(emb, emb)),
    radioButtons(ns("mode"), "Colour by", c("Gene", "Metadata"), inline = TRUE),
    conditionalPanel("input['mode'] == 'Gene'", ns = ns,
      selectizeInput(ns("gene"), "Gene", choices = NULL,
                     options = list(placeholder = "type a gene", maxOptions = 50))),
    conditionalPanel("input['mode'] == 'Metadata'", ns = ns,
      selectInput(ns("meta"), "Metadata", meta_cols)),
    selectInput(ns("cpalette"), "Palette", .scroll_continuous_palettes, selected = "viridis"),
    conditionalPanel("input['mode'] == 'Metadata'", ns = ns,
      selectInput(ns("dpalette"), "Palette (categorical)", names(.scroll_discrete_palettes))),
    sliderInput(ns("size"), "Spot size", 0.2, 4, 1.4, 0.2),
    if (has_img) checkboxInput(ns("image"), "Show tissue image", TRUE),
    tags$p(class = "scroll-desc", "Drag to zoom \u00b7 double-click to reset."))
  bslib::layout_columns(
    col_widths = c(3, 9), class = "scroll-panel",
    div(class = "scroll-controls", do.call(.scroll_group, c(list("Controls"), ctl))),
    div(class = "scroll-plot",
        div(class = "scroll-plot-bar",
            .scroll_size_slider(),
            .scroll_dl_button(ns("png"), "PNG"), .scroll_dl_button(ns("pdf"), "PDF"),
            .scroll_dl_button(ns("csv"), "CSV")),
        .scroll_spin(plotOutput(ns("plot"), height = "560px",
            brush = shiny::brushOpts(ns("brush"), resetOnNew = TRUE),
            dblclick = ns("dblclick")))))
}

spatial_server <- function(id, data, cells_r = shiny::reactive(data$cells),
                           view_r = shiny::reactive(NULL), theme_r = shiny::reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    m <- data$manifest
    updateSelectizeInput(session, "gene", server = TRUE,
                         choices = .scroll_features_of(m, m$default_assay))
    # scope the metadata choices to the active subset view (scoped module scores etc.)
    observeEvent(view_r(), {
      cols <- c(.scroll_cat_cols(m, view_r()), .scroll_num_cols(m, view_r()))
      cur <- input$meta
      updateSelectInput(session, "meta", choices = cols,
                        selected = if (!is.null(cur) && cur %in% cols) cur else cols[1])
    }, ignoreNULL = FALSE)
    zoom <- reactiveVal(NULL)
    # brush selection -> zoom to that window; double-click -> reset to full extent
    observeEvent(input$brush, {
      b <- input$brush
      if (!is.null(b)) zoom(list(x = c(b$xmin, b$xmax), y = c(b$ymin, b$ymax)))
    })
    observeEvent(input$dblclick, zoom(NULL))
    # reset zoom when the colour target or tissue map changes (a fresh view)
    observeEvent(list(input$mode, input$gene, input$meta, input$embedding), zoom(NULL))

    frame_r <- reactive(.scroll_spatial_frame(cells_r(), input, data))
    plot_r <- reactive({
      fr <- frame_r()
      show_img <- is.null(input$image) || isTRUE(input$image)
      img <- if (show_img) .scroll_spatial_image(data, fr$emb) else NULL
      .scroll_spatial_plot(fr$df, img, fr$values, fr$is_num,
                           input$size %||% 1.4, fr$lab, zoom(),
                           input$cpalette %||% "viridis", input$dpalette %||% "Tableau 10") +
        .scroll_ggtheme(theme_r())
    })
    output$plot <- renderPlot(plot_r())
    .scroll_plot_downloads(output, plot_r, id)
    output$csv <- .scroll_csv_handler(reactive({
      fr <- frame_r()
      d <- data.frame(x = fr$df$.x, y = fr$df$.y); d[[fr$lab]] <- fr$values
      if (!is.null(fr$df$cell)) d <- cbind(cell = fr$df$cell, d)
      d
    }), paste0("scroll_", id, ".csv"))
  })
}

.scroll_spatial_panels <- function() {
  gate <- function(m) length(.scroll_spatial_embeddings(m)) > 0
  list(list(id = "spatial", label = "Spatial", title = "Tissue map",
            desc = "Cells/spots in tissue space, coloured by a gene or metadata. Drag to zoom.",
            when = gate, ui = spatial_ui, server = spatial_server))
}
