# Spatial (10x Visium / imaging-based) support. Spot/cell coordinates are just a
# 2-D embedding, so DimPlot and FeaturePlot already work on spatial data once the
# tissue coordinates are exported as a reduction. What spatial adds on top is the
# **tissue image**: at build time `scroll_build(..., spatial = spatial_spec(...))`
# extracts `GetTissueCoordinates()` into a `spatial` embedding (in image-pixel
# space, y-flipped for ggplot) and bakes the H&E array to a small raster asset;
# the manifest records an `images` block and marks the embedding `kind: spatial`.
# A built-in Spatial panel then draws points over the tissue image with a fixed
# aspect ratio. RNA-only / non-spatial projects are unaffected (the panel is gated
# on the manifest images block, like the VDJ panels).

#' Describe a dataset's spatial image + coordinates for `scroll_build()`
#'
#' Names the tissue image (Seurat FOV / `@images` slot) whose spot coordinates and
#' H&E array should be exported. Pass the result as
#' `scroll_build(..., spatial = spatial_spec(...))` to add a `spatial` embedding
#' (from `GetTissueCoordinates()`), bake the tissue image, and enable the Spatial
#' panel. Defaults auto-detect the first image, so a standard Visium object needs
#' only `spatial_spec()`.
#'
#' @param image Name of the image / FOV to export (an entry of
#'   `SeuratObject::Images(object)`). `NULL` (default) uses the first one.
#' @param name Embedding name to give the exported coordinates (default
#'   `"spatial"`). Appears as a reduction in DimPlot/FeaturePlot and is the
#'   Spatial panel's coordinate system.
#' @param image_asset Whether to bake the tissue image as an underlay asset
#'   (default `TRUE`). `FALSE` exports coordinates only (the Spatial panel then
#'   draws points on a blank background).
#' @return A `scroll_spatial_spec` list.
#' @seealso [scroll_build()]
#' @export
spatial_spec <- function(image = NULL, name = "spatial", image_asset = TRUE) {
  structure(list(image = image, name = name, image_asset = isTRUE(image_asset)),
            class = "scroll_spatial_spec")
}

# ---- build-time: coordinates -> reduction + baked image ---------------------

# Pull the (cell x 2) tissue coordinates for `image` as a plain matrix with
# columns (col, row) in full-resolution pixel space, robust across Visium v1/v2.
.scroll_tissue_xy <- function(object, image) {
  co <- SeuratObject::GetTissueCoordinates(object[[image]])
  if (is.null(co) || !nrow(co)) stop("Image '", image, "' has no tissue coordinates.", call. = FALSE)
  rn <- if (!is.null(co$cell)) as.character(co$cell) else rownames(co)
  # v2: columns x (col), y (row); v1: imagecol / imagerow
  col <- if (!is.null(co$x)) co$x else co$imagecol
  row <- if (!is.null(co$y)) co$y else co$imagerow
  if (is.null(col) || is.null(row))
    stop("Could not find x/y (or imagecol/imagerow) tissue coordinates for '", image, "'.",
         call. = FALSE)
  m <- cbind(col = as.numeric(col), row = as.numeric(row))
  rownames(m) <- rn
  m
}

# The lowres scale factor mapping full-res pixels -> stored image pixels (1 if none).
# Reads the `scale.factors` slot directly (a Visium image object carries it) rather
# than depending on a `ScaleFactors` generic that varies across Seurat versions.
.scroll_image_scalefactor <- function(img) {
  lr <- NULL
  if (methods::.hasSlot(img, "scale.factors")) {
    sf <- methods::slot(img, "scale.factors")
    lr <- if (is.list(sf)) sf$lowres else sf["lowres"]
  }
  if (is.null(lr) || !is.finite(lr) || lr <= 0) 1 else as.numeric(lr)
}

# Prepare a spatial spec: returns the coordinate DimReduc (image-pixel space,
# y-flipped so tissue is upright in ggplot), the embedding name, and — when an
# image is baked — the raster + its [0,W] x [0,H] extent. Coordinates are aligned
# to the baked image so points overlay it directly.
.scroll_prepare_spatial <- function(object, spec) {
  images <- SeuratObject::Images(object)
  if (!length(images)) stop("spatial: object has no images (SeuratObject::Images() is empty).",
                            call. = FALSE)
  image <- spec$image %||% images[[1]]
  if (!image %in% images)
    stop("spatial: image '", image, "' not found; available: ",
         paste(images, collapse = ", "), call. = FALSE)

  xy <- .scroll_tissue_xy(object, image)
  img_obj <- object[[image]]
  arr <- tryCatch(SeuratObject::GetImage(img_obj, mode = "raw"), error = function(e) NULL)
  if (is.null(arr) && methods::.hasSlot(img_obj, "image")) arr <- img_obj@image

  raster <- NULL; W <- NULL; H <- NULL
  if (isTRUE(spec$image_asset) && !is.null(arr) && length(dim(arr)) == 3) {
    lr <- .scroll_image_scalefactor(img_obj)
    H <- dim(arr)[1]; W <- dim(arr)[2]
    px <- xy[, "col"] * lr                 # image column (x, left->right)
    py <- xy[, "row"] * lr                 # image row (top->bottom)
    coord <- cbind(px, H - py)             # flip row so ggplot y is upright
    raster <- grDevices::as.raster(arr)    # H x W character matrix (row 1 = top)
  } else {
    # coordinates-only: use raw pixels, flip y so the tissue isn't upside down
    coord <- cbind(xy[, "col"], max(xy[, "row"]) - xy[, "row"])
  }
  colnames(coord) <- c(paste0(spec$name, "_1"), paste0(spec$name, "_2"))
  rownames(coord) <- rownames(xy)

  reduction <- suppressWarnings(SeuratObject::CreateDimReducObject(
    embeddings = coord, key = paste0(gsub("[^A-Za-z0-9]", "", spec$name), "_"),
    assay = SeuratObject::DefaultAssay(object)))

  list(name = spec$name, image = image, reduction = reduction,
       raster = raster, width = W, height = H)
}

# Write the baked tissue raster (+ extent) and return the manifest `images` block.
# NULL when no image was baked (coordinates-only).
.scroll_write_spatial <- function(outdir, prep) {
  if (is.null(prep$raster)) return(NULL)
  dir.create(file.path(outdir, "spatial"), showWarnings = FALSE, recursive = TRUE)
  f <- file.path("spatial", paste0(prep$name, ".rds"))
  saveRDS(list(raster = prep$raster, width = prep$width, height = prep$height),
          file.path(outdir, f))
  block <- list(list(embedding = prep$name, image = prep$image, file = f,
                     width = prep$width, height = prep$height))
  names(block) <- prep$name
  block
}

# ---- runtime: read the baked image ------------------------------------------

# Load the baked raster + extent for a spatial embedding, or NULL.
.scroll_spatial_image <- function(data, embedding) {
  blk <- data$manifest$images[[embedding]]
  if (is.null(blk) || is.null(blk$file)) return(NULL)
  p <- file.path(data$dir, blk$file)
  if (!file.exists(p)) return(NULL)
  tryCatch(readRDS(p), error = function(e) NULL)
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
.scroll_spatial_plot <- function(df, img, values, is_num, size, lab, zoom = NULL) {
  df$.col <- values
  p <- ggplot2::ggplot(df, ggplot2::aes(.data$.x, .data$.y))
  if (!is.null(img))
    p <- p + ggplot2::annotation_raster(img$raster, xmin = 0, xmax = img$width,
                                        ymin = 0, ymax = img$height, interpolate = TRUE)
  p <- p + .scroll_point_layer(ggplot2::aes(color = .data$.col), size = size,
                               raster = isTRUE(nrow(df) > .scroll_raster_threshold))
  p <- p + if (is_num) .scroll_continuous_scale("viridis", name = lab)
           else ggplot2::scale_color_manual(values = .scroll_discrete_colors(df$.col), name = lab)
  xl <- if (!is.null(zoom)) zoom$x else if (!is.null(img)) c(0, img$width) else range(df$.x)
  yl <- if (!is.null(zoom)) zoom$y else if (!is.null(img)) c(0, img$height) else range(df$.y)
  p + ggplot2::coord_fixed(xlim = xl, ylim = yl, expand = FALSE) +
    ggplot2::labs(title = lab) + ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                   legend.position = "right")
}

# The spatial coordinates + colour vector for the current controls, or a stop().
.scroll_spatial_frame <- function(cells, input, data) {
  emb <- .scroll_spatial_embeddings(data$manifest)
  if (!length(emb)) stop("No spatial embedding in this project.")
  df <- .scroll_embedding_xy(cells, emb[[1]])
  if (!nrow(df)) stop("No cells with spatial coordinates.")
  if (identical(input$mode %||% "Gene", "Metadata")) {
    col <- input$meta
    if (is.null(col) || !col %in% names(df)) stop("Pick a metadata column.")
    list(df = df, values = df[[col]], is_num = is.numeric(df[[col]]), lab = col, emb = emb[[1]])
  } else {
    feat <- input$gene
    if (is.null(feat) || !nzchar(feat)) stop("Type a gene to colour by.")
    v <- .scroll_expr_vector(df, data$query1(data$manifest$default_assay, feat))
    list(df = df, values = v, is_num = TRUE, lab = feat, emb = emb[[1]])
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
  emb <- emb[[1]]
  has_img <- !is.null(m$images[[emb]])
  meta_cols <- c(.scroll_cat_cols(m), .scroll_num_cols(m))
  ctl <- list(
    radioButtons(ns("mode"), "Colour by", c("Gene", "Metadata"), inline = TRUE),
    conditionalPanel(sprintf("input['%s'] == 'Gene'", ns("mode")), ns = ns,
      selectizeInput(ns("gene"), "Gene", choices = NULL,
                     options = list(placeholder = "type a gene", maxOptions = 50))),
    conditionalPanel(sprintf("input['%s'] == 'Metadata'", ns("mode")), ns = ns,
      selectInput(ns("meta"), "Metadata", meta_cols)),
    sliderInput(ns("size"), "Spot size", 0.2, 4, 1.4, 0.2),
    if (has_img) checkboxInput(ns("image"), "Show tissue image", TRUE),
    tags$p(class = "scroll-desc", "Drag to zoom \u00b7 double-click to reset."))
  bslib::layout_columns(
    col_widths = c(3, 9), class = "scroll-panel",
    div(class = "scroll-controls", do.call(.scroll_group, c(list("Controls"), ctl))),
    div(class = "scroll-plot",
        div(class = "scroll-plot-bar",
            .scroll_dl_button(ns("png"), "PNG"), .scroll_dl_button(ns("pdf"), "PDF")),
        .scroll_spin(plotOutput(ns("plot"), height = "560px",
            brush = shiny::brushOpts(ns("brush"), resetOnNew = TRUE),
            dblclick = ns("dblclick")))))
}

spatial_server <- function(id, data, cells_r = shiny::reactive(data$cells),
                           view_r = shiny::reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    m <- data$manifest
    emb <- .scroll_spatial_embeddings(m)[[1]]
    updateSelectizeInput(session, "gene", server = TRUE,
                         choices = .scroll_features_of(m, m$default_assay))
    zoom <- reactiveVal(NULL)
    # brush selection -> zoom to that window; double-click -> reset to full extent
    observeEvent(input$brush, {
      b <- input$brush
      if (!is.null(b)) zoom(list(x = c(b$xmin, b$xmax), y = c(b$ymin, b$ymax)))
    })
    observeEvent(input$dblclick, zoom(NULL))
    # reset zoom when the colour target changes (a fresh view, not a stale window)
    observeEvent(list(input$mode, input$gene, input$meta), zoom(NULL))

    frame_r <- reactive(.scroll_spatial_frame(cells_r(), input, data))
    plot_r <- reactive({
      fr <- frame_r()
      show_img <- is.null(input$image) || isTRUE(input$image)
      img <- if (show_img) .scroll_spatial_image(data, fr$emb) else NULL
      .scroll_spatial_plot(fr$df, img, fr$values, fr$is_num,
                           input$size %||% 1.4, fr$lab, zoom())
    })
    output$plot <- renderPlot(plot_r())
    .scroll_plot_downloads(output, plot_r, id)
  })
}

.scroll_spatial_panels <- function() {
  gate <- function(m) length(.scroll_spatial_embeddings(m)) > 0
  list(list(id = "spatial", label = "Spatial", title = "Tissue map",
            desc = "Cells/spots in tissue space, coloured by a gene or metadata. Drag to zoom.",
            when = gate, ui = spatial_ui, server = spatial_server))
}
