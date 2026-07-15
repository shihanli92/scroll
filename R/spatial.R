# Spatial (10x Visium / imaging-based) support. Spot/cell coordinates are just a
# 2-D embedding, so DimPlot and FeaturePlot already work on spatial data once the
# tissue coordinates are exported as a reduction. What spatial adds on top is the
# **tissue image**: at build time `scroll_build(..., spatial = spatial_spec(...))`
# extracts `GetTissueCoordinates()` into a `spatial` embedding (in image-pixel
# space, y-flipped for ggplot) and bakes the H&E array to a small raster asset;
# the manifest records an `images` block and marks the embedding `kind: spatial`.
# A built-in Spatial panel then draws points over the tissue image with a fixed
# aspect ratio. It is gated on a spatial *embedding* existing (not the image), so it
# also serves imaging platforms (Xenium/CosMx) that have centroids but no H&E raster;
# RNA-only / non-spatial projects never show it.

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

