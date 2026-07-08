#' Build a scroll project from a Seurat object
#'
#' Reads a processed Seurat object exactly once and emits a self-contained
#' project directory of lightweight on-disk artifacts: a small Parquet table of
#' cell metadata + embeddings (`cells.parquet`), a feature-partitioned Parquet
#' store of expression (`expr/<assay>/`), a `manifest.yaml`, and authored
#' scaffolding (`config.yaml`, `story.qmd`). The running app never touches the
#' object again.
#'
#' This is the Phase 1 walking skeleton: it wires a single assay end to end.
#'
#' @param object A processed Seurat object, or a path to an `.rds` file.
#' @param outdir Output project directory (created if missing).
#' @param assays Assays to export. `NULL` exports the default assay only.
#' @param embeddings Reductions to export. `NULL` exports all reductions.
#' @param meta_cols Metadata columns to export. `NULL` infers a sensible set.
#' @param quantize If `TRUE`, expression is quantized to `uint8` (256 levels),
#'   with a per-assay `max` recorded in the manifest for dequantization.
#' @param overwrite If `TRUE`, an existing `outdir` is removed first.
#'
#' @return `outdir`, invisibly.
#' @export
scroll_build <- function(object, outdir,
                         assays = NULL, embeddings = NULL, meta_cols = NULL,
                         quantize = TRUE, overwrite = FALSE) {
  object <- .scroll_load_object(object)

  if (is.null(assays)) assays <- SeuratObject::DefaultAssay(object)
  if (is.null(embeddings)) embeddings <- SeuratObject::Reductions(object)
  md <- object[[]]
  if (is.null(meta_cols)) meta_cols <- .scroll_infer_meta_cols(md)

  .scroll_check_inputs(object, assays, embeddings, meta_cols)

  if (dir.exists(outdir)) {
    if (overwrite) unlink(outdir, recursive = TRUE)
    else if (length(list.files(outdir)) > 0)
      stop("outdir '", outdir, "' exists and is not empty; use overwrite = TRUE.",
           call. = FALSE)
  }
  dir.create(file.path(outdir, "expr"), recursive = TRUE, showWarnings = FALSE)

  # --- cells.parquet: metadata + all embeddings, the only globally-loaded file
  cells <- .scroll_extract_cells(object, md, meta_cols, embeddings)
  arrow::write_parquet(cells, file.path(outdir, "cells.parquet"))

  # --- expr/<assay>/: long, feature-partitioned expression
  assay_info <- list()
  for (a in assays) {
    assay_info[[a]] <- .scroll_export_assay(object, a, outdir, quantize)
  }

  # --- manifest + authored scaffolding
  .scroll_write_manifest(outdir, object, assay_info, embeddings, md, meta_cols,
                         n_cells = nrow(cells), quantize = quantize)
  scroll_scaffold(outdir, assay_info, embeddings, meta_cols, md)

  message("scroll project built at: ", normalizePath(outdir))
  invisible(outdir)
}

# Load + normalize the object to a Seurat v5 object.
.scroll_load_object <- function(object) {
  if (is.character(object)) {
    if (!file.exists(object)) stop("No such file: ", object, call. = FALSE)
    object <- readRDS(object)
  }
  if (!inherits(object, "Seurat"))
    stop("`object` must be a Seurat object or a path to an .rds of one.",
         call. = FALSE)
  # Upgrade v3/v4 assays to v5 so GetAssayData(layer=) is well defined.
  default <- SeuratObject::DefaultAssay(object)
  if (!methods::is(object[[default]], "Assay5")) {
    object <- SeuratObject::UpdateSeuratObject(object)
  }
  object
}

# Pick categorical annotations + standard QC numerics.
.scroll_infer_meta_cols <- function(md) {
  keep <- character(0)
  for (col in colnames(md)) {
    v <- md[[col]]
    is_cat <- is.factor(v) || is.character(v) || is.logical(v)
    n_lvl <- length(unique(v))
    if (is_cat && n_lvl > 1 && n_lvl <= 200) keep <- c(keep, col)
  }
  qc <- grep("^(nCount|nFeature|percent)", colnames(md), value = TRUE)
  unique(c(keep, qc))
}

.scroll_check_inputs <- function(object, assays, embeddings, meta_cols) {
  bad_a <- setdiff(assays, SeuratObject::Assays(object))
  if (length(bad_a)) stop("Unknown assay(s): ", paste(bad_a, collapse = ", "),
                          call. = FALSE)
  bad_r <- setdiff(embeddings, SeuratObject::Reductions(object))
  if (length(bad_r)) stop("Unknown reduction(s): ", paste(bad_r, collapse = ", "),
                          call. = FALSE)
  bad_m <- setdiff(meta_cols, colnames(object[[]]))
  if (length(bad_m)) stop("Unknown metadata column(s): ",
                          paste(bad_m, collapse = ", "), call. = FALSE)
  if (length(embeddings) == 0)
    stop("Object has no dimensional reductions to export.", call. = FALSE)
}

# One row per cell: id + metadata + <reduction>_<dim> coordinate columns.
.scroll_extract_cells <- function(object, md, meta_cols, embeddings) {
  cell_ids <- rownames(md)
  tab <- data.frame(cell = cell_ids, stringsAsFactors = FALSE)
  for (col in meta_cols) {
    v <- md[[col]]
    tab[[col]] <- if (is.factor(v)) as.character(v) else v
  }
  for (r in embeddings) {
    emb <- SeuratObject::Embeddings(object, reduction = r)
    idx <- match(cell_ids, rownames(emb))
    for (d in seq_len(ncol(emb))) {
      tab[[sprintf("%s_%d", r, d)]] <- emb[idx, d]
    }
  }
  tab
}

# Export one assay's `data` layer as long, feature-partitioned Parquet.
# Returns list(features=, max=, n_features=) for the manifest.
.scroll_export_assay <- function(object, assay, outdir, quantize) {
  mat <- SeuratObject::GetAssayData(object, assay = assay, layer = "data")
  if (is.null(mat) || nrow(mat) == 0 || length(mat@x) == 0)
    stop("Assay '", assay, "' has an empty `data` layer; normalize the object ",
         "before building.", call. = FALSE)

  feats <- rownames(mat)
  cell_names <- colnames(mat)
  .scroll_warn_unsafe_features(feats)

  trip <- Matrix::summary(mat)             # i (feature), j (cell), x (value)
  max_a <- max(trip$x)
  value <- trip$x
  if (quantize) {
    value <- as.integer(pmin(pmax(round(value / max_a * 255), 0L), 255L))
  }
  df <- data.frame(
    feature = feats[trip$i],
    cell    = cell_names[trip$j],
    value   = value,
    stringsAsFactors = FALSE
  )
  # Sort by feature so each partition is written contiguously: arrow can close a
  # partition's file before opening the next, keeping open handles bounded even
  # with tens of thousands of features.
  df <- df[order(df$feature), , drop = FALSE]

  tbl <- arrow::as_arrow_table(df)
  if (quantize) {
    tbl <- tbl$cast(arrow::schema(
      feature = arrow::utf8(), cell = arrow::utf8(), value = arrow::uint8()))
  }
  arrow::write_dataset(
    tbl, file.path(outdir, "expr", assay),
    partitioning = "feature", format = "parquet",
    max_partitions = length(feats) + 1L
  )
  list(features = feats, max = max_a, n_features = length(feats))
}

# Prefer a 2-D visualization embedding over PCA as the plotting default.
.scroll_default_embedding <- function(embeddings) {
  prefer <- c("wnn.umap", "umap", "tsne", "wnnUMAP", "UMAP", "TSNE")
  hit <- intersect(prefer, embeddings)
  if (length(hit)) hit[[1]] else embeddings[[1]]
}

.scroll_warn_unsafe_features <- function(feats) {
  unsafe <- grep("[/\\\\]", feats, value = TRUE)
  if (length(unsafe))
    warning("Some feature names contain path separators and may break ",
            "partitioning: ", paste(utils::head(unsafe, 5), collapse = ", "),
            call. = FALSE)
}
