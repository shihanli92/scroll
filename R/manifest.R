# Manifest read/write. The manifest is the app's index of what the build
# produced: assays + feature namespaces, embeddings + dims, metadata columns +
# levels, and the defaults a fresh dataset needs to render itself.

.scroll_write_manifest <- function(outdir, object, assay_info, embeddings, md,
                                   meta_cols, n_cells, quantize, has_counts = FALSE,
                                   cells = NULL, subsets = NULL, vdj = NULL,
                                   images = NULL, spatial_embeddings = character(),
                                   atac = NULL, peaks_assay = character()) {
  assays <- lapply(names(assay_info), function(a) {
    info <- assay_info[[a]]
    l <- list(max = info$max, n_features = info$n_features,
              features = as.list(info$features))
    if (a %in% peaks_assay) l$kind <- "peaks"        # accessibility, not expression
    l
  })
  names(assays) <- names(assay_info)

  # Per-embedding coverage: how many cells carry non-NA coordinates. A full
  # embedding covers every cell; a reprocessed-subset embedding covers fewer.
  emb <- lapply(embeddings, function(r) {
    cov <- if (!is.null(cells)) sum(!is.na(cells[[sprintf("%s_1", r)]])) else n_cells
    e <- list(dims = ncol(SeuratObject::Embeddings(object, reduction = r)),
              n_covered = as.integer(cov))
    if (r %in% spatial_embeddings) e$kind <- "spatial"   # tissue-map coordinates
    e
  })
  names(emb) <- embeddings

  # meta col -> owning subset, and per-subset membership (non-NA primary coords),
  # so subset-scoped columns get their levels/ranges over subset cells only.
  scope_map <- if (!is.null(subsets)) attr(subsets, "scope_map") else character(0)
  member <- list()
  if (!is.null(subsets) && !is.null(cells))
    for (nm in names(subsets))
      member[[nm]] <- !is.na(cells[[sprintf("%s_1", subsets[[nm]]$primary_embedding)]])

  meta <- lapply(meta_cols, function(col) {
    v <- md[[col]]
    sc <- if (col %in% names(scope_map)) scope_map[[col]] else NULL
    if (!is.null(sc) && !is.null(member[[sc]])) v <- v[member[[sc]]]
    if (is.factor(v) || is.character(v) || is.logical(v)) {
      entry <- list(type = "categorical",
                    levels = as.list(sort(unique(as.character(v)))))
    } else {
      entry <- list(type = "numeric",
                    range = list(min = min(v, na.rm = TRUE), max = max(v, na.rm = TRUE)))
    }
    if (!is.null(sc)) entry$scope <- sc
    entry
  })
  names(meta) <- meta_cols

  manifest <- list(
    scroll_version = as.character(utils::packageVersion("scroll")),
    store_version = 2L,          # v2: expr `cell` is an int32 global row-index
    cell_index = TRUE,           # (v1 stores lack this -> string-barcode read path)
    n_cells = n_cells,
    quantize = quantize,
    has_counts = has_counts,
    default_assay = SeuratObject::DefaultAssay(object),
    default_embedding = .scroll_default_embedding(embeddings),
    assays = assays,
    embeddings = emb,
    meta = meta
  )
  if (!is.null(subsets)) {
    sm <- lapply(names(subsets), function(nm) {
      s <- subsets[[nm]]
      list(label = s$label, embeddings = as.list(s$embeddings),
           primary_embedding = s$primary_embedding, meta = as.list(s$meta),
           n_cells = if (!is.null(member[[nm]])) as.integer(sum(member[[nm]])) else NA_integer_)
    })
    names(sm) <- names(subsets)
    manifest$subsets <- sm
  }
  if (!is.null(vdj)) manifest$vdj <- vdj    # repertoire block -> gates the VDJ panels
  if (!is.null(images)) manifest$images <- images   # tissue image(s) -> gates the Spatial panel
  if (!is.null(atac)) manifest$atac <- atac         # peak table -> gates the Peaks panel
  yaml::write_yaml(manifest, file.path(outdir, "manifest.yaml"))
  invisible(manifest)
}

#' Read a scroll project manifest
#'
#' @param dir A scroll project directory.
#' @return The parsed manifest as a list.
#' @export
scroll_manifest <- function(dir) {
  path <- file.path(dir, "manifest.yaml")
  if (!file.exists(path)) stop("No manifest.yaml in '", dir, "'.", call. = FALSE)
  yaml::read_yaml(path)
}
