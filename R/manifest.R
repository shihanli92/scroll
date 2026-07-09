# Manifest read/write. The manifest is the app's index of what the build
# produced: assays + feature namespaces, embeddings + dims, metadata columns +
# levels, and the defaults a fresh dataset needs to render itself.

.scroll_write_manifest <- function(outdir, object, assay_info, embeddings, md,
                                   meta_cols, n_cells, quantize, has_counts = FALSE) {
  assays <- lapply(names(assay_info), function(a) {
    info <- assay_info[[a]]
    list(max = info$max, n_features = info$n_features,
         features = as.list(info$features))
  })
  names(assays) <- names(assay_info)

  emb <- lapply(embeddings, function(r) {
    list(dims = ncol(SeuratObject::Embeddings(object, reduction = r)))
  })
  names(emb) <- embeddings

  meta <- lapply(meta_cols, function(col) {
    v <- md[[col]]
    if (is.factor(v) || is.character(v) || is.logical(v)) {
      list(type = "categorical", levels = as.list(sort(unique(as.character(v)))))
    } else {
      list(type = "numeric",
           range = list(min = min(v, na.rm = TRUE), max = max(v, na.rm = TRUE)))
    }
  })
  names(meta) <- meta_cols

  manifest <- list(
    scroll_version = as.character(utils::packageVersion("scroll")),
    n_cells = n_cells,
    quantize = quantize,
    has_counts = has_counts,
    default_assay = SeuratObject::DefaultAssay(object),
    default_embedding = .scroll_default_embedding(embeddings),
    assays = assays,
    embeddings = emb,
    meta = meta
  )
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
