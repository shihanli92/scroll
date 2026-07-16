# Streaming/incremental build: assemble one scroll project from many sources
# without ever holding them all in RAM. Reads one source at a time, exports its
# expression with a running GLOBAL cell offset, moves the part-files into the
# shared store (arrow concatenates parts within a bucket= partition), accumulates
# the small cells frame, and writes the manifest once at the end. Idempotent and
# append-safe: re-run as more sources arrive; sources already added are skipped.

# Move a source's exported feature partitions into the combined store, tagging the
# part-files so they never collide across sources.
.scroll_stream_move <- function(src_assay_dir, dst_assay_dir, tag) {
  for (fd in list.files(src_assay_dir)) {                 # bucket=<char> partition dirs
    dst <- file.path(dst_assay_dir, fd)
    if (!dir.exists(dst)) dir.create(dst, recursive = TRUE, showWarnings = FALSE)
    for (pp in list.files(file.path(src_assay_dir, fd)))
      file.rename(file.path(src_assay_dir, fd, pp),
                  file.path(dst, paste0(tag, "-", pp)))
  }
}

.scroll_safe_tag <- function(id) gsub("[^A-Za-z0-9]+", "_", as.character(id))

#' Build a scroll project by streaming many sources (large-dataset build)
#'
#' Produces the same flat-RAM store as [scroll_build()], but assembled from many
#' sources one at a time, so peak memory stays ~per-source instead of loading the
#' whole (possibly multi-million-cell) dataset at once. Each source is read by
#' `reader`, its expression appended to a shared v2 store under a running global
#' cell index, and the manifest written once at the end. The store is **compacted**
#' (one part-file per feature) on completion.
#'
#' Idempotent and **append-safe**: re-running with more `sources` adds only the new
#' ones (tracked in `outdir/.stream_sources.txt`), so a build can grow as data
#' arrives. Values are stored unquantized (`float32`); a global quantization max is
#' unknown while streaming, so `quantize = TRUE` is not supported here.
#'
#' @param outdir Project directory to create or append to.
#' @param sources A list/vector of sources (file paths, ids, …) passed one at a
#'   time to `reader`.
#' @param reader `function(source) -> Seurat object`, already processed (a
#'   log-normalized `data` layer, the reductions, and the metadata columns).
#' @param assays,embeddings Names to export; `NULL` takes the default assay / all
#'   reductions of the **first** source (then held constant across sources).
#' @param meta_cols Metadata columns to expose (required; the reader's objects may
#'   be heterogeneous, so the columns are declared explicitly).
#' @param quantize Must be `FALSE` (streaming stores float32; see Details).
#' @param id_of `function(source) -> character` id used for append-tracking and
#'   part-file tags (default `as.character`).
#' @param overwrite If `TRUE`, wipe `outdir` first (a fresh build).
#' @param verbose Print per-source progress.
#' @return Invisibly, `outdir`.
#' @seealso [scroll_build()]
#' @export
scroll_build_stream <- function(outdir, sources, reader, assays = NULL,
                                embeddings = NULL, meta_cols = NULL, quantize = FALSE,
                                id_of = as.character, overwrite = FALSE,
                                verbose = interactive()) {
  .scroll_need_seurat()
  if (isTRUE(quantize))
    stop("scroll_build_stream() writes an unquantized float32 store (a global ",
         "quantization max is unknown while streaming); use quantize = FALSE.", call. = FALSE)
  if (!is.function(reader)) stop("`reader` must be function(source) -> Seurat object.", call. = FALSE)
  if (is.null(meta_cols)) stop("`meta_cols` must be given for a streaming build.", call. = FALSE)
  sources <- as.list(sources)

  if (dir.exists(outdir) && isTRUE(overwrite)) unlink(outdir, recursive = TRUE)
  dir.create(file.path(outdir, "expr"), recursive = TRUE, showWarnings = FALSE)
  tmp_root <- file.path(outdir, ".stream_tmp")
  srcfile  <- file.path(outdir, ".stream_sources.txt")

  ids  <- vapply(sources, id_of, character(1))
  done <- if (file.exists(srcfile)) readLines(srcfile) else character(0)
  todo <- which(!(ids %in% done))
  if (!length(todo)) { message("scroll_build_stream: nothing new to add."); return(invisible(outdir)) }

  cells_path <- file.path(outdir, "cells.parquet")
  old_cells  <- if (file.exists(cells_path)) as.data.frame(arrow::read_parquet(cells_path)) else NULL
  offset <- if (is.null(old_cells)) 0L else nrow(old_cells)
  gmax   <- if (file.exists(file.path(outdir, "manifest.yaml")))
              tryCatch(lapply(scroll_manifest(outdir)$assays, function(a) a$max),
                       error = function(e) list()) else list()

  lock <- NULL; new_cells <- list(); last_seu <- NULL
  for (i in todo) {
    id <- ids[[i]]
    .scroll_step(verbose, sprintf("[stream] %s ...", id))
    seu <- reader(sources[[i]])
    seu <- tryCatch(SeuratObject::UpdateSeuratObject(seu), error = function(e) seu)
    md  <- as.data.frame(seu@meta.data)

    if (is.null(lock)) {
      a  <- if (is.null(assays)) SeuratObject::DefaultAssay(seu) else assays
      em <- if (is.null(embeddings)) names(seu@reductions) else embeddings
      lock <- list(assays = a, embeddings = em, meta_cols = meta_cols)
    }
    mc <- intersect(lock$meta_cols, colnames(md))
    cframe <- .scroll_extract_cells(seu, md, mc, lock$embeddings)

    for (assay in lock$assays) {
      ai <- .scroll_export_assay(seu, assay, tmp_root, quantize = FALSE, offset = offset)
      .scroll_stream_move(file.path(tmp_root, "expr", assay),
                          file.path(outdir, "expr", assay), tag = .scroll_safe_tag(id))
      gmax[[assay]] <- max(gmax[[assay]] %||% 0, ai$max)
    }
    unlink(file.path(tmp_root, "expr"), recursive = TRUE)

    new_cells[[id]] <- cframe
    offset <- offset + nrow(cframe)
    last_seu <- seu
    cat(id, "\n", sep = "", file = srcfile, append = TRUE)
    rm(seu); gc(FALSE)
  }
  unlink(tmp_root, recursive = TRUE)

  # merge cells (old + new, in the same order the global offsets were assigned)
  cells <- do.call(rbind, c(list(old_cells), unname(new_cells)))
  arrow::write_parquet(cells, cells_path, compression = "zstd")
  .scroll_compact_store(outdir, quantize = FALSE)

  SeuratObject::DefaultAssay(last_seu) <- lock$assays[1]
  assay_info <- lapply(lock$assays, function(assay) {
    feats <- rownames(SeuratObject::GetAssayData(last_seu, assay = assay, layer = "data"))
    list(features = feats, max = gmax[[assay]] %||% 1, n_features = length(feats))
  })
  names(assay_info) <- lock$assays
  .scroll_write_manifest(outdir, last_seu, assay_info, lock$embeddings, cells,
                         lock$meta_cols, n_cells = nrow(cells), quantize = FALSE,
                         has_counts = FALSE, cells = cells, subsets = NULL)
  scroll_scaffold_app(outdir)
  message("scroll (streamed) project built at: ", normalizePath(outdir),
          "  (", format(nrow(cells), big.mark = ","), " cells)")
  invisible(outdir)
}
