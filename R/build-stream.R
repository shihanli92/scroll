# Streaming/incremental build: assemble one scroll project from many sources
# without ever holding them all in RAM. Reads one source at a time, exports its
# expression with a running GLOBAL cell offset, moves the part-file into the shared
# store (one part per source; arrow concatenates them), and commits the source.
# Idempotent and append-safe: re-run as more sources arrive; sources already added
# are skipped. At the end of a run the project files (cells.parquet, stats/,
# manifest) are rebuilt from the commits and the per-source parts are merged into
# one feature-sorted file per assay.
#
# Commits live in outdir/.stream/: `<tag>.cells.parquet` (the source's cells) and
# `<tag>.rds` (its offset, cell count, settings, and per-assay features, gene sums,
# max and cell stats). The .rds is renamed into place last, so it is the commit
# point: a source either has one, and its expr / counts parts are complete, or it
# has not and whatever it left behind is removed on the next run. That keeps the
# global cell index of every stored value matched to its row of cells.parquet even
# after a run dies part-way. A run that finishes deletes .stream/: the project files
# then hold everything, and the next append adopts them as its base commit.

# Move a source's exported assay part-file(s) into the combined store, tagging them
# so they never collide across sources.
.scroll_stream_move <- function(src_assay_dir, dst_assay_dir, tag) {
  if (!dir.exists(dst_assay_dir)) dir.create(dst_assay_dir, recursive = TRUE, showWarnings = FALSE)
  for (pp in list.files(src_assay_dir, pattern = "\\.parquet$"))
    file.rename(file.path(src_assay_dir, pp),
                file.path(dst_assay_dir, paste0(tag, "-", pp)))
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
#' ones, so a build can grow as data arrives. Each source is committed on its own
#' (under the hidden `outdir/.stream/` folder), so a run that stops part-way --
#' an error, a full disk, a killed session -- keeps the sources it finished; run
#' it again to carry on. `.stream/` is removed once a run finishes. Values are stored unquantized (`float32`); a global quantization max is
#' unknown while streaming, so `quantize = TRUE` is not supported here.
#'
#' @param outdir Project directory to create or append to.
#' @param sources A list/vector of sources (file paths, ids, …) passed one at a
#'   time to `reader`.
#' @param reader `function(source) -> Seurat object`, already processed (a
#'   log-normalized `data` layer, the reductions, and the metadata columns).
#' @param assays,embeddings Names to export; `NULL` takes the default assay / all
#'   reductions of the **first** source (then held constant across sources).
#' @param meta_cols Metadata columns to expose (required; declared explicitly
#'   rather than inferred). The declared assays, embeddings, and meta columns must
#'   be present in **every** source — a source missing any is an error naming it
#'   (extra columns a source happens to carry are simply ignored).
#' @param quantize Must be `FALSE` (streaming stores float32; see Details).
#' @param counts If `TRUE`, also export each source's raw `counts` layer (one part
#'   per source under `counts/<assay>/`), enabling the Pseudobulk DE panel — same
#'   as [scroll_build()]'s `counts = TRUE`. An append run must match the setting
#'   the project was first built with.
#' @param ranks Store within-cell ranks for the rank-based signature scores (see
#'   [scroll_build()]). `NULL` (default) keeps an existing project's setting, and is
#'   `FALSE` for a new one; an append run must match the original setting.
#' @param id_of `function(source) -> character` id used for append-tracking and
#'   part-file tags (default `as.character`).
#' @param overwrite If `TRUE`, wipe `outdir` first (a fresh build).
#' @param verbose Print per-source progress.
#' @return Invisibly, `outdir`.
#' @seealso [scroll_build()]
#' @export
scroll_build_stream <- function(outdir, sources, reader, assays = NULL,
                                embeddings = NULL, meta_cols = NULL, quantize = FALSE,
                                counts = FALSE, id_of = as.character, overwrite = FALSE,
                                ranks = NULL, verbose = interactive()) {
  .scroll_need_seurat()
  .scroll_check_zstd()
  if (isTRUE(quantize))
    stop("scroll_build_stream() writes an unquantized float32 store (a global ",
         "quantization max is unknown while streaming); use quantize = FALSE.", call. = FALSE)
  if (!is.function(reader)) stop("`reader` must be function(source) -> Seurat object.", call. = FALSE)
  if (is.null(meta_cols)) stop("`meta_cols` must be given for a streaming build.", call. = FALSE)
  sources <- as.list(sources)
  ids  <- vapply(sources, id_of, character(1))
  tags <- .scroll_safe_tag(ids)
  if (anyDuplicated(ids))
    stop("scroll_build_stream: duplicate source id '", ids[anyDuplicated(ids)], "'.", call. = FALSE)
  if (anyDuplicated(tags)) {
    k <- anyDuplicated(tags)
    stop(sprintf("scroll_build_stream: source ids '%s' and '%s' map to the same file tag; ",
                 ids[match(tags[k], tags)], ids[k]),
         "give them ids that differ in letters or digits (`id_of`).", call. = FALSE)
  }

  kept_style <- if (dir.exists(outdir) && isTRUE(overwrite)) .scroll_keep_style(outdir)
  if (dir.exists(outdir) && isTRUE(overwrite)) unlink(outdir, recursive = TRUE)
  dir.create(file.path(outdir, "expr"), recursive = TRUE, showWarnings = FALSE)
  .scroll_restore_style(outdir, kept_style)       # saved plot styles survive a rebuild
  tmp_root <- file.path(outdir, ".stream_tmp")
  srcfile  <- file.path(outdir, ".stream_sources.txt")
  sdir <- .scroll_stream_dir(outdir)

  # recover from an interrupted run: finish / roll back a compaction, adopt a
  # project streamed by an older scroll, then drop anything no commit owns
  for (adir in list.dirs(file.path(outdir, "expr"), recursive = FALSE)) .scroll_compact_recover(adir)
  .scroll_stream_adopt(outdir)
  recs <- .scroll_stream_records(outdir)
  .scroll_stream_cleanup(outdir, recs)
  unlink(tmp_root, recursive = TRUE)

  done <- unlist(lapply(recs, `[[`, "ids"))
  todo <- which(!(ids %in% done))
  finalized <- length(recs) && file.exists(srcfile) &&
    setequal(readLines(srcfile, warn = FALSE), done)
  if (!length(todo) && finalized) {
    unlink(sdir, recursive = TRUE)                 # the project files hold it all
    message("scroll_build_stream: nothing new to add.")
    return(invisible(outdir))
  }

  set <- if (length(recs)) recs[[1]]$settings
  # A counts store must cover every cell or pseudobulk would silently drop the
  # sources without one, so an append run must match the original setting.
  if (!is.null(set) && !identical(isTRUE(set$counts), isTRUE(counts)))
    stop(sprintf(paste0("scroll_build_stream: this project was built with counts = %s; ",
                        "append with the same setting, or rebuild with overwrite = TRUE."),
                 isTRUE(set$counts)), call. = FALSE)
  if (is.null(ranks)) ranks <- isTRUE(set$ranks)
  if (!is.null(set) && !identical(isTRUE(ranks), isTRUE(set$ranks)))
    stop(sprintf(paste0("scroll_build_stream: this project was built with ranks = %s; ",
                        "append with the same setting, or rebuild with overwrite = TRUE."),
                 isTRUE(set$ranks)), call. = FALSE)
  if (!is.null(set) && !setequal(meta_cols, set$meta_cols))
    stop("scroll_build_stream: this project was built with meta_cols = c(",
         paste0('"', set$meta_cols, '"', collapse = ", "), "); append with the same ",
         "columns, or rebuild with overwrite = TRUE.", call. = FALSE)
  lock <- if (!is.null(set)) set[c("assays", "embeddings", "meta_cols")]
  offset <- sum(vapply(recs, function(r) as.numeric(r$n), numeric(1)))

  if (length(todo)) dir.create(sdir, showWarnings = FALSE)
  for (i in todo) {
    id <- ids[[i]]; tag <- tags[[i]]
    .scroll_step(verbose, sprintf("[stream] %s ...", id))
    seu <- reader(sources[[i]])
    seu <- tryCatch(SeuratObject::UpdateSeuratObject(seu), error = function(e) seu)
    md  <- as.data.frame(seu@meta.data)

    if (is.null(lock)) {
      a  <- if (is.null(assays)) SeuratObject::DefaultAssay(seu) else assays
      em <- if (is.null(embeddings)) names(seu@reductions) else embeddings
      lock <- list(assays = a, embeddings = em, meta_cols = meta_cols)
    }
    # Every source must carry the declared assays / embeddings / meta columns:
    # the cells frames are rbind-ed at the end, so a source missing a column would
    # otherwise fail deep in the merge. Fail here with the offending source named.
    miss_a <- setdiff(lock$assays, SeuratObject::Assays(seu))
    if (length(miss_a))
      stop(sprintf("[stream] source '%s' is missing declared assay(s): %s",
                   id, paste(miss_a, collapse = ", ")), call. = FALSE)
    miss_e <- setdiff(lock$embeddings, names(seu@reductions))
    if (length(miss_e))
      stop(sprintf("[stream] source '%s' is missing declared reduction(s): %s",
                   id, paste(miss_e, collapse = ", ")), call. = FALSE)
    miss_m <- setdiff(lock$meta_cols, colnames(md))
    if (length(miss_m))
      stop(sprintf("[stream] source '%s' is missing declared meta column(s): %s",
                   id, paste(miss_m, collapse = ", ")), call. = FALSE)
    cframe <- .scroll_extract_cells(seu, md, lock$meta_cols, lock$embeddings)

    info <- list()
    for (assay in lock$assays) {
      ai <- .scroll_export_assay(seu, assay, tmp_root, quantize = FALSE, offset = offset,
                                 cell_ids = cframe$cell, ranks = ranks)
      info[[assay]] <- list(features = ai$features, sums = unname(ai$stats$genes$sum),
                            max = ai$max, cstats = ai$stats$cells)
      .scroll_stream_move(file.path(tmp_root, "expr", assay),
                          file.path(outdir, "expr", assay), tag = tag)
      if (isTRUE(counts)) {
        cdir <- file.path(outdir, "counts", assay)
        dir.create(cdir, recursive = TRUE, showWarnings = FALSE)
        .scroll_export_counts(seu, assay, outdir, cell_ids = cframe$cell, offset = offset,
                              dest = file.path(cdir, paste0(tag, "-part.parquet")))
      }
    }
    unlink(file.path(tmp_root, "expr"), recursive = TRUE)

    # commit: the cells first, the record (renamed into place) last
    arrow::write_parquet(cframe, file.path(sdir, paste0(tag, ".cells.parquet")),
                         compression = "zstd")
    rec <- list(ids = id, tags = tag, key = tag, offset = offset, n = nrow(cframe),
                settings = c(lock, list(counts = isTRUE(counts), ranks = isTRUE(ranks))),
                emb_dims = vapply(lock$embeddings, function(r)
                  ncol(SeuratObject::Embeddings(seu, reduction = r)), integer(1)),
                assays = info)
    .scroll_save_rds_atomic(rec, file.path(sdir, paste0(tag, ".rds")))
    recs <- c(recs, list(rec))
    offset <- offset + nrow(cframe)
    rm(seu, cframe); gc(FALSE)
  }
  unlink(tmp_root, recursive = TRUE)

  cells <- .scroll_stream_finalize(outdir, recs)
  writeLines(unlist(lapply(recs, `[[`, "ids")), srcfile)   # what the project files now hold
  # the project files now carry everything an append needs (the next run adopts
  # them as its base), so the per-source copies can go
  unlink(sdir, recursive = TRUE)
  message("scroll (streamed) project built at: ", normalizePath(outdir),
          "  (", format(nrow(cells), big.mark = ","), " cells)")
  invisible(outdir)
}

.scroll_stream_dir <- function(outdir) file.path(outdir, ".stream")

# saveRDS to a temp name, then rename: the file is either absent or complete.
.scroll_save_rds_atomic <- function(x, path) {
  tmp <- paste0(path, ".tmp")
  saveRDS(x, tmp)
  if (!file.rename(tmp, path)) stop("Could not write '", path, "'.", call. = FALSE)
  invisible(path)
}

# The committed sources, in cell-index order. Their cell counts must tile the
# global index from 1 with no gap or overlap.
.scroll_stream_records <- function(outdir) {
  files <- list.files(.scroll_stream_dir(outdir), pattern = "\\.rds$", full.names = TRUE)
  recs <- lapply(files, readRDS)
  if (!length(recs)) return(list())
  recs <- recs[order(vapply(recs, function(r) as.numeric(r$offset), numeric(1)))]
  n <- vapply(recs, function(r) as.numeric(r$n), numeric(1))
  off <- vapply(recs, function(r) as.numeric(r$offset), numeric(1))
  if (!isTRUE(all.equal(off, c(0, cumsum(n)[-length(n)]))))
    stop("scroll_build_stream: the saved sources in '", .scroll_stream_dir(outdir),
         "' do not line up (cell offsets ", paste(off, collapse = ", "), "). Rebuild with ",
         "overwrite = TRUE.", call. = FALSE)
  recs
}

# Remove whatever no committed source owns: expr / counts parts of a source that
# never committed, and its half-written files in .stream/.
.scroll_stream_cleanup <- function(outdir, recs) {
  keep <- unlist(lapply(recs, `[[`, "tags"))
  owned <- function(files, rx) {
    tag <- sub(rx, "\\1", basename(files))
    !grepl(rx, basename(files)) | tag %in% keep
  }
  for (adir in list.dirs(file.path(outdir, "expr"), recursive = FALSE)) {
    f <- list.files(adir, pattern = "\\.parquet$", full.names = TRUE)
    unlink(f[!owned(f, "^([A-Za-z0-9_]+)-part-[0-9]+\\.parquet$")])
  }
  for (cdir in list.dirs(file.path(outdir, "counts"), recursive = FALSE)) {
    f <- list.files(cdir, pattern = "\\.parquet$", full.names = TRUE)
    unlink(f[!owned(f, "^([A-Za-z0-9_]+)-part\\.parquet$")])
  }
  sdir <- .scroll_stream_dir(outdir)
  keys <- vapply(recs, function(r) r$key, character(1))
  f <- list.files(sdir, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  cf <- f[grepl("\\.cells\\.parquet$", f)]
  unlink(c(f[grepl("\\.tmp$", f)],
           cf[!sub("\\.cells\\.parquet$", "", basename(cf)) %in% keys]))
  invisible(NULL)
}

# A finished project (or one streamed by scroll < 0.3.16) has no .stream/ commits:
# turn what it has (cells.parquet, stats, manifest, .stream_sources.txt) into one
# base commit.
# Sources whose expr parts were never merged belong to a run that died before
# writing the project files, so they are left out (and re-added on this run).
.scroll_stream_adopt <- function(outdir) {
  sdir <- .scroll_stream_dir(outdir)
  if (dir.exists(sdir)) return(invisible(NULL))
  cells_path <- file.path(outdir, "cells.parquet")
  mpath <- file.path(outdir, "manifest.yaml")
  if (!file.exists(cells_path) || !file.exists(mpath)) return(invisible(NULL))
  man <- scroll_manifest(outdir)
  cells <- as.data.frame(arrow::read_parquet(cells_path, mmap = FALSE))
  if (!identical(as.numeric(man$n_cells), as.numeric(nrow(cells))))
    stop("scroll_build_stream: this project's cells.parquet and manifest disagree on the ",
         "cell count (an earlier run stopped part-way). Rebuild with overwrite = TRUE.",
         call. = FALSE)
  srcfile <- file.path(outdir, ".stream_sources.txt")
  ids <- if (file.exists(srcfile)) unique(readLines(srcfile, warn = FALSE)) else character(0)
  tags <- .scroll_safe_tag(ids)
  loose <- unique(sub("^([A-Za-z0-9_]+)-part-[0-9]+\\.parquet$", "\\1",
                      grep("^[A-Za-z0-9_]+-part-[0-9]+\\.parquet$",
                           basename(list.files(file.path(outdir, "expr"), recursive = TRUE)),
                           value = TRUE)))
  keep <- !tags %in% loose
  info <- lapply(names(man$assays), function(a) {
    ma <- man$assays[[a]]
    sd <- file.path(outdir, "stats", a)
    st <- isTRUE(ma$stats) && file.exists(file.path(sd, "genes.parquet"))
    g <- if (st) as.data.frame(arrow::read_parquet(file.path(sd, "genes.parquet"), mmap = FALSE))
    list(features = if (st) g$feature else unlist(ma$features),
         sums = if (st) g$sum,
         max = ma$max %||% 1,
         cstats = if (st) as.data.frame(arrow::read_parquet(file.path(sd, "cells.parquet"), mmap = FALSE)))
  })
  names(info) <- names(man$assays)
  dir.create(sdir, showWarnings = FALSE)
  arrow::write_parquet(cells, file.path(sdir, "_base.cells.parquet"), compression = "zstd")
  rec <- list(ids = ids[keep], tags = tags[keep], key = "_base", offset = 0, n = nrow(cells),
              settings = list(assays = names(man$assays), embeddings = names(man$embeddings),
                              meta_cols = names(man$meta), counts = isTRUE(man$has_counts),
                              ranks = isTRUE(man$assays[[1]]$ranks)),
              emb_dims = vapply(man$embeddings, function(e) as.integer(e$dims), integer(1)),
              assays = info)
  .scroll_save_rds_atomic(rec, file.path(sdir, "_base.rds"))
  invisible(NULL)
}

# Rebuild the project files from the commits: cells.parquet, stats/, the merged
# expr store, the manifest (features = the union over sources) and the scaffold.
.scroll_stream_finalize <- function(outdir, recs) {
  sdir <- .scroll_stream_dir(outdir)
  cells <- do.call(rbind, lapply(recs, function(r) as.data.frame(
    arrow::read_parquet(file.path(sdir, paste0(r$key, ".cells.parquet")), mmap = FALSE))))
  rownames(cells) <- NULL
  .scroll_rewrite_parquet(cells, file.path(outdir, "cells.parquet"))
  set <- recs[[1]]$settings
  assay_info <- list()
  for (a in set$assays) {
    parts <- lapply(recs, function(r) r$assays[[a]])
    feats <- unique(unlist(lapply(parts, `[[`, "features")))
    info <- list(features = feats, n_features = length(feats), ranks = isTRUE(set$ranks),
                 max = max(vapply(parts, function(p) as.numeric(p$max %||% 0), numeric(1))))
    # gene / cell statistics need every source: a project first built without
    # them (an older scroll) stays without them
    if (all(vapply(parts, function(p) !is.null(p$sums), logical(1)))) {
      sm <- stats::setNames(numeric(length(feats)), feats)
      for (p in parts) sm[p$features] <- sm[p$features] + p$sums
      st <- list(genes = data.frame(feature = feats, sum = unname(sm),
                                    mean = unname(sm) / nrow(cells), stringsAsFactors = FALSE),
                 cells = do.call(rbind, lapply(parts, `[[`, "cstats")))
      .scroll_write_stats(outdir, a, st)
      info$stats <- st
    }
    assay_info[[a]] <- info
  }
  .scroll_compact_store(outdir, quantize = FALSE)
  .scroll_write_manifest(outdir, NULL, assay_info, set$embeddings, cells, set$meta_cols,
                         n_cells = nrow(cells), quantize = FALSE,
                         has_counts = isTRUE(set$counts), cells = cells, subsets = NULL,
                         default_assay = set$assays[[1]],
                         emb_dims = recs[[length(recs)]]$emb_dims)
  scroll_scaffold_app(outdir)
  cells
}
