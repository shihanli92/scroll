# Runtime data handle: open a built project once per app process and expose the
# cells table + manifest + config plus LRU-cached, dequantizing feature queries.
# The running app never loads the Seurat object (flat RAM).

# Read a baked side-car asset (repertoire table, tissue raster, peak table) from a
# project, returning NULL if it is absent or unreadable. Shared by the modality
# accessors (.scroll_vdj_read / .scroll_spatial_image / .scroll_atac_read).
.scroll_read_asset <- function(path, format = c("parquet", "rds")) {
  format <- match.arg(format)
  if (!file.exists(path)) return(NULL)
  tryCatch(
    if (format == "rds") readRDS(path)
    else as.data.frame(arrow::read_parquet(path)),
    error = function(e) NULL)
}

# Load the artifacts once (shared across sessions of one app process).
.scroll_load <- function(dir) {
  con <- scroll_connect(dir)
  manifest <- scroll_manifest(dir)
  # projects built before 0.3.16 could name an assay that was not exported
  if (length(manifest$assays))
    manifest$default_assay <- .scroll_pick_default_assay(manifest$default_assay,
                                                         names(manifest$assays))
  cells <- as.data.frame(arrow::read_parquet(file.path(dir, "cells.parquet")))
  # `.gidx` is the canonical global row-index a v2 expr store joins on; it rides
  # along on subset/filtered cells so joins stay correct under the app-bar filter.
  cells$.gidx <- seq_len(nrow(cells))
  d <- list(
    dir = normalizePath(dir),
    cells = cells,
    manifest = manifest,
    cell_index = isTRUE(manifest$cell_index),   # TRUE for v2 stores (int cell key)
    config = scroll_config(dir),
    con = con
  )
  # bound query helpers that dequantize to normalized units, memoized per feature
  # by one LRU (bounded in bytes) so cosmetic re-renders, the vector-export path,
  # re-selecting a recent gene, and a gene list edited one gene at a time never
  # re-read what they already have.
  cache <- .scroll_lru(max_bytes = .scroll_query_cache_bytes())
  d$query1 <- function(assay, feature) {
    .scroll_cached_features(cache, manifest, con, assay, feature)[[1]]
  }
  d$queryN <- function(assay, features) {
    parts <- .scroll_cached_features(cache, manifest, con, assay, features)
    n <- vapply(parts, nrow, integer(1))
    data.frame(feature = rep(names(parts), n),
               cell = unlist(lapply(parts, `[[`, "cell"), use.names = FALSE) %||% integer(0),
               value = unlist(lapply(parts, `[[`, "value"), use.names = FALSE) %||% numeric(0),
               stringsAsFactors = FALSE)
  }
  d
}

# The query cache's size, in MB: option `scroll.query_cache_mb` (default 512).
.scroll_query_cache_bytes <- function()
  as.numeric(getOption("scroll.query_cache_mb", 512)) * 2^20

# Each requested feature's (cell, value) table, dequantized, from the cache; the
# features it lacks are read in one store query and cached one by one. Returns a
# named list in `features` order (unique), each a data.frame (0 rows if unexpressed).
.scroll_cached_features <- function(cache, manifest, con, assay, features) {
  features <- unique(as.character(features))
  keys <- paste0(assay, "|", features)
  out <- stats::setNames(lapply(keys, cache$get), features)
  miss <- features[vapply(out, is.null, logical(1))]
  if (length(miss)) {
    long <- if (length(miss) == 1L) {                  # `==` prunes row groups best
              one <- scroll_query_feature(con, assay, miss)
              data.frame(feature = rep(miss, nrow(one)), cell = one$cell, value = one$value,
                         stringsAsFactors = FALSE)
            } else scroll_query_features(con, assay, miss)
    if (nrow(long)) long$value <- scroll_dequantize(long$value, manifest, assay)
    by <- split(long[, c("cell", "value"), drop = FALSE],
                factor(as.character(long$feature), levels = miss))
    for (f in miss) {
      v <- by[[f]]; rownames(v) <- NULL
      cache$set(paste0(assay, "|", f), v, size = 16 * nrow(v) + 1000)
      out[[f]] <- v
    }
  }
  out
}

# A bounded LRU over query results keyed by a string. Bounded by the entries'
# estimated bytes (and their number), so a long session -- or many sessions of one
# app process, which share it -- cannot grow memory without limit.
.scroll_lru <- function(max = 4096L, max_bytes = Inf) {
  store <- new.env(parent = emptyenv())
  order <- character(0); sizes <- numeric(0); total <- 0
  list(
    get = function(key) {
      if (!exists(key, store, inherits = FALSE)) return(NULL)
      order <<- c(setdiff(order, key), key)                # accessing bumps to MRU
      base::get(key, store)
    },
    set = function(key, val, size = 0) {
      if (!is.na(sizes[key])) total <<- total - sizes[[key]]
      order <<- c(setdiff(order, key), key)                # most-recently-used last
      assign(key, val, store); sizes[key] <<- size; total <<- total + size
      # evict least-recently-used (always keeping the entry just set)
      while (length(order) > 1L && (length(order) > max || total > max_bytes)) {
        k <- order[[1]]
        rm(list = k, envir = store); order <<- order[-1L]
        total <<- total - sizes[[k]]; sizes <<- sizes[names(sizes) != k]
      }
    },
    size = function() total
  )
}
