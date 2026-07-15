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
  # bound query helpers that dequantize to normalized units, memoized by one LRU
  # so cosmetic re-renders, the vector-export path, and re-selecting a recent gene
  # never re-hit the store.
  cache <- .scroll_lru(256L)
  d$query1 <- .scroll_cached_query(cache, manifest, "1",
                                   function(a, f) scroll_query_feature(con, a, f))
  d$queryN <- .scroll_cached_query(cache, manifest, "N",
                                   function(a, f) scroll_query_features(con, a, f))
  d
}

# One dequantizing, LRU-cached query closure. `prefix` ("1"/"N") disambiguates the
# single- vs multi-feature cache keys; `run(assay, features)` does the store read.
# For a single feature `sort(unique(.))` is a no-op, so both closures share one key
# scheme.
.scroll_cached_query <- function(cache, manifest, prefix, run) {
  function(assay, features) {
    key <- paste0(prefix, "|", assay, "|", paste(sort(unique(features)), collapse = ","))
    hit <- cache$get(key)
    if (is.null(hit)) {
      hit <- run(assay, features)
      if (nrow(hit)) hit$value <- scroll_dequantize(hit$value, manifest, assay)
      cache$set(key, hit)
    }
    hit
  }
}

# A bounded LRU over query results keyed by a string. Bounded so a long session
# browsing many genes cannot grow memory without limit; each entry is one
# feature's small sparse long table, so a generous bound is cheap.
.scroll_lru <- function(max = 256L) {
  store <- new.env(parent = emptyenv())
  order <- character(0)
  list(
    get = function(key) {
      if (!exists(key, store, inherits = FALSE)) return(NULL)
      order <<- c(setdiff(order, key), key)                # accessing bumps to MRU
      base::get(key, store)
    },
    set = function(key, val) {
      order <<- c(setdiff(order, key), key)                # most-recently-used last
      assign(key, val, store)
      while (length(order) > max) {                        # evict least-recently-used
        rm(list = order[[1]], envir = store); order <<- order[-1L]
      }
    }
  )
}
