# Runtime query layer. A single feature lookup reads only that feature's
# partition off disk via arrow's dataset partition pruning, so runtime RAM stays
# flat regardless of matrix size. The `arrow::Dataset` for each assay is built
# once and cached in the connection handle, so repeated lookups never re-crawl
# the (potentially tens of thousands of) feature partitions.

#' Open a scroll project for querying
#'
#' Returns a lightweight handle over a built project's Parquet store. The handle
#' caches one `arrow::Dataset` per assay (built lazily on first use), so feature
#' lookups avoid re-scanning the partition directory. Release it with
#' [scroll_disconnect()].
#'
#' @param dir A scroll project directory (must contain `expr/`).
#' @return A `scroll_con` handle. Close it with [scroll_disconnect()].
#' @export
scroll_connect <- function(dir) {
  if (!dir.exists(file.path(dir, "expr")))
    stop("No expr/ store in '", dir, "'; is this a scroll project?", call. = FALSE)
  con <- new.env(parent = emptyenv())
  con$dir <- normalizePath(dir)
  con$datasets <- new.env(parent = emptyenv())   # assay -> arrow Dataset (cached)
  class(con) <- "scroll_con"
  con
}

#' Close a scroll project handle
#'
#' Drops the handle's cached datasets. There is no live database connection to
#' close (queries read Parquet directly via arrow), so this is a safe no-op that
#' may be called more than once.
#'
#' @param con A handle from [scroll_connect()].
#' @return `invisible(NULL)`.
#' @export
scroll_disconnect <- function(con) {
  if (inherits(con, "scroll_con") && is.environment(con$datasets))
    rm(list = ls(con$datasets, all.names = TRUE), envir = con$datasets)
  invisible(NULL)
}

# Validated `expr/<assay>` directory, rejecting any assay name that is not a
# plain path segment. `assay` comes from the manifest (not runtime user input),
# so this is defense-in-depth that keeps path construction injection-proof.
.scroll_assay_dir <- function(dir, assay) {
  if (length(assay) != 1L || is.na(assay) || !nzchar(assay) ||
      basename(assay) != assay || grepl("[/\\\\]", assay))
    stop("Invalid assay name: ", assay, call. = FALSE)
  file.path(dir, "expr", assay)
}

# The cached arrow Dataset for an assay's expr/ subtree (built once per handle).
.scroll_dataset <- function(con, assay) {
  path <- .scroll_assay_dir(con$dir, assay)   # validates assay; stops on invalid
  if (!is.null(con$datasets[[assay]])) return(con$datasets[[assay]])
  ds <- arrow::open_dataset(path)
  con$datasets[[assay]] <- ds
  ds
}

#' Query one feature's expression from the Parquet store
#'
#' Reads only the `feature=<feature>` partition of the assay subtree. Values are
#' returned as stored (quantized `uint8` when the project was built with
#' `quantize = TRUE`); use [scroll_dequantize()] to map back to normalized units.
#'
#' @param con A handle from [scroll_connect()].
#' @param assay Assay name (subtree under `expr/`).
#' @param feature Feature name to look up.
#' @return A data.frame with columns `cell` and `value`. Zero-valued cells are
#'   absent (the matrix is sparse); callers treat missing cells as 0.
#' @export
scroll_query_feature <- function(con, assay, feature) {
  ds <- .scroll_dataset(con, assay)
  out <- dplyr::collect(dplyr::select(
    dplyr::filter(ds, .data$feature == !!feature), "cell", "value"))
  as.data.frame(out, stringsAsFactors = FALSE)
}

#' Query several features at once
#'
#' Reads only the named features' partitions (partition pruning). Used by
#' aggregate views (e.g. dotplot) that draw a bounded marker panel; still touches
#' only those partitions, never the whole store.
#'
#' @param con A handle from [scroll_connect()].
#' @param assay Assay name.
#' @param features Character vector of feature names.
#' @return A data.frame with columns `feature`, `cell`, `value` (zero-valued
#'   cells absent).
#' @export
scroll_query_features <- function(con, assay, features) {
  features <- unique(features)
  if (length(features) == 0)
    return(data.frame(feature = character(), cell = character(), value = numeric()))
  ds <- .scroll_dataset(con, assay)
  out <- dplyr::collect(dplyr::select(
    dplyr::filter(ds, .data$feature %in% !!features), "feature", "cell", "value"))
  as.data.frame(out, stringsAsFactors = FALSE)
}

#' Query all features for a set of cells
#'
#' Returns the long (feature, cell, value) records for the given cells across
#' every feature — used to reconstruct a contrast's expression matrix for live
#' DE. Unlike a feature lookup this scans all partitions (the cell filter does
#' not prune the feature-partitioned store), so it is the one query that touches
#' the whole assay; call it on demand, not per interaction.
#'
#' @param con A handle from [scroll_connect()].
#' @param assay Assay name.
#' @param cells Character vector of cell ids.
#' @param dict If `TRUE`, return the `feature` column dictionary-encoded (an R
#'   factor) instead of a character vector. arrow builds the dictionary during the
#'   scan, so downstream `unique()`/`match()` over the tens of millions of repeated
#'   gene names become cheap integer ops — a large speedup when reconstructing a
#'   DE contrast's matrix. The factor's levels are the store's own feature names.
#' @return A data.frame with columns `feature`, `cell`, `value`.
#' @export
scroll_query_cells <- function(con, assay, cells, dict = FALSE) {
  if (length(cells) == 0)
    return(data.frame(feature = if (dict) factor() else character(),
                      cell = character(), value = numeric()))
  ds <- .scroll_dataset(con, assay)
  q  <- dplyr::filter(ds, .data$cell %in% !!cells)
  if (dict)
    q <- dplyr::mutate(q, feature = arrow::cast(
      .data$feature, arrow::dictionary(index_type = arrow::int32(),
                                       value_type = arrow::utf8())))
  out <- dplyr::collect(dplyr::select(q, "feature", "cell", "value"))
  as.data.frame(out, stringsAsFactors = FALSE)
}

# Path to an assay's raw-counts Parquet (counts/<assay>.parquet), rejecting an
# assay name that is not a plain path segment (defense-in-depth, mirrors
# .scroll_assay_dir).
.scroll_counts_path <- function(dir, assay) {
  if (length(assay) != 1L || is.na(assay) || !nzchar(assay) ||
      basename(assay) != assay || grepl("[/\\\\]", assay))
    stop("Invalid assay name: ", assay, call. = FALSE)
  file.path(dir, "counts", paste0(assay, ".parquet"))
}

#' Aggregate raw counts into pseudobulk samples
#'
#' Sums an assay's raw counts (the `counts/<assay>.parquet` store, written by
#' [scroll_build()] with `counts = TRUE`) over a cell-to-sample mapping. The join
#' and summation run in arrow's engine, so the full genes x cells matrix is never
#' materialized in R. Used by [scroll_pseudobulk_de()].
#'
#' @param con A handle from [scroll_connect()].
#' @param assay Assay name.
#' @param mapping A data.frame with columns `cell` and `psample` (a cell may map
#'   to several pseudobulk samples, e.g. overlapping pseudo-replicates).
#' @return A data.frame with columns `feature`, `psample`, `count` (summed).
#' @export
scroll_aggregate_counts <- function(con, assay, mapping) {
  stopifnot(all(c("cell", "psample") %in% names(mapping)))
  path <- .scroll_counts_path(con$dir, assay)
  if (!file.exists(path))
    stop("No counts store for assay '", assay, "'; rebuild with counts = TRUE.",
         call. = FALSE)
  if (!nrow(mapping))
    return(data.frame(feature = character(), psample = character(), count = numeric()))
  # keep `cell` in its native type (v2 int32 index / v1 barcode string) so the
  # join matches the counts store's `cell` column type.
  map <- data.frame(cell = mapping$cell,
                    psample = as.character(mapping$psample),
                    stringsAsFactors = FALSE)
  ds <- arrow::open_dataset(path)
  agg <- tryCatch({
    joined <- dplyr::inner_join(ds, arrow::arrow_table(map), by = "cell")
    grouped <- dplyr::summarise(
      dplyr::group_by(joined, .data$feature, .data$psample),
      count = sum(.data$value), .groups = "drop")
    dplyr::collect(grouped)
  }, error = function(e) {
    # Fallback if the arrow engine declines the many-to-many join: collect the
    # (bounded) selected cells, then join + sum in R.
    sel <- dplyr::collect(dplyr::filter(ds, .data$cell %in% !!map$cell))
    dplyr::summarise(
      dplyr::group_by(dplyr::inner_join(sel, map, by = "cell",
                                        relationship = "many-to-many"),
                      .data$feature, .data$psample),
      count = sum(.data$value), .groups = "drop")
  })
  as.data.frame(agg, stringsAsFactors = FALSE)
}

#' Map stored (possibly quantized) values back to normalized expression
#'
#' @param values Numeric/integer vector of stored values.
#' @param manifest A manifest list (from [scroll_manifest()]).
#' @param assay Assay the values came from.
#' @return Numeric vector in normalized units.
#' @export
scroll_dequantize <- function(values, manifest, assay) {
  if (!isTRUE(manifest$quantize)) return(as.numeric(values))
  max_a <- manifest$assays[[assay]]$max
  as.numeric(values) / 255 * max_a
}
