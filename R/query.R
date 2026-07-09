# Runtime query layer. A single feature lookup reads only that feature's
# partition off disk via duckdb predicate pushdown, so runtime RAM stays flat
# regardless of matrix size.

#' Open a duckdb connection for a scroll project
#'
#' @param dir A scroll project directory (must contain `expr/`).
#' @return A DBIConnection. Close it with [scroll_disconnect()].
#' @export
scroll_connect <- function(dir) {
  if (!dir.exists(file.path(dir, "expr")))
    stop("No expr/ store in '", dir, "'; is this a scroll project?", call. = FALSE)
  con <- DBI::dbConnect(duckdb::duckdb())
  attr(con, "scroll_dir") <- normalizePath(dir)
  con
}

#' Close a scroll duckdb connection
#' @param con A connection from [scroll_connect()].
#' @export
scroll_disconnect <- function(con) {
  DBI::dbDisconnect(con, shutdown = TRUE)
}

#' Query one feature's expression from the Parquet store
#'
#' Reads only the `feature=<feature>` partition of the assay subtree. Values are
#' returned as stored (quantized `uint8` when the project was built with
#' `quantize = TRUE`); use [scroll_dequantize()] to map back to normalized units.
#'
#' @param con A connection from [scroll_connect()].
#' @param assay Assay name (subtree under `expr/`).
#' @param feature Feature name to look up.
#' @return A data.frame with columns `cell` and `value`. Zero-valued cells are
#'   absent (the matrix is sparse); callers treat missing cells as 0.
#' @export
scroll_query_feature <- function(con, assay, feature) {
  dir <- attr(con, "scroll_dir")
  glob <- file.path(dir, "expr", assay, "**", "*.parquet")
  sql <- paste0(
    "SELECT cell, value FROM read_parquet(?, hive_partitioning = true) ",
    "WHERE feature = ?"
  )
  DBI::dbGetQuery(con, sql, params = list(glob, feature))
}

#' Query several features at once
#'
#' Reads only the named features' partitions (via `feature IN (...)` predicate
#' pushdown). Used by aggregate views (e.g. dotplot) that draw a bounded marker
#' panel; still touches only those partitions, never the whole store.
#'
#' @param con A connection from [scroll_connect()].
#' @param assay Assay name.
#' @param features Character vector of feature names.
#' @return A data.frame with columns `feature`, `cell`, `value` (zero-valued
#'   cells absent).
#' @export
scroll_query_features <- function(con, assay, features) {
  features <- unique(features)
  if (length(features) == 0)
    return(data.frame(feature = character(), cell = character(), value = numeric()))
  dir <- attr(con, "scroll_dir")
  glob <- file.path(dir, "expr", assay, "**", "*.parquet")
  placeholders <- paste(rep("?", length(features)), collapse = ", ")
  sql <- paste0(
    "SELECT feature, cell, value FROM read_parquet(?, hive_partitioning = true) ",
    "WHERE feature IN (", placeholders, ")"
  )
  DBI::dbGetQuery(con, sql, params = c(list(glob), as.list(features)))
}

#' Query all features for a set of cells
#'
#' Returns the long (feature, cell, value) records for the given cells across
#' every feature — used to reconstruct a contrast's expression matrix for live
#' DE. Unlike a feature lookup this scans all partitions (the cell filter does
#' not prune the feature-partitioned store), so it is the one query that touches
#' the whole assay; call it on demand, not per interaction.
#'
#' @param con A connection from [scroll_connect()].
#' @param assay Assay name.
#' @param cells Character vector of cell ids.
#' @return A data.frame with columns `feature`, `cell`, `value`.
#' @export
scroll_query_cells <- function(con, assay, cells) {
  if (length(cells) == 0)
    return(data.frame(feature = character(), cell = character(), value = numeric()))
  dir <- attr(con, "scroll_dir")
  glob <- file.path(dir, "expr", assay, "**", "*.parquet")
  duckdb::duckdb_register(con, "scroll_cellsel", data.frame(cell = cells, stringsAsFactors = FALSE))
  on.exit(duckdb::duckdb_unregister(con, "scroll_cellsel"), add = TRUE)
  DBI::dbGetQuery(con, paste0(
    "SELECT feature, cell, value FROM read_parquet(?, hive_partitioning = true) ",
    "WHERE cell IN (SELECT cell FROM scroll_cellsel)"), params = list(glob))
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
