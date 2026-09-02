# Incremental project update: add or replace reductions / metadata / subset views
# on an already-built project WITHOUT re-exporting the expression store. A
# reduction or metadata column lives only in `cells.parquet` (a few MB) + the
# manifest; the heavy `expr/` and `counts/` stores are untouched. `scroll_build()`
# is all-or-nothing (it wipes `outdir` and re-writes everything), so use this when
# you only want to append coordinate/metadata columns to an existing app.

#' Add or update reductions / metadata on a built project (no expr rebuild)
#'
#' Rewrites only `cells.parquet` and `manifest.yaml`, aligning new columns to the
#' existing cells **by barcode** (cells absent from `object` get `NA`, i.e. a
#' partial-coverage subset embedding). The feature-partitioned `expr/` and
#' `counts/` stores are left byte-identical, so this is seconds rather than a full
#' [scroll_build()] re-export.
#'
#' The arguments mirror [scroll_build()], so the *same* object you would rebuild
#' from can instead be passed here to append just the new pieces. Columns that
#' already exist are overwritten (e.g. to relabel a metadata column); new ones are
#' added. Subset views are registered exactly as in [scroll_build()] /
#' [scroll_add_subset()].
#'
#' @param dir A built scroll project directory (must contain `cells.parquet` and
#'   `manifest.yaml`).
#' @param object A Seurat object (or path to `.rds`) whose barcodes overlap the
#'   project. Only the requested reductions / metadata columns are read from it.
#' @param embeddings Character vector of reduction names in `object` to add as
#'   embedding coordinate columns (e.g. `c("umap_tcell")`). `NULL` adds none.
#' @param meta_cols Character vector of metadata columns in `object` to add /
#'   replace. `NULL` adds none.
#' @param subsets Optional named list of subset-view specs (as in [scroll_build()]);
#'   defaults to any spec attached by [scroll_add_subset()]. Their embeddings /
#'   meta must be present after this update (existing or newly added).
#' @param max_levels Cap on cached categorical `levels:` per column, as in
#'   [scroll_build()] (default 200). Columns above it record only `n_levels` and
#'   the app recomputes the level set at runtime.
#' @param verbose If `TRUE`, report what was written.
#' @return `dir`, invisibly.
#' @seealso [scroll_build()], [scroll_add_subset()]
#' @export
scroll_update <- function(dir, object, embeddings = NULL, meta_cols = NULL,
                          subsets = NULL, max_levels = .SCROLL_MAX_LEVELS,
                          verbose = interactive()) {
  .scroll_need_seurat()
  cells_path <- file.path(dir, "cells.parquet")
  if (!file.exists(cells_path))
    stop("'", dir, "' has no cells.parquet; is it a built scroll project?", call. = FALSE)
  man <- scroll_manifest(dir)
  object <- .scroll_load_object(object)
  md <- object[[]]

  embeddings <- as.character(embeddings %||% character(0))
  meta_cols  <- as.character(meta_cols  %||% character(0))
  if (is.null(subsets)) subsets <- SeuratObject::Misc(object, "scroll_subsets")

  bad_e <- setdiff(embeddings, SeuratObject::Reductions(object))
  if (length(bad_e))
    stop("object has no reduction(s): ", paste(bad_e, collapse = ", "), call. = FALSE)
  bad_m <- setdiff(meta_cols, colnames(md))
  if (length(bad_m))
    stop("object has no metadata column(s): ", paste(bad_m, collapse = ", "), call. = FALSE)
  if (!length(embeddings) && !length(meta_cols) && is.null(subsets)) {
    .scroll_step(verbose, "Nothing to update (no embeddings / meta_cols / subsets).")
    return(invisible(dir))
  }

  cells <- as.data.frame(arrow::read_parquet(cells_path))

  # extract the new columns from `object` and align to the existing cells by
  # barcode; a cell missing from `object` gets NA (partial coverage).
  new <- .scroll_extract_cells(object, md, meta_cols, embeddings)
  idx <- match(cells$cell, new$cell)
  add_cols <- setdiff(names(new), "cell")
  for (col in add_cols) cells[[col]] <- new[[col]][idx]

  arrow::write_parquet(cells, cells_path, compression = "zstd")

  # --- merge manifest -------------------------------------------------------
  for (r in embeddings) {
    man$embeddings[[r]] <- list(
      dims = ncol(SeuratObject::Embeddings(object, reduction = r)),
      n_covered = as.integer(sum(!is.na(cells[[sprintf("%s_1", r)]]))))
  }

  # normalize subset specs against the POST-update column universe so a view can
  # reference embeddings / meta added in this same call.
  subs_norm <- if (!is.null(subsets))
    .scroll_normalize_subsets(subsets,
                              union(names(man$embeddings), embeddings),
                              union(names(man$meta), meta_cols)) else NULL
  scope_map <- if (!is.null(subs_norm)) attr(subs_norm, "scope_map") else character(0)
  member <- list()
  if (!is.null(subs_norm))
    for (nm in names(subs_norm))
      member[[nm]] <- !is.na(cells[[sprintf("%s_1", subs_norm[[nm]]$primary_embedding)]])

  for (col in meta_cols) {
    sc <- if (col %in% names(scope_map)) scope_map[[col]] else NULL
    v <- if (!is.null(sc) && !is.null(member[[sc]])) cells[[col]][member[[sc]]] else cells[[col]]
    entry <- .scroll_meta_entry(v, max_levels)
    if (!is.null(sc)) entry$scope <- sc
    man$meta[[col]] <- entry
  }
  if (length(meta_cols)) .scroll_report_trimmed(man$meta[meta_cols])

  if (!is.null(subs_norm)) {
    for (nm in names(subs_norm)) {
      s <- subs_norm[[nm]]
      man$subsets[[nm]] <- list(
        label = s$label, embeddings = as.list(s$embeddings),
        primary_embedding = s$primary_embedding, meta = as.list(s$meta),
        n_cells = as.integer(sum(member[[nm]])))
    }
  }

  # unicode = TRUE to match .scroll_write_manifest(): otherwise re-writing the
  # manifest here would re-escape non-ASCII feature names to `<U+XXXX>` and break
  # their match against the store's feature column.
  yaml::write_yaml(man, file.path(dir, "manifest.yaml"), unicode = TRUE)
  .scroll_step(verbose, sprintf(
    "Updated %s: +%d embedding(s), +%d meta col(s)%s (expr store untouched).",
    dir, length(embeddings), length(meta_cols),
    if (!is.null(subs_norm)) sprintf(", %d subset view(s)", length(subs_norm)) else ""))
  invisible(dir)
}

#' Add (or replace) a metadata column in a built project -- no Seurat object needed
#'
#' Writes a new column into a built project's `cells.parquet` and registers it in
#' `manifest.yaml`, so it becomes selectable in the app on the next restart
#' (Color-by / group-by and DotPlot/Violin/Proportions for a **categorical** column;
#' FeaturePlot / Biaxial for a **numeric** one). Unlike [scroll_update()] this needs
#' only the on-disk project, not the original Seurat object -- handy for editing a
#' deployed app. The column must be something you can compute from what the project
#' already holds: other metadata columns, or expression queried with
#' [scroll_query_feature()]. The expression store is untouched.
#'
#' @param dir A built scroll project directory (with `cells.parquet` + `manifest.yaml`).
#' @param name Name of the column to add. Cannot be `"cell"` (the barcode column); an
#'   existing metadata column is replaced only when `overwrite = TRUE`.
#' @param values The column, one of: a vector of length `n_cells` in `cells.parquet`
#'   **row order**; a **named** vector keyed by cell barcode (unmatched cells become
#'   `NA`); a `data.frame` with a `cell` column and one value column; or a
#'   `function(cells)` returning one of those (it receives the loaded `cells`
#'   data.frame, so you can derive the column from existing columns). Character /
#'   factor / logical values register as categorical, numeric as numeric.
#' @param scope Optional name of a subset view (a key of the manifest `subsets`); when
#'   set, the column is scoped to that view -- its levels are computed over the view's
#'   members only, so it surfaces only inside that View.
#' @param overwrite Replace an existing column of the same name (default `FALSE`).
#' @param max_levels Cap on cached categorical levels written to the manifest; above it
#'   the level list is omitted and recomputed at runtime (default 200).
#' @return `invisible(dir)`.
#' @seealso [scroll_update()] to add columns/embeddings from a Seurat object.
#' @export
scroll_add_meta <- function(dir, name, values, scope = NULL, overwrite = FALSE,
                            max_levels = .SCROLL_MAX_LEVELS) {
  cells_path <- file.path(dir, "cells.parquet")
  if (!file.exists(cells_path) || !file.exists(file.path(dir, "manifest.yaml")))
    stop("'", dir, "' is not a built scroll project (need cells.parquet + manifest.yaml).",
         call. = FALSE)
  if (length(name) != 1L || is.na(name) || !nzchar(name))
    stop("`name` must be a single non-empty string.", call. = FALSE)
  if (identical(name, "cell"))
    stop("`name` cannot be 'cell' (the barcode column).", call. = FALSE)

  cells <- as.data.frame(arrow::read_parquet(cells_path), stringsAsFactors = FALSE)
  n <- nrow(cells)
  if (name %in% names(cells) && !isTRUE(overwrite))
    stop("Column '", name, "' already exists; pass overwrite = TRUE to replace it.",
         call. = FALSE)

  if (is.function(values)) values <- values(cells)
  if (is.data.frame(values)) {
    if (!("cell" %in% names(values)) || ncol(values) < 2L)
      stop("A data.frame `values` needs a `cell` column and one value column.", call. = FALSE)
    vcol <- setdiff(names(values), "cell")[1]
    nm <- values[[vcol]]; names(nm) <- as.character(values$cell); values <- nm
  }
  if (!is.null(names(values))) {
    v <- unname(values[match(cells$cell, names(values))])   # align by barcode; NA if absent
  } else {
    if (length(values) != n)
      stop("`values` has length ", length(values), " but the project has ", n,
           " cells; pass a length-", n, " vector (row order) or a barcode-named vector.",
           call. = FALSE)
    v <- values
  }
  if (is.factor(v)) v <- as.character(v)

  # scope: compute the manifest entry over the subset's members only (levels don't
  # leak into whole-dataset selectors); the stored column keeps the caller's values.
  man <- scroll_manifest(dir)
  if (!is.null(scope)) {
    if (is.null(man$subsets[[scope]]))
      stop("No subset view '", scope, "' in the manifest.", call. = FALSE)
    emb   <- as.character(man$subsets[[scope]]$embeddings)[1]
    coord <- sprintf("%s_1", emb)
    if (!coord %in% names(cells))
      stop("Subset '", scope, "' primary embedding '", emb, "' not in cells.parquet.",
           call. = FALSE)
    member <- !is.na(cells[[coord]])
    entry  <- .scroll_meta_entry(v[member], max_levels)
    entry$scope <- scope
  } else {
    entry <- .scroll_meta_entry(v, max_levels)
  }

  cells[[name]] <- v
  arrow::write_parquet(cells, cells_path, compression = "zstd")
  man$meta[[name]] <- entry
  # unicode = TRUE to match .scroll_write_manifest() (keeps non-ASCII feature names).
  yaml::write_yaml(man, file.path(dir, "manifest.yaml"), unicode = TRUE)
  invisible(dir)
}
