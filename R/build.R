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
#'   with a per-assay `max` recorded in the manifest for dequantization. Values
#'   below ~`max/510` round to zero, so fraction-expressing statistics (dotplot
#'   dot size, DE `pct.1`/`pct.2`) slightly under-count low-expression cells; use
#'   `quantize = FALSE` for exact values.
#' @param counts If `TRUE`, also export each assay's raw `counts` layer to a
#'   `counts/<assay>.parquet` store (integer, unquantized). This roughly doubles
#'   expression storage but enables pseudobulk differential expression
#'   ([scroll_pseudobulk_de()], which needs counts for edgeR/limma-voom).
#' @param subsets Optional named list declaring reprocessed **subset views** — a
#'   slice of the data re-embedded in isolation (e.g. a T-cell-only UMAP) plus any
#'   subset-only metadata (subclusters). Each element is `list(label =, embeddings
#'   =, meta =)`: `embeddings` names one or more exported reductions that cover only
#'   the subset's cells (the first is the view's primary embedding; membership =
#'   cells with non-`NA` coordinates in it); `meta` names metadata columns that are
#'   meaningful only within the subset. In the app a subset becomes a **View** that
#'   restricts every panel to those cells and exposes its embedding + metadata.
#'   Defaults to any spec attached by [scroll_add_subset()].
#' @param overwrite If `TRUE`, an existing `outdir` is removed first.
#' @param verbose If `TRUE`, report progress: a step message per phase and a
#'   progress bar over each assay's feature-partition export (the slow step).
#'   Defaults to [interactive()], so it is quiet in scripts/CI; pass `TRUE` to
#'   force progress in a non-interactive session.
#'
#' @return `outdir`, invisibly.
#' @export
scroll_build <- function(object, outdir,
                         assays = NULL, embeddings = NULL, meta_cols = NULL,
                         quantize = TRUE, counts = FALSE, subsets = NULL,
                         overwrite = FALSE, verbose = interactive()) {
  .scroll_need_seurat()
  object <- .scroll_load_object(object)

  if (is.null(assays)) assays <- SeuratObject::DefaultAssay(object)
  if (is.null(embeddings)) embeddings <- SeuratObject::Reductions(object)
  md <- object[[]]
  if (is.null(meta_cols)) meta_cols <- .scroll_infer_meta_cols(md)
  # a spec attached by scroll_add_subset() lives in @misc (survives Seurat ops)
  if (is.null(subsets)) subsets <- SeuratObject::Misc(object, "scroll_subsets")

  .scroll_check_inputs(object, assays, embeddings, meta_cols)
  subsets <- .scroll_normalize_subsets(subsets, embeddings, meta_cols)

  if (dir.exists(outdir)) {
    if (overwrite) unlink(outdir, recursive = TRUE)
    else if (length(list.files(outdir)) > 0)
      stop("outdir '", outdir, "' exists and is not empty; use overwrite = TRUE.",
           call. = FALSE)
  }
  dir.create(file.path(outdir, "expr"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "de"), showWarnings = FALSE)
  if (counts) dir.create(file.path(outdir, "counts"), showWarnings = FALSE)

  # --- cells.parquet: metadata + all embeddings, the only globally-loaded file
  .scroll_step(verbose, "Extracting cell metadata + embeddings...")
  cells <- .scroll_extract_cells(object, md, meta_cols, embeddings)
  .scroll_check_subset_coverage(subsets, cells)
  arrow::write_parquet(cells, file.path(outdir, "cells.parquet"))

  # --- expr/<assay>/: long, feature-partitioned expression (+ optional counts)
  assay_info <- list()
  for (a in assays) {
    assay_info[[a]] <- .scroll_export_assay(object, a, outdir, quantize, verbose = verbose)
    if (counts) {
      .scroll_step(verbose, sprintf("Exporting raw counts for '%s'...", a))
      .scroll_export_counts(object, a, outdir)
    }
  }

  # --- manifest + authored scaffolding
  .scroll_step(verbose, "Writing manifest + app scaffold...")
  .scroll_write_manifest(outdir, object, assay_info, embeddings, md, meta_cols,
                         n_cells = nrow(cells), quantize = quantize,
                         has_counts = counts, cells = cells, subsets = subsets)
  scroll_scaffold_app(outdir)

  message("scroll project built at: ", normalizePath(outdir))
  invisible(outdir)
}

# SeuratObject is only needed for the offline build (reading the Seurat object);
# the runtime app never uses it, so it lives in Suggests. Fail early and clearly
# if a build is attempted without it.
.scroll_need_seurat <- function() {
  if (!requireNamespace("SeuratObject", quietly = TRUE))
    stop("scroll's build step needs the 'SeuratObject' package. Install it with ",
         "install.packages('SeuratObject'). (The runtime app does not need it.)",
         call. = FALSE)
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

# Emit a progress step message when verbose.
.scroll_step <- function(verbose, ...) if (isTRUE(verbose)) message(...)

# Export one assay's `data` layer as long, feature-partitioned Parquet.
# Returns list(features=, max=, n_features=) for the manifest. The write is done
# in feature batches so a progress bar can advance over the (potentially tens of
# thousands of) partitions; batching also keeps peak memory and open file handles
# bounded. Each feature lands in exactly one batch, so every partition holds one
# part file, identical in content to a single write.
.scroll_export_assay <- function(object, assay, outdir, quantize, verbose = FALSE) {
  mat <- SeuratObject::GetAssayData(object, assay = assay, layer = "data")
  if (is.null(mat) || nrow(mat) == 0 || length(mat@x) == 0)
    stop("Assay '", assay, "' has an empty `data` layer; normalize the object ",
         "before building.", call. = FALSE)

  feats <- rownames(mat)
  cell_names <- colnames(mat)
  .scroll_warn_unsafe_features(feats)
  nfeat <- length(feats)
  .scroll_step(verbose, sprintf("Exporting assay '%s' (%s features)...",
                                assay, format(nfeat, big.mark = ",")))

  trip <- Matrix::summary(mat)             # i (feature), j (cell), x (value)
  max_a <- max(trip$x)
  value <- trip$x
  if (quantize) {
    value <- as.integer(pmin(pmax(round(value / max_a * 255), 0L), 255L))
  }
  # Group nonzeros by feature (contiguous blocks), so a feature batch is a slice.
  o <- order(trip$i)
  fi <- trip$i[o]; cj <- trip$j[o]; vv <- value[o]
  bnd <- c(0L, cumsum(tabulate(fi, nbins = nfeat)))   # feature f rows: (bnd[f]+1):bnd[f+1]

  dest <- file.path(outdir, "expr", assay)
  sch  <- arrow::schema(feature = arrow::utf8(), cell = arrow::utf8(),
                        value = arrow::uint8())
  bsize  <- max(256L, ceiling(nfeat / 40))            # ~40 bar ticks at most
  starts <- seq.int(1L, nfeat, by = bsize)
  pb <- if (isTRUE(verbose)) utils::txtProgressBar(min = 0, max = length(starts), style = 3)
  on.exit(if (!is.null(pb)) close(pb), add = TRUE)

  for (bi in seq_along(starts)) {
    lo <- starts[bi]; hi <- min(lo + bsize - 1L, nfeat)
    r0 <- bnd[lo] + 1L; r1 <- bnd[hi + 1L]
    if (r1 >= r0) {                                    # some features here are nonzero
      idx <- r0:r1
      tbl <- arrow::as_arrow_table(data.frame(
        feature = feats[fi[idx]], cell = cell_names[cj[idx]], value = vv[idx],
        stringsAsFactors = FALSE))
      if (quantize) tbl <- tbl$cast(sch)
      arrow::write_dataset(
        tbl, dest, partitioning = "feature", format = "parquet",
        basename_template = sprintf("part-%d-{i}.parquet", bi),
        max_partitions = (hi - lo + 1L) + 1L)
    }
    if (!is.null(pb)) utils::setTxtProgressBar(pb, bi)
  }
  list(features = feats, max = max_a, n_features = nfeat)
}

# Export one assay's raw `counts` layer as a single long Parquet
# (counts/<assay>.parquet, columns feature/cell/value, integer, unquantized).
# Unlike expr/ this is NOT feature-partitioned: pseudobulk always full-scans all
# genes to aggregate, and one file scans far faster than tens of thousands of
# partitions. Only the nonzero entries are stored (sparse).
.scroll_export_counts <- function(object, assay, outdir) {
  mat <- SeuratObject::GetAssayData(object, assay = assay, layer = "counts")
  if (is.null(mat) || nrow(mat) == 0 || length(mat@x) == 0)
    stop("Assay '", assay, "' has an empty `counts` layer; cannot export counts ",
         "for pseudobulk.", call. = FALSE)
  feats <- rownames(mat); cell_names <- colnames(mat)
  trip <- Matrix::summary(mat)             # i (feature), j (cell), x (count)
  df <- data.frame(feature = feats[trip$i], cell = cell_names[trip$j],
                   value = as.integer(round(trip$x)), stringsAsFactors = FALSE)
  tbl <- arrow::as_arrow_table(df)$cast(arrow::schema(
    feature = arrow::utf8(), cell = arrow::utf8(), value = arrow::int32()))
  arrow::write_parquet(tbl, file.path(outdir, "counts", paste0(assay, ".parquet")))
  invisible(NULL)
}

# Prefer a 2-D visualization embedding over PCA as the plotting default.
.scroll_default_embedding <- function(embeddings) {
  prefer <- c("wnn.umap", "umap", "tsne", "wnnUMAP", "UMAP", "TSNE")
  hit <- intersect(prefer, embeddings)
  if (length(hit)) hit[[1]] else embeddings[[1]]
}

# --- subset views ------------------------------------------------------------

# Validate + normalize the `subsets` spec into a named list of
# list(label, embeddings, primary_embedding, meta). Attaches a `scope_map`
# attribute (meta column -> owning subset). Errors on unknown/duplicate refs.
.scroll_normalize_subsets <- function(subsets, embeddings, meta_cols) {
  if (is.null(subsets) || !length(subsets)) return(NULL)
  if (is.null(names(subsets)) || any(!nzchar(names(subsets))))
    stop("`subsets` must be a named list (subset name -> spec).", call. = FALSE)
  if (anyDuplicated(names(subsets)))
    stop("`subsets` has duplicate names.", call. = FALSE)
  scope_map <- character(0)                    # meta col -> subset name
  out <- lapply(names(subsets), function(nm) {
    s <- subsets[[nm]]
    emb <- as.character(s$embeddings %||% character(0))
    if (!length(emb))
      stop("subset '", nm, "' must name at least one embedding.", call. = FALSE)
    bad <- setdiff(emb, embeddings)
    if (length(bad))
      stop("subset '", nm, "' names unexported embedding(s): ",
           paste(bad, collapse = ", "), "; add them to `embeddings=`.",
           call. = FALSE)
    meta <- as.character(s$meta %||% character(0))
    badm <- setdiff(meta, meta_cols)
    if (length(badm))
      stop("subset '", nm, "' names unexported meta column(s): ",
           paste(badm, collapse = ", "), "; add them to `meta_cols=`.",
           call. = FALSE)
    for (mc in meta) {
      if (mc %in% names(scope_map))
        stop("meta column '", mc, "' is claimed by more than one subset ('",
             scope_map[[mc]], "' and '", nm, "').", call. = FALSE)
      scope_map[[mc]] <<- nm
    }
    list(label = as.character(s$label %||% nm), embeddings = emb,
         primary_embedding = emb[[1]], meta = meta)
  })
  names(out) <- names(subsets)
  attr(out, "scope_map") <- scope_map
  out
}

# A declared subset embedding should be partial (a reprocessed slice). Error on
# zero coverage (barcodes almost certainly don't align); warn on full coverage.
.scroll_check_subset_coverage <- function(subsets, cells) {
  if (is.null(subsets)) return(invisible())
  n <- nrow(cells)
  for (nm in names(subsets)) {
    pe <- subsets[[nm]]$primary_embedding
    cov <- sum(!is.na(cells[[sprintf("%s_1", pe)]]))
    if (cov == 0)
      stop("subset '", nm, "' embedding '", pe, "' covers no cells; check that ",
           "its barcodes match the parent object.", call. = FALSE)
    if (cov == n)
      warning("subset '", nm, "' embedding '", pe, "' covers all ", n,
              " cells; a subset view expects a partial (reprocessed) embedding.",
              call. = FALSE)
  }
  invisible()
}

#' Attach a reprocessed subset onto a parent object
#'
#' Convenience assembly for [scroll_build()]'s `subsets`: copies one or more
#' dimensional reductions (and optional metadata columns) from a separately
#' reprocessed child object onto `object`, aligned **by cell barcode**, and
#' records the subset spec so a later `scroll_build(object, ...)` exposes it as a
#' linked **View**. This is plumbing only — it performs no normalization,
#' clustering, or embedding; the child object must already be reprocessed.
#'
#' @param object The parent Seurat object (all cells).
#' @param sub_object A reprocessed child object (a subset of `object`'s cells,
#'   re-embedded in isolation). Its cell barcodes must match `object`'s.
#' @param name Short id for the subset (e.g. `"tcell"`).
#' @param embeddings Named character vector `new_name = source_reduction` naming
#'   reduction(s) in `sub_object` to copy (e.g. `c(umap_tcell = "umap")`). An
#'   unnamed value reuses the source name. Only the child's cells get coordinates;
#'   parent cells outside the subset are left uncovered (`NA`).
#' @param label Human-readable view label (defaults to `name`).
#' @param meta Metadata columns in `sub_object` to copy onto `object` as
#'   subset-scoped columns (`NA` for non-members).
#' @return `object` with the reduction(s)/metadata added and a `scroll_subsets`
#'   attribute recording the spec.
#' @export
scroll_add_subset <- function(object, sub_object, name, embeddings,
                              label = name, meta = NULL) {
  .scroll_need_seurat()
  if (!inherits(object, "Seurat") || !inherits(sub_object, "Seurat"))
    stop("`object` and `sub_object` must be Seurat objects.", call. = FALSE)
  new_names <- names(embeddings) %||% as.character(embeddings)
  new_names[!nzchar(new_names)] <- as.character(embeddings)[!nzchar(new_names)]
  parent_cells <- colnames(object)
  for (i in seq_along(embeddings)) {
    src <- as.character(embeddings)[[i]]; newn <- new_names[[i]]
    if (!src %in% SeuratObject::Reductions(sub_object))
      stop("sub_object has no reduction '", src, "'.", call. = FALSE)
    em <- SeuratObject::Embeddings(sub_object, reduction = src)
    keep <- rownames(em) %in% parent_cells
    if (!any(keep))
      stop("no barcodes of sub_object reduction '", src,
           "' match the parent object.", call. = FALSE)
    em <- em[keep, , drop = FALSE]
    colnames(em) <- sprintf("%s_%d", newn, seq_len(ncol(em)))
    object[[newn]] <- SeuratObject::CreateDimReducObject(
      embeddings = em, key = sprintf("%s_", gsub("[^A-Za-z0-9]", "", newn)),
      assay = SeuratObject::DefaultAssay(object))
  }
  for (mc in meta) {
    src_vals <- sub_object[[mc]][[1]]
    col <- rep(NA_character_, length(parent_cells)); names(col) <- parent_cells
    idx <- match(colnames(sub_object), parent_cells)
    ok <- !is.na(idx)
    col[idx[ok]] <- as.character(src_vals)[ok]
    object[[mc]] <- unname(col)
  }
  spec <- list(label = label, embeddings = unname(new_names),
               meta = as.character(meta %||% character(0)))
  # Store in @misc (not a bare attribute): downstream Seurat operations
  # (SetAssayData, subset, ...) preserve @misc but drop object attributes.
  prev <- SeuratObject::Misc(object, "scroll_subsets") %||% list()
  prev[[name]] <- spec
  SeuratObject::Misc(object, "scroll_subsets") <- prev
  object
}

.scroll_warn_unsafe_features <- function(feats) {
  unsafe <- grep("[/\\\\]", feats, value = TRUE)
  if (length(unsafe))
    warning("Some feature names contain path separators and may break ",
            "partitioning: ", paste(utils::head(unsafe, 5), collapse = ", "),
            call. = FALSE)
}
