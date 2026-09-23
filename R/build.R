#' Build a scroll project from a Seurat object
#'
#' Reads a processed Seurat object exactly once and emits a self-contained
#' project directory of lightweight on-disk artifacts: a small Parquet table of
#' cell metadata + embeddings (`cells.parquet`), a feature-partitioned Parquet
#' store of expression (`expr/<assay>/`), a `manifest.yaml`, and app scaffolding
#' (`config.yaml`, `app.R`). The running app never touches the object again.
#'
#' @param object A processed Seurat object, or a path to an `.rds` file.
#' @param outdir Output project directory (created if missing).
#' @param assays Assays to export. `NULL` exports the default assay only.
#' @param embeddings Reductions to export. `NULL` exports all reductions.
#' @param meta_cols Metadata columns to export. `NULL` infers a sensible set.
#' @param max_levels Cap on how many distinct values of a categorical column are
#'   cached as a `levels:` list in `manifest.yaml`. Columns above the cap (e.g. a
#'   clone id or barcode with thousands of values) still export normally and stay
#'   fully usable — the manifest just records their `n_levels` count and the app
#'   recomputes the level set from `cells.parquet` when a control needs it. Keeps
#'   the manifest small and hand-editable. Defaults to 200.
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
#' @param vdj Optional [vdj_spec()] describing per-cell TCR/BCR columns. When
#'   supplied, a compact `repertoire/` store (clone table + diversity / gene-usage
#'   residuals / tissue correlation) is baked and the VDJ panels are enabled.
#' @param spatial Optional [spatial_spec()] (or `TRUE` for defaults) describing a
#'   tissue image, or a **list** of specs for multi-FOV / multi-slide objects. For
#'   each spec `GetTissueCoordinates()` is exported as a spatial embedding, the H&E
#'   image is baked as a raster asset, and the Spatial panel is enabled. Multiple
#'   specs must have distinct `name`s.
#' @param atac Optional [atac_spec()] (or `TRUE` for defaults) naming a peaks assay.
#'   When supplied, the assay is marked `kind: peaks`, each peak's coordinates are
#'   parsed from its name, a peak-annotation table (with nearest gene when a Signac
#'   annotation is available) is baked, and the Peaks panel is enabled.
#' @param panels Optional character vector of panel ids to show **exactly**, in this
#'   order — an explicit override of the automatic data-driven gating (e.g.
#'   `c("dimplot", "featureplot", "spatial")`). Written to `config.yaml` as `panels:`
#'   and also hand-editable there. `NULL` (default) keeps the automatic behaviour.
#' @param exclude_panels Optional character vector of panel ids to hide (e.g.
#'   `c("de", "pseudobulk")`); applied after gating / the `panels` allowlist.
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
                         vdj = NULL, spatial = NULL, atac = NULL,
                         panels = NULL, exclude_panels = NULL,
                         max_levels = .SCROLL_MAX_LEVELS,
                         overwrite = FALSE, verbose = interactive()) {
  if (!is.null(panels) && !is.character(panels))
    stop("`panels` must be a character vector of panel ids, or NULL.", call. = FALSE)
  if (!is.null(exclude_panels) && !is.character(exclude_panels))
    stop("`exclude_panels` must be a character vector of panel ids, or NULL.", call. = FALSE)
  .scroll_need_seurat()
  .scroll_check_zstd()
  object <- .scroll_load_object(object)

  if (is.null(assays)) assays <- SeuratObject::DefaultAssay(object)
  if (is.null(embeddings)) embeddings <- SeuratObject::Reductions(object)
  md <- object[[]]
  # a spec attached by scroll_add_subset() lives in @misc (survives Seurat ops)
  if (is.null(subsets)) subsets <- SeuratObject::Misc(object, "scroll_subsets")
  if (is.null(meta_cols)) {
    # inferred: always keep the subsets' declared columns, which inference could
    # otherwise drop (e.g. a high-cardinality scoped label) and then fail validation
    meta_cols <- union(.scroll_infer_meta_cols(md),
                       unlist(lapply(subsets, `[[`, "meta"), use.names = FALSE))
  }

  # Spatial: extract tissue coordinates into a `spatial` reduction (image-pixel
  # space) so DimPlot/FeaturePlot/the Spatial panel can plot on it. Injected
  # before validation so it counts as an exported embedding.
  # Accept a single spatial_spec (or TRUE for defaults) or a LIST of specs — the
  # latter bakes several tissue maps (multi-FOV / multi-slide) as distinct spatial
  # embeddings, each with its own image asset.
  spatial_preps <- list()
  if (!is.null(spatial)) {
    if (isTRUE(spatial)) spatial <- spatial_spec()
    specs <- if (inherits(spatial, "scroll_spatial_spec")) list(spatial) else spatial
    for (sp in specs) {
      if (!inherits(sp, "scroll_spatial_spec"))
        stop("spatial: each entry must be a spatial_spec().", call. = FALSE)
      prep <- .scroll_prepare_spatial(object, sp)
      if (!is.null(spatial_preps[[prep$name]]))
        stop("spatial: duplicate embedding name '", prep$name,
             "'; give each spatial_spec() a distinct `name`.", call. = FALSE)
      # the reduction key may collide case-insensitively with the assay key
      # (e.g. `spatial_` vs the `Spatial` assay); harmless — coords are keyed by
      # the reduction *name*, so silence the Seurat key-rename notice.
      suppressWarnings(object[[prep$name]] <- prep$reduction)
      embeddings <- union(embeddings, prep$name)
      spatial_preps[[prep$name]] <- prep
    }
  }

  .scroll_check_inputs(object, assays, embeddings, meta_cols)
  subsets <- .scroll_normalize_subsets(subsets, embeddings, meta_cols)

  if (dir.exists(outdir)) {
    if (overwrite) unlink(outdir, recursive = TRUE)
    else if (length(list.files(outdir)) > 0)
      stop("outdir '", outdir, "' exists and is not empty; use overwrite = TRUE.",
           call. = FALSE)
  }
  dir.create(file.path(outdir, "expr"), recursive = TRUE, showWarnings = FALSE)
  if (counts) dir.create(file.path(outdir, "counts"), showWarnings = FALSE)

  # --- cells.parquet: metadata + all embeddings, the only globally-loaded file
  .scroll_step(verbose, "Extracting cell metadata + embeddings...")
  cells <- .scroll_extract_cells(object, md, meta_cols, embeddings)
  .scroll_check_subset_coverage(subsets, cells)
  arrow::write_parquet(cells, file.path(outdir, "cells.parquet"), compression = "zstd")

  # --- expr/<assay>/: long, feature-partitioned expression (+ optional counts)
  assay_info <- list()
  for (a in assays) {
    assay_info[[a]] <- .scroll_export_assay(object, a, outdir, quantize, verbose = verbose,
                                            cell_ids = cells$cell)
    if (counts) {
      .scroll_step(verbose, sprintf("Exporting raw counts for '%s'...", a))
      .scroll_export_counts(object, a, outdir, cell_ids = cells$cell)
    }
  }

  # --- compact each assay to a single part-file (a single-object build already
  # writes one pre-sorted file; this matters for streaming/appended stores, whose
  # per-source parts are merged + feature-sorted here).
  .scroll_compact_store(outdir, quantize)

  # --- VDJ repertoire store (optional): bake the clone table + diversity/usage
  vdj_block <- NULL
  if (!is.null(vdj)) {
    .scroll_step(verbose, "Baking VDJ repertoire store...")
    vdj_block <- .scroll_bake_vdj(md, outdir, vdj)
  }

  # --- Spatial tissue image(s) (optional): bake the H&E raster asset(s)
  images_block <- NULL
  if (length(spatial_preps)) {
    .scroll_step(verbose, "Baking spatial tissue image(s)...")
    for (prep in spatial_preps)
      images_block <- c(images_block, .scroll_write_spatial(outdir, prep))
  }

  # --- scATAC peak annotation (optional): parse coords + nearest gene
  atac_block <- NULL
  if (!is.null(atac)) {
    if (isTRUE(atac)) atac <- atac_spec()
    if (!atac$assay %in% assays)
      stop("atac: peaks assay '", atac$assay, "' must be listed in `assays=`.", call. = FALSE)
    .scroll_step(verbose, "Baking scATAC peak table...")
    atac_block <- .scroll_bake_atac(object, outdir, atac)
  }

  # --- manifest + authored scaffolding
  .scroll_step(verbose, "Writing manifest + app scaffold...")
  .scroll_write_manifest(outdir, object, assay_info, embeddings, md, meta_cols,
                         n_cells = nrow(cells), quantize = quantize,
                         has_counts = counts, cells = cells, subsets = subsets,
                         vdj = vdj_block, images = images_block,
                         spatial_embeddings = if (length(spatial_preps)) names(spatial_preps) else character(),
                         atac = atac_block,
                         peaks_assay = if (!is.null(atac_block)) atac_block$assay else character(),
                         max_levels = max_levels)
  scroll_scaffold_app(outdir, panels = panels, exclude_panels = exclude_panels)

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

# scroll writes AND reads every store with compression = "zstd". A "minimal"
# arrow build (common on HPC when install.packages() can't fetch the prebuilt
# libarrow) ships without codecs, so a store cannot be built or even read —
# otherwise arrow fails deep inside with a cryptic "NotImplemented: Support for
# codec 'zstd' not built". Fail early at the entry points with a fix instead.
.scroll_check_zstd <- function() {
  ok <- tryCatch(isTRUE(arrow::arrow_info()$capabilities[["zstd"]]),
                 error = function(e) FALSE)
  if (!ok)
    stop("The installed 'arrow' was built without the zstd codec, which scroll ",
         "needs to read and write its stores.\n",
         "Reinstall a full arrow build in a FRESH R session:\n",
         "  Sys.setenv(NOT_CRAN = \"true\", LIBARROW_MINIMAL = \"false\", ",
         "LIBARROW_BINARY = \"true\"); install.packages(\"arrow\")\n",
         "then confirm: arrow::arrow_info()$capabilities[[\"zstd\"]]  # must be TRUE\n",
         "(On HPC without internet, `conda install -c conda-forge r-arrow` also works.)",
         call. = FALSE)
  invisible(TRUE)
}

# GetAssayData(layer=) can hand back a DENSE base matrix (some v5 layers, small or
# ADT assays), which has no @x slot and no sparse triplet — so scroll's export
# (mat@x, Matrix::summary) fails with "no applicable method for `@`". Coerce any
# non-Csparse input (dense base matrix, dgeMatrix, dgTMatrix) to a dgCMatrix.
.scroll_as_sparse <- function(mat) {
  if (methods::is(mat, "CsparseMatrix")) return(mat)
  methods::as(Matrix::Matrix(mat, sparse = TRUE), "CsparseMatrix")
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
    if (is_cat && n_lvl > 1 && n_lvl <= .SCROLL_MAX_LEVELS) keep <- c(keep, col)
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

# Guard the load-bearing invariant behind the v2 integer cell-index: the store's
# `cell` column is the assay matrix COLUMN position, which is only a valid row index
# into cells.parquet if the matrix columns are in the SAME order as the cell metadata
# rows (`cell_ids`, i.e. rownames(md) / cells$cell). Seurat maintains this, but a
# mismatch would silently mis-join every expression value, so fail the build loudly.
.scroll_check_cell_order <- function(mat, cell_ids, assay) {
  if (is.null(cell_ids)) return(invisible())
  if (!identical(colnames(mat), as.character(cell_ids)))
    stop(sprintf(
      "Assay '%s': the %d assay-matrix columns are not identical (in order) to the %s cell-metadata rows, so the store's integer cell-index would be mis-aligned. Reorder the assay to match colnames == rownames(metadata) before building.",
      assay, ncol(mat), format(length(cell_ids), big.mark = ",")),
      call. = FALSE)
}

# Export one assay's `data` layer as long, feature-partitioned Parquet.
# Returns list(features=, max=, n_features=) for the manifest. The write is done
# in feature batches so a progress bar can advance over the (potentially tens of
# thousands of) partitions; batching also keeps peak memory and open file handles
# bounded. Each feature lands in exactly one batch, so every partition holds one
# part file, identical in content to a single write.
.scroll_export_assay <- function(object, assay, outdir, quantize, verbose = FALSE,
                                 offset = 0L, cell_ids = NULL) {
  # `offset` shifts the (1-based) local cell index onto a GLOBAL row index so a
  # streaming builder can append one source at a time into a shared store
  # (scroll_build uses offset = 0: local index == global index).
  mat <- SeuratObject::GetAssayData(object, assay = assay, layer = "data")
  if (!is.null(mat) && !methods::is(mat, "CsparseMatrix")) mat <- .scroll_as_sparse(mat)
  if (is.null(mat) || nrow(mat) == 0 || length(mat@x) == 0)
    stop("Assay '", assay, "' has an empty `data` layer; normalize the object ",
         "before building.", call. = FALSE)
  .scroll_check_cell_order(mat, cell_ids, assay)   # int cell-index <-> cells.parquet row order

  feats <- rownames(mat)
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

  dest  <- file.path(outdir, "expr", assay)
  # v2 store: one Parquet file per assay (expr/<assay>/part-<offset>.parquet) with
  # `feature` as a real column. Rows are written sorted by feature NAME so Parquet
  # row-group min/max statistics prune a single-gene query to just that gene's row
  # group(s) -- the pruning, not any directory split, is what keeps a lookup cheap
  # (benchmarked equivalent to the former first-letter bucketing). The read path
  # (`arrow::open_dataset` on expr/<assay> + a `feature` filter) is layout-agnostic,
  # so it reads this single file exactly as it read the old bucketed subdirs.
  vtype <- if (quantize) arrow::uint8() else arrow::float32()
  # `feature` is large_utf8 (int64 offsets): one gene name per nonzero can exceed
  # 2 GB of string bytes on a large assay, which overflows plain utf8's int32
  # offsets ("Failed casting from large_string to string: input array too large").
  sch   <- arrow::schema(feature = arrow::large_utf8(), cell = arrow::int32(), value = vtype)

  ford <- order(feats)                                 # features alphabetical by name
  # gather each feature's contiguous nonzero rows, in sorted-feature order
  idx <- unlist(lapply(ford, function(f)
    if (bnd[f + 1L] > bnd[f]) (bnd[f] + 1L):bnd[f + 1L] else integer(0)),
    use.names = FALSE)
  if (length(idx)) {
    tbl <- arrow::as_arrow_table(data.frame(
      feature = feats[fi[idx]], cell = offset + cj[idx], value = vv[idx],   # global cell index
      stringsAsFactors = FALSE))$cast(sch)
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    # part name carries the source offset so a streaming builder's per-source parts
    # coexist in one dir (compacted into a single file at the end of the run).
    arrow::write_parquet(tbl, file.path(dest, sprintf("part-%d.parquet", offset)),
                         compression = "zstd")
  }
  list(features = feats, max = max_a, n_features = nfeat)
}

# Export one assay's raw `counts` layer as a single long Parquet
# (counts/<assay>.parquet, columns feature/cell/value, integer, unquantized).
# Unlike expr/ this is NOT feature-partitioned: pseudobulk always full-scans all
# genes to aggregate, and one file scans far faster than tens of thousands of
# partitions. Only the nonzero entries are stored (sparse).
.scroll_export_counts <- function(object, assay, outdir, cell_ids = NULL) {
  mat <- SeuratObject::GetAssayData(object, assay = assay, layer = "counts")
  if (!is.null(mat) && !methods::is(mat, "CsparseMatrix")) mat <- .scroll_as_sparse(mat)
  if (is.null(mat) || nrow(mat) == 0 || length(mat@x) == 0)
    stop("Assay '", assay, "' has an empty `counts` layer; cannot export counts ",
         "for pseudobulk.", call. = FALSE)
  .scroll_check_cell_order(mat, cell_ids, assay)   # int cell-index <-> cells.parquet row order
  feats <- rownames(mat)
  trip <- Matrix::summary(mat)             # i (feature), j (cell), x (count)
  df <- data.frame(feature = feats[trip$i], cell = trip$j,          # int32 cell-index (v2)
                   value = as.integer(round(trip$x)), stringsAsFactors = FALSE)
  tbl <- arrow::as_arrow_table(df)$cast(arrow::schema(
    feature = arrow::large_utf8(), cell = arrow::int32(), value = arrow::int32()))
  arrow::write_parquet(tbl, file.path(outdir, "counts", paste0(assay, ".parquet")),
                       compression = "zstd")
  invisible(NULL)
}

# Compact each expr/<assay>/ that has >1 part-file (a streaming build appends one
# part per source) into a single zstd part-file, re-sorted by feature so row-group
# statistics prune. Preserves the int32/float32(uint8) schema and reads one assay at
# a time (flat memory). No-op for assays already holding a single file (the common
# single-object build, which .scroll_export_assay writes pre-sorted).
.scroll_compact_store <- function(outdir, quantize) {
  vtype <- if (quantize) arrow::uint8() else arrow::float32()
  sch <- arrow::schema(feature = arrow::large_utf8(), cell = arrow::int32(), value = vtype)
  root <- file.path(outdir, "expr")
  for (adir in list.dirs(root, recursive = FALSE)) {     # expr/<assay> dirs
    parts <- list.files(adir, pattern = "\\.parquet$", full.names = TRUE)
    if (length(parts) <= 1L) next
    tab <- dplyr::collect(arrow::open_dataset(adir))      # (feature, cell, value) in RAM
    tab <- tab[order(tab$feature), , drop = FALSE]        # re-sort so row-group stats prune
    unlink(parts)
    arrow::write_parquet(arrow::as_arrow_table(tab)$cast(sch),
                         file.path(adir, "part-0.parquet"), compression = "zstd")
  }
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
#'   parent cells outside the subset are left uncovered (`NA`). The first is the
#'   view's primary embedding (membership = cells with coordinates in it).
#'   `"auto"` (default) copies every reduction that is new in `sub_object` or whose
#'   coordinates differ from `object`'s on the shared cells (i.e. it was re-run on
#'   the subset), as `<name>_<reduction>` (e.g. `tcell_umap`), ordered UMAP, then
#'   t-SNE, then the rest so a 2-D map is the primary. Reductions the child merely
#'   inherited unchanged are skipped. Every copied dim becomes a `cells.parquet`
#'   column, so a 50-dim PCA/Harmony adds 50 (mostly `NA`) columns -- pass an
#'   explicit vector to keep only the maps you want.
#' @param label Human-readable view label (defaults to `name`).
#' @param meta Metadata columns in `sub_object` to copy onto `object` as
#'   subset-scoped columns (`NA` for non-members). `"auto"` detects them: every
#'   column that is new in `sub_object`, or whose values differ from `object`'s on
#'   the shared cells (e.g. re-run clusters or module scores). Unchanged inherited
#'   columns (`sample`, QC metrics, ...) are skipped, and the chosen columns are
#'   reported with a message -- check them before building.
#' @param prefix Prefix for the copied column names. Defaults to `"<name>_"` with
#'   `meta = "auto"` (so a subset's `seurat_clusters` becomes
#'   `tcell_seurat_clusters` rather than overwriting the parent's) and to `""`
#'   (copy under the same name) for an explicit `meta`.
#' @return `object` with the reduction(s)/metadata added and a `scroll_subsets`
#'   attribute recording the spec.
#' @export
scroll_add_subset <- function(object, sub_object, name, embeddings = "auto",
                              label = name, meta = NULL,
                              prefix = if (identical(meta, "auto")) paste0(name, "_") else "") {
  .scroll_need_seurat()
  if (!inherits(object, "Seurat") || !inherits(sub_object, "Seurat"))
    stop("`object` and `sub_object` must be Seurat objects.", call. = FALSE)
  force(prefix)   # evaluate against the caller's `meta` before it is resolved below
  if (identical(embeddings, "auto")) {
    reds <- .scroll_subset_reductions(object, sub_object)
    if (!length(reds))
      stop("scroll_add_subset('", name, "'): sub_object has no reduction that is new ",
           "or differs from the parent's -- re-embed the subset (e.g. RunUMAP), or ",
           "name one explicitly with `embeddings = c(new_name = \"reduction\")`.",
           call. = FALSE)
    embeddings <- stats::setNames(reds, paste0(name, "_", reds))
    dims <- vapply(reds, function(r) ncol(SeuratObject::Embeddings(sub_object, r)), integer(1))
    message(sprintf("scroll_add_subset('%s'): copying %d embedding(s) as %s", name,
                    length(reds), paste0(names(embeddings), " (", dims, "d)", collapse = ", ")))
  }
  if (identical(meta, "auto")) {
    meta <- .scroll_subset_meta_cols(object, sub_object)
    message(sprintf("scroll_add_subset('%s'): %s", name,
      if (length(meta)) paste0("copying ", length(meta), " subset column(s) as ",
                               paste0(prefix, meta, collapse = ", "))
      else "no new or changed metadata columns found"))
  }
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
    # Preserve the source column's type so numeric scores (e.g. module scores) stay
    # numeric rather than being stringified into a high-cardinality categorical. A
    # factor is flattened to character (matching .scroll_extract_cells); non-members
    # get the type-appropriate NA. (The `as.character`-everything path here used to
    # misclassify continuous scoped columns as categorical.)
    fill <- if (is.numeric(src_vals)) NA_real_
            else if (is.logical(src_vals)) NA
            else NA_character_
    src <- if (is.factor(src_vals)) as.character(src_vals) else src_vals
    col <- rep(fill, length(parent_cells)); names(col) <- parent_cells
    idx <- match(colnames(sub_object), parent_cells)
    ok <- !is.na(idx)
    col[idx[ok]] <- src[ok]
    object[[paste0(prefix, mc)]] <- unname(col)
  }
  meta_out <- if (length(meta)) paste0(prefix, meta) else character(0)
  spec <- list(label = label, embeddings = unname(new_names),
               meta = as.character(meta_out))
  # Store in @misc (not a bare attribute): downstream Seurat operations
  # (SetAssayData, subset, ...) preserve @misc but drop object attributes.
  prev <- SeuratObject::Misc(object, "scroll_subsets") %||% list()
  prev[[name]] <- spec
  SeuratObject::Misc(object, "scroll_subsets") <- prev
  object
}

# Reductions a reprocessed child adds or re-ran relative to its parent: new names,
# plus same-named ones whose coordinates differ on the shared cells (a different
# dim count counts as different). Ordered UMAP > t-SNE > rest, so the view's
# primary embedding (the first) is a 2-D map when there is one.
.scroll_subset_reductions <- function(parent, child) {
  reds <- SeuratObject::Reductions(child)
  if (!length(reds)) return(character(0))
  have <- SeuratObject::Reductions(parent)
  changed <- vapply(reds, function(r) {
    if (!r %in% have) return(TRUE)
    ce <- SeuratObject::Embeddings(child, r); pe <- SeuratObject::Embeddings(parent, r)
    cells <- intersect(rownames(ce), rownames(pe))
    if (!length(cells) || ncol(ce) != ncol(pe)) return(TRUE)
    !isTRUE(all.equal(unname(ce[cells, , drop = FALSE]), unname(pe[cells, , drop = FALSE])))
  }, logical(1))
  reds <- reds[changed]
  rank <- ifelse(grepl("umap", reds, ignore.case = TRUE), 1L,
                 ifelse(grepl("tsne", reds, ignore.case = TRUE), 2L, 3L))
  reds[order(rank, seq_along(reds))]
}

# Metadata columns a reprocessed child adds or changes relative to its parent,
# compared on the shared cells: new columns, plus shared ones whose values differ
# (re-run clusters, recomputed scores). Numerics compare with a tolerance so float
# noise is not a change; everything else compares as character, so a factor whose
# levels were merely re-ordered is not flagged.
.scroll_subset_meta_cols <- function(parent, child) {
  cells <- intersect(colnames(child), colnames(parent))
  if (!length(cells))
    stop("no barcodes of sub_object match the parent object.", call. = FALSE)
  pm <- parent[[]][cells, , drop = FALSE]
  cm <- child[[]][cells, , drop = FALSE]
  shared <- intersect(names(cm), names(pm))
  same <- vapply(shared, function(col) {
    a <- pm[[col]]; b <- cm[[col]]
    if (is.numeric(a) && is.numeric(b))
      isTRUE(all.equal(unname(a), unname(b), check.attributes = FALSE))
    else identical(as.character(a), as.character(b))
  }, logical(1))
  c(setdiff(names(cm), names(pm)), shared[!same])
}

.scroll_warn_unsafe_features <- function(feats) {
  unsafe <- grep("[/\\\\]", feats, value = TRUE)
  if (length(unsafe))
    warning("Some feature names contain path separators and may break ",
            "partitioning: ", paste(utils::head(unsafe, 5), collapse = ", "),
            call. = FALSE)
}
