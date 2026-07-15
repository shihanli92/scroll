# scATAC (chromatin accessibility) support. Peaks are just an assay, so once a
# peaks assay is exported, FeaturePlot / DotPlot / Violin / DE work on it for free
# (accessibility instead of expression). What ATAC adds is genomic context:
# `scroll_build(..., atac = atac_spec("peaks"))` marks the assay `kind: peaks`,
# parses each peak's chr/start/end from its "chr-start-end" name, and bakes a
# `peaks.parquet` annotation table (with the nearest gene when a Signac
# ChromatinAssay annotation is available). A built-in Peaks panel then lets you
# search accessibility by gene -> nearby peak. Coverage/genome-browser tracks and
# motif views are out of scope. Non-ATAC projects are unaffected (the panel is
# gated on the manifest `atac` block).

#' Describe a dataset's scATAC peaks assay for `scroll_build()`
#'
#' Marks which exported assay holds chromatin-accessibility **peaks** and bakes a
#' peak-annotation table. Pass the result as `scroll_build(..., atac = atac_spec())`.
#' The peaks assay must still be listed in `assays=` so its accessibility matrix is
#' exported like any other assay; `atac_spec()` additionally records genomic
#' coordinates (parsed from the peak names, e.g. `chr1-9776-10668`) and, when a
#' Signac `ChromatinAssay` annotation is present, each peak's nearest gene.
#'
#' @param assay Name of the peaks assay (default `"peaks"`). Must be one of the
#'   exported `assays`.
#' @param nearest_gene Whether to compute each peak's nearest gene via
#'   [Signac::ClosestFeature()] when the assay carries an annotation (default
#'   `TRUE`). Falls back to coordinates-only if unavailable.
#' @return A `scroll_atac_spec` list.
#' @seealso [scroll_build()]
#' @export
atac_spec <- function(assay = "peaks", nearest_gene = TRUE) {
  structure(list(assay = assay, nearest_gene = isTRUE(nearest_gene)),
            class = "scroll_atac_spec")
}

# ---- build-time: parse coordinates + bake the peak table --------------------

# Parse "chr-start-end" / "chr:start-end" peak names into a coordinate frame.
# Rows that don't match stay NA (kept, so the accessibility store is unaffected).
.scroll_parse_peaks <- function(feats) {
  m <- regmatches(feats, regexec("^(.+?)[:-](\\d+)-(\\d+)$", feats))
  ok <- lengths(m) == 4L
  chr <- rep(NA_character_, length(feats))
  start <- rep(NA_real_, length(feats)); end <- start
  chr[ok]   <- vapply(m[ok], `[`, "", 2L)
  start[ok] <- as.numeric(vapply(m[ok], `[`, "", 3L))
  end[ok]   <- as.numeric(vapply(m[ok], `[`, "", 4L))
  data.frame(feature = feats, chr = chr, start = start, end = end,
             width = end - start + 1, stringsAsFactors = FALSE)
}

# Nearest gene per peak from a Signac ChromatinAssay annotation, or NULL. Returns a
# data.frame(feature, nearest_gene, distance) aligned to `feats`.
.scroll_peak_nearest_gene <- function(object, assay, feats) {
  if (!requireNamespace("Signac", quietly = TRUE)) return(NULL)
  cf <- tryCatch(Signac::ClosestFeature(object, regions = feats, assay = assay),
                 error = function(e) NULL)
  if (is.null(cf) || !nrow(cf)) return(NULL)
  gene <- cf$gene_name %||% cf$gene_id %||% cf$tx_id
  if (is.null(gene)) return(NULL)
  # ClosestFeature returns rows in `regions` order, keyed by query_region
  key <- cf$query_region %||% cf$feature %||% feats
  idx <- match(feats, key)
  data.frame(feature = feats,
             nearest_gene = as.character(gene)[idx],
             distance = (cf$distance %||% rep(NA_real_, nrow(cf)))[idx],
             stringsAsFactors = FALSE)
}

# Write peaks.parquet (coords + optional nearest gene) and return the manifest
# `atac` block. Called by scroll_build() when `atac=` is supplied.
.scroll_bake_atac <- function(object, outdir, spec) {
  assay <- spec$assay
  mat <- SeuratObject::GetAssayData(object, assay = assay, layer = "data")
  feats <- rownames(mat)
  if (is.null(feats) || !length(feats))
    stop("atac: peaks assay '", assay, "' has no features.", call. = FALSE)

  tab <- .scroll_parse_peaks(feats)
  ng <- if (isTRUE(spec$nearest_gene)) .scroll_peak_nearest_gene(object, assay, feats) else NULL
  if (!is.null(ng)) {
    tab$nearest_gene <- ng$nearest_gene
    tab$distance <- ng$distance
  } else {
    tab$nearest_gene <- NA_character_
    tab$distance <- NA_real_
  }

  dir.create(file.path(outdir, "atac"), showWarnings = FALSE, recursive = TRUE)
  arrow::write_parquet(tab, file.path(outdir, "atac", "peaks.parquet"), compression = "zstd")

  list(assay = assay, n_peaks = length(feats),
       n_parsed = sum(!is.na(tab$chr)),
       has_nearest_gene = any(!is.na(tab$nearest_gene)),
       file = file.path("atac", "peaks.parquet"))
}

