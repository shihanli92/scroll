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

# ---- runtime: read the baked table ------------------------------------------

.scroll_atac_read <- function(data) {
  blk <- data$manifest$atac
  if (is.null(blk) || is.null(blk$file)) return(NULL)
  .scroll_read_asset(file.path(data$dir, blk$file), "parquet")
}

# ---- the Peaks panel --------------------------------------------------------

# Gene-centric accessibility search: pick a gene -> a nearby peak -> its
# accessibility on the embedding. Gated on the manifest `atac` block.
.scroll_atac_panels <- function() {
  gate <- function(m) !is.null(m$atac)

  genes_of <- function(data) {
    pk <- .scroll_atac_read(data)
    if (is.null(pk) || !any(!is.na(pk$nearest_gene))) return(character())
    sort(unique(stats::na.omit(pk$nearest_gene)))
  }
  # peaks matching the current search: by nearest gene, or by a text/region filter
  # over the peak names (the fallback when no gene annotation was baked).
  has_genes <- function(data) {
    pk <- .scroll_atac_read(data); !is.null(pk) && any(!is.na(pk$nearest_gene))
  }
  # Region search leads when there's no gene annotation, so the panel is usable out
  # of the box; the peak list pre-fills (even before typing) so something is shown.
  mode_choices <- function(data) if (has_genes(data)) c("Gene", "Region") else c("Region", "Gene")

  peaks_for <- function(input, data) {
    pk <- .scroll_atac_read(data)
    if (is.null(pk)) return(character())
    mode <- input$mode %||% mode_choices(data)[[1]]
    if (identical(mode, "Region")) {
      q <- input$region
      hits <- if (is.null(q) || !nzchar(q)) pk$feature       # empty -> show the first peaks
              else pk$feature[grepl(q, pk$feature, fixed = TRUE)]
    } else {
      g <- input$gene
      hits <- if (!is.null(g) && nzchar(g) && any(!is.na(pk$nearest_gene)))
        pk$feature[!is.na(pk$nearest_gene) & pk$nearest_gene == g] else character()
    }
    utils::head(sort(hits), 200L)                    # bound the client-side list
  }

  controls <- list(
    scroll_input_choice("mode", "Find peak by", choices = mode_choices),
    scroll_show_when(scroll_input_choice("gene", "Gene", choices = genes_of, widget = "select"),
                     control = "mode", equals = "Gene"),
    scroll_show_when(scroll_input_text("region", "Region contains (e.g. chr1-807)"),
                     control = "mode", equals = "Region"),
    scroll_input_levels("peak", "Peak", choices = peaks_for, watch = c("mode", "gene", "region")),
    scroll_input_embedding("embedding", "Embedding"))

  plot <- function(cells, input, data) {
    peak <- input$peak
    if (is.null(peak) || !length(peak) || !nzchar(peak[[1]]))
      stop("Pick a gene, then a nearby peak.")
    peak <- peak[[1]]
    assay <- data$manifest$atac$assay
    emb <- input$embedding %||% .scroll_view_embeddings(data$manifest, NULL)[[1]]
    df <- .scroll_embedding_xy(cells, emb)
    if (!nrow(df)) stop("No cells with coordinates for this embedding.")
    df$.acc <- .scroll_expr_vector(df, data$query1(assay, peak))
    pk <- .scroll_atac_read(data); row <- pk[match(peak, pk$feature), ]
    sub <- if (!is.null(row) && !is.na(row$nearest_gene))
      sprintf("nearest gene: %s", row$nearest_gene) else "chromatin accessibility"

    ggplot2::ggplot(df[order(df$.acc), ], ggplot2::aes(.data$.x, .data$.y, color = .data$.acc)) +
      .scroll_point_layer(size = 0.6, raster = isTRUE(nrow(df) > .scroll_raster_threshold)) +
      .scroll_continuous_scale("viridis", name = "accessibility") +
      ggplot2::labs(title = peak, subtitle = sub, x = paste0(emb, "_1"), y = paste0(emb, "_2")) +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                     panel.grid.minor = ggplot2::element_blank())
  }

  list(.scroll_gated_panel("peaks", "Peaks", "Chromatin accessibility",
    "Find a peak near a gene and colour the embedding by its accessibility.",
    gate, controls, plot))
}
