# The scATAC Peaks panel + its runtime accessor. Build-time peak parsing + bake
# logic lives in R/atac.R; this reads the baked peaks table and renders the
# gene/region-searchable accessibility panel.

# ---- runtime: read the baked table ------------------------------------------

.scroll_atac_read <- function(data) {
  blk <- data$manifest$atac
  if (is.null(blk) || is.null(blk$file)) return(NULL)
  .scroll_read_asset(file.path(data$dir, blk$file), "parquet")
}

# Longest peak list shipped to the client selectize (guards a broad region / a gene
# with thousands of peaks). Truncation is surfaced, not silent (see the plot subtitle).
.SCROLL_ATAC_MAX_PEAKS <- 500L

# Peaks matching a region query. A `chr:start-end` (or `chr-start-end`) query is a
# genomic-interval OVERLAP against the baked coordinates; anything else falls back to
# a substring match over the peak names (so a bare chr / coordinate prefix still works).
.scroll_atac_region_hits <- function(pk, q) {
  q <- trimws(q %||% "")
  if (!nzchar(q)) return(pk$feature)                       # empty -> show the first peaks
  iv <- regmatches(q, regexec("^(.+?)[:-](\\d+)-(\\d+)$", q))[[1]]
  if (length(iv) == 4L) {                                  # chr:start-end -> overlap
    chr <- iv[2]; lo <- min(as.numeric(iv[3:4])); hi <- max(as.numeric(iv[3:4]))
    ok <- !is.na(pk$chr) & pk$chr == chr &
          !is.na(pk$start) & !is.na(pk$end) & pk$start <= hi & pk$end >= lo
    return(pk$feature[ok])
  }
  pk$feature[grepl(q, pk$feature, fixed = TRUE)]           # substring fallback
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

  # All peaks matching the current search (uncapped) -- the choice list caps it, and
  # the plot uses the full count to report truncation.
  hits_for <- function(input, data) {
    pk <- .scroll_atac_read(data)
    if (is.null(pk)) return(character())
    mode <- input$mode %||% mode_choices(data)[[1]]
    if (identical(mode, "Region")) {
      .scroll_atac_region_hits(pk, input$region)
    } else {
      g <- input$gene
      if (!is.null(g) && nzchar(g) && any(!is.na(pk$nearest_gene)))
        pk$feature[!is.na(pk$nearest_gene) & pk$nearest_gene == g] else character()
    }
  }
  peaks_for <- function(input, data)
    utils::head(sort(hits_for(input, data)), .SCROLL_ATAC_MAX_PEAKS)  # bound the client list

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
    # Surface truncation instead of silently dropping hits past the cap.
    n_hits <- length(hits_for(input, data))
    if (n_hits > .SCROLL_ATAC_MAX_PEAKS)
      sub <- sprintf("%s \u00b7 showing first %d of %d matching peaks",
                     sub, .SCROLL_ATAC_MAX_PEAKS, n_hits)

    p <- ggplot2::ggplot(df[order(df$.acc), ], ggplot2::aes(.data$.x, .data$.y, color = .data$.acc)) +
      .scroll_point_layer(size = 0.6, raster = isTRUE(nrow(df) > .scroll_raster_threshold)) +
      .scroll_continuous_scale("viridis", name = "accessibility") +
      ggplot2::labs(title = peak, subtitle = sub, x = paste0(emb, "_1"), y = paste0(emb, "_2")) +
      .scroll_base_theme(axis_text = FALSE) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
    # CSV exports the peaks matching the current search (bounded), with their coords.
    hit_tab <- pk[pk$feature %in% peaks_for(input, data), , drop = FALSE]
    attr(p, "scroll_source") <- as.data.frame(hit_tab); p
  }

  list(.scroll_gated_panel("peaks", "Peaks", "Chromatin accessibility",
    "Find a peak near a gene and colour the embedding by its accessibility.",
    gate, controls, plot, csv = TRUE))
}
