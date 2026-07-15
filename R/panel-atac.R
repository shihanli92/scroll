# The scATAC Peaks panel + its runtime accessor. Build-time peak parsing + bake
# logic lives in R/atac.R; this reads the baked peaks table and renders the
# gene/region-searchable accessibility panel.

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
