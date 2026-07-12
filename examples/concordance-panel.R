# Concordance panel, written with scroll's declarative builder: run ONE contrast
# in TWO cell populations and scatter the per-gene logFCs. Compare this to a
# hand-written Shiny module - there is no moduleServer / observeEvent / renderPlot
# / download-handler boilerplate, and empty or degenerate selections surface a
# message instead of a blank plot.
#
#   library(scroll); library(ggplot2)
#   source("concordance-panel.R")          # registers the panel
#   scroll_serve("seu_app")
#
# Usage note: pick a Contrast column whose levels vary WITHIN each population
# (e.g. a condition/sample column), and a Population column of clusters/celltypes.
# For the Wilcoxon engine no counts are needed; Pseudobulk needs a build with
# counts = TRUE.

`%||%` <- function(a, b) if (is.null(a) || !length(a)) b else a

# The plot: join two DE tables by gene and scatter logFC vs logFC, coloured by
# concordance, with a y=x line, a Pearson r, and optional gene labels.
delta_de_map <- function(de_a, de_b, side1, side2, xlab = "logFC - A", ylab = "logFC - B",
                         padj = 0.05, lfc = 0, genes = character(0),
                         colors = c(up1 = "#92c5de", up2 = "#ca0020", ns = "#bdbdbd")) {
  a <- de_a[, c("gene", "logFC", "p_val_adj")]; names(a) <- c("gene", "lfc_a", "padj_a")
  b <- de_b[, c("gene", "logFC", "p_val_adj")]; names(b) <- c("gene", "lfc_b", "padj_b")
  df <- merge(a, b, by = "gene")
  df <- df[stats::complete.cases(df[, c("lfc_a", "lfc_b")]), , drop = FALSE]
  if (nrow(df) < 3) stop("Too few shared genes to compare (", nrow(df), ").")
  both <- (df$padj_a < padj & abs(df$lfc_a) >= lfc) & (df$padj_b < padj & abs(df$lfc_b) >= lfc)
  up1 <- paste("Up in", side1); up2 <- paste("Up in", side2)
  df$reg <- ifelse(both & df$lfc_a > 0 & df$lfc_b > 0, up1,
             ifelse(both & df$lfc_a < 0 & df$lfc_b < 0, up2, "NS"))
  df$reg <- factor(df$reg, levels = c(up1, "NS", up2))
  r <- stats::cor(df$lfc_a, df$lfc_b, use = "complete.obs")
  pv <- stats::cor.test(df$lfc_a, df$lfc_b)$p.value
  p_lab <- if (pv < 2.2e-16) "p < 2.2e-16" else paste0("p = ", signif(pv, 3))
  pal <- stats::setNames(c(colors[["up1"]], colors[["ns"]], colors[["up2"]]), c(up1, "NS", up2))
  siz <- stats::setNames(c(3, 1, 3), c(up1, "NS", up2))
  lab <- if (length(genes))
    ggrepel::geom_text_repel(data = df[df$gene %in% genes, ], ggplot2::aes(label = gene),
                             min.segment.length = 0, size = 3.5)
  ggplot2::ggplot(df[order(df$reg != "NS"), ], ggplot2::aes(lfc_a, lfc_b)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted") +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::geom_point(ggplot2::aes(fill = reg, size = reg), shape = 21, colour = "grey30", stroke = 0.2) +
    ggplot2::scale_fill_manual(values = pal) + ggplot2::scale_size_manual(values = siz) + lab +
    ggplot2::annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -1.6,
                      label = paste0("r = ", round(r, 3)), fontface = "italic", size = 5) +
    ggplot2::annotate("text", x = Inf, y = -Inf, hjust = 1.06, vjust = -0.5,
                      label = p_lab, fontface = "italic", size = 5) +
    ggplot2::coord_fixed() + ggplot2::labs(x = xlab, y = ylab) + ggplot2::theme_minimal() +
    ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA, colour = "black", linewidth = 0.5),
                   legend.position = "none", aspect.ratio = 1, panel.grid = ggplot2::element_blank())
}

register_plot_panel("concordance",
  label = "Concordance", title = "logFC concordance across populations",
  desc  = "Run one contrast in two cell populations and compare per-gene logFC.",
  after = "de",
  controls = list(
    scroll_input_choice("engine", "DE engine", c("Wilcoxon", "Pseudobulk")),
    scroll_input_column("group",  "Contrast column", "categorical"),
    scroll_input_levels("ident1", "Group 1", from = "group", required = TRUE),
    scroll_input_levels("ident2", "vs.",     from = "group", none = TRUE),
    scroll_input_column("popcol", "Population column", "categorical", selected = 2),
    scroll_input_levels("popA",   "Population A (x)", from = "popcol", required = TRUE),
    scroll_input_levels("popB",   "Population B (y)", from = "popcol", required = TRUE, selected = 2),
    scroll_input_numeric("padj", "Adj. p cutoff", 0.05, 0, 1, 0.01),
    scroll_input_slider("lfc",   "|logFC| cutoff", 0, 3, 0, 0.1),
    scroll_input_text("genes",   "Label genes (comma/space)")),
  plot = function(cells, input, data) {
    assay <- data$manifest$default_assay
    pick  <- function(v) cells[as.character(cells[[input$popcol]]) %in% v, , drop = FALSE]
    id2   <- if (length(input$ident2)) input$ident2 else NULL
    run <- function(sub, tag) {
      n1 <- sum(as.character(sub[[input$group]]) %in% input$ident1, na.rm = TRUE)
      if (n1 < 3)
        stop(sprintf("Population %s has only %d cells in Group 1 (%s) - pick a contrast that varies within it.",
                     tag, n1, paste(input$ident1, collapse = "/")))
      if (identical(input$engine, "Pseudobulk")) {
        if (!isTRUE(data$manifest$has_counts))
          stop("Pseudobulk needs a counts store - rebuild with scroll_build(..., counts = TRUE).")
        scroll_pseudobulk_de(data, assay, aggregate_cols = input$group,
                             ident1 = input$ident1, ident2 = id2, replicate_col = NULL, cells = sub)
      } else {
        scroll_de(data, assay, input$group, input$ident1, id2, min_pct = 0, cells = sub)
      }
    }
    genes <- strsplit(input$genes %||% "", "[,[:space:]]+")[[1]]; genes <- genes[nzchar(genes)]
    delta_de_map(run(pick(input$popA), "A"), run(pick(input$popB), "B"),
                 side1 = paste(input$ident1, collapse = "/"),
                 side2 = if (length(input$ident2)) paste(input$ident2, collapse = "/") else "rest",
                 xlab  = paste0("logFC - ", paste(input$popA, collapse = "/")),
                 ylab  = paste0("logFC - ", paste(input$popB, collapse = "/")),
                 padj  = input$padj %||% 0.05, lfc = input$lfc %||% 0, genes = genes)
  })
