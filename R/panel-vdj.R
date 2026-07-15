# The VDJ (TCR/BCR) repertoire panels + their runtime accessors. The build-time
# spec + bake logic lives in R/vdj.R; these read the baked repertoire/ store and
# render the four gated panels (clone overview, gene usage, CDR3 length, diversity).

# ---- runtime: read the baked store ------------------------------------------

.scroll_vdj_read <- function(data, f)
  .scroll_read_asset(file.path(data$dir, "repertoire", f), "parquet")
# categorical rep_cells columns available to split/colour by
.scroll_vdj_split_cols <- function(data) {
  rc <- .scroll_vdj_read(data, "rep_cells.parquet"); if (is.null(rc)) return(character())
  cand <- setdiff(names(rc), c("clone_id", "clone_count", grep("_len$", names(rc), value = TRUE),
                               "cdr3_combined", names(data$manifest$vdj$segments)))
  cand <- setdiff(cand, unlist(data$manifest$vdj$segments))
  cand[vapply(cand, function(c) is.character(rc[[c]]) && any(!is.na(rc[[c]])), logical(1))]
}

# ---- the four panels --------------------------------------------------------

# Assembled as built-in panels gated on manifest$vdj (see .scroll_assemble_panels).
.scroll_vdj_panels <- function() {
  gate <- function(m) !is.null(m$vdj)
  mk <- function(id, label, title, desc, controls, plot)
    .scroll_gated_panel(id, label, title, desc, gate, controls, plot)

  clone_overview <- mk("clone_overview", "Clone overview", "Clonal expansion overview",
    "Rank-abundance of clone sizes and expansion-category composition (clone-level).",
    list(scroll_input_choice("view", "View", c("Rank-abundance", "Expansion composition")),
         scroll_input_choice("group", "Group by",
           choices = function(data) intersect(c("group", "tissue", "antigen"),
                                              names(.scroll_vdj_read(data, "rep_cells.parquet"))))),
    function(cells, input, data) {
      rc <- .scroll_vdj_read(data, "rep_cells.parquet")
      if (is.null(rc) || !nrow(rc)) stop("No repertoire store; rebuild with vdj =.")
      splitc <- if (!is.null(rc$antigen) && any(!is.na(rc$antigen))) "antigen" else "group"
      if (identical(input$view, "Rank-abundance")) {
        d <- rc |> dplyr::group_by(.data$clone_id, .data[[splitc]]) |>
          dplyr::summarise(count = dplyr::n(), .groups = "drop") |>
          dplyr::arrange(.data[[splitc]], dplyr::desc(.data$count)) |>
          dplyr::group_by(.data[[splitc]]) |>
          dplyr::mutate(rank = dplyr::row_number()) |> dplyr::ungroup()
        ggplot2::ggplot(d, ggplot2::aes(.data$rank, .data$count, colour = .data[[splitc]])) +
          ggplot2::geom_point(size = 0.5) + ggplot2::scale_y_log10() +
          ggplot2::facet_wrap(stats::as.formula(paste0("~`", splitc, "`")), scales = "free_x") +
          ggplot2::labs(x = "Clone rank", y = "Clone size (cells, log10)", colour = splitc,
                        title = "Clone rank-abundance") +
          ggplot2::theme_minimal() +
          ggplot2::theme(legend.position = "none", panel.grid.minor = ggplot2::element_blank())
      } else {
        grp <- if (input$group %in% names(rc) && any(!is.na(rc[[input$group]]))) input$group else "group"
        cl <- rc |> dplyr::group_by(.data$clone_id) |>
          dplyr::summarise(count = dplyr::n(),
            g = names(sort(table(.data[[grp]]), decreasing = TRUE))[1], .groups = "drop") |>
          dplyr::mutate(expansion = cut(.data$count, c(0, 1, 4, 19, 99, Inf),
            labels = c("Single", "Small", "Medium", "Large", "Hyperexpanded")))
        ggplot2::ggplot(cl, ggplot2::aes(.data$g, fill = .data$expansion)) +
          ggplot2::geom_bar(position = "fill", colour = "white", linewidth = 0.2) +
          ggplot2::scale_y_continuous(labels = scales::percent) +
          ggplot2::scale_fill_brewer(palette = "YlOrRd") +
          ggplot2::labs(x = grp, y = "Fraction of clones", fill = "Expansion",
                        title = "Expansion-category composition") +
          ggplot2::theme_minimal() +
          ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
      }
    })

  gene_usage <- mk("gene_usage", "V/J gene usage", "TCR/BCR V/J gene usage",
    "Per-group V/J segment frequency, or the chi-square standardized-residual heatmap.",
    list(scroll_input_choice("segment", "Segment",
           choices = function(data) unlist(data$manifest$vdj$segments)),
         scroll_input_choice("view", "View", c("Frequency", "Chi-square residuals"))),
    function(cells, input, data) {
      rc <- .scroll_vdj_read(data, "rep_cells.parquet"); seg <- input$segment
      if (is.null(rc) || !seg %in% names(rc)) stop("Segment '", seg, "' not baked.")
      d <- rc[!is.na(rc[[seg]]), , drop = FALSE]; d$gene <- as.character(d[[seg]])
      if (identical(input$view, "Chi-square residuals")) {
        ch <- .scroll_vdj_read(data, "chisq.parquet")
        ch <- if (!is.null(ch)) ch[ch$segment == seg, , drop = FALSE] else NULL
        if (is.null(ch) || !nrow(ch)) stop("No chi-square residuals for ", seg, ".")
        ch$residual <- pmax(pmin(ch$residual, 3), -3)
        ord <- names(sort(tapply(abs(ch$residual), ch$gene, max)))
        ch$gene <- factor(ch$gene, levels = ord)
        ggplot2::ggplot(ch, ggplot2::aes(.data$gene, .data$group, fill = .data$residual)) +
          ggplot2::geom_tile(colour = "grey90") +
          ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                        midpoint = 0, limits = c(-3, 3)) +
          ggplot2::labs(x = seg, y = NULL, fill = "std\nresidual",
                        title = paste(seg, "usage vs group")) +
          ggplot2::theme_minimal() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, size = 8),
                         panel.grid = ggplot2::element_blank())
      } else {
        f <- d |> dplyr::group_by(.data$gene, .data$group) |>
          dplyr::summarise(count = dplyr::n(), .groups = "drop") |>
          dplyr::group_by(.data$group) |>
          dplyr::mutate(freq = .data$count / sum(.data$count)) |> dplyr::ungroup()
        ggplot2::ggplot(f, ggplot2::aes(.data$gene, .data$freq)) +
          ggplot2::geom_line(ggplot2::aes(group = .data$gene), colour = "grey70", linewidth = 0.3) +
          ggplot2::geom_point(ggplot2::aes(fill = .data$group), shape = 21, size = 2.5) +
          ggplot2::labs(x = seg, y = "Frequency (within group)", fill = "group",
                        title = paste(seg, "usage")) +
          ggplot2::theme_minimal() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, size = 7),
                         panel.grid.major.x = ggplot2::element_blank())
      }
    })

  cdr3_length <- mk("cdr3_length", "CDR3 length", "CDR3 length distribution",
    "CDR3 amino-acid length per group (clone-deduplicated).",
    list(scroll_input_choice("chain", "Chain",
           choices = function(data) c(unlist(data$manifest$vdj$cdr3_chains), "Combined")),
         scroll_input_choice("style", "Style", c("Histogram", "Density")),
         scroll_input_choice("colorby", "Colour by",
           choices = function(data) .scroll_vdj_split_cols(data))),
    function(cells, input, data) {
      rc <- .scroll_vdj_read(data, "rep_cells.parquet")
      lc <- if (identical(input$chain, "Combined")) "cdr3_combined" else paste0("cdr3_", input$chain, "_len")
      if (is.null(rc) || !lc %in% names(rc)) stop("CDR3 length for '", input$chain, "' not baked.")
      col <- if (!is.null(input$colorby) && input$colorby %in% names(rc)) input$colorby else "group"
      d <- rc[!is.na(rc[[lc]]) & !is.na(rc[[col]]), , drop = FALSE]
      d$len <- d[[lc]]; d$col <- as.character(d[[col]])
      d <- d[!duplicated(paste(d$clone_id, d$col)), , drop = FALSE]
      if (!nrow(d)) stop("No clones with a length for this selection.")
      if (identical(input$style, "Density")) {
        ggplot2::ggplot(d, ggplot2::aes(.data$len, colour = .data$col)) +
          ggplot2::geom_density(adjust = 2, linewidth = 1) +
          ggplot2::labs(x = "CDR3 length", y = "Density", colour = col,
                        title = paste(input$chain, "CDR3 length")) +
          ggplot2::theme_minimal() + ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
      } else {
        h <- d |> dplyr::group_by(.data$len, .data$col) |>
          dplyr::summarise(count = dplyr::n(), .groups = "drop") |>
          dplyr::group_by(.data$col) |>
          dplyr::mutate(freq = .data$count / sum(.data$count)) |> dplyr::ungroup()
        ggplot2::ggplot(h, ggplot2::aes(.data$len, .data$freq, fill = .data$col)) +
          ggplot2::geom_col(position = ggplot2::position_dodge2(preserve = "single"),
                            colour = "black", linewidth = 0.2) +
          ggplot2::labs(x = "CDR3 length", y = "Frequency (within group)", fill = col,
                        title = paste(input$chain, "CDR3 length")) +
          ggplot2::theme_minimal() + ggplot2::theme(panel.grid = ggplot2::element_blank())
      }
    })

  diversity <- mk("diversity", "Diversity", "Repertoire diversity",
    "Shannon / Simpson / clonality / Gini per group (or cluster), and tissue correlation.",
    list(scroll_input_choice("view", "View", c("Diversity metric", "Tissue correlation")),
         scroll_input_choice("by", "Group by",
           choices = function(data) if (isTRUE(data$manifest$vdj$has_cluster)) c("group", "cluster") else "group"),
         scroll_input_choice("metric", "Metric",
           c("shannon", "simpson", "clonality", "gini", "top_clone_prop"))),
    function(cells, input, data) {
      if (identical(input$view, "Tissue correlation")) {
        tc <- .scroll_vdj_read(data, "tissue_corr.parquet")
        if (is.null(tc) || !nrow(tc)) stop("Tissue correlation not available (needs loc_col + broad_col).")
        lev <- unique(c(tc$g1, tc$g2))
        tc$g1 <- factor(tc$g1, levels = lev); tc$g2 <- factor(tc$g2, levels = rev(lev))
        return(ggplot2::ggplot(tc, ggplot2::aes(.data$g1, .data$g2, fill = .data$cor)) +
          ggplot2::geom_tile(colour = "white") +
          ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data$cor)), size = 3) +
          ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                        midpoint = 0, limits = c(-1, 1)) +
          ggplot2::labs(x = NULL, y = NULL, fill = "Pearson r", title = "Clone-frequency correlation") +
          ggplot2::coord_fixed() + ggplot2::theme_minimal() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                         panel.grid = ggplot2::element_blank()))
      }
      dv <- .scroll_vdj_read(data, "diversity.parquet")
      if (is.null(dv) || !nrow(dv)) stop("No diversity table baked.")
      d <- dv[dv$level_type == input$by, , drop = FALSE]
      if (!nrow(d) || !input$metric %in% names(d)) stop("No '", input$metric, "' for '", input$by, "'.")
      d$value <- d[[input$metric]]; d$level <- factor(d$level, levels = d$level[order(d$value)])
      ggplot2::ggplot(d, ggplot2::aes(.data$level, .data$value, fill = .data$value)) +
        ggplot2::geom_col(colour = "black", linewidth = 0.3) +
        ggplot2::scale_fill_viridis_c() +
        ggplot2::labs(x = input$by, y = input$metric, title = paste(input$metric, "by", input$by)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none",
                       axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
                       panel.grid.major.x = ggplot2::element_blank())
    })

  list(clone_overview, gene_usage, cdr3_length, diversity)
}
