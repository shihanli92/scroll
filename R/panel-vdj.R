# The VDJ (TCR/BCR) repertoire panels + their runtime accessors. The build-time
# spec + bake logic lives in R/vdj.R; these read the baked repertoire/ store and
# render the four gated panels (clone overview, gene usage, CDR3 length, diversity).

# ---- runtime: read the baked store ------------------------------------------

.scroll_vdj_read <- function(data, f)
  .scroll_read_asset(file.path(data$dir, "repertoire", f), "parquet")
# Distinct non-NA values of a vector.
.scroll_ndistinct <- function(v) length(unique(v[!is.na(v)]))

# Categorical rep_cells columns that are candidates for grouping or as a clone id
# (excludes the segment genes, CDR3 lengths, and the clone-size count).
.scroll_vdj_cat_cols <- function(rc, segments) {
  cand <- setdiff(names(rc), c("clone_id", "clone_count", "cdr3_combined",
                               grep("_len$", names(rc), value = TRUE), unlist(segments)))
  cand[vapply(cand, function(c) is.character(rc[[c]]) && any(!is.na(rc[[c]])), logical(1))]
}

# Split the categorical columns by cardinality, relative to the baked clone_id: a
# column near clone-scale cardinality is treated as a clone id (an alternate clonal
# definition), everything lower-cardinality is a grouping/colour axis. This keeps a
# 50k-level clone column out of the group-by menus and out of the colour scales.
.scroll_vdj_high_card <- function(rc, col, ref) .scroll_ndistinct(rc[[col]]) >= 0.5 * ref

# categorical rep_cells columns available to split/colour by (low-cardinality only —
# a high-cardinality clone column would mean tens of thousands of facets/colours).
.scroll_vdj_split_cols <- function(data) {
  rc <- .scroll_vdj_read(data, "rep_cells.parquet"); if (is.null(rc)) return(character())
  ref <- if ("clone_id" %in% names(rc)) .scroll_ndistinct(rc$clone_id) else nrow(rc)
  Filter(function(c) !.scroll_vdj_high_card(rc, c, ref),
         .scroll_vdj_cat_cols(rc, data$manifest$vdj$segments))
}

# Columns that can serve as a clone id: the baked `clone_id` plus any carried
# categorical column of comparable (high) cardinality — an alternate clonal
# definition, e.g. a nucleotide-vs-amino-acid CDR3 clonotype. Low-cardinality
# grouping columns (donor, antigen, cell type) are excluded.
.scroll_vdj_clone_cols <- function(data) {
  rc <- .scroll_vdj_read(data, "rep_cells.parquet")
  if (is.null(rc) || !"clone_id" %in% names(rc)) return("clone_id")
  ref <- .scroll_ndistinct(rc$clone_id)
  alt <- Filter(function(c) .scroll_vdj_high_card(rc, c, ref),
                .scroll_vdj_cat_cols(rc, data$manifest$vdj$segments))
  unique(c("clone_id", alt))
}

# The distinct levels of a group column named by the control `group_id` (for a group
# filter's choices). Falls back to the baked `group` column.
.scroll_vdj_col_levels <- function(input, data, group_id = "group") {
  rc <- .scroll_vdj_read(data, "rep_cells.parquet"); if (is.null(rc)) return(character())
  g <- if (!is.null(input[[group_id]]) && input[[group_id]] %in% names(rc)) input[[group_id]] else "group"
  if (!g %in% names(rc)) return(character())
  sort(unique(as.character(rc[[g]][!is.na(rc[[g]])])))
}

# The active group levels = the group column's levels narrowed by the group filter
# (`levels_id`). Shared by the colour control (to render one picker per level) and the
# plot (to build the colour scale), so the two stay in lock-step.
.scroll_vdj_level_set <- function(input, data, group_id = "group", levels_id = "group_levels") {
  lv <- .scroll_vdj_col_levels(input, data, group_id)
  f <- input[[levels_id]]
  if (length(f)) lv <- lv[lv %in% f]
  lv
}

# A palette + manual-per-group colour control for the VDJ panels, keyed to the active
# group levels. The Manual pickers are rendered dynamically (needs the builder to pass
# `output` to the bind), so they track the chosen group column + group filter.
.scroll_vdj_colour_control <- function(group_id = "group", levels_id = "group_levels") {
  scroll_input_custom("palette", "Colour",
    ui = function(ns, data)
      tagList(selectInput(ns("palette"), "Colour", .scroll_cat_palettes()),
              uiOutput(ns("palette_manual"))),
    bind = function(input, session, data, output) {
      ns <- session$ns
      output$palette_manual <- renderUI({
        if (!identical(input$palette, "Manual")) return(NULL)
        .scroll_manual_ui(ns, .scroll_vdj_level_set(input, data, group_id, levels_id))
      })
    })
}

# Named colour vector for the active group levels: the chosen palette, or the manual
# per-group pickers when palette == "Manual" (falling back until they populate).
.scroll_vdj_colours <- function(input, levels) {
  pal <- input$palette %||% "Tableau 10"
  if (identical(pal, "Manual")) {
    m <- .scroll_manual_colors(input, levels)
    if (!is.null(m) && all(levels %in% names(m))) return(m)
    pal <- "Tableau 10"
  }
  .scroll_discrete_colors(levels, pal)
}

# ---- the four panels --------------------------------------------------------

# Assembled as built-in panels gated on manifest$vdj (see .scroll_assemble_panels).
.scroll_vdj_panels <- function() {
  gate <- function(m) !is.null(m$vdj)
  mk <- function(id, label, title, desc, controls, plot, csv = TRUE)
    .scroll_gated_panel(id, label, title, desc, gate, controls, plot, csv = csv)

  clone_overview <- mk("clone_overview", "Clone overview", "Clonal expansion overview",
    "Rank-abundance of clone sizes and expansion-category composition (clone-level).",
    list(scroll_input_choice("view", "View", c("Rank-abundance", "Expansion composition")),
         # which column defines a clone (the baked clone_id, or an alternate clonal
         # definition carried in the store, e.g. a nucleotide vs amino-acid clonotype)
         scroll_input_choice("clone_col", "Clone ID",
           choices = function(data) .scroll_vdj_clone_cols(data), widget = "select"),
         scroll_input_choice("group", "Group by",
           choices = function(data) .scroll_vdj_split_cols(data), widget = "select"),
         # restrict to a subset of the Group-by column's levels (empty = all)
         scroll_input_levels("group_levels", "Groups", none = TRUE, watch = "group",
           choices = function(input, data) .scroll_vdj_col_levels(input, data, "group")),
         # palette + manual per-group colour picker (rank-abundance colours by group)
         scroll_show_when(.scroll_vdj_colour_control("group", "group_levels"),
                          control = "view", equals = "Rank-abundance")),
    function(cells, input, data) {
      rc <- .scroll_vdj_read(data, "rep_cells.parquet")
      if (is.null(rc) || !nrow(rc)) stop("No repertoire store; rebuild with vdj =.")
      # Clone id column drives what counts as a clone; Group-by drives the split/facet.
      cid <- if (!is.null(input$clone_col) && input$clone_col %in% names(rc)) input$clone_col else "clone_id"
      grp <- if (!is.null(input$group) && input$group %in% names(rc) &&
                 any(!is.na(rc[[input$group]]))) input$group else "group"
      rc <- rc[!is.na(rc[[cid]]) & !is.na(rc[[grp]]), , drop = FALSE]
      if (length(input$group_levels))                          # optional group filter
        rc <- rc[as.character(rc[[grp]]) %in% input$group_levels, , drop = FALSE]
      if (!nrow(rc)) stop("No clones for the selected groups.")
      if (identical(input$view, "Rank-abundance")) {
        d <- rc |> dplyr::group_by(.data[[cid]], .data[[grp]]) |>
          dplyr::summarise(count = dplyr::n(), .groups = "drop") |>
          dplyr::arrange(.data[[grp]], dplyr::desc(.data$count)) |>
          dplyr::group_by(.data[[grp]]) |>
          dplyr::mutate(rank = dplyr::row_number()) |> dplyr::ungroup()
        lv <- .scroll_vdj_level_set(input, data, "group", "group_levels")
        p <- ggplot2::ggplot(d, ggplot2::aes(.data$rank, .data$count, colour = .data[[grp]])) +
          ggplot2::geom_point(size = 0.5) + ggplot2::scale_y_log10() +
          ggplot2::scale_colour_manual(values = .scroll_vdj_colours(input, lv), name = grp) +
          ggplot2::facet_wrap(stats::as.formula(paste0("~`", grp, "`")), scales = "free_x") +
          ggplot2::labs(x = "Clone rank", y = "Clone size (cells, log10)", colour = grp,
                        title = "Clone rank-abundance") +
          ggplot2::theme_minimal() +
          ggplot2::theme(legend.position = "none", panel.grid.minor = ggplot2::element_blank())
        attr(p, "scroll_source") <- as.data.frame(d); p
      } else {
        # Clone size + each clone's dominant group. Vectorized: one grouped (clone x
        # group) count, then take the top-count group per clone via order + de-dup --
        # replaces a table() call per clone, which cost seconds over tens of thousands
        # of clones (52k here). Equivalent result, ~8x faster.
        cg  <- dplyr::count(rc, .data[[cid]], .data[[grp]], name = "n")
        cg  <- cg[order(cg[[cid]], -cg$n), , drop = FALSE]
        dom <- cg[!duplicated(cg[[cid]]), c(cid, grp), drop = FALSE]   # dominant group
        cl  <- dplyr::left_join(dplyr::count(rc, .data[[cid]], name = "count"),
                                dom, by = cid)
        names(cl)[names(cl) == grp] <- "g"
        cl$expansion <- cut(cl$count, c(0, 1, 4, 19, 99, Inf),
          labels = c("Single", "Small", "Medium", "Large", "Hyperexpanded"))
        p <- ggplot2::ggplot(cl, ggplot2::aes(.data$g, fill = .data$expansion)) +
          ggplot2::geom_bar(position = "fill", colour = "white", linewidth = 0.2) +
          ggplot2::scale_y_continuous(labels = scales::percent) +
          ggplot2::scale_fill_brewer(palette = "YlOrRd") +
          ggplot2::labs(x = grp, y = "Fraction of clones", fill = "Expansion",
                        title = "Expansion-category composition") +
          ggplot2::theme_minimal() +
          ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
        attr(p, "scroll_source") <- as.data.frame(cl); p
      }
    })

  gene_usage <- mk("gene_usage", "V/J gene usage", "TCR/BCR V/J gene usage",
    "Per-group V/J segment frequency, or the chi-square standardized-residual heatmap.",
    list(scroll_input_choice("segment", "Segment",
           choices = function(data) unlist(data$manifest$vdj$segments)),
         scroll_input_choice("view", "View", c("Frequency", "Chi-square residuals")),
         # Frequency can re-group by any baked categorical; the residual heatmap is
         # baked per the build-time group_col, so hide this control in that view.
         scroll_show_when(
           scroll_input_choice("group", "Group by",
             choices = function(data) .scroll_vdj_split_cols(data), widget = "select"),
           control = "view", equals = "Frequency"),
         scroll_show_when(
           scroll_input_levels("group_levels", "Groups", none = TRUE, watch = "group",
             choices = function(input, data) .scroll_vdj_col_levels(input, data, "group")),
           control = "view", equals = "Frequency"),
         scroll_show_when(.scroll_vdj_colour_control("group", "group_levels"),
                          control = "view", equals = "Frequency")),
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
        p <- ggplot2::ggplot(ch, ggplot2::aes(.data$gene, .data$group, fill = .data$residual)) +
          ggplot2::geom_tile(colour = "grey90") +
          ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                        midpoint = 0, limits = c(-3, 3)) +
          ggplot2::labs(x = seg, y = NULL, fill = "std\nresidual",
                        title = paste(seg, "usage vs group")) +
          ggplot2::theme_minimal() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, size = 8),
                         panel.grid = ggplot2::element_blank())
        attr(p, "scroll_source") <- as.data.frame(ch); p
      } else {
        grp <- if (!is.null(input$group) && input$group %in% names(d) &&
                   any(!is.na(d[[input$group]]))) input$group else "group"
        d <- d[!is.na(d[[grp]]), , drop = FALSE]
        if (length(input$group_levels))
          d <- d[as.character(d[[grp]]) %in% input$group_levels, , drop = FALSE]
        if (!nrow(d)) stop("No cells for the selected groups.")
        f <- d |> dplyr::group_by(.data$gene, .data[[grp]]) |>
          dplyr::summarise(count = dplyr::n(), .groups = "drop") |>
          dplyr::group_by(.data[[grp]]) |>
          dplyr::mutate(freq = .data$count / sum(.data$count)) |> dplyr::ungroup()
        lv <- .scroll_vdj_level_set(input, data, "group", "group_levels")
        p <- ggplot2::ggplot(f, ggplot2::aes(.data$gene, .data$freq)) +
          ggplot2::geom_line(ggplot2::aes(group = .data$gene), colour = "grey70", linewidth = 0.3) +
          ggplot2::geom_point(ggplot2::aes(fill = .data[[grp]]), shape = 21, size = 2.5) +
          ggplot2::scale_fill_manual(values = .scroll_vdj_colours(input, lv), name = grp) +
          ggplot2::labs(x = seg, y = "Frequency (within group)", fill = grp,
                        title = paste(seg, "usage")) +
          ggplot2::theme_minimal() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, size = 7),
                         panel.grid.major.x = ggplot2::element_blank())
        attr(p, "scroll_source") <- as.data.frame(f); p
      }
    })

  cdr3_length <- mk("cdr3_length", "CDR3 length", "CDR3 length distribution",
    "CDR3 amino-acid length per group (clone-deduplicated).",
    list(scroll_input_choice("chain", "Chain",
           choices = function(data) c(unlist(data$manifest$vdj$cdr3_chains), "Combined")),
         scroll_input_choice("style", "Style", c("Histogram", "Density")),
         scroll_input_choice("clone_col", "Clone ID",
           choices = function(data) .scroll_vdj_clone_cols(data), widget = "select"),
         scroll_input_choice("colorby", "Colour by",
           choices = function(data) .scroll_vdj_split_cols(data)),
         scroll_input_levels("group_levels", "Groups", none = TRUE, watch = "colorby",
           choices = function(input, data) .scroll_vdj_col_levels(input, data, "colorby")),
         .scroll_vdj_colour_control("colorby", "group_levels")),
    function(cells, input, data) {
      rc <- .scroll_vdj_read(data, "rep_cells.parquet")
      lc <- if (identical(input$chain, "Combined")) "cdr3_combined" else paste0("cdr3_", input$chain, "_len")
      if (is.null(rc) || !lc %in% names(rc)) stop("CDR3 length for '", input$chain, "' not baked.")
      col <- if (!is.null(input$colorby) && input$colorby %in% names(rc)) input$colorby else "group"
      cid <- if (!is.null(input$clone_col) && input$clone_col %in% names(rc)) input$clone_col else "clone_id"
      d <- rc[!is.na(rc[[lc]]) & !is.na(rc[[col]]), , drop = FALSE]
      d$len <- d[[lc]]; d$col <- as.character(d[[col]])
      if (length(input$group_levels)) d <- d[d$col %in% input$group_levels, , drop = FALSE]
      d <- d[!duplicated(paste(d[[cid]], d$col)), , drop = FALSE]     # one length per clone
      if (!nrow(d)) stop("No clones with a length for this selection.")
      lv <- .scroll_vdj_level_set(input, data, "colorby", "group_levels")
      cols <- .scroll_vdj_colours(input, lv)
      if (identical(input$style, "Density")) {
        p <- ggplot2::ggplot(d, ggplot2::aes(.data$len, colour = .data$col)) +
          ggplot2::geom_density(adjust = 2, linewidth = 1) +
          ggplot2::scale_colour_manual(values = cols, name = col) +
          ggplot2::labs(x = "CDR3 length", y = "Density", colour = col,
                        title = paste(input$chain, "CDR3 length")) +
          ggplot2::theme_minimal() + ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
        attr(p, "scroll_source") <- data.frame(length = d$len, group = d$col,
                                               clone = d[[cid]]); p
      } else {
        h <- d |> dplyr::group_by(.data$len, .data$col) |>
          dplyr::summarise(count = dplyr::n(), .groups = "drop") |>
          dplyr::group_by(.data$col) |>
          dplyr::mutate(freq = .data$count / sum(.data$count)) |> dplyr::ungroup()
        p <- ggplot2::ggplot(h, ggplot2::aes(.data$len, .data$freq, fill = .data$col)) +
          ggplot2::geom_col(position = ggplot2::position_dodge2(preserve = "single"),
                            colour = "black", linewidth = 0.2) +
          ggplot2::scale_fill_manual(values = cols, name = col) +
          ggplot2::labs(x = "CDR3 length", y = "Frequency (within group)", fill = col,
                        title = paste(input$chain, "CDR3 length")) +
          ggplot2::theme_minimal() + ggplot2::theme(panel.grid = ggplot2::element_blank())
        attr(p, "scroll_source") <- as.data.frame(h); p
      }
    })

  diversity <- mk("diversity", "Diversity", "Repertoire diversity",
    "Shannon / Simpson / clonality / Gini / Gini-Simpson per group, and tissue correlation.",
    list(scroll_input_choice("view", "View", c("Diversity metric", "Tissue correlation")),
         # Diversity is recomputed at runtime from the clone table, so it can group by
         # any baked categorical column (not just the build-time group/cluster) and
         # count clones under any clone definition.
         scroll_show_when(
           scroll_input_choice("clone_col", "Clone ID",
             choices = function(data) .scroll_vdj_clone_cols(data), widget = "select"),
           control = "view", equals = "Diversity metric"),
         scroll_show_when(
           scroll_input_choice("by", "Group by",
             choices = function(data) .scroll_vdj_split_cols(data), widget = "select"),
           control = "view", equals = "Diversity metric"),
         scroll_show_when(
           scroll_input_levels("group_levels", "Groups", none = TRUE, watch = "by",
             choices = function(input, data) .scroll_vdj_col_levels(input, data, "by")),
           control = "view", equals = "Diversity metric"),
         scroll_show_when(
           scroll_input_choice("metric", "Metric",
             c("shannon", "simpson", "clonality", "gini", "top_clone_prop", "paired_rate")),
           control = "view", equals = "Diversity metric"),
         scroll_show_when(.scroll_vdj_colour_control("by", "group_levels"),
                          control = "view", equals = "Diversity metric")),
    function(cells, input, data) {
      if (identical(input$view, "Tissue correlation")) {
        tc <- .scroll_vdj_read(data, "tissue_corr.parquet")
        if (is.null(tc) || !nrow(tc)) stop("Tissue correlation not available (needs loc_col + broad_col).")
        lev <- unique(c(tc$g1, tc$g2))
        tc$g1 <- factor(tc$g1, levels = lev); tc$g2 <- factor(tc$g2, levels = rev(lev))
        p <- ggplot2::ggplot(tc, ggplot2::aes(.data$g1, .data$g2, fill = .data$cor)) +
          ggplot2::geom_tile(colour = "white") +
          ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data$cor)), size = 3) +
          ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                        midpoint = 0, limits = c(-1, 1)) +
          ggplot2::labs(x = NULL, y = NULL, fill = "Pearson r", title = "Clone-frequency correlation") +
          ggplot2::coord_fixed() + ggplot2::theme_minimal() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                         panel.grid = ggplot2::element_blank())
        attr(p, "scroll_source") <- as.data.frame(tc); return(p)
      }
      rc <- .scroll_vdj_read(data, "rep_cells.parquet")
      if (is.null(rc) || !nrow(rc)) stop("No repertoire store; rebuild with vdj =.")
      by_col <- if (!is.null(input$by) && input$by %in% names(rc)) input$by else "group"
      cid <- if (!is.null(input$clone_col) && input$clone_col %in% names(rc)) input$clone_col else "clone_id"
      if (length(input$group_levels))                               # optional group filter
        rc <- rc[!is.na(rc[[by_col]]) & as.character(rc[[by_col]]) %in% input$group_levels, , drop = FALSE]
      len_cols <- paste0("cdr3_", unlist(data$manifest$vdj$cdr3_chains), "_len")
      len_cols <- len_cols[len_cols %in% names(rc)]
      d <- .scroll_vdj_diversity(rc, cid, by_col, len_cols)          # per-level metrics
      metric <- input$metric %||% "shannon"
      if (is.null(d) || !nrow(d) || !metric %in% names(d))
        stop("No '", metric, "' for '", by_col, "'.")
      if (all(is.na(d[[metric]])))
        stop("'", metric, "' is undefined for this dataset (needs two CDR3 chains).")
      d$value <- d[[metric]]; d$level <- factor(d$level, levels = d$level[order(d$value)])
      lv <- .scroll_vdj_level_set(input, data, "by", "group_levels")
      p <- ggplot2::ggplot(d, ggplot2::aes(.data$level, .data$value, fill = .data$level)) +
        ggplot2::geom_col(colour = "black", linewidth = 0.3) +
        ggplot2::scale_fill_manual(values = .scroll_vdj_colours(input, lv)) +
        ggplot2::labs(x = by_col, y = metric, title = paste(metric, "by", by_col)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none",
                       axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
                       panel.grid.major.x = ggplot2::element_blank())
      attr(p, "scroll_source") <- as.data.frame(d); p
    })

  list(clone_overview, gene_usage, cdr3_length, diversity)
}
