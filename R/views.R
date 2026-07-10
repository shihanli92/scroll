# Built-in views. Views are pure functions of data frames (no Shiny), so they
# unit-test without a running app. render_view() is the single dispatcher the
# app calls; it unpacks a context list (cells + whatever payload the view needs
# + toggle state) and routes to the right view. The full *_ui/*_server module
# contract + register_view() is a later phase; keeping render logic pure here is
# what lets that refactor stay mechanical.

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# --- controls / palettes ------------------------------------------------------

# Resolve a control value: state override, else params, else default.
.scroll_opt <- function(params, state, key, default) {
  v <- state[[key]]; if (!is.null(v)) return(v)
  v <- params[[key]]; if (!is.null(v)) return(v)
  default
}

# Colorblind-safe discrete palettes. Colors are assigned deterministically by
# sorted level name so a category (e.g. a cell type) keeps its color across
# every panel and every subset.
.scroll_discrete_palettes <- list(
  "Tableau 10" = c("#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#B07AA1",
                   "#76B7B2", "#EDC948", "#FF9DA7", "#9C755F", "#BAB0AC"),
  "Okabe-Ito"  = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                   "#D55E00", "#CC79A7", "#999999"),
  "Set2"       = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
                   "#FFD92F", "#E5C494", "#B3B3B3"),
  "Set1"       = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                   "#FFD92F", "#A65628", "#F781BF", "#999999"),
  "Dark2"      = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
                   "#E6AB02", "#A6761D", "#666666"),
  "Paired"     = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
                   "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A")
)
.scroll_brewer_seq <- c("Blues", "Reds", "Greens", "Purples", "YlOrRd", "YlGnBu", "OrRd")
.scroll_brewer_div <- c("RdBu", "RdYlBu", "Spectral", "PuOr", "BrBG")
.scroll_continuous_palettes <- c("viridis", "magma", "plasma", "inferno",
                                 "cividis", "turbo", "rocket", "mako",
                                 .scroll_brewer_seq, .scroll_brewer_div,
                                 "grey-purple", "grey-red", "grey-blue")

.scroll_discrete_colors <- function(values, palette = "Tableau 10") {
  pal <- .scroll_discrete_palettes[[palette]] %||% .scroll_discrete_palettes[["Tableau 10"]]
  lv <- sort(unique(as.character(values)))
  stats::setNames(pal[((seq_along(lv) - 1L) %% length(pal)) + 1L], lv)
}

# `limits` clip the scale (values beyond are squished to the end color, not
# dropped) — used by FeaturePlot's quantile caps.
.scroll_continuous_scale <- function(palette, name = NULL, limits = NULL) {
  palette <- palette %||% "viridis"
  # ColorBrewer palettes via scale_color_distiller (sequential low->high dark,
  # diverging reversed so warm = high).
  if (palette %in% .scroll_brewer_seq)
    return(ggplot2::scale_color_distiller(name = name, palette = palette,
             direction = 1, limits = limits, oob = scales::squish))
  if (palette %in% .scroll_brewer_div)
    return(ggplot2::scale_color_distiller(name = name, palette = palette,
             direction = -1, limits = limits, oob = scales::squish))
  vir <- function(opt) ggplot2::scale_color_viridis_c(
    name = name, option = opt, limits = limits, oob = scales::squish)
  grad <- function(hi) ggplot2::scale_color_gradient(
    low = "grey88", high = hi, name = name, limits = limits, oob = scales::squish)
  switch(palette,
    viridis = vir("viridis"), magma = vir("magma"), plasma = vir("plasma"),
    inferno = vir("inferno"), cividis = vir("cividis"), turbo = vir("turbo"),
    rocket = vir("rocket"), mako = vir("mako"),
    `grey-purple` = grad("#3b0f70"), `grey-red` = grad("#b2182b"),
    `grey-blue` = grad("#08519c"),
    vir("viridis")
  )
}

# Color-scale limits from min/max quantile fractions (or NULL = no clipping).
.scroll_expr_limits <- function(expr, clip) {
  if (is.null(clip) || length(clip) != 2) return(NULL)
  if (clip[1] <= 0 && clip[2] >= 1) return(NULL)
  lims <- unname(stats::quantile(expr, c(max(0, clip[1]), min(1, clip[2])), na.rm = TRUE))
  if (!all(is.finite(lims)) || lims[1] >= lims[2]) return(NULL)
  lims
}

# --- shared helpers -----------------------------------------------------------

# Effective embedding / grouping column honour toggles over config.
.scroll_eff_embedding <- function(params, state) state$embedding %||% params$embedding
.scroll_group_col <- function(params, state) state$split_by %||% params$group_by

# Attach an embedding's first two dims as .x / .y. Rows with no coordinates
# (NA) are dropped: a reprocessed-subset embedding only covers its own cells, so
# this cleanly plots just the subset instead of emitting ggplot NA warnings.
.scroll_embedding_xy <- function(cells, embedding) {
  xcol <- sprintf("%s_1", embedding)
  ycol <- sprintf("%s_2", embedding)
  if (!all(c(xcol, ycol) %in% names(cells)))
    stop("Embedding '", embedding, "' not found in cells table.", call. = FALSE)
  cells$.x <- cells[[xcol]]
  cells$.y <- cells[[ycol]]
  cells[!is.na(cells$.x) & !is.na(cells$.y), , drop = FALSE]
}

# Expression vector aligned to a cells table, zero-filled for absent cells.
.scroll_expr_vector <- function(cells, values) {
  expr <- stats::setNames(rep(0, nrow(cells)), cells$cell)
  if (!is.null(values) && nrow(values) > 0) expr[values$cell] <- values$value
  as.numeric(expr[cells$cell])
}

# A scatter point layer that rasterizes on screen for speed but stays a true
# geom_point for exports and tests. `raster = FALSE` (the default) always yields
# a vector geom_point, so direct view calls are unchanged. When `raster = TRUE`
# and scattermore is installed, points are drawn into a raster buffer (~10x
# faster for tens of thousands of points); absent scattermore it falls back to
# geom_point. scattermore takes a pixel `pointsize`, not the mm `size` aesthetic,
# so we approximate it from `size` (interactive-only; exports use vector points).
.scroll_point_layer <- function(mapping = NULL, size = 0.6, alpha = 1,
                                raster = FALSE, ...) {
  if (isTRUE(raster) && requireNamespace("scattermore", quietly = TRUE))
    scattermore::geom_scattermore(mapping = mapping,
                                  pointsize = max(1, round(size * 2)), alpha = alpha, ...)
  else
    ggplot2::geom_point(mapping = mapping, size = size, alpha = alpha, ...)
}

.scroll_base_scatter <- function(df, embedding, legend = TRUE) {
  ggplot2::ggplot(df, ggplot2::aes(x = .data$.x, y = .data$.y)) +
    ggplot2::labs(x = sprintf("%s 1", embedding), y = sprintf("%s 2", embedding)) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", linewidth = 0.7, fill = NA),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      legend.position = if (isTRUE(legend)) "right" else "none"
    )
}

# Optional facet by a categorical toggle column.
.scroll_maybe_facet <- function(p, cells, state) {
  sb <- state$split_by
  if (!is.null(sb) && nzchar(sb %||% "") && sb %in% names(cells))
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste0("~ `", sb, "`")))
  p
}

# Aspect ratio (panel height / width). Set via theme(aspect.ratio) so the panel
# reshapes WITHIN the fixed-size plot canvas (letterboxing) rather than growing
# it — a tall aspect never forces the page to scroll. Left off at 1 so the
# default fills the canvas as before.
.scroll_apply_aspect <- function(p, aspect = 1) {
  if (is.numeric(aspect) && length(aspect) == 1 && abs(aspect - 1) > 1e-6)
    p <- p + ggplot2::theme(aspect.ratio = aspect)
  p
}

# Scatter finish: facet, then aspect.
.scroll_finish_scatter <- function(p, df, state) {
  .scroll_apply_aspect(.scroll_maybe_facet(p, df, state), state$aspect %||% 1)
}

# Per-group centroid labels for a categorical scatter.
.scroll_group_labels <- function(df, col) {
  agg <- stats::aggregate(cbind(.x, .y) ~ .g, data = transform(df, .g = df[[col]]),
                          FUN = stats::median)
  ggplot2::geom_text(
    data = agg, ggplot2::aes(x = .data$.x, y = .data$.y, label = .data$.g),
    inherit.aes = FALSE, size = 4, fontface = "bold"
  )
}

# --- scatter views ------------------------------------------------------------

# Named colors for a categorical column: the chosen palette, with any manual
# per-group overrides layered on top (both keyed by level name).
.scroll_group_colors <- function(values, state) {
  cols <- .scroll_discrete_colors(values, .scroll_opt(list(), state, "palette", "Tableau 10"))
  manual <- state$manual_colors
  if (!is.null(manual) && length(manual)) {
    manual <- unlist(manual)
    keep <- names(manual) %in% names(cols) & nzchar(manual)
    cols[names(manual)[keep]] <- manual[keep]
  }
  cols
}

#' Embedding colored by a metadata column
#'
#' @param cells The globally-loaded cells data.frame.
#' @param params List with `embedding` and `color_by`.
#' @param state Optional toggle state (`embedding`, `split_by`, `show_labels`).
#' @return A ggplot.
#' @export
view_umap_colorby <- function(cells, params, state = list()) {
  embedding <- .scroll_eff_embedding(params, state)
  df <- .scroll_embedding_xy(cells, embedding)
  color_by <- params$color_by
  if (!color_by %in% names(df))
    stop("color_by column '", color_by, "' not in cells table.", call. = FALSE)
  df$.col <- df[[color_by]]
  size <- .scroll_opt(params, state, "point_size", 0.6)
  alpha <- .scroll_opt(params, state, "alpha", 0.85)
  raster <- isTRUE(state$raster)
  base <- .scroll_base_scatter(df, embedding, .scroll_opt(params, state, "legend", TRUE))

  # numeric color-by: a continuous gradient
  if (is.numeric(df$.col)) {
    p <- base +
      .scroll_point_layer(ggplot2::aes(color = .data$.col), size = size, alpha = alpha, raster = raster) +
      .scroll_continuous_scale(.scroll_opt(params, state, "palette", "viridis"), color_by) +
      ggplot2::labs(color = color_by)
    return(.scroll_finish_scatter(p, df, state))
  }

  # categorical: discrete colors, with optional group highlight + manual colors
  df$.col <- as.character(df$.col)
  cols <- .scroll_group_colors(df$.col, state)
  highlight <- intersect(state$highlight %||% character(0), df$.col)

  if (length(highlight)) {                       # selected groups keep color, rest grey
    bg <- df[!(df$.col %in% highlight), , drop = FALSE]
    fg <- df[df$.col %in% highlight, , drop = FALSE]
    p <- base +
      .scroll_point_layer(data = bg, color = "grey85", size = size, alpha = alpha, raster = raster) +
      .scroll_point_layer(ggplot2::aes(color = .data$.col), data = fg, size = size, alpha = alpha, raster = raster) +
      ggplot2::scale_color_manual(values = cols, limits = highlight) +
      ggplot2::labs(color = color_by)
    label_df <- fg
  } else {
    p <- base +
      .scroll_point_layer(ggplot2::aes(color = .data$.col), size = size, alpha = alpha, raster = raster) +
      ggplot2::scale_color_manual(values = cols) +
      ggplot2::labs(color = color_by)
    label_df <- df
  }
  # legend key glyphs, sized independently of the (small) plotted points
  p <- p + ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4)))
  if (isTRUE(.scroll_opt(params, state, "show_labels", FALSE)) && nrow(label_df)) {
    p <- p + .scroll_group_labels(label_df, ".col")
    # on-plot labels stand in for the legend only when the legend is switched off,
    # so the "Legend" toggle still works when cluster labels are shown
    if (!isTRUE(.scroll_opt(params, state, "legend", TRUE)))
      p <- p + ggplot2::guides(color = "none")
  }
  .scroll_finish_scatter(p, df, state)
}

#' Embedding colored by a feature's expression
#'
#' @param cells The globally-loaded cells data.frame.
#' @param params List with `embedding` (and `feature`, for the title).
#' @param values A data.frame(cell, value) in normalized units, or `NULL`.
#' @param state Optional toggle state.
#' @return A ggplot.
#' @export
view_feature_plot <- function(cells, params, values = NULL, state = list()) {
  embedding <- .scroll_eff_embedding(params, state)
  df <- .scroll_embedding_xy(cells, embedding)
  df$.expr <- .scroll_expr_vector(df, values)
  if (isTRUE(.scroll_opt(params, state, "order", TRUE)))
    df <- df[order(df$.expr), , drop = FALSE]   # expressing cells drawn on top
  size <- .scroll_opt(params, state, "point_size", 0.7)
  lims <- .scroll_expr_limits(df$.expr, .scroll_opt(params, state, "clip", NULL))
  p <- .scroll_base_scatter(df, embedding, .scroll_opt(params, state, "legend", TRUE)) +
    .scroll_point_layer(ggplot2::aes(color = .data$.expr), size = size, raster = isTRUE(state$raster)) +
    .scroll_continuous_scale(.scroll_opt(params, state, "palette", "grey-purple"),
                             params$feature %||% "expression", lims)
  .scroll_finish_scatter(p, df, state)
}

# --- aggregate views ----------------------------------------------------------

# The heavy dotplot work — aggregation, optional z-scoring, and hierarchical
# clustering — depends only on data inputs (features, group, scale, cluster), not
# cosmetics. Factored out so the Shiny server computes it once in a data reactive
# rather than on every palette/dot-size change.
.scroll_dotplot_assemble <- function(cells, features, group_by, expr_long,
                                     scale = FALSE, cluster = "off") {
  grp <- stats::setNames(as.character(cells[[group_by]]), cells$cell)
  ng <- table(as.character(cells[[group_by]]))          # cells per group
  groups <- names(ng)

  # complete (feature x group) grid so absent combinations render as empty dots
  grid <- expand.grid(feature = features, group = groups,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  if (!is.null(expr_long) && nrow(expr_long) > 0) {
    e <- expr_long
    e$group <- grp[e$cell]
    e <- e[!is.na(e$group), , drop = FALSE]
    s_sum <- stats::aggregate(value ~ feature + group, e, sum)
    s_pos <- stats::aggregate(value ~ feature + group, e,
                              function(x) sum(x > 0))
    names(s_sum)[3] <- "sum"; names(s_pos)[3] <- "npos"
    agg <- merge(merge(grid, s_sum, all.x = TRUE), s_pos, all.x = TRUE)
  } else {
    agg <- grid; agg$sum <- 0; agg$npos <- 0
  }
  agg$sum[is.na(agg$sum)] <- 0
  agg$npos[is.na(agg$npos)] <- 0
  agg$mean <- agg$sum / as.numeric(ng[agg$group])
  agg$frac <- agg$npos / as.numeric(ng[agg$group])

  if (isTRUE(scale)) {
    agg$mean <- stats::ave(agg$mean, agg$feature, FUN = function(x) {
      s <- stats::sd(x); if (is.na(s) || s == 0) x * 0 else (x - mean(x)) / s
    })
  }
  # matrix of the plotted color statistic, for hierarchical clustering
  M <- matrix(0, length(features), length(groups), dimnames = list(features, groups))
  M[cbind(match(agg$feature, features), match(agg$group, groups))] <- agg$mean

  hr <- if (cluster %in% c("rows", "both") && length(features) > 2) stats::hclust(stats::dist(M))
  hc <- if (cluster %in% c("columns", "both") && length(groups) > 2) stats::hclust(stats::dist(t(M)))
  feature_order <- if (!is.null(hr)) rownames(M)[hr$order] else features
  group_order   <- if (!is.null(hc)) colnames(M)[hc$order] else groups

  agg$feature <- factor(agg$feature, levels = rev(feature_order))
  agg$group   <- factor(agg$group, levels = group_order)
  list(agg = agg, hr = hr, hc = hc, group_by = group_by, scaled = isTRUE(scale))
}

#' Dot plot: mean expression x fraction expressing, across groups
#'
#' @param cells The cells data.frame (its group sizes are the denominators).
#' @param params List with `features`, `group_by`, and optionally `assay`.
#' @param expr_long A data.frame(feature, cell, value) in normalized units for
#'   the requested features (zero-valued cells absent).
#' @param state Optional toggle state (`split_by` overrides `group_by`).
#' @param assembly Optional precomputed result of the internal aggregation +
#'   clustering step; when supplied (as the Shiny app does), `cells`/`expr_long`
#'   are ignored for assembly and only cosmetics are applied.
#' @return A ggplot.
#' @export
view_dotplot <- function(cells, params, expr_long, state = list(), assembly = NULL) {
  if (is.null(assembly)) {
    group_by <- .scroll_group_col(params, state)
    if (is.null(group_by) || !group_by %in% names(cells))
      stop("dotplot needs a valid `group_by` metadata column.", call. = FALSE)
    assembly <- .scroll_dotplot_assemble(
      cells, unlist(params$features), group_by, expr_long,
      scale = isTRUE(.scroll_opt(params, state, "scale", FALSE)),
      cluster = .scroll_opt(params, state, "cluster", "off"))
  }
  agg <- assembly$agg; hr <- assembly$hr; hc <- assembly$hc
  group_by <- assembly$group_by; scaled <- assembly$scaled

  dsize <- .scroll_opt(params, state, "dot_size", c(1, 6))
  pal <- .scroll_opt(params, state, "palette", "magma")
  # Plot shape is controlled by the rendered canvas height (set in the module),
  # not theme(aspect.ratio): the latter letterboxes the main panel and detaches
  # the aplot dendrograms, so it must not be used with the trees.
  p <- ggplot2::ggplot(agg, ggplot2::aes(x = .data$group, y = .data$feature)) +
    ggplot2::geom_point(ggplot2::aes(size = .data$frac, color = .data$mean)) +
    ggplot2::scale_size(range = dsize, limits = c(0, 1), labels = scales::percent,
                        name = "% expressing") +
    .scroll_continuous_scale(pal, if (scaled) "z-score" else "mean expr.") +
    ggplot2::labs(x = group_by, y = NULL) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(color = "black", linewidth = 0.7, fill = NA),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  # aspect.ratio letterboxes the panel, which would detach the aplot trees, so
  # apply it only to the plain (untreed) dot plot; with dendrograms the composite
  # fills the fixed canvas.
  if (is.null(hr) && is.null(hc)) p <- .scroll_apply_aspect(p, state$aspect %||% 1)
  .scroll_dotplot_trees(p, hr, hc)
}

# Attach row/column dendrograms (aplot + ggtree) when clustering is on and the
# packages are available; otherwise return the plain (reordered) dot plot.
.scroll_dotplot_trees <- function(p, hr, hc) {
  if (is.null(hr) && is.null(hc)) return(p)
  if (!requireNamespace("ggtree", quietly = TRUE) ||
      !requireNamespace("aplot", quietly = TRUE) ||
      !requireNamespace("ape", quietly = TRUE)) return(p)
  # as.phylo() carries the merge heights as branch lengths, so the tree is drawn
  # ultrametric (an hclust dendrogram): every leaf lands at the same position
  # against the plot, instead of a topology-only cladogram with ragged tips.
  # ggtree emits its own deprecation warnings (aes_(), ...) — suppress that noise.
  suppressWarnings({
    pp <- p
    if (!is.null(hr))
      pp <- aplot::insert_left(pp, ggtree::ggtree(ape::as.phylo(hr)), width = 0.16)
    if (!is.null(hc))
      pp <- aplot::insert_top(pp, ggtree::ggtree(ape::as.phylo(hc)) +
                                ggtree::layout_dendrogram(), height = 0.16)
  })
  pp
}

# Shared theme for the bar/violin panels: black box, no vertical gridlines, a
# faint horizontal guide for reading values, optional legend.
.scroll_box_theme <- function(legend = TRUE) {
  list(
    ggplot2::theme_bw(base_size = 13),
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", linewidth = 0.7, fill = NA),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = if (isTRUE(legend)) "right" else "none")
  )
}

#' Violin: a feature's per-group distribution
#'
#' @param cells The cells data.frame.
#' @param params List with `feature`, `group_by`, optionally `assay`.
#' @param values A data.frame(cell, value) for the feature, or `NULL`.
#' @param state Optional toggle state (`palette`, `jitter`, `legend`;
#'   `split_by` overrides `group_by`).
#' @return A ggplot.
#' @export
view_violin <- function(cells, params, values = NULL, state = list()) {
  group_by <- .scroll_group_col(params, state)
  if (is.null(group_by) || !group_by %in% names(cells))
    stop("violin needs a valid `group_by` metadata column.", call. = FALSE)
  df <- data.frame(cell = cells$cell, group = as.character(cells[[group_by]]),
                   stringsAsFactors = FALSE)
  df$expr <- .scroll_expr_vector(cells, values)
  cols <- .scroll_group_colors(df$group, state)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$group, y = .data$expr,
                                        fill = .data$group)) +
    ggplot2::geom_violin(scale = "width", trim = TRUE, linewidth = 0.3)
  if (isTRUE(.scroll_opt(params, state, "jitter", FALSE)))
    p <- p + ggplot2::geom_jitter(size = 0.2, alpha = 0.3, width = 0.2,
                                  show.legend = FALSE)
  p <- p +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(x = group_by, y = params$feature %||% "expression", fill = group_by) +
    .scroll_box_theme(.scroll_opt(params, state, "legend", FALSE))
  .scroll_apply_aspect(p, state$aspect %||% 1)
}

# Long data.frame of all column pairs (one facet per pair), NA rows dropped.
# Factored out so the Shiny server can build it once in a data reactive rather
# than re-expanding on every cosmetic change.
.scroll_biaxial_df <- function(cells, params) {
  feats <- params$features
  feats <- feats[feats %in% names(cells)]
  feats <- feats[!duplicated(feats)]
  if (length(feats) < 2)
    stop("biaxial needs at least two numeric metadata columns.", call. = FALSE)
  color_by <- params$color_by
  if (is.null(color_by) || !color_by %in% names(cells))
    stop("biaxial needs a valid `color_by` metadata column.", call. = FALSE)
  col_vals <- as.character(cells[[color_by]])
  pairs <- utils::combn(feats, 2, simplify = FALSE)
  df <- do.call(rbind, lapply(pairs, function(p) data.frame(
    .x = as.numeric(cells[[p[1]]]), .y = as.numeric(cells[[p[2]]]), .col = col_vals,
    pair = paste(p[1], "vs", p[2]), stringsAsFactors = FALSE)))
  df[!is.na(df$.x) & !is.na(df$.y), , drop = FALSE]
}

#' Pairwise biaxial scatter of numeric metadata columns
#'
#' All unique pairs of the chosen numeric columns as a faceted grid of scatter
#' plots, coloured by a categorical selection (e.g. hashtag CLR values coloured by
#' demux assignment, or ADT/QC pairs coloured by cell type). Rows missing either
#' coordinate are dropped.
#'
#' @param cells The cells data.frame.
#' @param params List with `features` (>= 2 numeric metadata column names) and
#'   `color_by` (a categorical metadata column).
#' @param state Optional (`palette`, `point_size`, `alpha`, `legend`, `aspect`).
#' @param df Optional precomputed pair data.frame (as the Shiny app supplies); when
#'   given, `cells`/`params` are not re-expanded.
#' @return A ggplot (facet per column pair).
#' @export
view_biaxial <- function(cells, params, state = list(), df = NULL) {
  if (is.null(df)) df <- .scroll_biaxial_df(cells, params)
  cols <- .scroll_group_colors(df$.col, state)
  p <- ggplot2::ggplot(df, ggplot2::aes(.data$.x, .data$.y, color = .data$.col)) +
    .scroll_point_layer(size = state$point_size %||% 0.5, alpha = state$alpha %||% 0.6,
                        raster = isTRUE(state$raster)) +
    ggplot2::facet_wrap(~ pair, scales = "free") +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::labs(x = NULL, y = NULL, color = params$color_by) +
    ggplot2::guides(color = ggplot2::guide_legend(
      override.aes = list(size = 2, alpha = 1))) +
    .scroll_box_theme(.scroll_opt(params, state, "legend", TRUE))
  .scroll_apply_aspect(p, state$aspect %||% 1)
}

#' Stacked composition of one categorical within another
#'
#' @param cells The cells data.frame.
#' @param params List with `group_by` (x axis) and `fill_by` (composition).
#' @param state Optional toggle state (`palette`, `normalize`, `legend`).
#' @return A ggplot.
#' @export
view_proportions <- function(cells, params, state = list()) {
  x <- params$group_by
  fill <- params$fill_by
  for (col in c(x, fill)) if (is.null(col) || !col %in% names(cells))
    stop("proportions needs `group_by` and `fill_by` metadata columns.",
         call. = FALSE)
  tab <- as.data.frame(table(x = as.character(cells[[x]]),
                             fill = as.character(cells[[fill]])),
                       stringsAsFactors = FALSE)
  normalize <- isTRUE(.scroll_opt(params, state, "normalize", TRUE))
  # no lower expansion so the bars sit flush on the x-axis; small headroom on top
  yexp <- ggplot2::expansion(mult = c(0, 0.05))
  if (normalize) {
    totals <- stats::aggregate(Freq ~ x, tab, sum)
    tab <- merge(tab, totals, by = "x", suffixes = c("", ".total"))
    tab$y <- ifelse(tab$Freq.total > 0, tab$Freq / tab$Freq.total, 0)
    yscale <- ggplot2::scale_y_continuous(labels = scales::percent, expand = yexp)
    ylab <- "composition"
  } else {
    tab$y <- tab$Freq
    yscale <- ggplot2::scale_y_continuous(expand = yexp)
    ylab <- "cells"
  }
  cols <- .scroll_group_colors(tab$fill, state)
  p <- ggplot2::ggplot(tab, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$fill)) +
    ggplot2::geom_col(width = 0.8, color = "black", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = cols) +
    yscale +
    ggplot2::labs(x = x, y = ylab, fill = fill)
  p <- p + .scroll_box_theme(.scroll_opt(params, state, "legend", TRUE))
  .scroll_apply_aspect(p, state$aspect %||% 1)
}

#' Differential-expression table for a contrast
#'
#' @param de_data A data.frame (a contrast table read from `de/`), or `NULL`.
#' @param params List with `contrast` and optionally `top` (rows to keep).
#' @return A data.frame ready for tabular rendering.
#' @export
view_de_table <- function(de_data, params) {
  if (is.null(de_data))
    return(data.frame(message = sprintf(
      "No DE table found for contrast '%s' in de/.", params$contrast %||% "?")))
  n <- params$top %||% 50
  utils::head(de_data, n)
}

#' Volcano plot of a differential-expression result
#'
#' @param de A data.frame from [scroll_de()] (needs `gene`, `logFC`,
#'   `p_val_adj`).
#' @param params List with `lfc` (fold-change cutoff), `padj` (adjusted-p
#'   cutoff), and `label_n` (number of genes to label).
#' @param state Optional toggle state (`aspect`).
#' @return A ggplot.
#' @export
view_volcano <- function(de, params = list(), state = list()) {
  if (is.null(de) || nrow(de) == 0)
    return(ggplot2::ggplot() +
             ggplot2::annotate("text", 0, 0, label = "No DE results.") +
             ggplot2::theme_void())
  lfc <- params$lfc %||% 1
  pcut <- params$padj %||% 0.05
  d <- de
  d$neglog <- -log10(pmax(d$p_val_adj, 1e-300))     # avoid Inf when padj underflows to 0
  d$sig <- factor("ns", levels = c("down", "ns", "up"))
  d$sig[d$p_val_adj < pcut & d$logFC >= lfc] <- "up"
  d$sig[d$p_val_adj < pcut & d$logFC <= -lfc] <- "down"
  cols <- c(down = "#2563A8", ns = "grey78", up = "#C4453B")

  p <- ggplot2::ggplot(d, ggplot2::aes(.data$logFC, .data$neglog, color = .data$sig)) +
    ggplot2::geom_point(size = 1, alpha = 0.75) +
    ggplot2::scale_color_manual(values = cols, guide = "none") +
    ggplot2::geom_vline(xintercept = c(-lfc, lfc), linetype = "dashed", color = "grey60") +
    ggplot2::geom_hline(yintercept = -log10(pcut), linetype = "dashed", color = "grey60") +
    ggplot2::labs(x = "logFC", y = "-log10 adjusted p") +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(color = "black", linewidth = 0.7, fill = NA))

  n <- params$label_n %||% 15
  lab <- d[d$sig != "ns", , drop = FALSE]
  lab <- utils::head(lab[order(-lab$neglog), , drop = FALSE], n)
  if (n > 0 && nrow(lab)) {
    aes_lab <- ggplot2::aes(x = .data$logFC, y = .data$neglog, label = .data$gene)
    p <- p + if (requireNamespace("ggrepel", quietly = TRUE))
      ggrepel::geom_text_repel(data = lab, mapping = aes_lab, inherit.aes = FALSE,
                               size = 3, color = "black", max.overlaps = 20)
    else
      ggplot2::geom_text(data = lab, mapping = aes_lab, inherit.aes = FALSE,
                         size = 3, color = "black", vjust = -0.6)
  }
  .scroll_apply_aspect(p, state$aspect %||% 1)
}

#' Stability plot for pseudobulk DE across random pseudo-replicate draws
#'
#' Median logFC vs selection frequency (how often a gene is a hit across runs);
#' genes above the consistency cutoff and the logFC threshold are highlighted.
#'
#' @param df A data.frame from [scroll_pseudobulk_stability()] (`gene`,
#'   `median_logFC`, `sel_freq`).
#' @param params List: `lfc` (effect threshold), `cut` (selection-frequency
#'   cutoff), `label_n` (genes to label).
#' @param state List: `aspect`.
#' @return A ggplot.
#' @export
view_stability <- function(df, params = list(), state = list()) {
  if (is.null(df) || !nrow(df))
    return(ggplot2::ggplot() +
             ggplot2::annotate("text", 0, 0, label = "No stability results.") +
             ggplot2::theme_void())
  lfc <- params$lfc %||% 1; cut <- params$cut %||% 0.8; n <- params$label_n %||% 15
  df$consistent <- df$sel_freq >= cut & abs(df$median_logFC) >= lfc
  p <- ggplot2::ggplot(df, ggplot2::aes(.data$median_logFC, .data$sel_freq,
                                        color = .data$consistent)) +
    ggplot2::geom_point(size = 1, alpha = 0.75) +
    ggplot2::scale_color_manual(values = c(`FALSE` = "grey78", `TRUE` = "#C4453B"),
                                guide = "none") +
    ggplot2::geom_vline(xintercept = c(-lfc, lfc), linetype = "dashed", color = "grey60") +
    ggplot2::geom_hline(yintercept = cut, linetype = "dashed", color = "grey60") +
    ggplot2::labs(x = "median logFC", y = "selection frequency") +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(color = "black", linewidth = 0.7, fill = NA))
  lab <- df[df$consistent, , drop = FALSE]
  lab <- utils::head(lab[order(-lab$sel_freq, -abs(lab$median_logFC)), , drop = FALSE], n)
  if (n > 0 && nrow(lab)) {
    aes_lab <- ggplot2::aes(x = .data$median_logFC, y = .data$sel_freq, label = .data$gene)
    p <- p + if (requireNamespace("ggrepel", quietly = TRUE))
      ggrepel::geom_text_repel(data = lab, mapping = aes_lab, inherit.aes = FALSE,
                               size = 3, color = "black", max.overlaps = 20)
    else
      ggplot2::geom_text(data = lab, mapping = aes_lab, inherit.aes = FALSE,
                         size = 3, color = "black", vjust = -0.6)
  }
  .scroll_apply_aspect(p, state$aspect %||% 1)
}

# --- dispatch -----------------------------------------------------------------

#' Which output kind a view produces
#'
#' @param view_type A view type name.
#' @return `"plot"` or `"table"`.
#' @export
scroll_view_kind <- function(view_type) {
  if (identical(view_type, "de_table")) "table" else "plot"
}

#' Dispatch a section's view to its renderer
#'
#' The single entry point the app calls. `ctx` carries `cells`, `params`,
#' `state` (toggles), and whatever payload the view needs (`feature_values`,
#' `feature_name`, `expr_long`, `de_data`). When a feature search is active it
#' retargets the focused view: scatter views recolor, violin switches feature,
#' dotplot includes the gene; proportions/de_table are unaffected.
#'
#' @param view_type One of the built-in view types.
#' @param ctx A context list (see Details).
#' @return A ggplot (for plot views) or a data.frame (for `de_table`).
#' @export
render_view <- function(view_type, ctx) {
  cells <- ctx$cells
  params <- ctx$params %||% list()
  state <- ctx$state %||% list()
  fname <- ctx$feature_name
  if (!is.null(fname)) params$feature <- fname

  switch(view_type,
    umap_colorby = if (!is.null(ctx$feature_values))
      view_feature_plot(cells, params, ctx$feature_values, state)
      else view_umap_colorby(cells, params, state),
    feature_plot = view_feature_plot(cells, params, ctx$feature_values, state),
    dotplot      = view_dotplot(cells, params, ctx$expr_long, state),
    violin       = view_violin(cells, params, ctx$feature_values, state),
    proportions  = view_proportions(cells, params, state),
    de_table     = view_de_table(ctx$de_data, params),
    stop("Unknown view type: ", view_type, call. = FALSE)
  )
}
