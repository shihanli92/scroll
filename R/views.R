# Built-in views. Views are pure functions of data frames (no Shiny), so they
# unit-test without a running app. render_view() is the single dispatcher the
# app calls; it unpacks a context list (cells + whatever payload the view needs
# + toggle state) and routes to the right view. The full *_ui/*_server module
# contract + register_view() is a later phase; keeping render logic pure here is
# what lets that refactor stay mechanical.

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# --- shared helpers -----------------------------------------------------------

# Effective embedding / grouping column honour toggles over config.
.scroll_eff_embedding <- function(params, state) state$embedding %||% params$embedding
.scroll_group_col <- function(params, state) state$split_by %||% params$group_by

# Attach an embedding's first two dims as .x / .y.
.scroll_embedding_xy <- function(cells, embedding) {
  xcol <- sprintf("%s_1", embedding)
  ycol <- sprintf("%s_2", embedding)
  if (!all(c(xcol, ycol) %in% names(cells)))
    stop("Embedding '", embedding, "' not found in cells table.", call. = FALSE)
  cells$.x <- cells[[xcol]]
  cells$.y <- cells[[ycol]]
  cells
}

# Expression vector aligned to a cells table, zero-filled for absent cells.
.scroll_expr_vector <- function(cells, values) {
  expr <- stats::setNames(rep(0, nrow(cells)), cells$cell)
  if (!is.null(values) && nrow(values) > 0) expr[values$cell] <- values$value
  as.numeric(expr[cells$cell])
}

.scroll_base_scatter <- function(df, embedding) {
  ggplot2::ggplot(df, ggplot2::aes(x = .data$.x, y = .data$.y)) +
    ggplot2::labs(x = sprintf("%s 1", embedding), y = sprintf("%s 2", embedding)) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      legend.position = "right"
    )
}

# Optional facet by a categorical toggle column.
.scroll_maybe_facet <- function(p, cells, state) {
  sb <- state$split_by
  if (!is.null(sb) && nzchar(sb %||% "") && sb %in% names(cells))
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste0("~ `", sb, "`")))
  p
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
  p <- .scroll_base_scatter(df, embedding) +
    ggplot2::geom_point(ggplot2::aes(color = .data$.col), size = 0.6, alpha = 0.8) +
    ggplot2::labs(color = color_by)
  if (is.numeric(df$.col)) {
    p <- p + ggplot2::scale_color_viridis_c()
  } else if (isTRUE(state$show_labels)) {
    p <- p + .scroll_group_labels(df, ".col") + ggplot2::guides(color = "none")
  }
  .scroll_maybe_facet(p, df, state)
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
  df <- df[order(df$.expr), , drop = FALSE]   # expressing cells drawn on top
  p <- .scroll_base_scatter(df, embedding) +
    ggplot2::geom_point(ggplot2::aes(color = .data$.expr), size = 0.6) +
    ggplot2::scale_color_gradient(low = "grey88", high = "#3b0f70") +
    ggplot2::labs(color = params$feature %||% "expression")
  .scroll_maybe_facet(p, df, state)
}

# --- aggregate views ----------------------------------------------------------

#' Dot plot: mean expression x fraction expressing, across groups
#'
#' @param cells The cells data.frame (its group sizes are the denominators).
#' @param params List with `features`, `group_by`, and optionally `assay`.
#' @param expr_long A data.frame(feature, cell, value) in normalized units for
#'   the requested features (zero-valued cells absent).
#' @param state Optional toggle state (`split_by` overrides `group_by`).
#' @return A ggplot.
#' @export
view_dotplot <- function(cells, params, expr_long, state = list()) {
  features <- unlist(params$features)
  group_by <- .scroll_group_col(params, state)
  if (is.null(group_by) || !group_by %in% names(cells))
    stop("dotplot needs a valid `group_by` metadata column.", call. = FALSE)

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
  agg$feature <- factor(agg$feature, levels = rev(features))

  ggplot2::ggplot(agg, ggplot2::aes(x = .data$group, y = .data$feature)) +
    ggplot2::geom_point(ggplot2::aes(size = .data$frac, color = .data$mean)) +
    ggplot2::scale_size_area(max_size = 8, labels = scales::percent,
                             name = "% expressing") +
    ggplot2::scale_color_viridis_c(name = "mean expr.") +
    ggplot2::labs(x = group_by, y = NULL) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(panel.grid.major = ggplot2::element_line(color = "grey92"),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Violin: a feature's per-group distribution
#'
#' @param cells The cells data.frame.
#' @param params List with `feature`, `group_by`, optionally `assay`.
#' @param values A data.frame(cell, value) for the feature, or `NULL`.
#' @param state Optional toggle state (`split_by` overrides `group_by`).
#' @return A ggplot.
#' @export
view_violin <- function(cells, params, values = NULL, state = list()) {
  group_by <- .scroll_group_col(params, state)
  if (is.null(group_by) || !group_by %in% names(cells))
    stop("violin needs a valid `group_by` metadata column.", call. = FALSE)
  df <- data.frame(cell = cells$cell, group = as.character(cells[[group_by]]),
                   stringsAsFactors = FALSE)
  df$expr <- .scroll_expr_vector(cells, values)
  ggplot2::ggplot(df, ggplot2::aes(x = .data$group, y = .data$expr,
                                   fill = .data$group)) +
    ggplot2::geom_violin(scale = "width", trim = TRUE, linewidth = 0.3) +
    ggplot2::labs(x = group_by, y = params$feature %||% "expression") +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Stacked composition of one categorical within another
#'
#' @param cells The cells data.frame.
#' @param params List with `group_by` (x axis) and `fill_by` (composition).
#' @param state Optional toggle state.
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
  totals <- stats::aggregate(Freq ~ x, tab, sum)
  tab <- merge(tab, totals, by = "x", suffixes = c("", ".total"))
  tab$frac <- ifelse(tab$Freq.total > 0, tab$Freq / tab$Freq.total, 0)
  ggplot2::ggplot(tab, ggplot2::aes(x = .data$x, y = .data$frac,
                                    fill = .data$fill)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::labs(x = x, y = "composition", fill = fill) +
    ggplot2::theme_minimal(base_size = 13)
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
