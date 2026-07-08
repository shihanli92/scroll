# Built-in views. Phase 1 keeps views as pure render functions of data frames
# (no Shiny), dispatched by type through render_view(). The full *_ui/*_server
# module contract + register_view() is a later phase; keeping the render logic
# pure here makes it unit-testable without a running app.

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

#' Embedding colored by a metadata column
#'
#' @param cells The globally-loaded cells data.frame.
#' @param params List with `embedding` and `color_by`.
#' @return A ggplot.
#' @export
view_umap_colorby <- function(cells, params) {
  df <- .scroll_embedding_xy(cells, params$embedding)
  color_by <- params$color_by
  if (!color_by %in% names(df))
    stop("color_by column '", color_by, "' not in cells table.", call. = FALSE)
  df$.col <- df[[color_by]]
  p <- .scroll_base_scatter(df, params$embedding) +
    ggplot2::geom_point(ggplot2::aes(color = .data$.col), size = 0.6, alpha = 0.8) +
    ggplot2::labs(color = color_by)
  if (is.numeric(df$.col))
    p <- p + ggplot2::scale_color_viridis_c()
  p
}

#' Embedding colored by a feature's expression
#'
#' @param cells The globally-loaded cells data.frame.
#' @param params List with `embedding` (and `feature`, for the title).
#' @param values A data.frame(cell, value) in normalized units, or `NULL`.
#'   Cells absent from `values` are treated as 0 (sparse).
#' @return A ggplot.
#' @export
view_feature_plot <- function(cells, params, values = NULL) {
  df <- .scroll_embedding_xy(cells, params$embedding)
  expr <- stats::setNames(rep(0, nrow(df)), df$cell)
  if (!is.null(values) && nrow(values) > 0) {
    expr[values$cell] <- values$value
  }
  df$.expr <- as.numeric(expr[df$cell])
  # Draw expressing cells last so they sit on top of the zero background.
  df <- df[order(df$.expr), , drop = FALSE]
  .scroll_base_scatter(df, params$embedding) +
    ggplot2::geom_point(ggplot2::aes(color = .data$.expr), size = 0.6) +
    ggplot2::scale_color_gradient(low = "grey88", high = "#3b0f70") +
    ggplot2::labs(color = params$feature %||% "expression")
}

#' Dispatch a section's view to its renderer
#'
#' The single entry point the app calls. When `feature_values` is supplied, the
#' active feature search overrides the section's native coloring (the persistent
#' search box retargets whatever view is in focus); otherwise the view renders
#' its configured coloring.
#'
#' @param view_type One of `"umap_colorby"`, `"feature_plot"`.
#' @param cells The cells data.frame.
#' @param params The section's params from config.yaml.
#' @param feature_values Optional data.frame(cell, value) for an active feature.
#' @param feature_name Optional active feature name (for labelling).
#' @return A ggplot.
#' @export
render_view <- function(view_type, cells, params,
                        feature_values = NULL, feature_name = NULL) {
  if (!is.null(feature_values)) {
    params$feature <- feature_name %||% params$feature
    return(view_feature_plot(cells, params, feature_values))
  }
  switch(view_type,
    umap_colorby = view_umap_colorby(cells, params),
    feature_plot = view_feature_plot(cells, params, NULL),
    stop("Unknown view type: ", view_type, call. = FALSE)
  )
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
