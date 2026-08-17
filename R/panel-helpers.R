# Helpers shared across the built-in panels (R/panel-*.R): manifest-derived
# control choices (categorical/numeric columns, embeddings, colour-by, view
# scoping), the manual palette picker, cosmetic-debounce + rasterization, and the
# per-plot toolbar / download / empty-state scaffolding.

.scroll_nz <- function(x) if (is.null(x) || !length(x) || !nzchar(x)) NULL else x

# Debounce a reactive snapshot of cosmetic inputs so a slider drag coalesces into
# one redraw (~120 ms after the last move) instead of one render per tick. Pass a
# `reactive({...})` (constructed at the call site so input dependencies are
# tracked). Data inputs are kept out of this reactive so a gene change stays
# immediate.
.scroll_cosmetic <- function(r, millis = 120) shiny::debounce(r, millis)

# Above this many plotted points, rasterize the scatter on screen (kept below the
# small test project's n so tests/direct view calls stay vector).
.scroll_raster_threshold <- 15000L

# Whether to rasterize on screen: the per-panel toggle (default on) AND the point
# count exceeding the threshold. With the toggle off, points stay vector even on
# large datasets (slower, but exact/zoomable).
.scroll_use_raster <- function(on, n) isTRUE(on %||% TRUE) && n > .scroll_raster_threshold

# --- manifest-derived control choices -----------------------------------------

# Metadata columns split by whether they belong to a subset view. A column's
# `$scope` (set at build time) names its owning subset; unscoped columns are
# global. The plain `.scroll_cat_cols`/`.scroll_num_cols` return global columns
# only, so subset-scoped columns never leak into whole-dataset selectors; the
# `_view_` variants add back the active subset's scoped columns.
.scroll_scope_of <- function(m, col) m$meta[[col]]$scope
.scroll_col_in_view <- function(m, col, view) {
  sc <- .scroll_scope_of(m, col)
  is.null(sc) || identical(sc, view)
}
.scroll_cat_cols <- function(m, view = NULL) {
  cols <- names(Filter(function(x) identical(x$type, "categorical"), m$meta))
  cols[vapply(cols, function(c) .scroll_col_in_view(m, c, view), logical(1))]
}
.scroll_num_cols <- function(m, view = NULL) {
  cols <- names(Filter(function(x) identical(x$type, "numeric"), m$meta))
  cols[vapply(cols, function(c) .scroll_col_in_view(m, c, view), logical(1))]
}
.scroll_reductions <- function(m) names(m$embeddings)
.scroll_assays_of <- function(m) names(m$assays)
.scroll_features_of <- function(m, assay) unlist(m$assays[[assay]]$features)

# Distinct-value count of a categorical column. `n_levels` is always written by
# recent builds; fall back to the cached `levels` length for older manifests.
.scroll_n_levels <- function(m, col) {
  m$meta[[col]]$n_levels %||% length(m$meta[[col]]$levels)
}

# The level set of a categorical column. The manifest caches `levels` only up to
# `max_levels` (see `.scroll_meta_entry`); above that the cache is absent and we
# recompute from the in-RAM cells table (cheap; NA -> a missing category, which
# for a subset-scoped column drops non-members, matching the build-time scoping).
.scroll_meta_levels <- function(data, col) {
  lv <- data$manifest$meta[[col]]$levels
  if (!is.null(lv)) return(unlist(lv, use.names = FALSE))
  v <- data$cells[[col]]
  if (is.null(v)) return(character(0))
  sort(unique(as.character(v[!is.na(v)])))
}

# Is there a categorical column with >=2 levels, i.e. a grouping a contrast can
# be built from? DE/Pseudobulk are pointless without one (e.g. a single-region
# Visium slide has `region` but only one level), so they gate on this. Uses the
# level *count* so a high-cardinality column (levels uncached) still qualifies.
.scroll_has_contrast <- function(m) {
  cats <- .scroll_cat_cols(m)
  any(vapply(cats, function(c) .scroll_n_levels(m, c) >= 2, logical(1)))
}
.scroll_subset_names <- function(m) names(m$subsets)

# Embeddings usable in a view: whole-dataset shows only full-coverage reductions;
# a subset view shows that subset's (partial) embeddings. Old manifests without
# `n_covered` treat every embedding as full (back-compatible).
.scroll_global_embeddings <- function(m) {
  rn <- names(m$embeddings)
  rn[vapply(rn, function(r) {
    nc <- m$embeddings[[r]]$n_covered
    is.null(nc) || nc >= m$n_cells
  }, logical(1))]
}
.scroll_view_embeddings <- function(m, view = NULL) {
  if (is.null(view) || is.null(m$subsets[[view]])) return(.scroll_global_embeddings(m))
  as.character(unlist(m$subsets[[view]]$embeddings))
}

# Embedding a live contrast preview (DE / Pseudobulk) should draw on: the config
# `default_embedding` for the whole dataset (so it matches DimPlot, not just the
# first global reduction), else the active subset view's primary embedding. Mirrors
# the reduction-default logic in dimplot_server.
.scroll_preview_embedding <- function(data, view) {
  reds <- .scroll_view_embeddings(data$manifest, view)
  sel <- if (is.null(view)) .scroll_default(data, "default_embedding", reds[[1]]) else reds[[1]]
  if (!sel %in% reds) sel <- reds[[1]]
  sel
}

# Restrict cells to a subset view (members = non-NA primary-embedding coords).
.scroll_view_cells <- function(cells, m, view) {
  if (is.null(view) || is.null(m$subsets[[view]])) return(cells)
  col <- sprintf("%s_1", m$subsets[[view]]$primary_embedding)
  if (!col %in% names(cells)) return(cells)
  cells[!is.na(cells[[col]]), , drop = FALSE]
}

.scroll_colorby_choices <- function(m, view = NULL) {
  ch <- list()
  cats <- .scroll_cat_cols(m, view); nums <- .scroll_num_cols(m, view)
  if (length(cats)) ch[["Cell annotations"]] <- as.list(stats::setNames(cats, cats))
  if (length(nums)) ch[["Numeric / QC"]] <- as.list(stats::setNames(nums, nums))
  ch
}

# Refresh a panel's categorical `selectInput`s when the active view changes so
# subset-scoped columns appear only inside their subset view. `prepend` is a
# leading choice map (e.g. c("None" = "")); the current selection is preserved
# when still valid, else falls back to `default` (or the prepend value).
# Keep a `selectInput`'s choices in sync with the active subset view, using
# `cols_fn(view)` to list the columns valid in that view (so scoped columns appear
# only in their view). Preserves the current pick when still valid.
.scroll_bind_view_cols <- function(input, session, view_r, cols_fn, ids,
                                   prepend = NULL, default = NULL) {
  observeEvent(view_r(), {
    cols <- cols_fn(view_r())
    for (id in ids) {
      cur <- input[[id]]
      sel <- if (!is.null(cur) && cur %in% cols) cur
             else if (!is.null(default) && default %in% cols) default
             else if (!is.null(prepend)) unname(prepend)[[1]]
             else if (length(cols)) cols[[1]] else NULL
      updateSelectInput(session, id,
                        choices = c(prepend, stats::setNames(cols, cols)),
                        selected = sel)
    }
  }, ignoreNULL = FALSE)
}

# Categorical / numeric view-aware selector binders (scoped columns surface only
# in their own view).
.scroll_bind_view_cats <- function(input, session, view_r, m, ids,
                                   prepend = NULL, default = NULL)
  .scroll_bind_view_cols(input, session, view_r,
                         function(v) .scroll_cat_cols(m, v), ids, prepend, default)

.scroll_bind_view_nums <- function(input, session, view_r, m, ids,
                                   prepend = NULL, default = NULL)
  .scroll_bind_view_cols(input, session, view_r,
                         function(v) .scroll_num_cols(m, v), ids, prepend, default)

.scroll_default <- function(data, key, fallback) {
  data$config[[key]] %||% data$manifest[[key]] %||% fallback
}

# Per-level colour pickers for the "Manual" palette, shared by the categorical
# panels. `.scroll_manual_ui` builds colourInputs (ids col_1..col_n) seeded from
# the Tableau-10 defaults; `.scroll_manual_colors` reads them back into a named
# vector (NULL until set). A panel wires them by: adding "Manual" to its palette
# choices, a `uiOutput(ns("manual"))`, and passing `manual_colors` into the view
# state (which `.scroll_group_colors` honours when palette == "Manual").
# Above this many levels, per-level colour pickers become an unusable wall of
# widgets (and a DOM blow-up), so fall back to a note pointing at a palette.
.SCROLL_MANUAL_CAP <- 30L
.scroll_manual_ui <- function(ns, levels) {
  if (length(levels) > .SCROLL_MANUAL_CAP)
    return(tags$p(class = "scroll-desc",
                  sprintf("Manual colours are unavailable for %d-level columns - pick a palette instead.",
                          length(levels))))
  defaults <- .scroll_discrete_colors(levels, "Tableau 10")
  # compact circular swatches that wrap into rows (rather than tall full-width
  # inputs), so many levels stay manageable; the level name is a caption + title.
  div(class = "scroll-manual-grid",
    lapply(seq_along(levels), function(i)
      div(class = "scroll-swatch", title = levels[i],
        colourpicker::colourInput(ns(paste0("col_", i)), label = NULL,
                                  value = defaults[[levels[i]]],
                                  showColour = "both", closeOnClick = TRUE),
        span(class = "scroll-swatch-label", levels[i]))))
}
.scroll_manual_colors <- function(input, levels) {
  vals <- lapply(seq_along(levels), function(i) input[[paste0("col_", i)]])
  names(vals) <- levels
  vals <- vals[!vapply(vals, is.null, logical(1))]
  if (length(vals)) unlist(vals) else NULL
}
# Categorical palette choices including the Manual option.
.scroll_cat_palettes <- function() c(names(.scroll_discrete_palettes), "Manual")

.scroll_group <- function(title, ...) {
  div(class = "scroll-cgroup", div(class = "scroll-cgroup-h", title), ...)
}

# --- per-plot export ----------------------------------------------------------

# A compact download button for a plot/table toolbar.
.scroll_dl_button <- function(id, label)
  downloadButton(id, label, class = "btn-sm scroll-dl", icon = shiny::icon("download"))

# Wrap an output in a loading spinner when shinycssloaders is available (a
# guarded Suggests dep); otherwise return the output unchanged.
.scroll_spin <- function(tag) {
  if (requireNamespace("shinycssloaders", quietly = TRUE))
    shinycssloaders::withSpinner(tag, type = 6, color = "#2563A8", size = 0.6,
                                 proxy.height = "260px")
  else tag
}

# Standard plot area: a small toolbar (PNG + PDF, and optionally a CSV of the
# plot's source data) above the plot output. Pass csv = TRUE to add the CSV button
# (the server must then wire output$csv, e.g. via .scroll_plot_downloads(csv_r=)).
.scroll_plot_area <- function(ns, height = "460px", csv = FALSE)
  div(class = "scroll-plot",
      div(class = "scroll-plot-bar",
          .scroll_dl_button(ns("png"), "PNG"),
          .scroll_dl_button(ns("pdf"), "PDF"),
          if (isTRUE(csv)) .scroll_dl_button(ns("csv"), "CSV")),
      .scroll_spin(plotOutput(ns("plot"), height = height)))

# Image downloadHandler for a plot reactive, raster (PNG) or vector (PDF). Uses
# a device + print() (not ggsave) so it renders both bare ggplots and the
# DotPlot's aplot composite, which ggsave() rejects as a non-ggplot.
.scroll_img_handler <- function(plot_r, name, format = c("png", "pdf")) {
  format <- match.arg(format)
  downloadHandler(
    filename = function() name,
    content = function(file) {
      if (format == "pdf") grDevices::pdf(file, width = 8, height = 6, bg = "white")
      else grDevices::png(file, width = 8, height = 6, units = "in", res = 150, bg = "white")
      on.exit(grDevices::dev.off())
      print(plot_r())
    })
}

# Register the standard PNG + PDF download outputs (ids "png"/"pdf") for a plot
# reactive on a module's `output`. Filenames stem from the panel id. When `csv_r`
# (a data.frame reactive of the plot's source data) is supplied, also wire output$csv.
.scroll_plot_downloads <- function(output, plot_r, id, csv_r = NULL) {
  output$png <- .scroll_img_handler(plot_r, paste0("scroll_", id, ".png"), "png")
  output$pdf <- .scroll_img_handler(plot_r, paste0("scroll_", id, ".pdf"), "pdf")
  if (!is.null(csv_r)) output$csv <- .scroll_csv_handler(csv_r, paste0("scroll_", id, ".csv"))
  invisible()
}

# CSV downloadHandler for a data.frame reactive.
.scroll_csv_handler <- function(df_r, name)
  downloadHandler(
    filename = function() name,
    content = function(file) utils::write.csv(df_r(), file, row.names = FALSE))

# Shared aspect-ratio control, applied uniformly as a rendered-height multiplier
# (works for every panel including the aplot dendrogram composite, where
# theme(aspect.ratio) would detach the trees). Every panel adds the input and
# sets its renderPlot height via .scroll_plot_height(); future panels get it for
# free by doing the same.
.scroll_aspect_input <- function(ns) sliderInput(ns("aspect"), "Aspect ratio", 0.4, 3, 1, 0.1)

# Placeholder body for a panel that can't run on this dataset (e.g. no categorical
# grouping column) — avoids a hard `cats[[1]]` crash at UI build. This is the second
# of two gating layers: a panel's `when` predicate drops it from the whole app when
# NO mounted dataset supports it, while this inline guard covers the per-dataset case
# in scroll_multi_app (a panel kept because another tab supports it must still render
# a friendly message for a tab that doesn't).
.scroll_empty_panel <- function(msg)
  div(class = "scroll-panel", tags$p(class = "scroll-desc", msg))

