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

# On-screen scatter point cap: above this, the interactive (rasterized) draw uses a
# deterministic subsample so a large embedding renders ~2x faster; exports (raster =
# FALSE) always draw every point. Override with options(scroll.onscreen_cap = N);
# NA disables. Chosen so a rasterized cloud looks the same as the full set.
.SCROLL_ONSCREEN_CAP <- 50000L
.scroll_onscreen_cap <- function() getOption("scroll.onscreen_cap", .SCROLL_ONSCREEN_CAP)

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
                                   prepend = NULL, default = NULL,
                                   multiple = FALSE, allow_empty = FALSE) {
  observeEvent(view_r(), {
    cols <- cols_fn(view_r())
    for (id in ids) {
      cur <- input[[id]]
      if (multiple) {                                  # selectize interaction picker
        sel <- intersect(cur, cols)
        if (!length(sel) && !allow_empty && length(cols)) sel <- cols[[1]]
        updateSelectizeInput(session, id,
                             choices = stats::setNames(cols, cols), selected = sel)
      } else {                                         # single selectInput (may prepend "None")
        sel <- if (length(cur) == 1 && cur %in% cols) cur
               else if (!is.null(default) && default %in% cols) default
               else if (!is.null(prepend)) unname(prepend)[[1]]
               else if (length(cols)) cols[[1]] else NULL
        updateSelectInput(session, id,
                          choices = c(prepend, stats::setNames(cols, cols)), selected = sel)
      }
    }
  }, ignoreNULL = FALSE)
}

# Categorical / numeric view-aware selector binders (scoped columns surface only
# in their own view). `multiple = TRUE` drives a selectize interaction picker
# (keeping the still-valid subset of the pick); `allow_empty` lets it stay empty
# rather than falling back to the first column.
.scroll_bind_view_cats <- function(input, session, view_r, m, ids,
                                   prepend = NULL, default = NULL, ...)
  .scroll_bind_view_cols(input, session, view_r,
                         function(v) .scroll_cat_cols(m, v), ids, prepend, default, ...)

.scroll_bind_view_nums <- function(input, session, view_r, m, ids,
                                   prepend = NULL, default = NULL, ...)
  .scroll_bind_view_cols(input, session, view_r,
                         function(v) .scroll_num_cols(m, v), ids, prepend, default, ...)

.scroll_default <- function(data, key, fallback) {
  data$config[[key]] %||% data$manifest[[key]] %||% fallback
}

# Wire a server-side, assay-aware gene selectize with the behaviours every gene box
# shares, so the panels don't each hand-roll them:
#  - repopulate for the active assay, keeping the still-valid current pick (or
#    `defaults` when nothing is selected yet);
#  - accept a pasted, delimited gene list when `paste_id` is given (case-insensitive
#    match, unknown tokens reported);
#  - when `warn_n`/`mark_id` are given, warn past `warn_n` genes (via the `warn_msg`
#    sprintf template) and mirror the current pick into a `mark_id` "Label genes"
#    selectize.
# `assay_r` is the active-assay reactive; `m` the manifest.
.scroll_bind_gene_box <- function(input, session, m, assay_r, id,
                                  paste_id = NULL, defaults = NULL,
                                  warn_n = NULL, warn_msg = NULL, mark_id = NULL) {
  observeEvent(assay_r(), {
    feats <- .scroll_features_of(m, assay_r())
    cur <- isolate(input[[id]]) %||% defaults
    updateSelectizeInput(session, id, choices = feats, server = TRUE,
                         selected = intersect(cur, feats))
  })
  if (!is.null(paste_id))
    observeEvent(input[[paste_id]], {
      feats <- .scroll_features_of(m, assay_r())
      parsed <- .scroll_parse_gene_list(input[[paste_id]], feats)
      req(length(parsed$ok) > 0 || length(parsed$missing) > 0)
      sel <- unique(c(isolate(input[[id]]), parsed$ok))
      updateSelectizeInput(session, id, choices = feats, server = TRUE, selected = sel)
      if (length(parsed$missing))
        showNotification(paste("Not in this assay:", paste(parsed$missing, collapse = ", ")),
                         type = "warning", duration = 6)
    })
  if (!is.null(warn_n) || !is.null(mark_id))
    observeEvent(input[[id]], {
      if (!is.null(warn_n) && length(input[[id]]) > warn_n)
        showNotification(sprintf(warn_msg %||% "%d genes selected - this may be slow to compute.",
                                 length(input[[id]])), type = "warning", duration = 5)
      if (!is.null(mark_id))                            # mirror the pick into "Label genes"
        updateSelectizeInput(session, mark_id, choices = input[[id]],
                             selected = intersect(isolate(input[[mark_id]]), input[[id]]))
    }, ignoreNULL = FALSE)
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
# `prefix` names the picker inputs (`<prefix>_<i>`); pass a distinct prefix when a
# panel has more than one manual-colour control (e.g. group vs expansion), so their
# inputs don't collide. `defaults` optionally seeds the initial swatch colours.
.scroll_manual_ui <- function(ns, levels, prefix = "col", defaults = NULL) {
  if (length(levels) > .SCROLL_MANUAL_CAP)
    return(tags$p(class = "scroll-desc",
                  sprintf("Manual colours are unavailable for %d-level columns - pick a palette instead.",
                          length(levels))))
  defaults <- defaults %||% .scroll_discrete_colors(levels, "Tableau 10")
  # compact circular swatches that wrap into rows (rather than tall full-width
  # inputs), so many levels stay manageable; the level name is a caption + title.
  div(class = "scroll-manual-grid",
    lapply(seq_along(levels), function(i)
      div(class = "scroll-swatch", title = levels[i],
        colourpicker::colourInput(ns(paste0(prefix, "_", i)), label = NULL,
                                  value = defaults[[levels[i]]],
                                  showColour = "both", closeOnClick = TRUE),
        span(class = "scroll-swatch-label", levels[i]))))
}
.scroll_manual_colors <- function(input, levels, prefix = "col") {
  vals <- lapply(seq_along(levels), function(i) input[[paste0(prefix, "_", i)]])
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
    # hide.ui = FALSE keeps the previous render visible under the spinner while a plot
    # recomputes, so a control change doesn't blank the panel before the new image lands
    shinycssloaders::withSpinner(tag, type = 6, color = "#2563A8", size = 0.6,
                                 proxy.height = "260px", hide.ui = FALSE)
  else tag
}

# Lazy on-screen rendering: gate a plot so it recomputes only while its card is on
# screen, and otherwise keeps its last render. A View/filter change then redraws
# only the panels in view; an off-screen panel refreshes when scrolled to. The
# `onscreen` input is driven by the IntersectionObserver in .scroll_lazy_js and
# defaults on, so a panel renders before/without JS (and under testServer).
#
# `build` is a zero-arg function returning the plot. It is wrapped in a memoized
# reactive so that scrolling a panel off screen and back does NOT rebuild it: the
# rebuild only happens if `build`'s own dependencies (data/controls) changed while
# it was away. While off screen the gate takes no dependency on `built`, so a
# View/filter change does not invalidate an off-screen panel until it returns.
#
# `style_r` (optional) is the panel's Style-sheet value; its labels are applied to
# whatever `build` returns (.scroll_apply_style), inside the memoized reactive.
.scroll_lazy_plot <- function(input, build, style_r = NULL) {
  built <- reactive(.scroll_style_plot(build(), style_r))   # memoized; reports the layers
  last <- NULL
  reactive({
    if (isTRUE(input$onscreen %||% TRUE)) last <<- built()
    shiny::req(last)
    last
  })
}

# Wire `output$plot` from a panel's (lazy) plot reactive, optionally caching the
# rendered IMAGE. When `cache` is non-NULL (the app passes "app" -> Shiny's shared
# app-level cache), `bindCache` memoizes the drawn PNG keyed on `key_r()`, so
# returning to a (view + controls) state serves the image without re-drawing
# (~ms vs hundreds of ms). `cache = NULL` (the default, and every testServer path)
# is the unchanged plain `renderPlot`. The key MUST capture every input `build()`
# reads -- including `input$onscreen` -- or a cache hit could serve a stale image.
.scroll_render_cached <- function(output, plot_r, key_r, cache = NULL) {
  r <- shiny::renderPlot(plot_r())
  if (!is.null(cache)) r <- shiny::bindCache(r, key_r(), cache = cache)
  output$plot <- r
}

# Standard plot area: a small toolbar (PNG + PDF, and optionally a CSV of the
# plot's source data) above the plot output. Pass csv = TRUE to add the CSV button
# (the server must then wire output$csv, e.g. via .scroll_plot_downloads(csv_r=)).
# A compact download-scale slider for a plot toolbar: sets the module's `dl_scale`
# input (1-5), which the PNG/PDF handlers pass to ggsave(scale=) so the *exported*
# figure is larger/smaller. It does not change the on-screen plot. A live "Nx" readout
# (the next sibling) shows the current value.
.scroll_size_slider <- function(ns)
  tags$span(class = "scroll-size-wrap", title = "Export (download) scale",
    tags$input(type = "range", class = "scroll-size", min = "1", max = "5", step = "0.5",
               value = "1", `data-input` = ns("dl_scale"),
               oninput = paste0("Shiny.setInputValue(this.dataset.input, parseFloat(this.value));",
                                "this.nextElementSibling.textContent=this.value+'\u00d7';")),
    tags$span(class = "scroll-size-val", "1\u00d7"))

# Default plot height: fills the card but bounded so a tall/portrait monitor does
# not stretch a plot to ~1700px. min(viewport-height, 90% of the container width)
# keeps it roughly landscape; clamped to [320, 1100]px and snapped to a whole CSS
# pixel (mirrors the width snap) to avoid HiDPI sub-pixel re-render loops. `90cqw`
# resolves against .scroll-content (container-type:inline-size), falling back to
# viewport units where no container ancestor exists.
.SCROLL_PLOT_H <- "calc(round(down, clamp(320px, min(100vh - 190px, 90cqw), 1100px), 1px))"

.scroll_plot_area <- function(ns, height = .SCROLL_PLOT_H, csv = FALSE)
  div(class = "scroll-plot",
      div(class = "scroll-plot-bar",
          .scroll_style_button(ns),
          .scroll_size_slider(ns),
          .scroll_dl_button(ns("png"), "PNG"),
          .scroll_dl_button(ns("pdf"), "PDF"),
          if (isTRUE(csv)) .scroll_dl_button(ns("csv"), "CSV")),
      # the hold reserves the plot's height so that when an off-screen card is
      # suspended (its output display:none'd for lazy rendering) the card does not
      # collapse and jump the scroll position.
      div(class = "scroll-plot-hold", style = sprintf("min-height:%s;", height),
          .scroll_spin(plotOutput(ns("plot"), height = height))))

# The download-scale slider value (the ggsave `scale`), read from a module `session`
# (or reactive domain); 1 when absent, clamped to (0, 5].
.scroll_dl_scale <- function(dom = shiny::getDefaultReactiveDomain()) {
  v <- if (!is.null(dom)) tryCatch(shiny::isolate(dom$input[["dl_scale"]]),
                                   error = function(e) NULL) else NULL
  s <- suppressWarnings(as.numeric(v))
  if (isTRUE(is.finite(s)) && s > 0) min(s, 5) else 1
}

# Write a plot to a raster (PNG) or vector (PDF) file at a ggsave `scale`. A bare ggplot
# goes through ggsave; the DotPlot's aplot composite (which ggsave rejects) uses a device
# + print with the canvas grown by scale -- the same effect ggsave's `scale` has (it
# multiplies the output dimensions in inches).
.scroll_write_plot <- function(file, p, format = c("png", "pdf"), scale = 1) {
  format <- match.arg(format)
  if (inherits(p, "ggplot")) {
    ggplot2::ggsave(file, plot = p, device = format, width = 8, height = 6,
                    units = "in", dpi = 150, scale = scale, bg = "white")
  } else {
    if (format == "pdf") grDevices::pdf(file, width = 8 * scale, height = 6 * scale, bg = "white")
    else grDevices::png(file, width = 8 * scale, height = 6 * scale, units = "in", res = 150, bg = "white")
    on.exit(grDevices::dev.off())
    print(p)
  }
}

# Image downloadHandler for a plot reactive. `scale` is a 0-arg fn returning the
# download-scale slider value, so a larger value exports a bigger figure.
.scroll_img_handler <- function(plot_r, name, format = c("png", "pdf"), scale = function() 1) {
  format <- match.arg(format)
  downloadHandler(
    filename = function() name,
    content = function(file) .scroll_write_plot(file, plot_r(), format, scale()))
}

# Register the standard PNG + PDF download outputs (ids "png"/"pdf") for a plot
# reactive on a module's `output`. Filenames stem from the panel id. When `csv_r`
# (a data.frame reactive of the plot's source data) is supplied, also wire output$csv.
.scroll_plot_downloads <- function(output, plot_r, id, csv_r = NULL) {
  # capture the module session now (during server setup) so the handler reads the
  # right (namespaced) dl_scale input when it later fires.
  dom <- shiny::getDefaultReactiveDomain()
  scale <- function() .scroll_dl_scale(dom)
  output$png <- .scroll_img_handler(plot_r, paste0("scroll_", id, ".png"), "png", scale)
  output$pdf <- .scroll_img_handler(plot_r, paste0("scroll_", id, ".pdf"), "pdf", scale)
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

