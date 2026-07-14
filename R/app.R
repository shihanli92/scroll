# The runtime: a polished bslib explorer with one scrolling section per analysis
# type, each with its own controls. The app never loads the Seurat object -- it
# reads the built artifacts (cells.parquet once globally; expression one feature
# at a time via arrow), so runtime RAM stays flat.

# --- data handle --------------------------------------------------------------

# Read a baked side-car asset (repertoire table, tissue raster, peak table) from a
# project, returning NULL if it is absent or unreadable. Shared by the modality
# accessors (.scroll_vdj_read / .scroll_spatial_image / .scroll_atac_read).
.scroll_read_asset <- function(path, format = c("parquet", "rds")) {
  format <- match.arg(format)
  if (!file.exists(path)) return(NULL)
  tryCatch(
    if (format == "rds") readRDS(path)
    else as.data.frame(arrow::read_parquet(path)),
    error = function(e) NULL)
}

# Load the artifacts once (shared across sessions of one app process).
.scroll_load <- function(dir) {
  con <- scroll_connect(dir)
  manifest <- scroll_manifest(dir)
  cells <- as.data.frame(arrow::read_parquet(file.path(dir, "cells.parquet")))
  # `.gidx` is the canonical global row-index a v2 expr store joins on; it rides
  # along on subset/filtered cells so joins stay correct under the app-bar filter.
  cells$.gidx <- seq_len(nrow(cells))
  d <- list(
    dir = normalizePath(dir),
    cells = cells,
    manifest = manifest,
    cell_index = isTRUE(manifest$cell_index),   # TRUE for v2 stores (int cell key)
    config = scroll_config(dir),
    con = con
  )
  # bound query helpers that dequantize to normalized units, memoized by an LRU
  # so cosmetic re-renders, the vector-export path, and re-selecting a recent gene
  # never re-hit the store.
  cache <- .scroll_lru(256L)
  d$query1 <- function(assay, feature) {
    key <- paste0("1|", assay, "|", feature)
    hit <- cache$get(key)
    if (is.null(hit)) {
      hit <- scroll_query_feature(con, assay, feature)
      if (nrow(hit)) hit$value <- scroll_dequantize(hit$value, manifest, assay)
      cache$set(key, hit)
    }
    hit
  }
  d$queryN <- function(assay, features) {
    key <- paste0("N|", assay, "|", paste(sort(unique(features)), collapse = ","))
    hit <- cache$get(key)
    if (is.null(hit)) {
      hit <- scroll_query_features(con, assay, features)
      if (nrow(hit)) hit$value <- scroll_dequantize(hit$value, manifest, assay)
      cache$set(key, hit)
    }
    hit
  }
  d
}

# A bounded LRU over query results keyed by a string. Bounded so a long session
# browsing many genes cannot grow memory without limit; each entry is one
# feature's small sparse long table, so a generous bound is cheap.
.scroll_lru <- function(max = 256L) {
  store <- new.env(parent = emptyenv())
  order <- character(0)
  list(
    get = function(key) {
      if (!exists(key, store, inherits = FALSE)) return(NULL)
      order <<- c(setdiff(order, key), key)                # accessing bumps to MRU
      base::get(key, store)
    },
    set = function(key, val) {
      order <<- c(setdiff(order, key), key)                # most-recently-used last
      assign(key, val, store)
      while (length(order) > max) {                        # evict least-recently-used
        rm(list = order[[1]], envir = store); order <<- order[-1L]
      }
    }
  )
}

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

# Is there a categorical column with >=2 levels, i.e. a grouping a contrast can
# be built from? DE/Pseudobulk are pointless without one (e.g. a single-region
# Visium slide has `region` but only one level), so they gate on this.
.scroll_has_contrast <- function(m) {
  cats <- .scroll_cat_cols(m)
  any(vapply(cats, function(c) length(m$meta[[c]]$levels) >= 2, logical(1)))
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
.scroll_bind_view_cats <- function(input, session, view_r, m, ids,
                                    prepend = NULL, default = NULL) {
  observeEvent(view_r(), {
    cats <- .scroll_cat_cols(m, view_r())
    for (id in ids) {
      cur <- input[[id]]
      sel <- if (!is.null(cur) && cur %in% cats) cur
             else if (!is.null(default) && default %in% cats) default
             else if (!is.null(prepend)) unname(prepend)[[1]]
             else if (length(cats)) cats[[1]] else NULL
      updateSelectInput(session, id,
                        choices = c(prepend, stats::setNames(cats, cats)),
                        selected = sel)
    }
  }, ignoreNULL = FALSE)
}

.scroll_default <- function(data, key, fallback) {
  data$config[[key]] %||% data$manifest[[key]] %||% fallback
}

# Per-level colour pickers for the "Manual" palette, shared by the categorical
# panels. `.scroll_manual_ui` builds colourInputs (ids col_1..col_n) seeded from
# the Tableau-10 defaults; `.scroll_manual_colors` reads them back into a named
# vector (NULL until set). A panel wires them by: adding "Manual" to its palette
# choices, a `uiOutput(ns("manual"))`, and passing `manual_colors` into the view
# state (which `.scroll_group_colors` honours when palette == "Manual").
.scroll_manual_ui <- function(ns, levels) {
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

# --- DimPlot panel ------------------------------------------------------------

dimplot_ui <- function(id, data) {
  ns <- NS(id)
  m <- data$manifest
  reductions <- .scroll_view_embeddings(m, NULL)   # whole-dataset (full) embeddings
  cats <- .scroll_cat_cols(m)
  first_cat <- if (length(cats)) cats[[1]] else .scroll_num_cols(m)[[1]]
  bslib::layout_columns(
    col_widths = c(3, 9), class = "scroll-panel",
    div(
      class = "scroll-controls",
      .scroll_group("Embedding",
        selectInput(ns("reduction"), "Reduction", reductions,
                    selected = .scroll_default(data, "default_embedding", reductions[[1]])),
        selectInput(ns("colorby"), "Color by", .scroll_colorby_choices(m), selected = first_cat)),
      .scroll_group("Groups",
        selectizeInput(ns("highlight"), "Highlight", choices = NULL, multiple = TRUE,
                       options = list(placeholder = "All groups"))),
      .scroll_group("Appearance",
        selectInput(ns("palette"), "Palette", names(.scroll_discrete_palettes)),
        uiOutput(ns("manual")),
        sliderInput(ns("size"), "Point size", 0.1, 5, 0.6, 0.1),
        sliderInput(ns("alpha"), "Opacity", 0.1, 1, 0.85, 0.05),
        bslib::input_switch(ns("labels"), "Cluster labels", TRUE),
        bslib::input_switch(ns("legend"), "Legend", TRUE),
        bslib::input_switch(ns("raster"), "Rasterize (fast)", TRUE)),
      .scroll_group("Layout",
        selectInput(ns("split"), "Split by",
                    c("None" = "", stats::setNames(cats, cats))),
        .scroll_aspect_input(ns))
    ),
    .scroll_plot_area(ns, "460px")
  )
}

dimplot_server <- function(id, data, cells_r = reactive(data$cells),
                           view_r = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    m <- data$manifest
    is_cat <- reactive(identical(m$meta[[input$colorby]]$type, "categorical"))
    levels_of <- reactive(if (is_cat()) unlist(m$meta[[input$colorby]]$levels) else character(0))

    # When the active view changes: restrict reductions to that view's embeddings
    # (subset views default to the primary sub-embedding) and add the subset's
    # scoped columns to Color-by. Split-by follows the same scoped-column rule.
    observeEvent(view_r(), {
      reds <- .scroll_view_embeddings(m, view_r())
      sel_red <- if (is.null(view_r()))
        .scroll_default(data, "default_embedding", reds[[1]]) else reds[[1]]
      if (!sel_red %in% reds) sel_red <- reds[[1]]
      updateSelectInput(session, "reduction", choices = reds, selected = sel_red)
      cb <- .scroll_colorby_choices(m, view_r())
      flat <- unlist(cb, use.names = FALSE)
      cur <- input$colorby
      updateSelectInput(session, "colorby", choices = cb,
                        selected = if (!is.null(cur) && cur %in% flat) cur else flat[[1]])
    }, ignoreNULL = FALSE)
    .scroll_bind_view_cats(input, session, view_r, m, "split", prepend = c("None" = ""))

    observeEvent(input$colorby, {
      if (is_cat()) {
        updateSelectInput(session, "palette",
                          choices = c(names(.scroll_discrete_palettes), "Manual"), selected = "Tableau 10")
        updateSelectizeInput(session, "highlight", choices = levels_of(), selected = character(0))
      } else {
        updateSelectInput(session, "palette",
                          choices = .scroll_continuous_palettes, selected = "viridis")
        updateSelectizeInput(session, "highlight", choices = character(0), selected = character(0))
      }
    })

    # per-group color pickers, shown only when the palette is "Manual"
    output$manual <- renderUI({
      if (!is_cat() || !identical(input$palette, "Manual")) return(NULL)
      .scroll_manual_ui(session$ns, levels_of())
    })
    manual_colors <- reactive({
      if (!is_cat() || !identical(input$palette, "Manual")) return(NULL)
      .scroll_manual_colors(input, levels_of())
    })

    # DATA reactive (cells + params) invalidates only on data-input changes.
    data_r <- reactive({
      req(input$reduction, input$colorby)
      cells <- cells_r()
      list(cells = cells, embedding = input$reduction, color_by = input$colorby,
           n = nrow(cells))
    })
    # COSMETIC reactive, debounced; restyle-only inputs (incl. highlight/manual).
    cosmetic_r <- .scroll_cosmetic(reactive(
      list(palette = input$palette, point_size = input$size, alpha = input$alpha,
           show_labels = isTRUE(input$labels), legend = isTRUE(input$legend),
           split_by = .scroll_nz(input$split), aspect = input$aspect,
           highlight = input$highlight, manual_colors = manual_colors())))
    build <- function(raster) {
      d <- data_r(); st <- cosmetic_r(); st$raster <- raster
      view_umap_colorby(d$cells, list(embedding = d$embedding, color_by = d$color_by), st)
    }
    plot_r   <- reactive(build(.scroll_use_raster(input$raster, data_r()$n)))  # rasterized on screen
    export_r <- reactive(build(FALSE))                                         # vector for downloads
    output$plot <- renderPlot(plot_r())
    .scroll_plot_downloads(output, export_r, id)
  })
}

# --- FeaturePlot panel --------------------------------------------------------

featureplot_ui <- function(id, data) {
  ns <- NS(id)
  m <- data$manifest
  assays <- .scroll_assays_of(m); cats <- .scroll_cat_cols(m)
  bslib::layout_columns(
    col_widths = c(3, 9), class = "scroll-panel",
    div(
      class = "scroll-controls",
      .scroll_group("Feature",
        selectizeInput(ns("feature"), "Gene", choices = NULL, multiple = FALSE,
                       options = list(placeholder = "Search a gene...", maxOptions = 50)),
        if (length(.scroll_num_cols(m)))
          selectInput(ns("metacol"), "...or numeric column",
                      c("(use gene)" = "", stats::setNames(.scroll_num_cols(m), .scroll_num_cols(m)))),
        if (length(assays) > 1)
          selectInput(ns("assay"), "Assay", assays, selected = m$default_assay)),
      .scroll_group("Embedding",
        selectInput(ns("reduction"), "Reduction", .scroll_view_embeddings(m, NULL),
                    selected = .scroll_default(data, "default_embedding", .scroll_view_embeddings(m, NULL)[[1]]))),
      .scroll_group("Appearance",
        selectInput(ns("palette"), "Palette", .scroll_continuous_palettes, selected = "grey-purple"),
        sliderInput(ns("size"), "Point size", 0.1, 5, 0.7, 0.1),
        sliderInput(ns("clip"), "Color quantiles (%)", 0, 100, c(0, 100), 1),
        bslib::input_switch(ns("order"), "Expressing cells on top", TRUE),
        bslib::input_switch(ns("legend"), "Legend", TRUE),
        bslib::input_switch(ns("raster"), "Rasterize (fast)", TRUE)),
      .scroll_group("Layout",
        selectInput(ns("split"), "Split by", c("None" = "", stats::setNames(cats, cats))),
        .scroll_aspect_input(ns))
    ),
    .scroll_plot_area(ns, "460px")
  )
}

featureplot_server <- function(id, data, cells_r = reactive(data$cells),
                               view_r = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    m <- data$manifest
    assay <- reactive(input$assay %||% m$default_assay)
    observeEvent(view_r(), {
      reds <- .scroll_view_embeddings(m, view_r())
      sel_red <- if (is.null(view_r()))
        .scroll_default(data, "default_embedding", reds[[1]]) else reds[[1]]
      if (!sel_red %in% reds) sel_red <- reds[[1]]
      updateSelectInput(session, "reduction", choices = reds, selected = sel_red)
      # numeric-column choices follow the view (scoped numerics appear in-view)
      nums <- .scroll_num_cols(m, view_r()); cur <- input$metacol
      updateSelectInput(session, "metacol",
                        choices = c("(use gene)" = "", stats::setNames(nums, nums)),
                        selected = if (!is.null(cur) && cur %in% nums) cur else "")
    }, ignoreNULL = FALSE)
    .scroll_bind_view_cats(input, session, view_r, m, "split", prepend = c("None" = ""))
    # repopulate the gene list for the active assay; drop a selection that does
    # not exist in the newly chosen assay (else it silently queries empty)
    observeEvent(assay(), {
      feats <- .scroll_features_of(m, assay())
      cur <- isolate(input$feature)
      keep <- if (!is.null(cur) && cur %in% feats) cur else character(0)
      updateSelectizeInput(session, "feature", choices = feats, server = TRUE, selected = keep)
    })
    # DATA reactive: cells + the (cached) expression query. Cosmetic drags do not
    # invalidate it, so they never re-hit the store.
    data_r <- reactive({
      req(input$reduction)
      cells <- cells_r()
      metacol <- .scroll_nz(input$metacol)
      # a numeric metadata column is coloured like expression (same continuous
      # controls: quantile clip, order-on-top, palette) but read from `cells`,
      # not queried; picking one takes precedence over the gene selector.
      if (!is.null(metacol) && metacol %in% names(cells)) {
        vals <- data.frame(cell = cells$cell,
                           value = suppressWarnings(as.numeric(cells[[metacol]])),
                           stringsAsFactors = FALSE)
        return(list(cells = cells, embedding = input$reduction, feature = metacol,
                    values = vals, n = nrow(cells)))
      }
      feat <- .scroll_nz(input$feature)
      validate(need(!is.null(feat),
                    "Search for a gene, or pick a numeric column, to colour the embedding."))
      list(cells = cells, embedding = input$reduction, feature = feat,
           values = data$query1(assay(), feat), n = nrow(cells))
    })
    cosmetic_r <- .scroll_cosmetic(reactive(
      list(palette = input$palette, point_size = input$size,
           order = isTRUE(input$order), legend = isTRUE(input$legend),
           clip = input$clip / 100, split_by = .scroll_nz(input$split),
           aspect = input$aspect)))
    build <- function(raster) {
      d <- data_r(); st <- cosmetic_r(); st$raster <- raster
      view_feature_plot(d$cells, list(embedding = d$embedding, feature = d$feature),
                        d$values, st)
    }
    plot_r   <- reactive(build(.scroll_use_raster(input$raster, data_r()$n)))
    export_r <- reactive(build(FALSE))
    output$plot <- renderPlot(plot_r())
    .scroll_plot_downloads(output, export_r, id)
  })
}

# --- DotPlot panel ------------------------------------------------------------

dotplot_ui <- function(id, data) {
  ns <- NS(id)
  m <- data$manifest
  assays <- .scroll_assays_of(m); cats <- .scroll_cat_cols(m)
  if (!length(cats)) return(.scroll_empty_panel("Needs a categorical metadata column (none found)."))
  bslib::layout_columns(
    col_widths = c(3, 9), class = "scroll-panel",
    div(
      class = "scroll-controls",
      .scroll_group("Genes",
        selectizeInput(ns("markers"), "Marker genes", choices = NULL, multiple = TRUE,
                       options = list(placeholder = "Add genes...", maxOptions = 50,
                                      plugins = list("remove_button")))),
      .scroll_group("Grouping",
        selectInput(ns("group"), "Group by", stats::setNames(cats, cats), selected = cats[[1]]),
        if (length(assays) > 1)
          selectInput(ns("assay"), "Assay", assays, selected = m$default_assay)),
      .scroll_group("Appearance",
        bslib::input_switch(ns("scale"), "Scale expression (z-score)", TRUE),
        selectInput(ns("palette"), "Palette", .scroll_continuous_palettes, selected = "magma"),
        sliderInput(ns("dotrange"), "Dot size", 0, 10, c(1, 6), 0.5)),
      .scroll_group("Layout",
        selectInput(ns("cluster"), "Cluster (hclust)",
                    c("Off" = "off", "Rows" = "rows", "Columns" = "columns", "Both" = "both")),
        .scroll_aspect_input(ns))
    ),
    .scroll_plot_area(ns, "520px")
  )
}

dotplot_server <- function(id, data, cells_r = reactive(data$cells),
                           view_r = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    m <- data$manifest
    .scroll_bind_view_cats(input, session, view_r, m, "group")
    assay <- reactive(input$assay %||% m$default_assay)
    defaults <- intersect(unlist(data$config$markers), .scroll_features_of(m, m$default_assay))
    # repopulate markers for the active assay, keeping only those present in it
    observeEvent(assay(), {
      feats <- .scroll_features_of(m, assay())
      cur <- isolate(input$markers) %||% defaults
      updateSelectizeInput(session, "markers", choices = feats, server = TRUE,
                           selected = intersect(cur, feats))
    })
    # DATA reactive: query + aggregation + hclust. `scale` and `cluster` change
    # the aggregation/clustering, so they are DATA inputs (not cosmetic); palette
    # and dot size are cosmetic. No rasterization (dots = features x groups).
    data_r <- reactive({
      req(input$group)
      feats <- input$markers
      validate(need(length(feats) > 0, "Add one or more marker genes to build the panel."))
      cells <- cells_r()
      list(cells = cells, group_by = input$group, features = feats,
           assembly = .scroll_dotplot_assemble(
             cells, feats, input$group, data$queryN(assay(), feats),
             scale = isTRUE(input$scale), cluster = input$cluster))
    })
    cosmetic_r <- .scroll_cosmetic(reactive(
      list(palette = input$palette, dot_size = input$dotrange, aspect = input$aspect)))
    plot_r <- reactive({
      d <- data_r()
      view_dotplot(d$cells, list(group_by = d$group_by, features = d$features),
                   NULL, cosmetic_r(), assembly = d$assembly)
    })
    output$plot <- renderPlot(plot_r())
    .scroll_plot_downloads(output, plot_r, id)
  })
}

# --- Violin panel -------------------------------------------------------------

violin_ui <- function(id, data) {
  ns <- NS(id)
  m <- data$manifest
  assays <- .scroll_assays_of(m); cats <- .scroll_cat_cols(m)
  if (!length(cats)) return(.scroll_empty_panel("Needs a categorical metadata column (none found)."))
  bslib::layout_columns(
    col_widths = c(3, 9), class = "scroll-panel",
    div(
      class = "scroll-controls",
      .scroll_group("Feature",
        selectizeInput(ns("feature"), "Gene", choices = NULL, multiple = FALSE,
                       options = list(placeholder = "Search a gene...", maxOptions = 50)),
        if (length(assays) > 1)
          selectInput(ns("assay"), "Assay", assays, selected = m$default_assay)),
      .scroll_group("Grouping",
        selectInput(ns("group"), "Group by", stats::setNames(cats, cats), selected = cats[[1]])),
      .scroll_group("Appearance",
        selectInput(ns("palette"), "Palette", .scroll_cat_palettes()),
        uiOutput(ns("manual")),
        bslib::input_switch(ns("jitter"), "Show points", FALSE),
        bslib::input_switch(ns("legend"), "Legend", FALSE)),
      .scroll_group("Layout", .scroll_aspect_input(ns))
    ),
    .scroll_plot_area(ns, "460px")
  )
}

violin_server <- function(id, data, cells_r = reactive(data$cells),
                          view_r = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    m <- data$manifest
    .scroll_bind_view_cats(input, session, view_r, m, "group")
    assay <- reactive(input$assay %||% m$default_assay)
    # repopulate the gene list for the active assay; drop a selection that does
    # not exist in the newly chosen assay (else it silently queries empty)
    observeEvent(assay(), {
      feats <- .scroll_features_of(m, assay())
      cur <- isolate(input$feature)
      keep <- if (!is.null(cur) && cur %in% feats) cur else character(0)
      updateSelectizeInput(session, "feature", choices = feats, server = TRUE, selected = keep)
    })
    # per-group color pickers when the palette is "Manual" (levels of Group by)
    lvl_r <- reactive({ req(input$group); unlist(m$meta[[input$group]]$levels) })
    output$manual <- renderUI(
      if (identical(input$palette, "Manual")) .scroll_manual_ui(session$ns, lvl_r()))
    manual_colors <- reactive(
      if (identical(input$palette, "Manual")) .scroll_manual_colors(input, lvl_r()))
    # DATA reactive: cells + query (cosmetic changes no longer re-hit the store).
    data_r <- reactive({
      req(input$group)
      feat <- .scroll_nz(input$feature)
      validate(need(!is.null(feat), "Search for a gene to plot its distribution."))
      cells <- cells_r()
      list(cells = cells, feature = feat, group_by = input$group,
           values = data$query1(assay(), feat))
    })
    cosmetic_r <- .scroll_cosmetic(reactive(
      list(palette = input$palette, jitter = isTRUE(input$jitter),
           legend = isTRUE(input$legend), aspect = input$aspect,
           manual_colors = manual_colors())))
    plot_r <- reactive({
      d <- data_r()
      view_violin(d$cells, list(feature = d$feature, group_by = d$group_by),
                  d$values, cosmetic_r())
    })
    output$plot <- renderPlot(plot_r())
    .scroll_plot_downloads(output, plot_r, id)
  })
}

# --- Proportions panel --------------------------------------------------------

proportions_ui <- function(id, data) {
  ns <- NS(id)
  cats <- .scroll_cat_cols(data$manifest)
  if (!length(cats)) return(.scroll_empty_panel("Needs categorical metadata columns (none found)."))
  x_default <- if (length(cats) >= 2) cats[[2]] else cats[[1]]
  bslib::layout_columns(
    col_widths = c(3, 9), class = "scroll-panel",
    div(
      class = "scroll-controls",
      .scroll_group("Composition",
        selectInput(ns("group"), "Group by (x)", stats::setNames(cats, cats), selected = x_default),
        selectInput(ns("fill"), "Fill by", stats::setNames(cats, cats), selected = cats[[1]])),
      .scroll_group("Appearance",
        selectInput(ns("palette"), "Palette", .scroll_cat_palettes()),
        uiOutput(ns("manual")),
        bslib::input_switch(ns("normalize"), "Normalize to 100%", TRUE),
        bslib::input_switch(ns("legend"), "Legend", TRUE)),
      .scroll_group("Layout", .scroll_aspect_input(ns))
    ),
    .scroll_plot_area(ns, "460px")
  )
}

proportions_server <- function(id, data, cells_r = reactive(data$cells),
                               view_r = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    m <- data$manifest
    .scroll_bind_view_cats(input, session, view_r, m, c("group", "fill"))
    # per-fill color pickers when the palette is "Manual" (levels of Fill by)
    lvl_r <- reactive({ req(input$fill); unlist(m$meta[[input$fill]]$levels) })
    output$manual <- renderUI(
      if (identical(input$palette, "Manual")) .scroll_manual_ui(session$ns, lvl_r()))
    manual_colors <- reactive(
      if (identical(input$palette, "Manual")) .scroll_manual_colors(input, lvl_r()))
    data_r <- reactive({
      req(input$group, input$fill)
      list(cells = cells_r(), group_by = input$group, fill_by = input$fill)
    })
    cosmetic_r <- .scroll_cosmetic(reactive(
      list(palette = input$palette, normalize = isTRUE(input$normalize),
           legend = isTRUE(input$legend), aspect = input$aspect,
           manual_colors = manual_colors())))
    plot_r <- reactive({
      d <- data_r()
      view_proportions(d$cells, list(group_by = d$group_by, fill_by = d$fill_by), cosmetic_r())
    })
    output$plot <- renderPlot(plot_r())
    .scroll_plot_downloads(output, plot_r, id)
  })
}

# --- Biaxial panel ------------------------------------------------------------

# Prefer hashtag / antibody CLR columns for the default axes, else the first few
# numerics — so an HTO or CITE-seq dataset opens on its biaxial signal plots.
.scroll_default_biaxial <- function(nums) {
  hit <- grep("^hto_|hashtag|adt_|_adt$|^ab_", nums, value = TRUE, ignore.case = TRUE)
  if (length(hit) >= 2) hit else utils::head(nums, 3)
}

biaxial_ui <- function(id, data) {
  ns <- NS(id)
  m <- data$manifest
  nums <- .scroll_num_cols(m); cats <- .scroll_cat_cols(m)
  if (length(nums) < 2)
    return(.scroll_empty_panel("Needs at least two numeric metadata columns (none found)."))
  if (!length(cats))
    return(.scroll_empty_panel("Needs a categorical column to colour by (none found)."))
  color_default <- if ("hto" %in% cats) "hto" else cats[[1]]
  bslib::layout_columns(
    col_widths = c(3, 9), class = "scroll-panel",
    div(
      class = "scroll-controls",
      .scroll_group("Axes",
        selectizeInput(ns("features"), "Numeric columns", choices = nums,
                       selected = .scroll_default_biaxial(nums), multiple = TRUE,
                       options = list(placeholder = "Pick 2+ numeric columns"))),
      .scroll_group("Colour",
        selectInput(ns("colorby"), "Colour by", stats::setNames(cats, cats),
                    selected = color_default)),
      .scroll_group("Appearance",
        selectInput(ns("palette"), "Palette", .scroll_cat_palettes()),
        uiOutput(ns("manual")),
        sliderInput(ns("size"), "Point size", 0.1, 3, 0.5, 0.1),
        sliderInput(ns("alpha"), "Opacity", 0.1, 1, 0.6, 0.05),
        bslib::input_switch(ns("legend"), "Legend", TRUE),
        bslib::input_switch(ns("raster"), "Rasterize (fast)", TRUE)),
      .scroll_group("Layout", .scroll_aspect_input(ns))
    ),
    .scroll_plot_area(ns, "460px")
  )
}

biaxial_server <- function(id, data, cells_r = reactive(data$cells),
                           view_r = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    m <- data$manifest
    .scroll_bind_view_cats(input, session, view_r, m, "colorby")
    # per-level color pickers when the palette is "Manual" (levels of Colour by)
    lvl_r <- reactive({ req(input$colorby); unlist(m$meta[[input$colorby]]$levels) })
    output$manual <- renderUI(
      if (identical(input$palette, "Manual")) .scroll_manual_ui(session$ns, lvl_r()))
    manual_colors <- reactive(
      if (identical(input$palette, "Manual")) .scroll_manual_colors(input, lvl_r()))
    # DATA reactive: build the (potentially large) pairwise long df once; cosmetic
    # drags no longer re-expand it. `params` is carried for the colour label.
    data_r <- reactive({
      req(input$colorby)
      feats <- input$features
      validate(need(length(feats) >= 2, "Pick at least two numeric columns."))
      params <- list(features = feats, color_by = input$colorby)
      list(params = params, df = .scroll_biaxial_df(cells_r(), params))
    })
    cosmetic_r <- .scroll_cosmetic(reactive(
      list(palette = input$palette, point_size = input$size,
           alpha = input$alpha, legend = isTRUE(input$legend), aspect = input$aspect,
           manual_colors = manual_colors())))
    build <- function(raster) {
      d <- data_r(); st <- cosmetic_r(); st$raster <- raster
      view_biaxial(NULL, d$params, st, df = d$df)
    }
    plot_r   <- reactive(build(.scroll_use_raster(input$raster, nrow(data_r()$df))))
    export_r <- reactive(build(FALSE))
    output$plot <- renderPlot(plot_r())
    .scroll_plot_downloads(output, export_r, id)
  })
}

# --- DE panel -----------------------------------------------------------------

# One contrast, computed once via presto, shown two ways: a ranked marker
# Table and a Volcano. A single `result` reactive feeds both tabs; the table's
# min-%-expressing and top-N are display filters (so the volcano keeps every
# gene), and the volcano thresholds only restyle the plot.
de_ui <- function(id, data) {
  ns <- NS(id)
  m <- data$manifest
  cats <- .scroll_cat_cols(m); assays <- .scroll_assays_of(m)
  if (!length(cats)) return(.scroll_empty_panel("Needs a categorical metadata column (none found)."))
  bslib::layout_columns(
    col_widths = c(3, 9), class = "scroll-panel",
    div(
      class = "scroll-controls",
      .scroll_group("Contrast",
        selectInput(ns("group"), "Group by", stats::setNames(cats, cats), selected = cats[[1]]),
        selectizeInput(ns("ident1"), "Group 1", choices = NULL, multiple = TRUE,
                       options = list(plugins = list("remove_button"))),
        selectizeInput(ns("ident2"), "vs.", choices = NULL, multiple = TRUE,
                       options = list(plugins = list("remove_button"),
                                      placeholder = "rest (all other cells)")),
        if (length(assays) > 1) selectInput(ns("assay"), "Assay", assays, selected = m$default_assay)),
      .scroll_group("Table",
        sliderInput(ns("minpct"), "Min % expressing", 0, 50, 10, 1),
        sliderInput(ns("topn"), "Show top", 10, 300, 50, 10),
        bslib::input_switch(ns("export_all"), "Export all genes (CSV)", FALSE)),
      .scroll_group("Volcano",
        sliderInput(ns("lfc"), "logFC cutoff", 0, 3, 1, 0.1),
        numericInput(ns("padj"), "Adj. p cutoff", 0.05, min = 0, max = 1, step = 0.01),
        sliderInput(ns("labeln"), "Label top", 0, 40, 15, 1),
        .scroll_aspect_input(ns)),
      actionButton(ns("compute"), "Compute DE", class = "btn-primary", width = "100%")
    ),
    div(class = "scroll-plot",
        bslib::navset_tab(
          bslib::nav_panel("Table",
            div(class = "scroll-plot-bar", .scroll_dl_button(ns("csv"), "CSV")),
            div(class = "scroll-table",
                .scroll_spin(if (.scroll_has_dt()) DT::dataTableOutput(ns("table"))
                             else tableOutput(ns("table"))))),
          bslib::nav_panel("Volcano",
            div(class = "scroll-plot-bar",
                .scroll_dl_button(ns("png"), "PNG"),
                .scroll_dl_button(ns("pdf"), "PDF")),
            .scroll_spin(plotOutput(ns("plot"), height = "520px")))))
  )
}

.scroll_has_dt <- function() requireNamespace("DT", quietly = TRUE)

de_server <- function(id, data, cells_r = reactive(data$cells),
                      view_r = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    m <- data$manifest
    assay <- reactive(input$assay %||% m$default_assay)
    .scroll_bind_view_cats(input, session, view_r, m, "group")

    observeEvent(input$group, {
      lv <- unlist(m$meta[[input$group]]$levels)
      updateSelectizeInput(session, "ident1", choices = lv, selected = lv[[1]])
      updateSelectizeInput(session, "ident2", choices = lv, selected = character(0))
    })

    # min_pct = 0 so the volcano keeps every gene; the table applies its own
    # %-expressing filter below. A failed presto run (e.g. a contrast with too
    # few cells) is caught and surfaced as a friendly inline message, not a raw
    # Shiny error.
    result <- eventReactive(input$compute, {
      req(input$ident1)
      tryCatch(
        list(ok = scroll_de(data, assay(), input$group, input$ident1,
                            if (length(input$ident2)) input$ident2 else NULL,
                            min_pct = 0, cells = cells_r())),
        error = function(e) list(err = conditionMessage(e)))
    })

    de_df <- reactive({
      validate(need(input$compute > 0, "Pick a contrast and click Compute DE."))
      r <- result()
      validate(need(is.null(r$err), r$err))
      r$ok
    })

    table_rows <- function() {
      res <- de_df()
      res <- res[pmax(res$pct.1, res$pct.2) >= input$minpct / 100, , drop = FALSE]
      utils::head(res, input$topn)
    }
    if (.scroll_has_dt())
      output$table <- DT::renderDataTable(
        DT::formatSignif(
          DT::datatable(table_rows(), rownames = FALSE, options = list(pageLength = 15, dom = "tip")),
          columns = c("p_val", "p_val_adj"), digits = 3))
    else
      output$table <- renderTable(table_rows())

    volcano_r <- reactive({
      view_volcano(de_df(),
                   params = list(lfc = input$lfc, padj = input$padj, label_n = input$labeln),
                   state = list(aspect = input$aspect))
    })
    output$plot <- renderPlot(volcano_r())
    .scroll_plot_downloads(output, volcano_r, id)
    output$csv <- .scroll_csv_handler(
      reactive(if (isTRUE(input$export_all)) de_df() else table_rows()),
      paste0("scroll_", id, ".csv"))
  })
}

# --- Pseudobulk DE panel ------------------------------------------------------

# Combined-interaction levels observed in `cells` for the given columns.
.scroll_combo_choices <- function(cells, cols) {
  if (!length(cols)) return(character(0))
  combo <- .scroll_combo_levels(cells, cols)
  sort(unique(combo[!is.na(combo)]))
}

pseudobulk_de_ui <- function(id, data) {
  ns <- NS(id)
  m <- data$manifest
  cats <- .scroll_cat_cols(m); assays <- .scroll_assays_of(m)
  if (!isTRUE(m$has_counts))
    return(.scroll_empty_panel(paste(
      "Pseudobulk DE needs a counts store. Rebuild the project with",
      "scroll_build(..., counts = TRUE).")))
  if (!length(cats)) return(.scroll_empty_panel("Needs a categorical metadata column (none found)."))
  bslib::layout_columns(
    col_widths = c(3, 9), class = "scroll-panel",
    div(
      class = "scroll-controls",
      .scroll_group("Contrast",
        selectizeInput(ns("aggregate_by"), "Aggregate by", stats::setNames(cats, cats),
                       multiple = TRUE, selected = cats[[1]],
                       options = list(plugins = list("remove_button"))),
        selectizeInput(ns("ident1"), "Group 1", choices = NULL, multiple = TRUE,
                       options = list(plugins = list("remove_button"))),
        selectizeInput(ns("ident2"), "vs.", choices = NULL, multiple = TRUE,
                       options = list(plugins = list("remove_button"),
                                      placeholder = "rest (all other cells)")),
        selectInput(ns("replicate"), "Replicate",
                    c("No replicate (pseudo)" = "no_replicate", stats::setNames(cats, cats))),
        if (length(assays) > 1) selectInput(ns("assay"), "Assay", assays, selected = m$default_assay)),
      .scroll_group("Pseudobulk",
        sliderInput(ns("mincells"), "Min cells / sample", 3, 200, 10, 1),
        numericInput(ns("npseudo"), "Pseudo-reps (if no replicate)", 3, min = 2, max = 10),
        numericInput(ns("cellsper"), "Cells / pseudo-rep", 50, min = 5, max = 2000)),
      .scroll_group("Stability",
        numericInput(ns("runs"), "Runs (re-sample; pseudo only)", 1, min = 1, max = 100),
        sliderInput(ns("stabcut"), "Consistency cutoff", 0, 1, 0.8, 0.05)),
      .scroll_group("Table",
        sliderInput(ns("topn"), "Show top", 10, 300, 50, 10),
        bslib::input_switch(ns("export_all"), "Export all genes (CSV)", FALSE)),
      .scroll_group("Volcano / stability",
        sliderInput(ns("lfc"), "logFC cutoff", 0, 3, 1, 0.1),
        numericInput(ns("padj"), "Adj. p cutoff", 0.05, min = 0, max = 1, step = 0.01),
        sliderInput(ns("labeln"), "Label top", 0, 40, 15, 1),
        .scroll_aspect_input(ns)),
      actionButton(ns("compute"), "Compute pseudobulk DE", class = "btn-primary", width = "100%")
    ),
    div(class = "scroll-plot",
        uiOutput(ns("note")),
        bslib::navset_tab(
          bslib::nav_panel("Table",
            div(class = "scroll-plot-bar", .scroll_dl_button(ns("csv"), "CSV")),
            div(class = "scroll-table",
                .scroll_spin(if (.scroll_has_dt()) DT::dataTableOutput(ns("table"))
                             else tableOutput(ns("table"))))),
          bslib::nav_panel("Plot",
            div(class = "scroll-plot-bar",
                .scroll_dl_button(ns("png"), "PNG"), .scroll_dl_button(ns("pdf"), "PDF")),
            .scroll_spin(plotOutput(ns("plot"), height = "520px")))))
  )
}

pseudobulk_de_server <- function(id, data, cells_r = reactive(data$cells)) {
  moduleServer(id, function(input, output, session) {
    m <- data$manifest
    assay <- reactive(input$assay %||% m$default_assay)

    # combined levels depend on the chosen aggregate-by columns + active subset
    observeEvent(list(input$aggregate_by, cells_r()), {
      combos <- .scroll_combo_choices(cells_r(), input$aggregate_by)
      updateSelectizeInput(session, "ident1", choices = combos,
                           selected = intersect(isolate(input$ident1), combos))
      updateSelectizeInput(session, "ident2", choices = combos,
                           selected = intersect(isolate(input$ident2), combos))
    })

    # stability = re-run the pseudo-replication K times (pseudo mode only). The
    # K runs are computed once here; de_df() aggregates them cheaply when the
    # lfc/padj cutoffs change.
    result <- eventReactive(input$compute, {
      req(input$aggregate_by, input$ident1)
      args <- list(
        data, assay(), aggregate_cols = input$aggregate_by, ident1 = input$ident1,
        ident2 = if (length(input$ident2)) input$ident2 else NULL,
        replicate_col = input$replicate, min_cells = input$mincells,
        n_pseudo = input$npseudo, cells_per_pseudo = input$cellsper, cells = cells_r())
      stability <- isTRUE(input$runs > 1) && identical(input$replicate, "no_replicate")
      tryCatch(
        if (stability)
          list(runs = do.call(.scroll_pseudobulk_runs, c(args, list(runs = input$runs))))
        else list(ok = do.call(scroll_pseudobulk_de, args)),
        error = function(e) list(err = conditionMessage(e)))
    })

    de_df <- reactive({
      validate(need(input$compute > 0, "Pick a contrast and click Compute pseudobulk DE."))
      r <- result()
      validate(need(is.null(r$err), r$err))
      if (!is.null(r$runs)) .scroll_stability_aggregate(r$runs, input$lfc, input$padj)
      else r$ok
    })

    output$note <- renderUI({
      req(input$compute > 0); r <- result(); if (!is.null(r$err)) return(NULL)
      mk <- function(txt) div(class = "scroll-desc", style = "margin:0 0 8px; color:#B45309", txt)
      if (!is.null(r$runs))
        mk(sprintf(paste("Stability across %d pseudo-replicate runs - sel_freq is a",
                         "robustness heuristic, not a p-value; prefer real replicates."),
                   length(r$runs)))
      else if (isTRUE(attr(r$ok, "pseudo"))) {
        n <- attr(r$ok, "n_samples")
        mk(sprintf(paste("Pseudo-replicates used (no biological replicates): %d vs %d",
                         "samples - treat p-values with caution."), n[["group1"]], n[["group2"]]))
      } else NULL
    })

    # display columns differ between single-run and stability results
    table_rows <- function() {
      d <- de_df()
      if ("sel_freq" %in% names(d))
        d <- d[, c("gene", "sel_freq", "median_logFC", "sign_agree", "median_padj", "n_tested")]
      utils::head(d, input$topn)
    }
    if (.scroll_has_dt())
      output$table <- DT::renderDataTable({
        d <- table_rows()
        DT::formatSignif(
          DT::datatable(d, rownames = FALSE, options = list(pageLength = 15, dom = "tip")),
          columns = intersect(c("p_val", "p_val_adj", "median_padj"), names(d)), digits = 3)
      })
    else
      output$table <- renderTable(table_rows())

    plot_r <- reactive({
      d <- de_df()
      if ("sel_freq" %in% names(d))
        view_stability(d, params = list(lfc = input$lfc, cut = input$stabcut,
                                        label_n = input$labeln), state = list(aspect = input$aspect))
      else
        view_volcano(d, params = list(lfc = input$lfc, padj = input$padj,
                                      label_n = input$labeln), state = list(aspect = input$aspect))
    })
    output$plot <- renderPlot(plot_r())
    .scroll_plot_downloads(output, plot_r, id)
    output$csv <- .scroll_csv_handler(
      reactive(if (isTRUE(input$export_all)) de_df() else table_rows()),
      paste0("scroll_", id, ".csv"))
  })
}

# --- section registry + shell -------------------------------------------------

# The built-in panels, in scroll order. Each entry: display meta + its module's
# ui/server. Display numbers are assigned at assembly time (see
# .scroll_assemble_panels), so appended custom panels number correctly.
.scroll_builtin_panels <- function() c(list(
  list(id = "dimplot", label = "DimPlot",
       title = "Cells by annotation",
       desc = "The embedding coloured by any cell metadata column.",
       ui = dimplot_ui, server = dimplot_server),
  list(id = "featureplot", label = "FeaturePlot",
       title = "Gene expression",
       desc = "The embedding coloured by a gene's expression.",
       ui = featureplot_ui, server = featureplot_server),
  list(id = "biaxial", label = "Biaxial",
       title = "Biaxial signal",
       desc = "Pairwise scatters of numeric columns (e.g. hashtags) coloured by a selection.",
       when = function(m) length(.scroll_num_cols(m)) >= 2,
       ui = biaxial_ui, server = biaxial_server),
  list(id = "dotplot", label = "DotPlot",
       title = "Marker panel",
       desc = "Mean expression and fraction expressing across groups.",
       when = function(m) length(.scroll_cat_cols(m)) >= 1,
       ui = dotplot_ui, server = dotplot_server),
  list(id = "violin", label = "Violin",
       title = "Expression distribution",
       desc = "A gene's per-group expression distribution.",
       when = function(m) length(.scroll_cat_cols(m)) >= 1,
       ui = violin_ui, server = violin_server),
  list(id = "proportions", label = "Proportions",
       title = "Composition",
       desc = "Stacked composition of one annotation within another.",
       when = function(m) length(.scroll_cat_cols(m)) >= 2,
       ui = proportions_ui, server = proportions_server),
  list(id = "de", label = "DE",
       title = "Differential expression",
       desc = "Wilcoxon markers for a contrast (live via presto): ranked table + volcano.",
       when = .scroll_has_contrast,
       ui = de_ui, server = de_server),
  list(id = "pseudobulk", label = "Pseudobulk DE",
       title = "Pseudobulk differential expression",
       desc = "Aggregate cells into pseudobulk samples and test with edgeR/limma-voom.",
       when = function(m) isTRUE(m$has_counts) && .scroll_has_contrast(m),
       ui = pseudobulk_de_ui, server = pseudobulk_de_server)
),
  # Modality panels: each carries a `when(manifest)` predicate and only surfaces
  # for projects whose manifest has the matching block (see .scroll_assemble_panels).
  .scroll_vdj_panels(), .scroll_spatial_panels(), .scroll_atac_panels())

# Mutable registry of user-added panels (session-global, like knitr's engines).
.scroll_registry <- new.env(parent = emptyenv())
.scroll_registry$panels <- list()

#' Register a custom panel in the scroll explorer
#'
#' Adds a user-defined analysis panel to every [scroll_app()] built afterwards in
#' this session. A panel follows the same module contract as the built-ins — a
#' pair of functions:
#'
#' * `ui(id, data)` returns a UI tag; namespace its inputs with [shiny::NS()]`(id)`.
#' * `server(id, data, cells_r)` wires the module, typically via
#'   [shiny::moduleServer()]. `cells_r()` is a reactive of the *active* cells
#'   (already narrowed by the app-bar subset filter); use it, not `data$cells`,
#'   so the panel honours the global filter.
#'
#' `data` is the shared handle: `data$cells` (a data.frame of metadata +
#' embedding coordinates), `data$manifest`, `data$config`, and the query helpers
#' `data$query1(assay, feature)` / `data$queryN(assay, features)`, which return
#' dequantized long expression (`cell`, `value` / `feature`, `cell`, `value`).
#'
#' Registering an `id` that matches an existing panel — built-in or custom —
#' replaces it in place, so you can override a built-in. Clear all custom panels
#' with [scroll_reset_panels()].
#'
#' @param id Unique panel id: a letter followed by letters, digits, or
#'   underscores. Used as the section anchor and the Shiny module namespace.
#' @param ui,server The panel's UI and server functions (see contract above).
#' @param label Short rail label (defaults to `id`).
#' @param title,desc Section heading and one-line description.
#' @param after Id of the panel to insert this one after; `NULL` (default)
#'   appends at the end. Ignored when replacing an existing id.
#' @param before Id of the panel to insert this one before (e.g. `"dimplot"` to
#'   put it at the top). Takes precedence over `after`. Ignored when replacing an
#'   existing id.
#' @return Invisibly, `id`.
#' @examples
#' # A minimal custom panel: a cells-per-group bar chart.
#' count_ui <- function(id, data) {
#'   ns <- shiny::NS(id)
#'   cols <- names(Filter(function(x) identical(x$type, "categorical"),
#'                        data$manifest$meta))
#'   shiny::tagList(shiny::selectInput(ns("grp"), "Group", cols),
#'                  shiny::plotOutput(ns("plot")))
#' }
#' count_server <- function(id, data, cells_r = shiny::reactive(data$cells)) {
#'   shiny::moduleServer(id, function(input, output, session) {
#'     output$plot <- shiny::renderPlot({
#'       shiny::req(input$grp)
#'       barplot(table(cells_r()[[input$grp]]))
#'     })
#'   })
#' }
#' register_panel("counts", count_ui, count_server, label = "Counts",
#'                title = "Cells per group")
#' scroll_reset_panels()   # (undo, so the example leaves no state)
#' @export
register_panel <- function(id, ui, server, label = id, title = label,
                           desc = NULL, after = NULL, before = NULL) {
  .scroll_check_panel_id(id)
  if (!is.function(ui) || !is.function(server))
    stop("`ui` and `server` must be functions.", call. = FALSE)
  spec <- list(id = id, label = label, title = title, desc = desc %||% "",
               ui = ui, server = server, after = after, before = before)
  reg <- .scroll_registry$panels
  ids <- vapply(reg, `[[`, "", "id")
  reg[[if (id %in% ids) which(ids == id) else length(reg) + 1L]] <- spec
  .scroll_registry$panels <- reg
  invisible(id)
}

#' Clear all custom panels registered with [register_panel()]
#'
#' @return Invisibly, `NULL`.
#' @export
scroll_reset_panels <- function() {
  .scroll_registry$panels <- list()
  invisible(NULL)
}

.scroll_check_panel_id <- function(id) {
  if (!is.character(id) || length(id) != 1L || is.na(id) ||
      !grepl("^[A-Za-z][A-Za-z0-9_]*$", id))
    stop("`id` must be a single name: a letter, then letters/digits/underscores.",
         call. = FALSE)
  invisible(id)
}

# Built-ins + registered panels, in final scroll order, each stamped with a
# display number by position. A registered id matching an existing panel replaces
# it in place; otherwise it is inserted before `before`, else after `after`, else
# appended.
.scroll_assemble_panels <- function(manifests = NULL) {
  panels <- .scroll_builtin_panels()
  for (spec in .scroll_registry$panels) {
    ids <- vapply(panels, `[[`, "", "id")
    if (spec$id %in% ids) {
      panels[[which(ids == spec$id)]] <- spec
    } else if (!is.null(spec$before) && spec$before %in% ids) {
      panels <- append(panels, list(spec), after = which(ids == spec$before) - 1L)
    } else if (!is.null(spec$after) && spec$after %in% ids) {
      panels <- append(panels, list(spec), after = which(ids == spec$after))
    } else {
      panels <- c(panels, list(spec))
    }
  }
  # Modality gating: a panel may declare `when = function(manifest)`; keep it only
  # if some supplied manifest satisfies it (so e.g. VDJ panels appear only for
  # projects built with `vdj=`). `manifests` is one manifest or a list of them
  # (scroll_multi_app); NULL leaves every panel in (used by tests).
  if (!is.null(manifests)) {
    if (!is.null(manifests$scroll_version)) manifests <- list(manifests)
    keep <- vapply(panels, function(p) is.null(p$when) ||
      any(vapply(manifests, function(m) isTRUE(tryCatch(p$when(m), error = function(e) FALSE)),
                 logical(1))), logical(1))
    panels <- panels[keep]
  }
  for (i in seq_along(panels)) panels[[i]]$num <- sprintf("%02d", i)
  panels
}

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

# Standard plot area: a small toolbar (PNG + PDF download) above the plot output.
.scroll_plot_area <- function(ns, height = "460px")
  div(class = "scroll-plot",
      div(class = "scroll-plot-bar",
          .scroll_dl_button(ns("png"), "PNG"),
          .scroll_dl_button(ns("pdf"), "PDF")),
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
# reactive on a module's `output`. Filenames stem from the panel id.
.scroll_plot_downloads <- function(output, plot_r, id) {
  output$png <- .scroll_img_handler(plot_r, paste0("scroll_", id, ".png"), "png")
  output$pdf <- .scroll_img_handler(plot_r, paste0("scroll_", id, ".pdf"), "pdf")
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

# The global cell-subset filter: restrict `cells` to rows whose `col` value is in
# `vals`; no-op when the filter is inactive.
.scroll_subset_cells <- function(cells, col, vals) {
  if (is.null(col) || is.null(vals) || !length(vals) || !col %in% names(cells))
    return(cells)
  cells[as.character(cells[[col]]) %in% vals, , drop = FALSE]
}

.scroll_stat <- function(value, label)
  div(class = "scroll-stat", span(class = "scroll-stat-v", value),
      span(class = "scroll-stat-l", label))

.scroll_appbar <- function(data, title, ns = identity) {
  m <- data$manifest
  assay <- m$default_assay
  cats <- .scroll_cat_cols(m)
  # brand wordmark; `brand:` in config.yaml overrides the default "scroll"
  brand <- list(span(class = "scroll-logo", .scroll_nz(data$config$brand) %||% "scroll"))
  if (!is.null(title))
    brand <- c(brand, list(span(class = "scroll-slash", "/"),
                           span(class = "scroll-dataset", title)))
  # Reprocessed subset views: pick a linked view (restricts every panel + swaps
  # to the subset's embedding). Shown only when the build declared `subsets`.
  subs <- .scroll_subset_names(m)
  view_ui <- if (length(subs)) div(
    class = "scroll-subset",
    span(class = "scroll-subset-label", "View"),
    selectInput(ns("scroll_view"), NULL,
                c("Whole dataset" = "",
                  stats::setNames(subs, vapply(subs, function(s) m$subsets[[s]]$label, ""))),
                width = "160px"))
  # Ad-hoc categorical subset filter. Off by default (subset views cover the
  # common case); set `subset_filter: true` in config.yaml to re-enable.
  subset_ui <- if (isTRUE(data$config$subset_filter) && length(cats)) div(
    class = "scroll-subset",
    span(class = "scroll-subset-label", "Subset"),
    selectInput(ns("scroll_subset_col"), NULL,
                c("All cells" = "", stats::setNames(cats, cats)), width = "150px"),
    selectizeInput(ns("scroll_subset_val"), NULL, choices = NULL, multiple = TRUE,
                   width = "200px", options = list(placeholder = "all")))
  div(
    class = "scroll-appbar",
    div(class = "scroll-brand", brand),
    view_ui,
    subset_ui,
    div(class = "scroll-stats",
        .scroll_stat(textOutput(ns("scroll_ncells"), inline = TRUE), "cells"),
        .scroll_stat(format(m$assays[[assay]]$n_features, big.mark = ","), "genes"),
        .scroll_stat(paste(.scroll_assays_of(m), collapse = ", "), "assays"),
        .scroll_stat(paste(.scroll_reductions(m), collapse = ", "), "reductions"))
  )
}

.scroll_rail <- function(panels, ns = identity) {
  tags$nav(
    class = "scroll-rail",
    lapply(panels, function(s) tags$a(
      class = "scroll-rail-item", href = paste0("#", ns(s$id)),
      span(class = "scroll-rail-num", s$num), span(s$label))),
    div(class = "scroll-rail-foot", "auto-generated from manifest.yaml")
  )
}

.scroll_section_card <- function(sec, data, ns = identity) {
  bslib::card(
    id = ns(sec$id), class = "scroll-section", full_screen = FALSE,
    bslib::card_header(
      div(class = "scroll-eyebrow",
          span(class = "scroll-num", sec$num), span(class = "scroll-kicker", sec$label)),
      tags$h2(class = "scroll-title", sec$title),
      tags$p(class = "scroll-desc", sec$desc)),
    sec$ui(ns(sec$id), data)
  )
}

# The per-dataset body (app bar + rail + section cards), with every id passed
# through `ns` so multiple datasets can coexist in one page (see scroll_multi_app).
# Does NOT include page_fluid/theme/head -- those wrap it once at the page level.
.scroll_body <- function(data, title, panels, ns = identity) {
  tagList(
    .scroll_appbar(data, title, ns),
    div(
      class = "scroll-layout",
      .scroll_rail(panels, ns),
      div(class = "scroll-content",
          lapply(panels, function(s) .scroll_section_card(s, data, ns)))
    )
  )
}

.scroll_page <- function(data, title, panels) {
  bslib::page_fluid(
    theme = .scroll_theme(),
    tags$head(tags$style(HTML(.scroll_css())), tags$script(HTML(.scroll_spy_js()))),
    .scroll_body(data, title, panels)
  )
}

# The per-dataset server logic (active view/subset composition + panel wiring).
# Called flat by scroll_app() and inside a moduleServer() by scroll_multi_app(),
# so `input`/`output`/`session` carry the right namespace in both.
.scroll_wire <- function(input, output, session, data, panels) {
  m <- data$manifest
  active_view <- reactive(.scroll_nz(input$scroll_view))
  active_cells <- reactive({
    base <- .scroll_view_cells(data$cells, m, active_view())
    .scroll_subset_cells(base, .scroll_nz(input$scroll_subset_col),
                         input$scroll_subset_val)
  })
  observeEvent(input$scroll_subset_col, {
    col <- .scroll_nz(input$scroll_subset_col)
    lv <- if (is.null(col)) character(0) else unlist(m$meta[[col]]$levels)
    updateSelectizeInput(session, "scroll_subset_val", choices = lv,
                         selected = character(0), server = TRUE)
  })
  output$scroll_ncells <- renderText({
    n <- nrow(active_cells()); tot <- m$n_cells
    lab <- if (!is.null(active_view())) m$subsets[[active_view()]]$label
    base <- if (n < tot) sprintf("%s of %s", format(n, big.mark = ","),
                                 format(tot, big.mark = ",")) else format(tot, big.mark = ",")
    if (!is.null(lab)) paste0(base, " \u00b7 ", lab) else base
  })
  # Pass `active_view` only to panels that opt in (declare a `view_r` formal),
  # keeping register_panel()'s 3-arg server contract backward compatible.
  for (sec in panels) {
    args <- list(sec$id, data, active_cells)
    if ("view_r" %in% names(formals(sec$server))) args <- c(args, list(active_view))
    do.call(sec$server, args)
  }
}

.scroll_theme <- function() {
  bslib::bs_theme(
    version = 5,
    bg = "#FBFCFD", fg = "#111826", primary = "#2563A8",
    "border-color" = "#E4E8EE"
  )
}

#' Build the scroll explorer app
#'
#' Returns a `shiny::shinyApp` that reads a built project directory. Deploy the
#' scaffolded `app.R` (which calls this) to a Shiny Server, or run it locally
#' with [scroll_serve()].
#'
#' @param dir A built scroll project directory.
#' @return A `shiny.appobj`.
#' @export
scroll_app <- function(dir = ".") {
  data <- .scroll_load(dir)
  panels <- .scroll_assemble_panels(data$manifest)
  title <- .scroll_nz(data$config$title)

  ui <- .scroll_page(data, title, panels)
  server <- function(input, output, session)
    .scroll_wire(input, output, session, data, panels)
  # release the query handle's cached datasets when the app stops
  shiny::shinyApp(ui, server, onStart = function() {
    shiny::onStop(function() try(scroll_disconnect(data$con), silent = TRUE))
  })
}

#' Build a multi-dataset scroll explorer
#'
#' Mounts several built scroll projects behind one app: a tab per dataset, each
#' with its own app bar, manifest-driven panels, and query handle. Every dataset's
#' UI is namespaced (via [shiny::NS()]) so the projects coexist without id
#' collisions; the built-in and custom panels are the same registry for all.
#'
#' Custom panels that read baked assets should resolve them from the per-dataset
#' handle (`data$dir`) rather than a global option, so each tab reads its own
#' project's files.
#'
#' @param projects A named list/vector of built project directories. Names are used
#'   as tab labels (falling back to each project's `config.yaml` title, then a
#'   generic label).
#' @return A `shiny.appobj`.
#' @export
scroll_multi_app <- function(projects) {
  projects <- as.list(projects)
  if (!length(projects)) stop("scroll_multi_app(): supply at least one project directory.")
  labels <- names(projects)
  if (is.null(labels)) labels <- rep("", length(projects))
  ids <- paste0("ds", seq_along(projects))

  datas  <- lapply(projects, .scroll_load)
  panels <- .scroll_assemble_panels(lapply(datas, function(d) d$manifest))
  titles <- vapply(seq_along(datas), function(i) {
    .scroll_nz(datas[[i]]$config$title) %||%
      (if (nzchar(labels[[i]])) labels[[i]] else paste("Dataset", i))
  }, "")
  labels <- ifelse(nzchar(labels), labels, titles)

  tabs <- lapply(seq_along(datas), function(i)
    bslib::nav_panel(labels[[i]],
      .scroll_body(datas[[i]], titles[[i]], panels, shiny::NS(ids[[i]]))))

  ui <- bslib::page_fluid(
    theme = .scroll_theme(),
    tags$head(tags$style(HTML(.scroll_css())), tags$script(HTML(.scroll_spy_js()))),
    # .scroll-multi lets the CSS pin the dataset tab strip and drop each dataset's
    # app bar + rail below it, so the dataset selector stays visible while scrolling.
    div(class = "scroll-multi", do.call(bslib::navset_tab, tabs))
  )
  server <- function(input, output, session) {
    for (i in seq_along(datas)) local({
      ii <- i
      moduleServer(ids[[ii]], function(input, output, session)
        .scroll_wire(input, output, session, datas[[ii]], panels))
    })
  }
  shiny::shinyApp(ui, server, onStart = function() {
    shiny::onStop(function()
      for (d in datas) try(scroll_disconnect(d$con), silent = TRUE))
  })
}

#' Run the scroll explorer locally
#'
#' @param dir A built scroll project directory.
#' @param ... Passed to [shiny::runApp()] (e.g. `port`, `host`, `launch.browser`).
#' @export
scroll_serve <- function(dir = ".", ...) {
  shiny::runApp(scroll_app(dir), ...)
}
