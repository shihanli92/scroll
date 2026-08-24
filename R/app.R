# The runtime: a polished bslib explorer with one scrolling panel per analysis
# type, each with its own controls. The app never loads the Seurat object -- it
# reads the built artifacts (cells.parquet once globally; expression one feature
# at a time via arrow), so runtime RAM stays flat.


# --- panel registry + shell ---------------------------------------------------

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
#'   underscores. Used as the panel anchor and the Shiny module namespace.
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
#
# `include` / `exclude` are the config-driven panel override (config.yaml `panels:` /
# `exclude_panels:`): `include` shows exactly those ids in that order and *bypasses*
# the `when` gate (an explicit request wins); `exclude` drops ids afterwards. Unknown
# ids warn but don't error. When `include` is NULL the usual `when` gating applies.
.scroll_assemble_panels <- function(manifests = NULL, include = NULL, exclude = NULL) {
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
  ids_all <- vapply(panels, `[[`, "", "id")
  if (!is.null(include)) {
    # explicit allowlist: exactly these, in this order, gate bypassed
    include <- as.character(include)
    .scroll_warn_panel_ids(include, ids_all, "panels")
    panels <- panels[stats::na.omit(match(include[include %in% ids_all], ids_all))]
  } else if (!is.null(manifests)) {
    # Modality gating: a panel may declare `when = function(manifest)`; keep it only
    # if some supplied manifest satisfies it (so e.g. VDJ panels appear only for
    # projects built with `vdj=`). `manifests` is one manifest or a list of them
    # (scroll_multi_app); NULL leaves every panel in (used by tests).
    if (!is.null(manifests$scroll_version)) manifests <- list(manifests)
    keep <- vapply(panels, function(p) is.null(p$when) ||
      any(vapply(manifests, function(m) isTRUE(tryCatch(p$when(m), error = function(e) FALSE)),
                 logical(1))), logical(1))
    panels <- panels[keep]
  }
  if (!is.null(exclude)) {
    exclude <- as.character(exclude)
    .scroll_warn_panel_ids(exclude, ids_all, "exclude_panels")
    ids_now <- vapply(panels, `[[`, "", "id")
    panels <- panels[!ids_now %in% exclude]
  }
  for (i in seq_along(panels)) panels[[i]]$num <- sprintf("%02d", i)
  panels
}

# Warn (non-fatally) about panel ids in a config override that match no known panel.
.scroll_warn_panel_ids <- function(requested, known, key) {
  unknown <- setdiff(requested, known)
  if (length(unknown))
    warning(sprintf("scroll: unknown panel id(s) in `%s`: %s", key,
                    paste(unknown, collapse = ", ")), call. = FALSE)
}

# The global cell-subset filter: restrict `cells` to rows whose `col` value is in
# `vals`; no-op when the filter is inactive.
.scroll_subset_cells <- function(cells, col, vals) {
  if (is.null(col) || is.null(vals) || !length(vals) || !col %in% names(cells))
    return(cells)
  cells[as.character(cells[[col]]) %in% vals, , drop = FALSE]
}

# ---- global filter rail (right sidebar) -------------------------------------
# A set of per-column filters that narrow the cells EVERY panel sees, composed into
# the active-cells reactive. Categorical columns get a level multi-select, numeric
# columns a range slider. One deterministic spec list drives the UI, the apply step,
# and Reset, so their input ids line up. config `filters: false` turns the rail off;
# `filters: [col, col]` curates (and orders) which columns appear.

.scroll_filter_specs <- function(data) {
  m <- data$manifest; cfg <- data$config$filters
  if (isFALSE(cfg)) return(list())
  cols <- c(.scroll_cat_cols(m), .scroll_num_cols(m))
  if (is.character(cfg)) cols <- cfg[cfg %in% cols]            # curated subset, config order
  specs <- list()
  for (col in cols) {
    e <- m$meta[[col]]
    if (identical(e$type, "categorical")) {
      if (.scroll_n_levels(m, col) > .SCROLL_MAX_LEVELS) next  # skip barcode/clone-scale cols
      specs[[length(specs) + 1L]] <- list(col = col, type = "categorical",
                                          id = paste0("flt_", length(specs)))
    } else {
      lo <- e$range$min; hi <- e$range$max
      if (is.null(lo) || is.null(hi) || !is.finite(lo) || !is.finite(hi) || lo >= hi) next
      specs[[length(specs) + 1L]] <- list(col = col, type = "numeric",
                                          id = paste0("flt_", length(specs)), min = lo, max = hi)
    }
  }
  specs
}

# A collapsible <details> section (native, no JS). open = TRUE renders it expanded.
.scroll_details <- function(title, ..., open = TRUE)
  tags$details(class = "scroll-ctl-details", open = if (isTRUE(open)) NA,
    tags$summary(class = "scroll-ctl-summary", title),
    div(class = "scroll-ctl-body", ...))

.scroll_filters_ui <- function(data, ns = identity) {
  specs <- .scroll_filter_specs(data)
  if (!length(specs)) return(NULL)
  ctrls <- lapply(specs, function(s) {
    if (identical(s$type, "categorical"))
      div(class = "scroll-filter",
          selectizeInput(ns(s$id), s$col, choices = .scroll_meta_levels(data, s$col),
                         multiple = TRUE,
                         options = list(placeholder = "all", plugins = list("remove_button"))))
    else
      div(class = "scroll-filter",
          sliderInput(ns(s$id), s$col, min = s$min, max = s$max, value = c(s$min, s$max)))
  })
  .scroll_details("Filters", open = TRUE,
    div(class = "scroll-reset-row", actionLink(ns("scroll_filter_reset"), "Reset all")),
    ctrls)
}

# ---- global theme controls (right sidebar) ----------------------------------
# Font size / legend position / gridlines applied to every plot via .scroll_ggtheme().
# Each defaults to "Default" (a no-op), so rendering is unchanged until the user picks.
.scroll_theme_ui <- function(data, ns = identity) {
  if (isFALSE(data$config$theme_controls)) return(NULL)
  sel <- function(id, label, choices)
    div(class = "scroll-filter", selectInput(ns(id), label, choices))
  .scroll_details("Theme", open = TRUE,
    .scroll_details("Text & legend", open = TRUE,
      sel("scroll_theme_font", "Text size",
          c("Default" = "", "Small" = "10", "Medium" = "13", "Large" = "16")),
      sel("scroll_theme_legend", "Legend",
          c("Default" = "", "Right" = "right", "Bottom" = "bottom", "Hidden" = "none"))),
    .scroll_details("Axes", open = FALSE,
      sel("scroll_theme_axes", "Axis text",
          c("Default" = "", "Show" = "show", "Hide" = "hide")),
      sel("scroll_theme_axis_titles", "Axis titles",
          c("Default" = "", "Show" = "show", "Hide" = "hide")),
      sel("scroll_theme_angle", "X label angle",
          c("Default" = "", "0" = "0", "45" = "45", "90" = "90"))),
    .scroll_details("Panel", open = FALSE,
      sel("scroll_theme_grid", "Gridlines",
          c("Default" = "", "On" = "on", "Off" = "off")),
      sel("scroll_theme_border", "Panel border",
          c("Default" = "", "On" = "on", "Off" = "off")),
      sel("scroll_theme_bg", "Background",
          c("Default" = "", "White" = "white", "Grey" = "grey", "None" = "none"))))
}

# The right-hand control rail: a Theme section (always, unless disabled) plus the
# Filters section (when the project has filterable columns).
.scroll_controls_ui <- function(data, ns = identity) {
  th <- .scroll_theme_ui(data, ns); fl <- .scroll_filters_ui(data, ns)
  if (is.null(th) && is.null(fl)) return(NULL)
  tags$aside(class = "scroll-filters", th, fl)
}

# The global theme state as a reactive, read by panels that declare a `theme_r` formal.
.scroll_active_theme <- function(input)
  reactive(list(font   = .scroll_nz(input$scroll_theme_font),
                legend = .scroll_nz(input$scroll_theme_legend),
                grid   = .scroll_nz(input$scroll_theme_grid),
                border = .scroll_nz(input$scroll_theme_border),
                axes   = .scroll_nz(input$scroll_theme_axes),
                axis_titles = .scroll_nz(input$scroll_theme_axis_titles),
                angle  = .scroll_nz(input$scroll_theme_angle),
                bg     = .scroll_nz(input$scroll_theme_bg)))

# Narrow `cells` by every active filter (AND). An untouched control is a no-op: an
# empty categorical selection means "all", a full-range slider means "all". Numeric
# NAs are kept (a scoped column is simply absent for non-members, not filtered out).
.scroll_filter_cells <- function(cells, data, input) {
  specs <- .scroll_filter_specs(data)
  if (!length(specs) || !nrow(cells)) return(cells)
  keep <- rep(TRUE, nrow(cells))
  for (s in specs) {
    if (!s$col %in% names(cells)) next
    v <- cells[[s$col]]; val <- input[[s$id]]
    if (identical(s$type, "categorical")) {
      if (length(val)) keep <- keep & as.character(v) %in% val
    } else if (length(val) == 2L && (val[[1]] > s$min || val[[2]] < s$max)) {
      keep <- keep & (is.na(v) | (v >= val[[1]] & v <= val[[2]]))
    }
  }
  cells[keep, , drop = FALSE]
}

# Reset every filter control to its "all" default on the Reset link.
.scroll_bind_filters <- function(input, session, data) {
  specs <- .scroll_filter_specs(data)
  if (!length(specs)) return(invisible())
  observeEvent(input$scroll_filter_reset, {
    for (s in specs)
      if (identical(s$type, "categorical"))
        updateSelectizeInput(session, s$id, selected = character(0))
      else
        updateSliderInput(session, s$id, value = c(s$min, s$max))
  })
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

.scroll_panel_card <- function(sec, data, ns = identity) {
  bslib::card(
    id = ns(sec$id), class = "scroll-panel-card", full_screen = FALSE,
    bslib::card_header(
      div(class = "scroll-eyebrow",
          span(class = "scroll-num", sec$num), span(class = "scroll-kicker", sec$label)),
      tags$h2(class = "scroll-title", sec$title),
      tags$p(class = "scroll-desc", sec$desc)),
    sec$ui(ns(sec$id), data)
  )
}

# The per-dataset body (app bar + rail + panel cards), with every id passed
# through `ns` so multiple datasets can coexist in one page (see scroll_multi_app).
# Does NOT include page_fluid/theme/head -- those wrap it once at the page level.
.scroll_body <- function(data, title, panels, ns = identity) {
  controls <- .scroll_controls_ui(data, ns)
  tagList(
    .scroll_appbar(data, title, ns),
    div(
      class = paste0("scroll-layout", if (!is.null(controls)) " has-filters"),
      .scroll_rail(panels, ns),
      div(class = "scroll-content",
          lapply(panels, function(s) .scroll_panel_card(s, data, ns))),
      controls                                       # right-hand global control rail
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

# The per-dataset server logic, composed from four focused steps. Called flat by
# scroll_app() and inside a moduleServer() by scroll_multi_app(), so
# `input`/`output`/`session` carry the right namespace in both.
.scroll_wire <- function(input, output, session, data, panels) {
  active_view <- reactive(.scroll_nz(input$scroll_view))
  active_cells <- .scroll_active_cells(input, data, active_view)
  active_theme <- .scroll_active_theme(input)
  .scroll_bind_subset_control(input, session, data)
  .scroll_bind_filters(input, session, data)
  .scroll_render_ncells(output, data$manifest, active_cells, active_view)
  .scroll_mount_panels(data, panels, active_cells, active_view, active_theme)
}

# The active cell set as a reactive: the selected subset view, then narrowed by the
# app-bar categorical filter.
.scroll_active_cells <- function(input, data, active_view) {
  m <- data$manifest
  # debounced so dragging a filter slider coalesces into one redraw of the panels
  shiny::debounce(reactive({
    base <- .scroll_view_cells(data$cells, m, active_view())
    base <- .scroll_subset_cells(base, .scroll_nz(input$scroll_subset_col),
                                 input$scroll_subset_val)
    .scroll_filter_cells(base, data, input)          # right-rail global filters
  }), 150)
}

# Repopulate the app-bar subset-value selectize when its column changes.
# server = TRUE streams choices incrementally, so this stays fast even for a
# high-cardinality column (whose levels are recomputed from cells, not cached).
.scroll_bind_subset_control <- function(input, session, data) {
  observeEvent(input$scroll_subset_col, {
    col <- .scroll_nz(input$scroll_subset_col)
    lv <- if (is.null(col)) character(0) else .scroll_meta_levels(data, col)
    updateSelectizeInput(session, "scroll_subset_val", choices = lv,
                         selected = character(0), server = TRUE)
  })
}

# The app-bar "N of total [ \u00b7 view]" cell-count readout.
.scroll_render_ncells <- function(output, m, active_cells, active_view) {
  output$scroll_ncells <- renderText({
    n <- nrow(active_cells()); tot <- m$n_cells
    lab <- if (!is.null(active_view())) m$subsets[[active_view()]]$label
    base <- if (n < tot) sprintf("%s of %s", format(n, big.mark = ","),
                                 format(tot, big.mark = ",")) else format(tot, big.mark = ",")
    if (!is.null(lab)) paste0(base, " \u00b7 ", lab) else base
  })
}

# Mount each panel's server, threading `active_view` only to panels that declare a
# `view_r` formal (keeping register_panel()'s 3-arg server contract compatible).
.scroll_mount_panels <- function(data, panels, active_cells, active_view, active_theme = reactive(NULL)) {
  for (sec in panels) {
    fmls <- names(formals(sec$server))
    args <- list(sec$id, data, cells_r = active_cells)      # named so formal order can vary
    if ("view_r" %in% fmls)  args$view_r  <- active_view
    if ("theme_r" %in% fmls) args$theme_r <- active_theme
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
  # config.yaml `panels:` / `exclude_panels:` override the automatic gating
  panels <- .scroll_assemble_panels(data$manifest, include = data$config$panels,
                                    exclude = data$config$exclude_panels)
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
#' The `config.yaml` panel override (`panels:` / `exclude_panels:`) is honoured by
#' [scroll_app()]; here the tab panel list is shared, so the panels are the gated
#' union across all mounted projects.
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

# Find a panel spec by id: registered panels first (a custom panel replacing a
# built-in wins), then the built-ins.
.scroll_lookup_panel <- function(id) {
  all <- c(.scroll_registry$panels, .scroll_builtin_panels())
  ids <- vapply(all, `[[`, "", "id")
  if (id %in% ids) all[[which(ids == id)[1L]]] else NULL
}

#' Preview a single panel while developing it
#'
#' Mounts **one** panel in the full scroll harness — the app bar (with the
#' cell-subset filter and the subset-View selector) but none of the other panels —
#' so a custom panel can be driven for quick feedback without building or scrolling
#' through the whole app. Because the harness is the real one, the panel is tested
#' exactly as it will run: under the app-bar filter and any subset views.
#'
#' The development loop is:
#' ```r
#' source("panels/my_panel.R")               # (re-)registers the panel
#' scroll_preview_panel("my_panel", "projects/demo")
#' # edit the panel file -> re-source -> re-run
#' ```
#'
#' @param panel The id of a panel registered with [register_panel()] /
#'   [register_plot_panel()] (or a built-in id, e.g. `"dotplot"`).
#' @param dir A built scroll project directory to render the panel against.
#' @return A `shiny.appobj` — auto-prints/launches at the console, or pass to
#'   [shiny::runApp()].
#' @seealso [register_plot_panel()], [scroll_app()]
#' @export
scroll_preview_panel <- function(panel, dir = ".") {
  if (!is.character(panel) || length(panel) != 1L || is.na(panel))
    stop("`panel` must be a single registered or built-in panel id.", call. = FALSE)
  data <- .scroll_load(dir)
  spec <- .scroll_lookup_panel(panel)
  if (is.null(spec))
    stop("No panel with id '", panel, "'. Register it first with ",
         "register_panel() / register_plot_panel(), or check the id.", call. = FALSE)

  panels <- list(spec)
  ttl <- .scroll_nz(data$config$title) %||% "scroll"
  title <- paste0(ttl, " \u00b7 preview: ", spec$id)
  ui <- .scroll_page(data, title, panels)
  server <- function(input, output, session)
    .scroll_wire(input, output, session, data, panels)
  shiny::shinyApp(ui, server, onStart = function() {
    shiny::onStop(function() try(scroll_disconnect(data$con), silent = TRUE))
  })
}
