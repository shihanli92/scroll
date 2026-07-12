# Approachable custom-panel authoring. Two layers on top of register_panel():
#   - Layer 2 helpers (scroll_render_plot / scroll_bind_levels / scroll_columns)
#     remove the per-panel Shiny boilerplate and the silent-failure footguns.
#   - Layer 1 (register_plot_panel + scroll_input_* controls) is fully declarative:
#     the author writes only a plot function plus a list of declared controls.
# Both reuse the internal panel helpers in R/app.R (.scroll_plot_area,
# .scroll_plot_downloads, .scroll_group, .scroll_empty_panel, .scroll_cat_cols, ...).

# ---- Layer 2: power-user helpers --------------------------------------------

#' Metadata columns of a given type
#'
#' Convenience accessor over a scroll data handle's manifest, for populating
#' custom-panel controls.
#'
#' @param data A scroll data handle (as passed to a panel `ui`/`server`).
#' @param type `"categorical"`, `"numeric"`, or `"any"`.
#' @return A character vector of column names.
#' @export
scroll_columns <- function(data, type = c("categorical", "numeric", "any")) {
  type <- match.arg(type)
  m <- data$manifest
  switch(type,
         categorical = .scroll_cat_cols(m),
         numeric     = .scroll_num_cols(m),
         any         = names(m$meta))
}

#' Populate a selectize control from a metadata column's levels
#'
#' Wires the observe-and-update pattern custom panels need: whenever the column
#' named by `from` changes, the `id` selectize is repopulated with that column's
#' levels. `from` may be another control's input id (dynamic) or a literal
#' categorical column name (static, populated once).
#'
#' @param input,session The module's `input` and `session`.
#' @param id Input id of the selectize to populate (unnamespaced).
#' @param from A control id whose value is a column name, or a column name.
#' @param data The scroll data handle.
#' @param selected Which levels to preselect: an integer index/vector (e.g. `1`,
#'   `c(1, 2)`), or `"none"` for an empty default.
#' @return Invisibly, `NULL`.
#' @export
scroll_bind_levels <- function(input, session, id, from, data, selected = 1) {
  levs <- function(col) unlist(data$manifest$meta[[col]]$levels)
  pick <- function(lv) {
    if (identical(selected, "none") || !length(lv)) character(0)
    else if (is.numeric(selected)) lv[selected[selected <= length(lv)]]
    else intersect(as.character(selected), lv)
  }
  if (from %in% names(data$manifest$meta)) {                 # static column
    lv <- levs(from)
    updateSelectizeInput(session, id, choices = lv, selected = pick(lv))
  } else {                                                    # follows another control
    observeEvent(input[[from]], {
      req(input[[from]])
      lv <- levs(input[[from]])
      updateSelectizeInput(session, id, choices = lv, selected = pick(lv))
    })
  }
  invisible(NULL)
}

# A plot-like object our download path (print()-based) can render.
.scroll_is_plot <- function(x) inherits(x, c("ggplot", "gg", "patchwork", "aplot"))

#' Wire a plot output with error surfacing and PNG/PDF downloads
#'
#' One call replaces the usual `reactive`/`eventReactive` + `tryCatch` +
#' `renderPlot` + download-handler boilerplate. `fun` is evaluated inside a
#' `tryCatch`; any error (or a non-plot return) is shown **inline** in the plot
#' area rather than failing silently, and PNG/PDF handlers (`output$png` /
#' `output$pdf`, matching [scroll_render_plot()]'s companion UI
#' `.scroll_plot_area()`) are registered for the result.
#'
#' @param output The module's `output`.
#' @param id The panel id (used for download filenames).
#' @param fun A zero-argument function returning a ggplot (may `stop()` on bad
#'   input — the message is surfaced).
#' @param event `NULL` for a live plot, or a reactive gating recomputation (e.g.
#'   `shiny::reactive(input$go)` for a Compute button). When it is an action-button
#'   count, the plot waits for the first click.
#' @param placeholder Optional message shown before the first `event` fires.
#' @return Invisibly, the plot reactive.
#' @export
scroll_render_plot <- function(output, id, fun, event = NULL, placeholder = NULL) {
  safe <- function() {
    res <- tryCatch(fun(), error = function(e) conditionMessage(e))
    validate(need(.scroll_is_plot(res), res))
    res
  }
  plot_r <- if (is.null(event)) reactive(safe()) else eventReactive(event(), safe())
  output$plot <- renderPlot({
    if (!is.null(event) && !is.null(placeholder)) {
      ev <- tryCatch(event(), error = function(e) NULL)
      computed <- !is.null(ev) && (!is.numeric(ev) || ev > 0)
      validate(need(computed, placeholder))
    }
    plot_r()
  })
  .scroll_plot_downloads(output, plot_r, id)
  invisible(plot_r)
}

# ---- Layer 1: declarative control constructors ------------------------------
# Each returns a spec: list(kind, id, label, required, needs, ui(ns,data),
# bind(input,session,data)). `needs` names a required column type for the UI guard.

#' Declarative controls for [register_plot_panel()]
#'
#' Each constructor returns a control spec. `scroll_input_column` picks a metadata
#' column of a type; `scroll_input_levels` is a level selector bound (via
#' [scroll_bind_levels()]) to the column chosen by another control (`from`) or a
#' fixed column; `scroll_input_gene` is a server-side gene search over the default
#' assay; the rest are thin wrappers over the matching Shiny inputs.
#'
#' @param id Control id (also the name under which its value reaches `plot`'s
#'   `input`).
#' @param label Control label.
#' @param type For `scroll_input_column`: `"categorical"`, `"numeric"`, or `"any"`.
#' @param from For `scroll_input_levels`: a column-picking control's id, or a
#'   categorical column name.
#' @param multiple Allow multiple selections.
#' @param required If `TRUE`, an empty value shows a "Select <label>." message
#'   instead of computing.
#' @param none If `TRUE`, the level selector defaults to empty (e.g. a "vs. rest").
#' @param selected Default selection: integer index/vector, or `"none"`.
#' @param choices,value,min,max,step,placeholder,inline Passed to the underlying
#'   Shiny input.
#' @return A control spec (a list) for `register_plot_panel(controls = )`.
#' @name scroll_input
NULL

#' @rdname scroll_input
#' @export
scroll_input_column <- function(id, label, type = "categorical", selected = 1) {
  list(kind = "column", id = id, label = label, required = FALSE, needs = type,
       ui = function(ns, data) {
         cols <- scroll_columns(data, type)
         selectInput(ns(id), label, stats::setNames(cols, cols),
                     selected = cols[min(selected, length(cols))])
       },
       bind = function(input, session, data) NULL)
}

#' @rdname scroll_input
#' @export
scroll_input_levels <- function(id, label, from, multiple = TRUE, required = FALSE,
                                none = FALSE, selected = 1) {
  sel <- if (none) "none" else selected
  list(kind = "levels", id = id, label = label, required = required, needs = NULL,
       ui = function(ns, data)
         selectizeInput(ns(id), label, choices = NULL, multiple = multiple,
                        options = list(plugins = list("remove_button"),
                                       placeholder = if (none) "rest / none" else NULL)),
       bind = function(input, session, data)
         scroll_bind_levels(input, session, id, from, data, selected = sel))
}

#' @rdname scroll_input
#' @export
scroll_input_gene <- function(id, label, multiple = FALSE, required = FALSE) {
  list(kind = "gene", id = id, label = label, required = required, needs = NULL,
       ui = function(ns, data)
         selectizeInput(ns(id), label, choices = NULL, multiple = multiple,
                        options = list(placeholder = "type a gene", maxOptions = 50)),
       bind = function(input, session, data)
         updateSelectizeInput(session, id, server = TRUE,
                              choices = .scroll_features_of(data$manifest, data$manifest$default_assay)))
}

#' @rdname scroll_input
#' @export
scroll_input_numeric <- function(id, label, value, min = NA, max = NA, step = NA)
  list(kind = "numeric", id = id, label = label, required = FALSE, needs = NULL,
       ui = function(ns, data) numericInput(ns(id), label, value, min = min, max = max, step = step),
       bind = function(input, session, data) NULL)

#' @rdname scroll_input
#' @export
scroll_input_slider <- function(id, label, min, max, value, step = NULL)
  list(kind = "slider", id = id, label = label, required = FALSE, needs = NULL,
       ui = function(ns, data) sliderInput(ns(id), label, min, max, value, step = step),
       bind = function(input, session, data) NULL)

#' @rdname scroll_input
#' @export
scroll_input_choice <- function(id, label, choices, selected = NULL, inline = TRUE)
  list(kind = "choice", id = id, label = label, required = FALSE, needs = NULL,
       ui = function(ns, data) radioButtons(ns(id), label, choices,
                                            selected = selected %||% choices[[1]], inline = inline),
       bind = function(input, session, data) NULL)

#' @rdname scroll_input
#' @export
scroll_input_text <- function(id, label, placeholder = NULL)
  list(kind = "text", id = id, label = label, required = FALSE, needs = NULL,
       ui = function(ns, data) textInput(ns(id), label, placeholder = placeholder),
       bind = function(input, session, data) NULL)

# ---- Layer 1: the builder ----------------------------------------------------

# If a control needs a column type the dataset lacks, return a reason string.
.scroll_panel_missing <- function(controls, data) {
  for (ctl in controls) {
    if (!is.null(ctl$needs) && ctl$needs != "any" &&
        length(scroll_columns(data, ctl$needs)) == 0)
      return(paste0("This panel needs a ", ctl$needs, " metadata column (none found)."))
  }
  NULL
}

#' Register a plot panel declaratively (no Shiny boilerplate)
#'
#' A higher-level companion to [register_panel()]: describe the controls and a
#' single plot function, and scroll generates the whole module — UI, control
#' population, error surfacing, a Compute gate, downloads, spinner, and app-bar
#' subset/view awareness. The silent-failure footguns of a hand-written server
#' (a stray `req()`, an uncaught error) are handled for you.
#'
#' @param id Panel id (a single name).
#' @param plot `function(cells, input, data)` returning a ggplot. `cells` is the
#'   active (subset/view-filtered) cell table; `input` exposes each control by its
#'   id; `data` is the handle (`data$query1`, `data$con`, `data$manifest`, …).
#' @param controls A list of `scroll_input_*()` specs (see [scroll_input]).
#' @param label,title,desc,after,before As in [register_panel()].
#' @param compute If `TRUE` (default), the plot recomputes on a Compute button;
#'   `FALSE` makes it live (recompute on any control change).
#' @return Invisibly, `id`.
#' @examples
#' \dontrun{
#' register_plot_panel("counts",
#'   controls = list(scroll_input_column("grp", "Group", "categorical")),
#'   plot = function(cells, input, data) {
#'     ggplot2::ggplot(cells, ggplot2::aes(.data[[input$grp]])) + ggplot2::geom_bar()
#'   },
#'   label = "Counts", compute = FALSE)
#' }
#' @export
register_plot_panel <- function(id, plot, controls = list(), label = id, title = label,
                                desc = NULL, after = NULL, before = NULL, compute = TRUE) {
  if (!is.function(plot)) stop("`plot` must be a function(cells, input, data).", call. = FALSE)
  if (!is.list(controls) || (length(controls) &&
        !all(vapply(controls, function(c) is.list(c) && !is.null(c$ui), logical(1)))))
    stop("`controls` must be a list of scroll_input_*() specs.", call. = FALSE)

  ui <- function(id, data) {
    ns <- NS(id)
    miss <- .scroll_panel_missing(controls, data)
    if (!is.null(miss)) return(.scroll_empty_panel(miss))
    ctl_ui <- lapply(controls, function(ctl) ctl$ui(ns, data))
    go <- if (isTRUE(compute))
      actionButton(ns("scroll_compute"), "Compute", class = "btn-primary", width = "100%")
    bslib::layout_columns(
      col_widths = c(3, 9), class = "scroll-panel",
      div(class = "scroll-controls",
          do.call(.scroll_group, c(list("Controls"), ctl_ui)), go),
      .scroll_plot_area(ns))
  }

  server <- function(id, data, cells_r = reactive(data$cells), view_r = reactive(NULL)) {
    moduleServer(id, function(input, output, session) {
      for (ctl in controls) ctl$bind(input, session, data)
      req_ctls <- Filter(function(c) isTRUE(c$required), controls)
      body <- function() {
        for (c in req_ctls)
          validate(need(length(input[[c$id]]) > 0, paste0("Select ", c$label, ".")))
        plot(cells_r(), input, data)
      }
      event <- if (isTRUE(compute)) reactive(input$scroll_compute) else NULL
      placeholder <- if (isTRUE(compute)) "Set the controls, then click Compute." else NULL
      scroll_render_plot(output, id, body, event = event, placeholder = placeholder)
    })
  }

  register_panel(id, ui, server, label = label, title = title, desc = desc,
                 after = after, before = before)
}
