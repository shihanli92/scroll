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
#' @param view Optional subset-view name; when given, subset-scoped columns are
#'   restricted to that view (global columns always appear). `NULL` (default)
#'   returns the whole-dataset column set.
#' @return A character vector of column names.
#' @export
scroll_columns <- function(data, type = c("categorical", "numeric", "any"), view = NULL) {
  type <- match.arg(type)
  m <- data$manifest
  switch(type,
         categorical = .scroll_cat_cols(m, view),
         numeric     = .scroll_num_cols(m, view),
         any         = Filter(function(c) .scroll_col_in_view(m, c, view), names(m$meta)))
}

#' Embeddings (reductions) available in a scroll project
#'
#' @param data A scroll data handle.
#' @param view Optional subset-view name; a subset view returns its (possibly
#'   partial) embeddings, the whole dataset returns full-coverage reductions.
#' @return A character vector of embedding names.
#' @export
scroll_embeddings <- function(data, view = NULL)
  .scroll_view_embeddings(data$manifest, view)

#' Features (e.g. genes) of an assay in a scroll project
#'
#' @param data A scroll data handle.
#' @param assay Assay name; `NULL` (default) uses the manifest's default assay.
#' @return A character vector of feature names.
#' @export
scroll_features <- function(data, assay = NULL)
  .scroll_features_of(data$manifest, assay %||% data$manifest$default_assay)

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
#' @param control_ids Ids of the sibling controls, used only to catch a `from`
#'   that names neither a metadata column nor another control (a common typo);
#'   the builder passes this automatically.
#' @return Invisibly, `NULL`.
#' @export
scroll_bind_levels <- function(input, session, id, from, data, selected = 1,
                               control_ids = character()) {
  levs <- function(col) .scroll_meta_levels(data, col)
  pick <- function(lv) {
    if (identical(selected, "none") || !length(lv)) character(0)
    else if (is.numeric(selected)) lv[selected[selected <= length(lv)]]
    else intersect(as.character(selected), lv)
  }
  is_col <- from %in% names(data$manifest$meta)
  if (!is_col && !from %in% control_ids)
    warning(sprintf(paste0("scroll_input_levels('%s'): from = \"%s\" is neither a ",
                           "metadata column nor another control's id, so the level ",
                           "list will stay empty. Check the spelling of `from`."),
                    id, from), call. = FALSE)
  if (is_col) {                                              # static column
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
    # On error, `fun()` yields the message string (validation gates included);
    # a legitimate plot fn never returns a character, so is.character = error.
    res <- tryCatch(fun(), error = function(e) conditionMessage(e))
    if (is.character(res)) validate(need(FALSE, res))
    validate(need(.scroll_is_plot(res),
                  sprintf("Plot function returned %s, expected a ggplot.", class(res)[1])))
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
#' assay; `scroll_input_assay` / `scroll_input_embedding` pick an assay / reduction
#' from the manifest (the embedding narrows to the active subset view when
#' `view_aware`); `scroll_input_custom` wraps your own `ui`/`bind` as a control; the
#' rest are thin wrappers over the matching Shiny inputs. Wrap any of them in
#' [scroll_show_when()] for conditional visibility.
#'
#' @param id Control id (also the name under which its value reaches `plot`'s
#'   `input`).
#' @param label Control label.
#' @param type For `scroll_input_column`: `"categorical"`, `"numeric"`, or `"any"`;
#'   for `scroll_input_palette`: `"discrete"` or `"continuous"`.
#' @param from For `scroll_input_levels`: a column-picking control's id, or a
#'   categorical column name. Supply either `from` or `choices`.
#' @param multiple Allow multiple selections.
#' @param required If `TRUE`, an empty value shows a "Select <label>." message
#'   instead of computing.
#' @param none If `TRUE`, the level selector defaults to empty (e.g. a "vs. rest").
#' @param selected Default selection: integer index/vector, or `"none"`. For
#'   `scroll_input_palette`, a palette name.
#' @param view_aware For `scroll_input_column`/`scroll_input_embedding`: repopulate
#'   choices from the active subset view (requires the panel to be view-aware).
#' @param prefer For `scroll_input_column`: one or more column names to default-select
#'   when present (first match wins), falling back to `selected` otherwise.
#' @param widget For `scroll_input_choice`: `"radio"`, `"select"`, or `"auto"`
#'   (radio for a short static list, a selectize dropdown otherwise).
#' @param watch For `scroll_input_levels` with a `choices` function: control ids whose
#'   change repopulates the levels (empty = populate once at start-up).
#' @param ui For `scroll_input_custom`: a `function(ns, data)` returning the
#'   control's UI (namespace inputs with `ns()`).
#' @param bind For `scroll_input_custom`: an optional
#'   `function(input, session, data)` wiring server-side behaviour.
#' @param needs For `scroll_input_custom`: a required column type
#'   (`"categorical"`/`"numeric"`) that gates the panel's empty state, or `NULL`.
#' @param choices For most inputs a static vector; `scroll_input_choice` also accepts
#'   a `function(data)` and `scroll_input_levels` a `function(input, data)`, both
#'   evaluated to derive choices from the data handle (e.g. a baked asset under
#'   `data$dir`).
#' @param value,min,max,step,placeholder,inline Passed to the underlying Shiny input.
#' @return A control spec (a list) for `register_plot_panel(controls = )`.
#' @name scroll_input
NULL

#' @rdname scroll_input
#' @export
scroll_input_column <- function(id, label, type = "categorical", selected = 1,
                                view_aware = FALSE, prefer = NULL) {
  pick <- function(cols) {
    if (!is.null(prefer)) { p <- intersect(prefer, cols); if (length(p)) return(p[[1]]) }
    cols[min(selected, length(cols))]
  }
  list(kind = "column", id = id, label = label, required = FALSE, needs = type,
       ui = function(ns, data) {
         cols <- scroll_columns(data, type)
         selectInput(ns(id), label, stats::setNames(cols, cols), selected = pick(cols))
       },
       bind = if (isTRUE(view_aware))
         function(input, session, data, view_r = reactive(NULL))
           observeEvent(view_r(), {
             cols <- scroll_columns(data, type, view_r())
             cur <- input[[id]]
             updateSelectInput(session, id, choices = stats::setNames(cols, cols),
                               selected = if (!is.null(cur) && cur %in% cols) cur else pick(cols))
           }, ignoreNULL = FALSE)
       else function(input, session, data) NULL)
}

#' @rdname scroll_input
#' @export
scroll_input_levels <- function(id, label, from = NULL, choices = NULL, watch = NULL,
                                multiple = TRUE, required = FALSE, none = FALSE, selected = 1) {
  if (isTRUE(required) && isTRUE(none))
    warning(sprintf(paste0("scroll_input_levels('%s'): required = TRUE with none = TRUE ",
                           "defaults to an empty selection, so the panel stays blocked ",
                           "until the user picks a level."), id), call. = FALSE)
  if (is.null(from) && is.null(choices))
    stop("scroll_input_levels(): supply either `from` (a column/control) or `choices` (a function).",
         call. = FALSE)
  sel <- if (none) "none" else selected
  pick <- function(lv) {
    if (identical(sel, "none") || !length(lv)) character(0)
    else if (is.numeric(sel)) lv[sel[sel <= length(lv)]]
    else intersect(as.character(sel), lv)
  }
  list(kind = "levels", id = id, label = label, required = required, needs = NULL,
       from = from,
       ui = function(ns, data)
         selectizeInput(ns(id), label, choices = NULL, multiple = multiple,
                        options = list(plugins = list("remove_button"),
                                       placeholder = if (none) "rest / none" else "all")),
       # `choices` (a function(input, data)) takes precedence over `from`; `watch` names
       # the control ids whose change repopulates the levels (empty = populate once).
       bind = if (!is.null(choices))
         function(input, session, data, control_ids = character()) {
           repop <- function() {
             lv <- tryCatch(as.character(choices(input, data)), error = function(e) character(0))
             lv <- sort(unique(lv[!is.na(lv) & nzchar(lv)]))
             updateSelectizeInput(session, id, choices = lv, selected = pick(lv))
           }
           if (length(watch))
             observeEvent(lapply(watch, function(w) input[[w]]), repop(), ignoreNULL = FALSE)
           else repop()
         }
       else
         function(input, session, data, control_ids = character())
           scroll_bind_levels(input, session, id, from, data, selected = sel,
                              control_ids = control_ids))
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
scroll_input_choice <- function(id, label, choices, selected = NULL, inline = TRUE,
                                multiple = FALSE, widget = c("auto", "radio", "select")) {
  widget <- match.arg(widget)
  list(kind = "choice", id = id, label = label, required = FALSE, needs = NULL,
       ui = function(ns, data) {
         ch <- if (is.function(choices)) choices(data) else choices    # data-derived choices
         w  <- if (widget != "auto") widget
               else if (is.function(choices) || isTRUE(multiple) || length(ch) > 8) "select"
               else "radio"
         # `selected` may be a function(data) so the default can be data-derived
         # (e.g. the modal value) while `choices` stays in its own order.
         sel <- (if (is.function(selected)) selected(data) else selected) %||%
                (if (length(ch)) ch[[1]] else NULL)
         if (identical(w, "radio"))
           radioButtons(ns(id), label, ch, selected = sel, inline = inline)
         else
           selectizeInput(ns(id), label, ch, selected = sel, multiple = multiple,
                          options = if (isTRUE(multiple)) list(plugins = list("remove_button")) else list())
       },
       bind = function(input, session, data) NULL)
}

#' @rdname scroll_input
#' @export
scroll_input_text <- function(id, label, placeholder = NULL)
  list(kind = "text", id = id, label = label, required = FALSE, needs = NULL,
       ui = function(ns, data) textInput(ns(id), label, placeholder = placeholder),
       bind = function(input, session, data) NULL)

#' @rdname scroll_input
#' @export
scroll_input_palette <- function(id, label = "Palette",
                                 type = c("discrete", "continuous"), selected = NULL) {
  type <- match.arg(type)
  list(kind = "palette", id = id, label = label, required = FALSE, needs = NULL,
       ui = function(ns, data) {
         pals <- if (type == "discrete") names(.scroll_discrete_palettes)
                 else .scroll_continuous_palettes
         selectInput(ns(id), label, pals, selected = selected %||% pals[[1]])
       },
       bind = function(input, session, data) NULL)
}

#' @rdname scroll_input
#' @export
scroll_input_assay <- function(id, label = "Assay")
  list(kind = "assay", id = id, label = label, required = FALSE, needs = NULL,
       ui = function(ns, data) {
         assays <- .scroll_assays_of(data$manifest)
         sel <- data$manifest$default_assay %||% assays[1]
         inp <- selectInput(ns(id), label, stats::setNames(assays, assays), selected = sel)
         # One assay: keep the value available to `plot` but hide the useless control.
         if (length(assays) <= 1) div(style = "display:none", inp) else inp
       },
       bind = function(input, session, data) NULL)

#' @rdname scroll_input
#' @export
scroll_input_embedding <- function(id, label = "Embedding", view_aware = TRUE)
  list(kind = "embedding", id = id, label = label, required = FALSE, needs = NULL,
       ui = function(ns, data) {
         emb <- scroll_embeddings(data)
         selectInput(ns(id), label, stats::setNames(emb, emb),
                     selected = .scroll_default(data, "default_embedding", emb[1]))
       },
       bind = if (isTRUE(view_aware))
         function(input, session, data, view_r = reactive(NULL))
           observeEvent(view_r(), {
             emb <- scroll_embeddings(data, view_r())
             cur <- input[[id]]
             updateSelectInput(session, id, choices = stats::setNames(emb, emb),
                               selected = if (!is.null(cur) && cur %in% emb) cur else emb[1])
           }, ignoreNULL = FALSE)
       else function(input, session, data) NULL)

#' @rdname scroll_input
#' @export
scroll_input_custom <- function(id, label, ui, bind = NULL, needs = NULL,
                                required = FALSE) {
  if (!is.function(ui)) stop("`ui` must be a function(ns, data).", call. = FALSE)
  list(kind = "custom", id = id, label = label, required = required, needs = needs,
       ui = ui, bind = bind %||% function(input, session, data) NULL)
}

#' Show a control only when another control has a given value
#'
#' Decorates a `scroll_input_*()` spec so the builder wraps it in a
#' [shiny::conditionalPanel()] — the control appears only when the control named
#' `control` holds one of `equals`. Enables declarative cascading (e.g. show a
#' replicate picker only when an engine choice is `"Pseudobulk"`).
#'
#' @param control_spec A control spec from a `scroll_input_*()` constructor.
#' @param control Id of the control whose value gates visibility.
#' @param equals One or more values of `control` that reveal this control.
#' @return The decorated control spec.
#' @export
scroll_show_when <- function(control_spec, control, equals) {
  if (!is.list(control_spec) || is.null(control_spec$ui))
    stop("`control_spec` must be a scroll_input_*() spec.", call. = FALSE)
  control_spec$visible_when <- list(control = control, equals = as.character(equals))
  control_spec
}

# ---- Layer 1: the builder ----------------------------------------------------

# Call a control's bind, passing the optional `view_r` / `control_ids` / `output`
# only when its formals declare them (mirrors the panel-server view_r idiom in
# .scroll_wire). Keeps 3-arg author binds working while letting new controls opt into
# more — e.g. a control that renders a dynamic uiOutput needs `output`.
.scroll_call_bind <- function(bind, input, session, data, view_r, control_ids, output = NULL) {
  fmls <- names(formals(bind))
  extra <- list()
  if ("view_r" %in% fmls)      extra$view_r <- view_r
  if ("control_ids" %in% fmls) extra$control_ids <- control_ids
  if ("output" %in% fmls)      extra$output <- output
  do.call(bind, c(list(input, session, data), extra))
}

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
#' @param csv If `TRUE`, add a CSV button that exports the plot's source data. The
#'   `plot` function opts in by attaching the table to its result, e.g.
#'   `attr(p, "scroll_source") <- df; p`. Default `FALSE`.
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
                                desc = NULL, after = NULL, before = NULL, compute = TRUE,
                                csv = FALSE) {
  if (!is.function(plot)) stop("`plot` must be a function(cells, input, data).", call. = FALSE)
  if (!is.list(controls) || (length(controls) &&
        !all(vapply(controls, function(c) is.list(c) && !is.null(c$ui), logical(1)))))
    stop("`controls` must be a list of scroll_input_*() specs.", call. = FALSE)

  us <- .scroll_plot_panel_uiserver(plot, controls, compute, csv = csv)
  register_panel(id, us$ui, us$server, label = label, title = title, desc = desc,
                 after = after, before = before)
}

# Build the (ui, server) pair for a declarative plot panel from a plot function +
# control specs. Shared by register_plot_panel() and the built-in modality panels
# (VDJ/spatial/ATAC), so those get the same control population / error surfacing /
# Compute gate / downloads without re-implementing the module.
#
# `csv = TRUE` adds a CSV button to the plot toolbar and exports the plot's source
# data: the plot fn opts in by attaching it as attr(p, "scroll_source") <- df, which
# rides on the object the plot reactive already returns (no re-computation). A plot
# that attaches nothing yields an empty CSV.
.scroll_plot_panel_uiserver <- function(plot, controls, compute = TRUE, csv = FALSE) {
  ui <- function(id, data) {
    ns <- NS(id)
    miss <- .scroll_panel_missing(controls, data)
    if (!is.null(miss)) return(.scroll_empty_panel(miss))
    ctl_ui <- lapply(controls, function(ctl) {
      u <- ctl$ui(ns, data)
      vw <- ctl$visible_when
      if (is.null(vw)) return(u)
      cond <- if (length(vw$equals) == 1)
        sprintf("input['%s'] == '%s'", vw$control, vw$equals)
      else
        sprintf("[%s].indexOf(input['%s']) > -1",
                paste(sprintf("'%s'", vw$equals), collapse = ","), vw$control)
      shiny::conditionalPanel(cond, u, ns = ns)
    })
    go <- if (isTRUE(compute))
      actionButton(ns("scroll_compute"), "Compute", class = "btn-primary", width = "100%")
    bslib::layout_columns(
      col_widths = c(3, 9), class = "scroll-panel",
      div(class = "scroll-controls",
          do.call(.scroll_group, c(list("Controls"), ctl_ui)), go),
      .scroll_plot_area(ns, csv = isTRUE(csv)))
  }

  server <- function(id, data, cells_r = reactive(data$cells), view_r = reactive(NULL),
                     theme_r = reactive(NULL)) {
    moduleServer(id, function(input, output, session) {
      control_ids <- vapply(controls, function(c) c$id, character(1))
      for (ctl in controls)
        .scroll_call_bind(ctl$bind, input, session, data, view_r, control_ids, output)
      req_ctls <- Filter(function(c) isTRUE(c$required), controls)
      body <- function() {
        for (c in req_ctls)
          validate(need(length(input[[c$id]]) > 0, paste0("Select ", c$label, ".")))
        p <- plot(cells_r(), input, data)
        # apply the global theme controls to bare ggplots (not aplot/patchwork composites),
        # then an optional per-panel aspect ratio (any panel with an `aspect` slider; 1 =
        # unconstrained). `+` preserves the plot fn's scroll_source attr (used by CSV).
        if (inherits(p, "ggplot")) {
          p <- p + .scroll_ggtheme(theme_r())
          a <- input$aspect
          if (is.numeric(a) && length(a) == 1 && abs(a - 1) > 1e-6)
            p <- p + ggplot2::theme(aspect.ratio = a)
        }
        p
      }
      event <- if (isTRUE(compute)) reactive(input$scroll_compute) else NULL
      placeholder <- if (isTRUE(compute)) "Set the controls, then click Compute." else NULL
      plot_r <- scroll_render_plot(output, id, body, event = event, placeholder = placeholder)
      # Export the plot's source table (attached by the plot fn) when csv is on.
      if (isTRUE(csv))
        output$csv <- .scroll_csv_handler(
          reactive(attr(plot_r(), "scroll_source") %||% data.frame()),
          paste0("scroll_", id, ".csv"))
    })
  }

  list(ui = ui, server = server)
}

# A built-in modality panel spec built from the declarative builder and gated on a
# manifest predicate (see .scroll_assemble_panels). Used by the VDJ and ATAC panels;
# the spatial panel is hand-written instead because it owns brush-zoom state.
.scroll_gated_panel <- function(id, label, title, desc, when, controls, plot, csv = FALSE)
  c(list(id = id, label = label, title = title, desc = desc, when = when),
    .scroll_plot_panel_uiserver(plot, controls, compute = FALSE, csv = csv))

# ---- Reusable render helpers (parity with the built-in panels) ---------------

#' Reusable plotting helpers for custom views
#'
#' The pieces the built-in panels use, exposed so a custom `plot` (or a raw
#' [register_panel()] view) matches their look and performance:
#'
#' * `scroll_point_layer()` — a scatter layer that rasterizes on screen for speed
#'   (via scattermore, when installed and `raster = TRUE`) but stays a true
#'   `geom_point` otherwise, so exports and small plots are unchanged.
#' * `scroll_discrete_colors()` — a named colour vector for a categorical column,
#'   keyed by level so a group keeps its colour across panels.
#' * `scroll_continuous_scale()` — a ggplot2 continuous colour scale (viridis,
#'   ColorBrewer, or grey-to-hue) matching FeaturePlot.
#' * `scroll_group_labels()` — per-group median-centroid text labels for a scatter.
#'
#' @param mapping,size,alpha,raster,... For `scroll_point_layer()`: a ggplot2
#'   aesthetic mapping, point size, opacity, whether to rasterize on screen, and
#'   extra args passed to the geom.
#' @param values A vector of categorical values (its unique levels are coloured).
#' @param palette A palette name (discrete: e.g. `"Tableau 10"`; continuous: e.g.
#'   `"viridis"`, `"magma"`, a ColorBrewer name, or `"grey-red"`).
#' @param name,limits For `scroll_continuous_scale()`: legend title and optional
#'   `c(low, high)` clip (out-of-range values are squished, not dropped).
#' @param df A data.frame with the scatter coordinates and group column.
#' @param group,x,y For `scroll_group_labels()`: the grouping and coordinate
#'   column names in `df`.
#' @return A ggplot2 layer/scale (or, for `scroll_discrete_colors`, a named vector).
#' @name scroll_helpers
NULL

#' @rdname scroll_helpers
#' @export
scroll_point_layer <- function(mapping = NULL, size = 0.6, alpha = 1, raster = FALSE, ...)
  .scroll_point_layer(mapping = mapping, size = size, alpha = alpha, raster = raster, ...)

#' @rdname scroll_helpers
#' @export
scroll_discrete_colors <- function(values, palette = "Tableau 10")
  .scroll_discrete_colors(values, palette)

#' @rdname scroll_helpers
#' @export
scroll_continuous_scale <- function(palette = "viridis", name = NULL, limits = NULL)
  .scroll_continuous_scale(palette, name = name, limits = limits)

#' @rdname scroll_helpers
#' @export
scroll_group_labels <- function(df, group, x, y) {
  d <- data.frame(g = as.character(df[[group]]), x = df[[x]], y = df[[y]])
  agg <- stats::aggregate(cbind(x, y) ~ g, data = d, FUN = stats::median)
  ggplot2::geom_text(data = agg,
                     ggplot2::aes(x = .data$x, y = .data$y, label = .data$g),
                     inherit.aes = FALSE, size = 4, fontface = "bold")
}
