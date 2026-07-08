# The scroll x reactive coupling, centralized once so every view inherits it.
#
# The sticky visual is modelled as a SINGLE persistent output redrawn by a
# SINGLE reactive whose dependencies are (active_section, feature, toggles).
# closeread owns scroll position; Shiny owns the search box + toggles; they never
# fight over the DOM because there is one sticky, and a JS bridge
# (inst/app/cr-bridge.js) reports the active closeread trigger to Shiny as the
# global input `active_section`.
#
# scroll_app_ui()/scroll_app_server() are the one place the app calls; custom
# views added later plug into render_view() and inherit this coupling for free.

# Per-directory cache of the globally-loaded, session-shared artifacts.
.scroll_cache <- new.env(parent = emptyenv())

.scroll_data <- function(dir) {
  key <- normalizePath(dir)
  if (is.null(.scroll_cache[[key]])) {
    .scroll_cache[[key]] <- list(
      cells    = as.data.frame(arrow::read_parquet(file.path(dir, "cells.parquet"))),
      manifest = scroll_manifest(dir),
      config   = scroll_config(dir),
      con      = scroll_connect(dir)
    )
  }
  .scroll_cache[[key]]
}

# Index config sections by id for O(1) lookup of the active section.
.scroll_sections_by_id <- function(config) {
  stats::setNames(config$sections, vapply(config$sections, `[[`, "", "id"))
}

# Read a precomputed DE contrast table from de/ (parquet/csv/tsv), or NULL.
.scroll_read_de <- function(dir, contrast) {
  if (is.null(contrast)) return(NULL)
  for (ext in c("parquet", "csv", "tsv")) {
    path <- file.path(dir, "de", paste0(contrast, ".", ext))
    if (file.exists(path)) {
      return(switch(ext,
        parquet = as.data.frame(arrow::read_parquet(path)),
        csv     = utils::read.csv(path, check.names = FALSE),
        tsv     = utils::read.delim(path, check.names = FALSE)))
    }
  }
  NULL
}

# The coupling, as a pure function: given the active closeread section, the
# feature-search value, and the toggle state, assemble the context the single
# sticky should render. Kept free of Shiny so it is unit-testable; the Shiny
# wrapper only feeds it inputs.
.scroll_resolve <- function(active_section, feature, state, data) {
  state <- state %||% list()
  sections <- .scroll_sections_by_id(data$config)
  sec <- if (!is.null(active_section) && !is.null(sections[[active_section]]))
    sections[[active_section]] else data$config$sections[[1]]
  view <- sec$view
  params <- sec$params
  assay <- state$assay %||% params$assay %||% data$config$default_assay

  # subset toggle: restrict the cells all views see
  cells <- data$cells
  if (!is.null(state$subset_col) && !is.null(state$subset_val) &&
      nzchar(state$subset_val) && state$subset_col %in% names(cells)) {
    cells <- cells[as.character(cells[[state$subset_col]]) == state$subset_val, ,
                   drop = FALSE]
  }

  ctx <- list(view = view, params = params, state = state, cells = cells,
              feature_values = NULL, feature_name = NULL,
              expr_long = NULL, de_data = NULL, kind = scroll_view_kind(view))
  feature_active <- !is.null(feature) && nzchar(feature)

  if (view %in% c("umap_colorby", "feature_plot", "violin")) {
    feat <- if (feature_active) feature
            else if (!identical(view, "umap_colorby")) params$feature else NULL
    if (!is.null(feat)) {
      hit <- tryCatch(scroll_query_feature(data$con, assay, feat),
                      error = function(e) NULL)
      if (!is.null(hit) && nrow(hit) > 0) {
        hit$value <- scroll_dequantize(hit$value, data$manifest, assay)
        ctx$feature_values <- hit
        ctx$feature_name <- feat
      }
    }
  } else if (identical(view, "dotplot")) {
    feats <- unlist(params$features)
    if (feature_active) feats <- unique(c(feats, feature))
    el <- tryCatch(scroll_query_features(data$con, assay, feats),
                   error = function(e) NULL)
    if (!is.null(el) && nrow(el) > 0)
      el$value <- scroll_dequantize(el$value, data$manifest, assay)
    ctx$expr_long <- el
    ctx$params$features <- feats
  } else if (identical(view, "de_table")) {
    ctx$de_data <- .scroll_read_de(attr(data$con, "scroll_dir"), params$contrast)
  }
  ctx
}

# Toggle choices for the UI, derived from the manifest.
.scroll_toggle_choices <- function(manifest) {
  cats <- names(Filter(function(m) identical(m$type, "categorical"), manifest$meta))
  list(
    embeddings = names(manifest$embeddings),
    assays = names(manifest$assays),
    categoricals = cats
  )
}

#' UI for the shared scroll visual, search box, and toggles
#'
#' Emitted inside the closeread sticky in `story.qmd`. Returns the feature search
#' box, the toggle bar (embedding / split / subset / labels), and the single
#' sticky output (a plot, or a table for `de_table`), plus the JS bridge that
#' reports the active closeread trigger to Shiny.
#'
#' @param id An id prefix (namespaces this app's inputs/outputs).
#' @param dir The project directory (read for toggle choices).
#' @return A Shiny UI tag list.
#' @export
scroll_app_ui <- function(id = "main", dir = ".") {
  .scroll_require_shiny()
  ids <- .scroll_ids(id)
  manifest <- scroll_manifest(dir)
  ch <- .scroll_toggle_choices(manifest)
  subset_choices <- c("(all)" = "")
  if (length(ch$categoricals))
    subset_choices <- c(subset_choices,
                        stats::setNames(ch$categoricals, ch$categoricals))

  shiny::tagList(
    shiny::tags$style(shiny::HTML(.scroll_closeread_css(dir))),
    shiny::tags$script(shiny::HTML(.scroll_bridge_script())),
    shiny::div(
      class = "scroll-searchbox",
      shiny::textInput(ids$feature, label = NULL, placeholder = "Search a feature...")
    ),
    shiny::div(
      class = "scroll-toggles",
      shiny::selectInput(ids$embedding, "Embedding", choices = ch$embeddings,
                         selected = manifest$default_embedding),
      shiny::selectInput(ids$split, "Split by",
                         choices = c("(none)" = "", stats::setNames(ch$categoricals, ch$categoricals))),
      shiny::selectInput(ids$subset_col, "Subset column", choices = subset_choices),
      shiny::textInput(ids$subset_val, "Subset value", value = ""),
      shiny::checkboxInput(ids$labels, "Labels", value = FALSE)
    ),
    shiny::conditionalPanel(
      condition = sprintf("output['%s'] == 'plot'", ids$kind),
      shiny::plotOutput(ids$plot, height = "78vh")),
    shiny::conditionalPanel(
      condition = sprintf("output['%s'] == 'table'", ids$kind),
      shiny::div(class = "scroll-table", shiny::tableOutput(ids$table)))
  )
}

#' Server for the shared scroll visual
#'
#' Wires the single `(active_section, feature, toggles)` reactive to the single
#' sticky output.
#'
#' @param id The id prefix matching [scroll_app_ui()].
#' @param dir The scroll project directory.
#' @export
scroll_app_server <- function(id = "main", dir = ".") {
  .scroll_require_shiny()
  # A `server: shiny` deployment serves story_files/ but not _extensions/, so
  # closeread's linked CSS 404s and the sticky layout collapses. Serve the
  # extension assets ourselves so the story lays out the same on a Shiny Server
  # as under `quarto preview`.
  ext <- file.path(dir, "_extensions")
  if (dir.exists(ext)) shiny::addResourcePath("_extensions", normalizePath(ext))
  ids <- .scroll_ids(id)
  data <- .scroll_data(dir)
  domain <- shiny::getDefaultReactiveDomain()
  output <- domain$output
  input <- domain$input

  toggles <- shiny::reactive({
    list(
      embedding  = .scroll_blank_null(input[[ids$embedding]]),
      split_by   = .scroll_blank_null(input[[ids$split]]),
      subset_col = .scroll_blank_null(input[[ids$subset_col]]),
      subset_val = input[[ids$subset_val]],
      show_labels = isTRUE(input[[ids$labels]])
    )
  })

  view_ctx <- shiny::reactive({
    .scroll_resolve(input[["active_section"]], input[[ids$feature]], toggles(), data)
  })

  output[[ids$kind]] <- shiny::renderText(view_ctx()$kind)
  shiny::outputOptions(output, ids$kind, suspendWhenHidden = FALSE)

  output[[ids$plot]] <- shiny::renderPlot({
    v <- view_ctx()
    if (identical(v$kind, "plot")) render_view(v$view, v) else NULL
  })
  output[[ids$table]] <- shiny::renderTable({
    v <- view_ctx()
    if (identical(v$kind, "table")) render_view(v$view, v) else NULL
  })
  invisible(NULL)
}

.scroll_ids <- function(id) list(
  feature = paste0(id, "_feature"),
  embedding = paste0(id, "_embedding"),
  split = paste0(id, "_split"),
  subset_col = paste0(id, "_subset_col"),
  subset_val = paste0(id, "_subset_val"),
  labels = paste0(id, "_labels"),
  plot = paste0(id, "_plot"),
  table = paste0(id, "_table"),
  kind = paste0(id, "_kind")
)

.scroll_blank_null <- function(x) if (is.null(x) || !nzchar(x)) NULL else x

.scroll_require_shiny <- function() {
  if (!requireNamespace("shiny", quietly = TRUE))
    stop("The scroll app requires the 'shiny' package.", call. = FALSE)
}

# Inline closeread's stylesheet into the rendered HTML. A `server: shiny`
# deployment doesn't serve `_extensions/`, so closeread's linked CSS 404s and the
# sticky layout collapses; baking it into the page makes the layout robust in
# preview and on a Shiny Server alike.
.scroll_closeread_css <- function(dir) {
  hits <- Sys.glob(file.path(dir, "_extensions", "*", "closeread", "closeread.css"))
  if (!length(hits)) hits <- Sys.glob(file.path(dir, "_extensions", "closeread", "closeread.css"))
  if (!length(hits)) return("")
  paste(readLines(hits[[1]], warn = FALSE), collapse = "\n")
}

.scroll_bridge_script <- function() {
  path <- system.file("app", "cr-bridge.js", package = "scroll")
  if (nzchar(path) && file.exists(path)) return(paste(readLines(path), collapse = "\n"))
  local <- file.path("inst", "app", "cr-bridge.js")
  if (file.exists(local)) return(paste(readLines(local), collapse = "\n"))
  ""
}
