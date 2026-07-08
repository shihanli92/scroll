# The scroll x reactive coupling, centralized once so every view inherits it.
#
# The sticky visual is modelled as a SINGLE persistent Shiny output redrawn by a
# SINGLE reactive whose dependencies are (active_section, feature). closeread
# owns scroll position; Shiny owns the search box; they never fight over the DOM
# because there is exactly one output element, and a JS bridge
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

# The coupling, as a pure function: given the active closeread section and the
# current feature-search value, decide what the single sticky should draw. Kept
# free of Shiny so it is unit-testable; the Shiny wrapper only feeds it inputs.
.scroll_resolve <- function(active_section, feature, data) {
  sections <- .scroll_sections_by_id(data$config)
  sec <- if (!is.null(active_section) && !is.null(sections[[active_section]]))
    sections[[active_section]] else data$config$sections[[1]]

  feature_values <- NULL
  feature_name <- NULL
  if (!is.null(feature) && nzchar(feature)) {
    assay <- sec$params$assay %||% data$config$default_assay
    hit <- tryCatch(scroll_query_feature(data$con, assay, feature),
                    error = function(e) NULL)
    if (!is.null(hit) && nrow(hit) > 0) {
      hit$value <- scroll_dequantize(hit$value, data$manifest, assay)
      feature_values <- hit
      feature_name <- feature
    }
  }
  list(view = sec$view, params = sec$params,
       feature_values = feature_values, feature_name = feature_name)
}

#' UI for the shared scroll visual + persistent feature search
#'
#' Emitted inside the closeread sticky in `story.qmd`. Returns the search box and
#' the single sticky plot output, plus the JS bridge that reports the active
#' closeread trigger to Shiny.
#'
#' @param id An id prefix (namespaces this app's inputs/outputs).
#' @return A Shiny UI tag list.
#' @export
scroll_app_ui <- function(id = "main") {
  .scroll_require_shiny()
  feat <- paste0(id, "_feature")
  sticky <- paste0(id, "_sticky")
  shiny::tagList(
    shiny::tags$script(shiny::HTML(.scroll_bridge_script())),
    shiny::div(
      class = "scroll-searchbox",
      shiny::textInput(feat, label = NULL, placeholder = "Search a feature...")
    ),
    shiny::plotOutput(sticky, height = "80vh")
  )
}

#' Server for the shared scroll visual
#'
#' Wires the single `(active_section, feature)` reactive to the single sticky
#' output.
#'
#' @param id The id prefix matching [scroll_app_ui()].
#' @param dir The scroll project directory.
#' @export
scroll_app_server <- function(id = "main", dir = ".") {
  .scroll_require_shiny()
  feat <- paste0(id, "_feature")
  sticky <- paste0(id, "_sticky")

  data <- .scroll_data(dir)
  domain <- shiny::getDefaultReactiveDomain()
  output <- domain$output
  input <- domain$input

  current_view <- shiny::reactive({
    .scroll_resolve(input[["active_section"]], input[[feat]], data)
  })

  output[[sticky]] <- shiny::renderPlot({
    v <- current_view()
    render_view(v$view, data$cells, v$params, v$feature_values, v$feature_name)
  })
  invisible(NULL)
}

.scroll_require_shiny <- function() {
  if (!requireNamespace("shiny", quietly = TRUE))
    stop("The scroll app requires the 'shiny' package.", call. = FALSE)
}

.scroll_bridge_script <- function() {
  path <- system.file("app", "cr-bridge.js", package = "scroll")
  if (nzchar(path) && file.exists(path)) return(paste(readLines(path), collapse = "\n"))
  # Fallback when running from source (system.file empty before install).
  local <- file.path("inst", "app", "cr-bridge.js")
  if (file.exists(local)) return(paste(readLines(local), collapse = "\n"))
  ""
}
