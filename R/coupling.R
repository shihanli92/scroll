# The scroll app: per-section stickies driven by closeread's native scroll model.
#
# Each config section owns one closeread sticky (a Shiny output rendering that
# section's view). closeread pins the right sticky as the reader scrolls, so
# there is no scroll->Shiny bridge and no `active_section` input: scroll position
# is closeread's job, reactivity is Shiny's, and they never contend for a DOM
# element because each section has its own.
#
# A persistent header holds the feature search + toggles as global inputs; every
# section's sticky reactive reads them, so the search recolours whatever view is
# in focus. A section renders once and re-renders only when a global input
# changes -- scrolling itself triggers no recompute.

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

# Global input ids (persistent header) and the per-section output id.
.scroll_inputs <- list(
  feature = "scroll_feature", embedding = "scroll_embedding",
  split = "scroll_split", subset_col = "scroll_subset_col",
  subset_val = "scroll_subset_val", labels = "scroll_labels"
)
.scroll_out <- function(section_id) paste0("scroll_sticky_", section_id)

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

# Assemble the context one section's sticky should render, given the current
# global feature search and toggle state. Pure (no Shiny) -> unit-testable.
.scroll_section_ctx <- function(section, feature, state, data) {
  state <- state %||% list()
  view <- section$view
  params <- section$params
  assay <- state$assay %||% params$assay %||% data$config$default_assay

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

# Toggle choices for the header, derived from the manifest.
.scroll_toggle_choices <- function(manifest) {
  cats <- names(Filter(function(m) identical(m$type, "categorical"), manifest$meta))
  list(embeddings = names(manifest$embeddings),
       assays = names(manifest$assays),
       categoricals = cats)
}

#' Persistent header: feature search + toggles
#'
#' Placed once at the top of `story.qmd`. Holds the global inputs every section
#' reads (search box, embedding / split / subset / labels), and inlines
#' closeread's CSS so the sticky layout survives a `server: shiny` deployment
#' (which doesn't serve `_extensions/`).
#'
#' @param dir The project directory (read for toggle choices).
#' @return A Shiny UI tag list.
#' @export
scroll_app_ui <- function(dir = ".") {
  .scroll_require_shiny()
  manifest <- scroll_manifest(dir)
  ch <- .scroll_toggle_choices(manifest)
  ii <- .scroll_inputs
  cat_choices <- stats::setNames(ch$categoricals, ch$categoricals)

  shiny::tagList(
    shiny::tags$style(shiny::HTML(paste(.scroll_closeread_css(dir),
                                        .scroll_header_css(), sep = "\n"))),
    shiny::div(
      class = "scroll-header",
      shiny::textInput(ii$feature, NULL, placeholder = "Search a feature..."),
      shiny::selectInput(ii$embedding, "Embedding", ch$embeddings,
                         selected = manifest$default_embedding),
      shiny::selectInput(ii$split, "Split by", c("(none)" = "", cat_choices)),
      shiny::selectInput(ii$subset_col, "Subset", c("(all)" = "", cat_choices)),
      shiny::textInput(ii$subset_val, "= value", ""),
      shiny::checkboxInput(ii$labels, "Labels", FALSE)
    )
  )
}

#' A single section's sticky output
#'
#' Placed inside each `.cr-section`'s `.sticky` in `story.qmd`. Emits the right
#' output kind (plot, or table for `de_table`) for that section's view.
#'
#' @param id The section id (must exist in `config.yaml`).
#' @param dir The project directory.
#' @return A Shiny output tag.
#' @export
scroll_sticky_ui <- function(id, dir = ".") {
  .scroll_require_shiny()
  sec <- .scroll_sections_by_id(scroll_config(dir))[[id]]
  if (is.null(sec)) stop("No section '", id, "' in config.yaml.", call. = FALSE)
  out <- .scroll_out(id)
  if (identical(scroll_view_kind(sec$view), "table"))
    shiny::div(class = "scroll-table", shiny::tableOutput(out))
  else
    shiny::plotOutput(out, height = "80vh")
}

#' Server: render every section's sticky from the global inputs
#'
#' Registers one renderer per config section. Each reads the shared feature +
#' toggle inputs and its own params; closeread decides which is on screen.
#'
#' @param dir The scroll project directory.
#' @export
scroll_app_server <- function(dir = ".") {
  .scroll_require_shiny()
  # A `server: shiny` deployment serves story_files/ but not _extensions/; serve
  # the extension assets ourselves so url()-referenced assets resolve too.
  ext <- file.path(dir, "_extensions")
  if (dir.exists(ext)) shiny::addResourcePath("_extensions", normalizePath(ext))

  data <- .scroll_data(dir)
  ii <- .scroll_inputs
  domain <- shiny::getDefaultReactiveDomain()
  input <- domain$input
  output <- domain$output

  toggles <- shiny::reactive(list(
    embedding   = .scroll_blank_null(input[[ii$embedding]]),
    split_by    = .scroll_blank_null(input[[ii$split]]),
    subset_col  = .scroll_blank_null(input[[ii$subset_col]]),
    subset_val  = input[[ii$subset_val]],
    show_labels = isTRUE(input[[ii$labels]])
  ))
  feature <- shiny::reactive(input[[ii$feature]])

  for (section in data$config$sections) {
    local({
      s <- section
      draw <- function() {
        ctx <- .scroll_section_ctx(s, feature(), toggles(), data)
        render_view(ctx$view, ctx)
      }
      output[[.scroll_out(s$id)]] <- if (identical(scroll_view_kind(s$view), "table"))
        shiny::renderTable(draw()) else shiny::renderPlot(draw())
    })
  }
  invisible(NULL)
}

.scroll_header_css <- function() {
  paste(
    ".scroll-header{position:sticky;top:0;z-index:1000;display:flex;gap:14px;",
    "align-items:flex-end;flex-wrap:wrap;padding:10px 14px;margin-bottom:8px;",
    "background:rgba(255,255,255,0.96);border-bottom:1px solid #e5e5e5;}",
    ".scroll-header .form-group{margin-bottom:0;}",
    ".scroll-header .shiny-input-container{width:auto;min-width:150px;}",
    ".scroll-table{max-height:80vh;overflow:auto;font-size:0.85em;}",
    sep = "")
}

.scroll_blank_null <- function(x) if (is.null(x) || !nzchar(x)) NULL else x

.scroll_require_shiny <- function() {
  if (!requireNamespace("shiny", quietly = TRUE))
    stop("The scroll app requires the 'shiny' package.", call. = FALSE)
}

# Inline closeread's stylesheet into the page. A `server: shiny` deployment
# doesn't serve `_extensions/`, so closeread's linked CSS 404s and the sticky
# layout collapses; baking it in makes the layout robust everywhere.
.scroll_closeread_css <- function(dir) {
  hits <- Sys.glob(file.path(dir, "_extensions", "*", "closeread", "closeread.css"))
  if (!length(hits)) hits <- Sys.glob(file.path(dir, "_extensions", "closeread", "closeread.css"))
  if (!length(hits)) return("")
  paste(readLines(hits[[1]], warn = FALSE), collapse = "\n")
}
