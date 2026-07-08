# Scaffolding: emit an authorable config.yaml (ordered sections + view specs)
# and a starter story.qmd closeread narrative. A new dataset is *authored*
# (edit these two files), not *coded*.

#' Scaffold config.yaml and story.qmd for a project
#'
#' Called by [scroll_build()]; exported so authors can regenerate the
#' scaffolding. Existing files are not overwritten.
#'
#' @param outdir Project directory.
#' @param assay_info Per-assay info list from the build (features, max).
#' @param embeddings Character vector of exported reduction names.
#' @param meta_cols Character vector of exported metadata columns.
#' @param md The object metadata data.frame (used to pick sensible defaults).
#' @return `outdir`, invisibly.
#' @export
scroll_scaffold <- function(outdir, assay_info, embeddings, meta_cols, md) {
  default_embedding <- .scroll_default_embedding(embeddings)
  cats <- meta_cols[vapply(meta_cols, function(c)
    is.factor(md[[c]]) || is.character(md[[c]]) || is.logical(md[[c]]),
    logical(1))]
  default_assay <- names(assay_info)[[1]]
  features <- unlist(assay_info[[default_assay]]$features)

  .scroll_write_config(outdir, default_assay, default_embedding, cats, features)
  .scroll_write_story(outdir)
  invisible(outdir)
}

# Prefer recognisable marker genes for the demo panels, else fall back.
.scroll_pick_features <- function(features, n = 6) {
  markers <- c("CD3D", "CD8A", "IL7R", "CCR7", "MS4A1", "CD79A", "CD14", "LYZ",
               "FCGR3A", "NKG7", "GNLY", "PPBP")
  hit <- intersect(markers, features)
  if (length(hit) >= 3) utils::head(hit, n) else utils::head(features, n)
}

.scroll_write_config <- function(outdir, assay, embedding, cats, features) {
  color_by <- if (length(cats)) cats[[1]] else NULL
  feature <- features[[1]]
  panel <- .scroll_pick_features(features)

  sections <- list(
    list(id = "overview", view = "umap_colorby",
         params = list(embedding = embedding, color_by = color_by %||% "orig.ident")),
    list(id = "feature", view = "feature_plot",
         params = list(embedding = embedding, assay = assay, feature = feature))
  )
  if (length(cats) >= 1) {
    sections <- c(sections, list(
      list(id = "markers", view = "dotplot",
           params = list(assay = assay, group_by = cats[[1]], features = as.list(panel))),
      list(id = "distribution", view = "violin",
           params = list(assay = assay, group_by = cats[[1]], feature = feature))))
  }
  if (length(cats) >= 2) {
    sections <- c(sections, list(
      list(id = "composition", view = "proportions",
           params = list(group_by = cats[[2]], fill_by = cats[[1]]))))
  }

  config <- list(
    title = "A scroll story",
    default_assay = assay,
    default_embedding = embedding,
    sections = sections
  )
  path <- file.path(outdir, "config.yaml")
  if (!file.exists(path)) yaml::write_yaml(config, path)
}

.scroll_write_story <- function(outdir) {
  path <- file.path(outdir, "story.qmd")
  if (file.exists(path)) return(invisible())
  writeLines(.scroll_story_template(), path)
}

# A closeread + Shiny narrative. The sticky visual and the persistent search box
# come from scroll_app_ui(); server logic from scroll_app_server(). Each
# closeread trigger's id matches a section id in config.yaml; cr-bridge.js pushes
# the active trigger to Shiny as `active_section`.
.scroll_story_template <- function() {
  c(
    "---",
    "title: \"A scroll story\"",
    "format:",
    "  closeread-html:",
    "    cr-style:",
    "      narrative-background-color-overlay: \"#111111\"",
    "server: shiny",
    "---",
    "",
    "```{r}",
    "#| context: setup",
    "library(scroll)",
    "project_dir <- \".\"",
    "```",
    "",
    "::: {.cr-section}",
    "",
    "::: {#cr-sticky .sticky}",
    "```{r}",
    "scroll_app_ui(\"main\", dir = project_dir)",
    "```",
    ":::",
    "",
    "Welcome. Scroll to walk through the analysis; the visual on the",
    "right stays pinned and answers to both your scroll position and the",
    "feature search box above it. [@cr-overview]{#overview}",
    "",
    "Search any gene above - the pinned view recolors live to that",
    "feature's expression. [@cr-overview]{#feature}",
    "",
    "A marker panel: dot size is the fraction of cells expressing, colour is",
    "mean expression, across groups. [@cr-overview]{#markers}",
    "",
    "The same markers as per-group distributions. [@cr-overview]{#distribution}",
    "",
    "And how composition shifts across conditions. [@cr-overview]{#composition}",
    "",
    ":::",
    "",
    "```{r}",
    "#| context: server",
    "scroll_app_server(\"main\", dir = project_dir)",
    "```"
  )
}

#' Read a scroll project config
#'
#' @param dir A scroll project directory.
#' @return The parsed config as a list.
#' @export
scroll_config <- function(dir) {
  path <- file.path(dir, "config.yaml")
  if (!file.exists(path)) stop("No config.yaml in '", dir, "'.", call. = FALSE)
  yaml::read_yaml(path)
}
