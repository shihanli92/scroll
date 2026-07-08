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
  color_by <- if (length(cats)) cats[[1]] else meta_cols[[1]]
  default_assay <- names(assay_info)[[1]]
  default_feature <- assay_info[[default_assay]]$features[[1]]

  .scroll_write_config(outdir, default_assay, default_embedding, color_by,
                       cats, default_feature)
  .scroll_write_story(outdir)
  invisible(outdir)
}

.scroll_write_config <- function(outdir, assay, embedding, color_by, cats,
                                 feature) {
  sections <- list(
    list(id = "overview", view = "umap_colorby",
         params = list(embedding = embedding, color_by = color_by)),
    list(id = "feature", view = "feature_plot",
         params = list(embedding = embedding, assay = assay, feature = feature))
  )
  # Add a second colorby section if a second categorical exists.
  if (length(cats) >= 2) {
    sections <- append(sections, list(
      list(id = "second", view = "umap_colorby",
           params = list(embedding = embedding, color_by = cats[[2]]))),
      after = 1)
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
    "scroll_app_ui(\"main\")",
    "```",
    ":::",
    "",
    "Welcome. Scroll to walk through the analysis; the visual on the",
    "right stays pinned and answers to both your scroll position and the",
    "feature search box above it. [@cr-overview]{#overview}",
    "",
    "As we scroll, the same cells are recolored by the next annotation.",
    "[@cr-overview]{#second}",
    "",
    "Now search any gene above - the pinned view recolors live to that",
    "feature's expression. [@cr-overview]{#feature}",
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
