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

  sections <- .scroll_build_sections(default_assay, default_embedding, cats, features)
  .scroll_write_config(outdir, default_assay, default_embedding, sections)
  .scroll_write_story(outdir, sections)
  invisible(outdir)
}

# Prefer recognisable marker genes for the demo panels, else fall back.
.scroll_pick_features <- function(features, n = 6) {
  markers <- c("CD3D", "CD8A", "IL7R", "CCR7", "MS4A1", "CD79A", "CD14", "LYZ",
               "FCGR3A", "NKG7", "GNLY", "PPBP")
  hit <- intersect(markers, features)
  if (length(hit) >= 3) utils::head(hit, n) else utils::head(features, n)
}

# Build the ordered example sections from what the data supports.
.scroll_build_sections <- function(assay, embedding, cats, features) {
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
  sections
}

.scroll_write_config <- function(outdir, assay, embedding, sections) {
  config <- list(
    title = "A scroll story",
    default_assay = assay,
    default_embedding = embedding,
    sections = sections
  )
  path <- file.path(outdir, "config.yaml")
  if (!file.exists(path)) yaml::write_yaml(config, path)
}

.scroll_write_story <- function(outdir, sections) {
  path <- file.path(outdir, "story.qmd")
  if (file.exists(path)) return(invisible())
  writeLines(.scroll_story_template(sections), path)
}

# Default narrative copy per view type (authors edit these freely).
.scroll_trigger_prose <- function(view) switch(view,
  umap_colorby = "The cells, laid out by their embedding and coloured by annotation. Scroll on.",
  feature_plot = "Search any gene in the box above and the pinned view recolours live to its expression.",
  dotplot      = "A marker panel: dot size is the fraction of cells expressing, colour the mean expression, across groups.",
  violin       = "The same signal as a per-group distribution.",
  proportions  = "How composition shifts across groups.",
  de_table     = "Differential expression for this contrast.",
  "This section of the story."
)

# A closeread + Shiny narrative built from the config sections. One persistent
# sticky (`#cr-sticky`) holds scroll_app_ui(); every trigger focuses it via
# `@cr-sticky` and carries a `data-section` span that cr-bridge.js reports to
# Shiny as `active_section`. The section ids here match config.yaml.
.scroll_story_template <- function(sections) {
  head <- c(
    "---",
    "title: \"A scroll story\"",
    "format:",
    "  closeread-html:",
    "    cr-style:",
    "      narrative-background-color-overlay: \"#1a1a1a\"",
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
    ""
  )
  triggers <- unlist(lapply(sections, function(s) c(
    sprintf("%s <span data-section=\"%s\"></span> @cr-sticky",
            .scroll_trigger_prose(s$view), s$id),
    ""
  )))
  tail <- c(
    ":::",
    "",
    "```{r}",
    "#| context: server",
    "scroll_app_server(\"main\", dir = project_dir)",
    "```"
  )
  c(head, triggers, tail)
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
