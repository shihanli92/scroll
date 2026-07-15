# App scaffolding: emit the runtime `app.R` (a one-liner that launches the
# explorer) and a starter `config.yaml` of optional defaults. The app itself
# auto-populates every control from the manifest; config only sets defaults
# (title, default reduction/assay, a curated DotPlot marker list).

#' Scaffold the runtime app into a project
#'
#' Writes `app.R` (deployable to a Shiny Server) and a starter `config.yaml`
#' into a built project directory. Called by [scroll_build()]; exported so it can
#' be regenerated. Existing files are not overwritten.
#'
#' @param dir A built scroll project directory (must contain `manifest.yaml`).
#' @param panels Optional character vector of panel ids to show **exactly**, in this
#'   order, overriding the automatic data-driven gating (written to `config.yaml` as
#'   `panels:`). `NULL` (default) leaves the automatic behaviour.
#' @param exclude_panels Optional character vector of panel ids to hide (written as
#'   `exclude_panels:`); applied after gating / the `panels` allowlist.
#' @return `dir`, invisibly.
#' @export
scroll_scaffold_app <- function(dir, panels = NULL, exclude_panels = NULL) {
  man <- scroll_manifest(dir)
  assay <- man$default_assay
  features <- unlist(man$assays[[assay]]$features)

  # `title` is intentionally omitted -- set it in config.yaml to label the app
  # bar (e.g. title: "PBMC 3k"); without it the bar shows just the brand.
  config <- list(
    default_assay = assay,
    default_embedding = man$default_embedding,
    markers = as.list(.scroll_pick_features(features))
  )
  if (!is.null(panels)) config$panels <- as.list(as.character(panels))
  if (!is.null(exclude_panels)) config$exclude_panels <- as.list(as.character(exclude_panels))
  cfg_path <- file.path(dir, "config.yaml")
  if (!file.exists(cfg_path)) yaml::write_yaml(config, cfg_path)

  app_path <- file.path(dir, "app.R")
  if (!file.exists(app_path))
    writeLines(c(
      "library(scroll)",
      "",
      "# Add custom panels here with register_panel(id, ui, server, ...)",
      "# before scroll_app(); see ?register_panel.",
      "",
      "scroll_app(\".\")"), app_path)

  invisible(dir)
}

# Prefer recognisable marker genes for the default DotPlot panel, else fall back.
# Matching is case-insensitive so mouse symbols (Cd3d, Ms4a1, ...) are found too;
# the actual-case gene name is returned, ordered by the panel.
.scroll_pick_features <- function(features, n = 8) {
  markers <- c("CD3D", "CD8A", "IL7R", "CCR7", "MS4A1", "CD79A", "CD14", "LYZ",
               "FCGR3A", "NKG7", "GNLY", "PPBP", "CD19", "CD68", "ITGAX", "FOXP3")
  up <- toupper(features)
  hit <- features[up %in% markers]
  hit <- hit[!duplicated(toupper(hit))]
  hit <- hit[order(match(toupper(hit), markers))]
  if (length(hit) >= 3) utils::head(hit, n) else utils::head(features, n)
}

#' Read a scroll project config (optional defaults), or `NULL` if absent
#'
#' @param dir A scroll project directory.
#' @return The parsed config list, or `NULL` when there is no `config.yaml`.
#' @export
scroll_config <- function(dir) {
  path <- file.path(dir, "config.yaml")
  if (!file.exists(path)) return(NULL)
  yaml::read_yaml(path)
}
