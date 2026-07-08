# Thin wrappers over the Quarto CLI. The rendered story.qmd carries its Shiny
# dependencies, so the same command previews locally (the meeting machine) and
# produces the uploadable output for an open-source Shiny Server.

#' Preview a scroll project locally
#'
#' Runs `quarto preview` on the project's `story.qmd` with `server: shiny`.
#'
#' @param dir A scroll project directory.
#' @param ... Extra arguments passed to `quarto preview`.
#' @return Invisibly, the exit status.
#' @export
scroll_serve <- function(dir = ".", ...) {
  .scroll_run_quarto("preview", dir, ...)
}

#' Render a scroll project for deployment
#'
#' @param dir A scroll project directory.
#' @param ... Extra arguments passed to `quarto render`.
#' @return Invisibly, the exit status.
#' @export
scroll_render <- function(dir = ".", ...) {
  .scroll_run_quarto("render", dir, ...)
}

.scroll_run_quarto <- function(cmd, dir, ...) {
  quarto <- Sys.which("quarto")
  if (!nzchar(quarto))
    stop("Quarto CLI not found on PATH. Install it from https://quarto.org, ",
         "then run `quarto add qmd-lab/closeread` inside '", dir, "'.",
         call. = FALSE)
  story <- file.path(dir, "story.qmd")
  if (!file.exists(story))
    stop("No story.qmd in '", dir, "'. Did you run scroll_build()?", call. = FALSE)
  if (!.scroll_has_closeread(dir))
    stop("The closeread extension is not installed in '", dir, "'. Run ",
         "`quarto add qmd-lab/closeread` there first.", call. = FALSE)

  # Run from the project dir so relative paths in story.qmd resolve.
  old <- setwd(dir); on.exit(setwd(old), add = TRUE)
  status <- system2(quarto, args = c(cmd, "story.qmd", ...), wait = TRUE)
  invisible(status)
}

.scroll_has_closeread <- function(dir) {
  dir.exists(file.path(dir, "_extensions", "qmd-lab", "closeread")) ||
    dir.exists(file.path(dir, "_extensions", "closeread"))
}
