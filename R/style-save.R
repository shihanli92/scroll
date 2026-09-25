# Saved plot styles: <project>/style.yaml.
#
# The Style sheet's Save button writes every plot whose style differs from what is
# on disk into `style.yaml` next to config.yaml; each session reads the file when it
# starts, so a saved look persists across sessions (and for every viewer of a
# deployed app). One entry per panel id:
#
#   scroll_style: 1
#   panels:
#     dimplot:
#       theme:  {legend: bottom}            # the sheet families, as .scroll_style_families()
#       labels: {title: My UMAP}
#       scales: {xmin: '-5'}
#       layers: {Cells: {size: '1.2'}}      # per-layer overrides (R/style-layers.R)
#       plot:   {palette: Set2, aspect: 1.2} # the panel's own look controls (its style_ui)
#       manual: {col: {B: '#1F77B4'}}       # Manual-palette colours, by level name
#
# The file is hand-editable and read defensively: YAML expressions are never
# evaluated, unknown keys are dropped, and a malformed file reads as "no styles"
# (with a warning) rather than stopping the app. `style_save: false` in config.yaml
# hides the Save button (e.g. on a shared deployment).

.SCROLL_STYLE_FILE <- "style.yaml"
.SCROLL_STYLE_MAX_BYTES <- 512 * 1024
.SCROLL_STYLE_MAX_VALUES <- 50L        # longest vector kept for one plot control

.scroll_style_path <- function(dir) file.path(dir, .SCROLL_STYLE_FILE)

# Every scale key any panel can have (a panel keeps only those its caps expose).
.scroll_all_scale_keys <- function()
  .scroll_scale_keys(list(limits = c("x", "y"), trans = list(x = "log10", y = "log10"),
                          breaks = c("x", "y"), flip = TRUE, facet = TRUE))

# A YAML scalar as a non-empty string (labels truncated), else NULL.
.scroll_style_scalar <- function(v) {
  if (!is.atomic(v) || length(v) != 1L || is.na(v)) return(NULL)
  v <- substr(as.character(v), 1L, .SCROLL_LABEL_MAX)
  if (nzchar(v)) v
}

# One panel's entry, validated: only known families / keys, scalar values.
.scroll_style_clean <- function(p) {
  if (!is.list(p) || is.null(names(p))) return(list())
  out <- list()
  keys <- list(theme = .scroll_theme_keys(), labels = .scroll_label_keys(),
               scales = .scroll_all_scale_keys())
  for (f in names(keys)) {
    fam <- p[[f]]
    if (!is.list(fam) || is.null(names(fam))) next
    fam <- lapply(fam, .scroll_style_scalar)
    v <- .scroll_family_compact(fam[!vapply(fam, is.null, logical(1))], keys[[f]])
    if (length(v)) out[[f]] <- v
  }
  if (is.list(p$layers) && !is.null(names(p$layers))) {
    ly <- lapply(p$layers, function(l) {
      if (!is.list(l) || is.null(names(l))) return(list())
      l <- lapply(l, .scroll_style_scalar)
      l[!vapply(l, is.null, logical(1))]
    })
    ly <- .scroll_layers_norm(ly)
    if (length(ly)) out$layers <- ly
  }
  if (is.list(p$plot) && !is.null(names(p$plot))) {
    ok <- function(v) is.atomic(v) && length(v) >= 1L && length(v) <= .SCROLL_STYLE_MAX_VALUES &&
      !anyNA(v) && (!is.character(v) || all(nchar(v) <= .SCROLL_LABEL_MAX))
    pl <- p$plot[grepl("^[A-Za-z][A-Za-z0-9_.]*$", names(p$plot))]
    pl <- pl[vapply(pl, ok, logical(1))]
    if (length(pl)) out$plot <- pl
  }
  if (is.list(p$manual) && !is.null(names(p$manual))) {
    mn <- lapply(p$manual, function(m) {
      if (!is.list(m) || is.null(names(m))) return(NULL)
      m <- m[vapply(m, function(v) is.character(v) && length(v) == 1L && !is.na(v) &&
                                     .scroll_is_colour(v), logical(1))]
      if (length(m)) vapply(m, as.character, "")
    })
    mn <- mn[grepl("^[A-Za-z][A-Za-z0-9_]*$", names(mn)) & !vapply(mn, is.null, logical(1))]
    if (length(mn)) out$manual <- mn
  }
  out
}

# All saved styles for a project: a named list, panel id -> entry. Never errors.
.scroll_read_styles <- function(dir) {
  path <- .scroll_style_path(dir)
  if (is.null(dir) || !file.exists(path)) return(list())
  bad <- function(why) {
    warning("scroll: ignoring ", path, " (", why, ").", call. = FALSE)
    list()
  }
  if (isTRUE(file.size(path) > .SCROLL_STYLE_MAX_BYTES)) return(bad("file too large"))
  x <- tryCatch(yaml::read_yaml(path, eval.expr = FALSE), error = function(e) e)
  if (inherits(x, "error")) return(bad(conditionMessage(x)))
  if (!is.list(x) || !length(x$panels)) return(list())
  if (!is.list(x$panels) || is.null(names(x$panels))) return(bad("`panels:` is not a map"))
  out <- lapply(x$panels, .scroll_style_clean)
  out[vapply(out, length, integer(1)) > 0]
}

# Write all saved styles (atomically: a temp file renamed over the old one).
.scroll_write_styles <- function(dir, styles) {
  styles <- styles[vapply(styles, length, integer(1)) > 0]
  styles <- styles[order(names(styles))]
  # named colour vectors -> maps (as.yaml drops the names of an atomic vector)
  styles <- lapply(styles, function(p) {
    if (length(p$manual)) p$manual <- lapply(p$manual, as.list)
    p
  })
  body <-yaml::as.yaml(list(scroll_style = 1L, panels = styles), unicode = TRUE)
  tmp <- tempfile(".style-", tmpdir = dir, fileext = ".yaml")
  con <- file(tmp, open = "wb")
  writeLines(enc2utf8(c(
    "# scroll plot styles -- written by the Style sheet's Save button, read when the app starts.",
    "# Hand-editable. Delete a panel's entry (or this file) to go back to the default look.",
    body)), con, sep = "\n", useBytes = TRUE)
  close(con)
  if (!file.rename(tmp, .scroll_style_path(dir))) {
    unlink(tmp)
    stop("could not replace ", .scroll_style_path(dir), call. = FALSE)
  }
  invisible(.scroll_style_path(dir))
}

# Can this process write style.yaml? (A deployed app often runs as a user that can't.)
.scroll_style_writable <- function(dir) {
  path <- .scroll_style_path(dir)
  dir.exists(dir) && file.access(dir, 2L) == 0L &&
    (!file.exists(path) || file.access(path, 2L) == 0L)
}

.scroll_style_save_on <- function(data) !isFALSE(data$config$style_save) && !is.null(data$dir)

# Canonical comparison of two entries (YAML text), so an entry read back from disk
# and the same style snapshotted from live inputs compare equal despite types.
.scroll_style_same <- function(a, b)
  identical(yaml::as.yaml(a %||% list()), yaml::as.yaml(b %||% list()))

# The session-level save state for one dataset: what is on disk (the baseline a
# panel is "unsaved" against) and each panel's snapshot function.
.scroll_style_ctx <- function(data) {
  ctx <- new.env(parent = emptyenv())
  ctx$saved <- .scroll_read_styles(data$dir)
  ctx$snap <- list()
  ctx
}

# Save every panel whose style differs from its saved entry. Other panels' entries
# are re-read from disk first, so two sessions saving different plots don't clobber
# each other. Returns the ids written.
.scroll_save_styles <- function(data, ctx) {
  if (!.scroll_style_writable(data$dir))
    stop("the app can't write to the project folder, so the style can't be saved here.",
         call. = FALSE)
  cur <- lapply(ctx$snap, function(f) f())
  dirty <- names(cur)[!vapply(names(cur), function(k)
    .scroll_style_same(cur[[k]], ctx$saved[[k]]), logical(1))]
  if (!length(dirty)) return(character(0))
  disk <- .scroll_read_styles(data$dir)
  for (k in dirty) { disk[[k]] <- cur[[k]]; ctx$saved[[k]] <- cur[[k]] }
  .scroll_write_styles(data$dir, disk)
  dirty
}

# The input ids of a panel's look controls (its style_ui), found by rendering the UI
# once under a placeholder namespace. Wrapper ids (labels, uiOutputs) are harmless:
# they never carry an input value, so they are skipped when snapshotting.
.scroll_style_plot_ids <- function(sec, data) {
  if (!is.function(sec$style_ui)) return(character(0))
  tag <- tryCatch(sec$style_ui("scrollstyleid", data), error = function(e) NULL)
  if (is.null(tag)) return(character(0))
  html <- paste(as.character(tag), collapse = "")
  m <- regmatches(html, gregexpr('id="scrollstyleid-[^"]+"', html))[[1]]
  unique(sub('^id="scrollstyleid-(.*)"$', "\\1", m))
}

# Keep a project's style.yaml across an overwrite rebuild (which wipes the folder).
.scroll_keep_style <- function(dir) {
  path <- .scroll_style_path(dir)
  if (file.exists(path)) readBin(path, "raw", file.size(path))
}
.scroll_restore_style <- function(dir, bytes)
  if (length(bytes)) writeBin(bytes, .scroll_style_path(dir))
