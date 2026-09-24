# Per-layer controls for the Style sheet's "Layers" section. Generic: the sheet shows
# the layers of the plot actually drawn (reported by .scroll_style_plot), with a
# control for every constant aesthetic the geom supports that is NOT mapped in aes()
# (a constant would override a mapping), plus a curated set of geom / stat /
# position parameters. The overrides live in the panel's style value as
#   style$layers = list(<layer key> = list(<param> = "<value>"))
# and are applied by .scroll_apply_layers() to a CLONE of each matching layer
# (a child ggproto with new params), never mutating the plot it was given.

# ---- families + keys --------------------------------------------------------

# A stable family for a layer's geom. On-screen raster points (scattermore) and
# exported vector points share "points", so one setting drives both.
.scroll_layer_family <- function(l) {
  g <- class(l$geom)[1]
  switch(g,
    GeomPoint = , GeomScattermore = , GeomJitter = "points",
    GeomViolin = "violins", GeomBar = , GeomCol = "bars",
    GeomTile = , GeomRaster = , GeomRect = "tiles",
    GeomText = , GeomLabel = "text", GeomTextRepel = , GeomLabelRepel = "labels",
    GeomSegment = , GeomCurve = "segments", GeomPath = , GeomLine = , GeomStep = "lines",
    GeomDensity = , GeomArea = , GeomRibbon = "areas",
    GeomHline = , GeomVline = , GeomAbline = "reference lines",
    GeomBoxplot = "boxes", GeomPolygon = , GeomLogo = "shapes",
    tolower(sub("^Geom", "", g)))
}

.scroll_cap1 <- function(x) paste0(toupper(substr(x, 1, 1)), substring(x, 2))

# Tag a layer with a display role ("Cells", "Background cells", ...). Views tag the
# layers that would otherwise be ambiguous; the role becomes the layer's key and
# label, so a setting sticks to the same layer when others come and go.
.scroll_role <- function(l, role) { if (inherits(l, "LayerInstance")) l$scroll_role <- role; l }

# Keys for a plot's layers: the role when tagged, else "<Family> <n>" numbering the
# untagged layers of each family ("Points", "Points 2", ...). Layers sharing a role
# share a key (and so one setting) -- e.g. a volcano's two cutoff lines.
.scroll_layer_keys <- function(layers) {
  seen <- list()
  vapply(layers, function(l) {
    role <- tryCatch(l$scroll_role, error = function(e) NULL)
    if (is.character(role) && length(role) == 1L && nzchar(role)) return(role)
    fam <- .scroll_layer_family(l)
    n <- (seen[[fam]] %||% 0L) + 1L; seen[[fam]] <<- n
    if (n == 1L) .scroll_cap1(fam) else paste(.scroll_cap1(fam), n)
  }, character(1))
}

# ---- what can be adjusted ----------------------------------------------------

.SCROLL_SHAPES <- c("16 Filled circle" = "16", "19 Solid circle" = "19",
                    "21 Circle (fill + border)" = "21", "1 Open circle" = "1",
                    "15 Filled square" = "15", "22 Square (fill + border)" = "22",
                    "0 Open square" = "0", "17 Filled triangle" = "17",
                    "24 Triangle (fill + border)" = "24", "2 Open triangle" = "2",
                    "18 Filled diamond" = "18", "23 Diamond (fill + border)" = "23",
                    "3 Plus" = "3", "4 Cross" = "4", "8 Star" = "8")

# Constant aesthetics the Layers section can set, and their control types.
.SCROLL_LAYER_AES <- list(
  colour     = list(type = "colour", label = "Colour"),
  fill       = list(type = "colour", label = "Fill"),
  size       = list(type = "number", label = "Size", min = 0),
  linewidth  = list(type = "number", label = "Line width", min = 0),
  stroke     = list(type = "number", label = "Stroke", min = 0),
  alpha      = list(type = "number", label = "Opacity (0-1)", min = 0, max = 1),
  shape      = list(type = "select", label = "Shape", choices = .SCROLL_SHAPES),
  linetype   = list(type = "select", label = "Line type",
                    choices = c("Solid" = "solid", "Dashed" = "dashed", "Dotted" = "dotted",
                                "Dot-dash" = "dotdash", "Long dash" = "longdash",
                                "Two-dash" = "twodash", "None" = "blank")),
  fontface   = list(type = "select", label = "Font face",
                    choices = c("Plain" = "plain", "Bold" = "bold", "Italic" = "italic",
                                "Bold italic" = "bold.italic")),
  family     = list(type = "select", label = "Font",
                    choices = c("Sans" = "sans", "Serif" = "serif", "Mono" = "mono")),
  angle      = list(type = "number", label = "Angle"),
  hjust      = list(type = "number", label = "Horizontal justify"),
  vjust      = list(type = "number", label = "Vertical justify"),
  lineheight = list(type = "number", label = "Line height", min = 0),
  width      = list(type = "number", label = "Width", min = 0),
  height     = list(type = "number", label = "Height", min = 0))

.SCROLL_TF <- c("Yes" = "TRUE", "No" = "FALSE")

# Curated non-aesthetic parameters per family: slot = where the value lives.
.SCROLL_LAYER_PARAMS <- list(
  violins = list(trim  = list(slot = "stat", type = "select", label = "Trim tails", choices = .SCROLL_TF),
                 scale = list(slot = "stat", type = "select", label = "Scale widths by",
                              choices = c("Area" = "area", "Count" = "count", "Equal width" = "width"))),
  areas   = list(adjust = list(slot = "stat", type = "number", label = "Smoothing (adjust)", min = 0)),
  text    = list(check_overlap = list(slot = "geom", type = "select", label = "Hide overlapping",
                                      choices = .SCROLL_TF)),
  labels  = list(max.overlaps = list(slot = "geom", type = "number", label = "Max overlaps", min = 0),
                 box.padding = list(slot = "geom", type = "number", label = "Box padding", min = 0),
                 min.segment.length = list(slot = "geom", type = "number",
                                           label = "Min leader length", min = 0)))

# Position parameters (jittered points): name in the sheet -> field on the position.
.scroll_position_params <- function(l) {
  pos <- class(l$position)[1]
  if (identical(pos, "PositionJitter"))
    list(jitter_width = list(slot = "position", field = "width", type = "number",
                             label = "Jitter width", min = 0),
         jitter_height = list(slot = "position", field = "height", type = "number",
                              label = "Jitter height", min = 0))
  else if (identical(pos, "PositionJitterdodge"))
    list(jitter_width = list(slot = "position", field = "jitter.width", type = "number",
                             label = "Jitter width", min = 0))
  else list()
}

# Aesthetics mapped for a layer (its own mapping + the plot's, when inherited).
.scroll_mapped_aes <- function(l, p) {
  m <- names(l$mapping)
  if (isTRUE(l$inherit.aes) && !is.null(p$mapping)) m <- c(m, names(p$mapping))
  m[m == "color"] <- "colour"
  unique(m)
}

# A layer's current constant for a parameter, as a short string ("" if none set).
.scroll_layer_current <- function(l, param, spec) {
  v <- switch(spec$slot %||% "aes",
    aes = l$aes_params[[param]], geom = l$geom_params[[param]],
    stat = l$stat_params[[param]],
    position = tryCatch(l$position[[spec$field]], error = function(e) NULL))
  if (is.null(v) || !length(v) || (!is.character(v) && !is.numeric(v) && !is.logical(v))) return("")
  v <- v[[1]]
  if (is.na(v)) "" else if (is.numeric(v)) format(signif(v, 4)) else as.character(v)
}

# The controls a layer offers: list(param = spec + current).
.scroll_layer_controls <- function(l, p) {
  fam <- .scroll_layer_family(l)
  aes_ok <- intersect(names(.SCROLL_LAYER_AES), names(l$geom$default_aes))
  aes_ok <- setdiff(aes_ok, .scroll_mapped_aes(l, p))
  ctl <- lapply(stats::setNames(aes_ok, aes_ok), function(a) c(.SCROLL_LAYER_AES[[a]], list(slot = "aes")))
  ctl <- c(ctl, .SCROLL_LAYER_PARAMS[[fam]] %||% list(), .scroll_position_params(l))
  for (k in names(ctl)) ctl[[k]]$current <- .scroll_layer_current(l, k, ctl[[k]])
  ctl
}

# The main ggplot of a plot object (itself, a patchwork's last plot, an aplot's main
# panel), or NULL.
.scroll_main_ggplot <- function(p) {
  if (inherits(p, "ggplot")) return(p)
  if (inherits(p, "aplot"))
    return(tryCatch({ m <- p$plotlist[[p$layout[p$main_row, p$main_col]]]
                      if (inherits(m, "ggplot")) m else NULL }, error = function(e) NULL))
  NULL
}

# One entry per adjustable layer of the drawn plot: key, label, controls. Layers
# sharing a key collapse to one entry. Plain data (strings/lists) so two catalogs of
# the same plot compare identical().
.scroll_layer_catalog <- function(p) {
  g <- .scroll_main_ggplot(p)
  if (is.null(g) || !length(g$layers)) return(list())
  keys <- .scroll_layer_keys(g$layers)
  out <- list()
  for (i in seq_along(g$layers)) {
    k <- keys[[i]]
    if (!is.null(out[[k]])) next
    ctl <- .scroll_layer_controls(g$layers[[i]], g)
    if (length(ctl)) out[[k]] <- list(key = k, label = k, controls = ctl)
  }
  unname(out)
}

# ---- applying overrides ------------------------------------------------------

# A stored override string -> the value ggplot wants, or NULL if unusable.
.scroll_layer_value <- function(spec, v, param = "") {
  if (is.null(v) || !length(v) || !nzchar(v <- as.character(v)[1])) return(NULL)
  switch(spec$type,
    number = {
      x <- suppressWarnings(as.numeric(v))
      if (!is.finite(x) || (!is.null(spec$min) && x < spec$min) ||
          (!is.null(spec$max) && x > spec$max)) NULL else x
    },
    colour = if (isTRUE(tryCatch({ grDevices::col2rgb(v); TRUE }, error = function(e) FALSE))) v,
    select = {
      if (!v %in% spec$choices) NULL
      else if (v %in% c("TRUE", "FALSE")) as.logical(v)
      else if (identical(param, "shape")) as.integer(v)
      else v
    })
}

# A clone of layer `l` with overrides `ov` (param -> string) applied.
.scroll_override_layer <- function(l, ov, p) {
  ctl <- .scroll_layer_controls(l, p)
  aes <- list(); geom <- list(); stat <- list(); pos <- list()
  for (k in intersect(names(ov), names(ctl))) {
    spec <- ctl[[k]]; v <- .scroll_layer_value(spec, ov[[k]], k)
    if (is.null(v)) next
    switch(spec$slot,
      aes = { aes[[k]] <- v; if (k == "width" && "width" %in% names(l$stat_params)) stat$width <- v },
      geom = geom[[k]] <- v, stat = stat[[k]] <- v, position = pos[[spec$field]] <- v)
  }
  scatter <- identical(class(l$geom)[1], "GeomScattermore")
  if (scatter) {
    # raster points take a pixel size and have no shape/stroke; those two still apply
    # to the vector points drawn for exports and when rasterizing is off
    if (!is.null(aes$size)) { geom$pointsize <- max(1, round(2 * aes$size)); aes$size <- NULL }
    aes$shape <- NULL; aes$stroke <- NULL
  }
  if (!length(aes) && !length(geom) && !length(stat) && !length(pos)) return(l)
  parent <- l                                      # bind first: ggproto captures lazily
  args <- list()
  if (length(aes))  args$aes_params  <- utils::modifyList(parent$aes_params, aes)
  if (length(geom)) args$geom_params <- utils::modifyList(parent$geom_params, geom)
  if (length(stat)) args$stat_params <- utils::modifyList(parent$stat_params, stat)
  if (length(pos)) { ppos <- parent$position
                     args$position <- do.call(ggplot2::ggproto, c(list(NULL, ppos), pos)) }
  do.call(ggplot2::ggproto, c(list(NULL, parent), args))
}

# Apply `overrides` (key -> list(param -> string)) to a plot's layers. A patchwork
# gets them on every sub-plot; an aplot is returned as-is (its views apply the style
# to the main plot before attaching the dendrograms).
.scroll_apply_layers <- function(p, overrides) {
  if (!length(overrides)) return(p)
  if (inherits(p, "patchwork")) {
    p$patches$plots <- lapply(p$patches$plots, .scroll_apply_layers, overrides = overrides)
  }
  if (!inherits(p, "ggplot") || !length(p$layers)) return(p)
  keys <- .scroll_layer_keys(p$layers)
  if (!any(keys %in% names(overrides))) return(p)
  p$layers <- lapply(seq_along(p$layers), function(i) {
    ov <- overrides[[keys[[i]]]]
    if (length(ov)) .scroll_override_layer(p$layers[[i]], ov, p) else p$layers[[i]]
  })
  p
}

# Apply a panel's style to its ON-SCREEN plot and report the plot's layer catalog to
# the panel's Style sheet (the catalog reactiveVal rides on the style reactiveVal as
# an attribute; see .scroll_wire). Set only when it changed, so no reactive loops.
# Exports use plain .scroll_apply_style (nothing to report).
.scroll_style_plot <- function(p, style_r) {
  if (is.null(style_r)) return(p)
  cat_rv <- attr(style_r, "scroll_catalog")
  if (is.function(cat_rv)) {
    cat <- tryCatch(.scroll_layer_catalog(p), error = function(e) NULL)
    if (!is.null(cat)) shiny::isolate(if (!identical(cat_rv(), cat)) cat_rv(cat))
  }
  .scroll_apply_style(p, style_r())
}

# ---- sheet UI ---------------------------------------------------------------

.scroll_layer_input_id <- function(key, param)
  paste0("st_layer__", gsub("[^A-Za-z0-9]+", "_", key), "__", gsub("[^A-Za-z0-9]+", "_", param))

# The Layers section body for a catalog, seeded from stored overrides.
.scroll_layers_ui <- function(ns, catalog, overrides = list()) {
  if (!length(catalog)) return(tags$p(class = "scroll-desc", "This plot has no adjustable layers."))
  lapply(seq_along(catalog), function(i) {
    e <- catalog[[i]]; ov <- overrides[[e$key]] %||% list()
    ctls <- lapply(names(e$controls), function(k) {
      spec <- e$controls[[k]]; id <- ns(.scroll_layer_input_id(e$key, k))
      val <- ov[[k]] %||% ""; cur <- spec$current %||% ""
      switch(spec$type,
        number = div(class = "scroll-filter",
                     textInput(id, spec$label, value = val,
                               placeholder = if (nzchar(cur)) cur else "default")),
        select = div(class = "scroll-filter",
                     selectInput(id, spec$label,
                                 c(stats::setNames("", if (nzchar(cur)) paste0("Unchanged (", cur, ")")
                                                       else "Unchanged"), spec$choices),
                                 selected = val, selectize = FALSE)),
        colour = {
          mode_id <- paste0(id, "_mode")
          start <- if (nzchar(val)) val else if (nzchar(cur) && .scroll_is_colour(cur)) cur else "#000000"
          div(class = "scroll-filter scroll-layer-colour",
              selectInput(mode_id, spec$label,
                          c(stats::setNames("", if (nzchar(cur)) paste0("Unchanged (", cur, ")")
                                                else "Unchanged"), "Custom" = "custom"),
                          selected = if (nzchar(val)) "custom" else "", selectize = FALSE),
              conditionalPanel(sprintf("input['%s'] == 'custom'", mode_id),
                if (requireNamespace("colourpicker", quietly = TRUE))
                  colourpicker::colourInput(id, NULL, value = start, showColour = "background")
                else textInput(id, NULL, value = start)))
        })
    })
    .scroll_details(e$label, open = i == 1L, ctls)
  })
}

.scroll_is_colour <- function(v)
  isTRUE(tryCatch({ grDevices::col2rgb(v); TRUE }, error = function(e) FALSE))

# Read the Layers inputs for a catalog into overrides (key -> param -> string), or
# NULL while any of its controls has not bound yet.
.scroll_layers_snapshot <- function(input, catalog) {
  out <- list()
  for (e in catalog) {
    vals <- list()
    for (k in names(e$controls)) {
      id <- .scroll_layer_input_id(e$key, k)
      if (identical(e$controls[[k]]$type, "colour")) {
        mode <- input[[paste0(id, "_mode")]]
        if (is.null(mode)) return(NULL)
        if (identical(mode, "custom")) {
          v <- input[[id]]
          if (is.null(v)) return(NULL)
          vals[[k]] <- v
        }
      } else {
        v <- input[[id]]
        if (is.null(v)) return(NULL)
        if (nzchar(v)) vals[[k]] <- as.character(v)
      }
    }
    out[[e$key]] <- vals
  }
  out
}

# Canonical layer overrides: empty entries dropped, keys and params sorted, so two
# equivalent values compare identical() (no no-op redraws from the sync).
.scroll_layers_norm <- function(x) {
  if (!length(x)) return(list())
  x <- lapply(x, function(v) { v <- v[vapply(v, function(e) length(e) && nzchar(e[[1]]), logical(1))]
                               v[order(names(v))] })
  x <- x[vapply(x, length, integer(1)) > 0]
  x[order(names(x))]
}
