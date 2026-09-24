# Per-plot Style sheet: the "how it looks" options for one panel, opened from a
# paintbrush button in the plot toolbar and shown as a side sheet over the page.
#
# State model. Each panel owns one `reactiveVal` holding its style, created in
# .scroll_wire and handed to the panel as `style_r` (and, for older servers, as
# `theme_r = style_r()$theme`). The value is a plain list of families:
#   list(theme = list(<key> = <value>, ...), labels = list(...))
# (see .scroll_style_families()); an absent key means "no override". The sheet's
# inputs are inserted lazily the first time it opens, seeded from that value; from then on a debounced snapshot of the inputs feeds
# the value, and the value feeds back into the inputs (for "apply to all" and
# Reset), with an equality check on both sides so the two never ping-pong.
#
# Sheets live OUTSIDE the panel cards (see .scroll_body): the cards set
# container-type, which makes them the containing block for position:fixed
# descendants, so a sheet inside a card would be clipped to it.

# ---- theme family -----------------------------------------------------------

# The theme() override keys, in display order. Each maps to one block in
# .scroll_ggtheme(); every one defaults to "no override".
.scroll_theme_keys <- function()
  c("font", "font_family", "text_colour", "title_size", "title_style",
    "axis_title_size", "axis_text_size", "legend_text_size", "strip_text_size",
    "legend", "legend_dir", "legend_title", "legend_key",
    "axes", "axis_titles", "axis_ticks", "axis_line", "axis_colour", "angle", "yangle",
    "grid_major", "grid_minor", "grid_colour", "line_size", "line_colour",
    "border", "border_colour", "strip_bg", "margin")

# which theme keys are colour pickers (updated differently from the selects)
.scroll_theme_colour_keys <- function()
  c("text_colour", "axis_colour", "grid_colour", "line_colour", "border_colour",
    "strip_bg")

.SCROLL_THEME_PREFIX <- "st_theme_"

# Theme colour pickers start at the colour the plots already use (scroll's base
# theme: black text / lines / border, white gridlines / strip), and a picker left at
# its default means "no override" -- so an untouched sheet changes nothing. (A
# colourpicker cannot be empty: colourpicker >= 1.3 reports "" as "#FFFFFF", which
# once silently painted text white.) The pickers are opaque (no alpha slider).
.SCROLL_THEME_COLOUR_DEFAULTS <- c(text_colour = "#000000", axis_colour = "#000000",
                                   grid_colour = "#FFFFFF", line_colour = "#000000",
                                   border_colour = "#000000", strip_bg = "#FFFFFF")
# Legacy "no override" markers (fully-transparent values from 0.3.1's pickers).
.SCROLL_NO_COLOUR <- "#FFFFFF00"
.scroll_is_no_colour <- function(v)
  is.character(v) && length(v) == 1L && !is.na(v) &&
    toupper(v) %in% c(.SCROLL_NO_COLOUR, "TRANSPARENT")
# A picked colour as stored: NULL for a legacy no-override marker; a zero-alpha pick
# (0.3.1 pickers kept the transparent default's alpha) is read as the opaque colour.
.scroll_colour_value <- function(v) {
  if (.scroll_is_no_colour(v)) return(NULL)
  if (grepl("^#[0-9A-Fa-f]{6}00$", v)) substr(v, 1L, 7L) else v
}
# A theme colour key's value as stored: NULL when it is the key's default.
.scroll_theme_colour <- function(k, v) {
  v <- .scroll_colour_value(v)
  if (is.null(v) || identical(toupper(v), .SCROLL_THEME_COLOUR_DEFAULTS[[k]])) NULL else v
}

# Snapshot the theme inputs into a named list over every key (NULL = no override),
# the shape .scroll_ggtheme() reads.
.scroll_theme_values <- function(input, prefix = .SCROLL_THEME_PREFIX)
  stats::setNames(lapply(.scroll_theme_keys(),
    function(k) .scroll_nz(input[[paste0(prefix, k)]])), .scroll_theme_keys())

# Canonical form of a theme list: only set keys, in key order. Two themes that
# mean the same thing compare identical(), which is what stops the input <-> value
# sync from re-rendering a plot for a no-op.
.scroll_theme_compact <- function(x) .scroll_family_compact(x, .scroll_theme_keys())

# The Theme section of a sheet. `values` seeds each control (lazy insertion
# restores whatever the panel's style already holds, e.g. after "apply to all").
.scroll_theme_inputs <- function(ns, values = list(), prefix = .SCROLL_THEME_PREFIX) {
  v <- function(k) values[[k]] %||% ""
  id <- function(k) ns(paste0(prefix, k))
  sel <- function(k, label, choices)
    div(class = "scroll-filter",
        selectInput(id(k), label, choices, selected = v(k), selectize = FALSE))
  # colour picker (empty = no override); a hex text field without colourpicker
  col <- function(k, label)
    div(class = "scroll-filter scroll-colour",
        if (requireNamespace("colourpicker", quietly = TRUE))
          colourpicker::colourInput(id(k), label,
                                    value = .scroll_nz(v(k)) %||% .SCROLL_THEME_COLOUR_DEFAULTS[[k]],
                                    showColour = "background")
        else textInput(id(k), label, value = v(k), placeholder = "#hex or name"))
  sz <- function(s, m, l) c("Default" = "", "Small" = s, "Medium" = m, "Large" = l)
  showhide <- c("Default" = "", "Show" = "show", "Hide" = "hide")
  onoff    <- c("Default" = "", "On" = "on", "Off" = "off")
  .scroll_details("Theme", open = TRUE,
    .scroll_details("Text & fonts", open = TRUE,
      sel("font", "Text size", sz("11", "13", "16")),
      sel("font_family", "Font",
          c("Default" = "", "Sans" = "sans", "Serif" = "serif", "Mono" = "mono")),
      col("text_colour", "Text colour"),
      sel("title_size", "Title size", sz("14", "18", "22")),
      sel("title_style", "Title style",
          c("Default" = "", "Plain" = "plain", "Bold" = "bold", "Italic" = "italic")),
      sel("axis_title_size", "Axis-title size", sz("11", "13", "15")),
      sel("axis_text_size", "Axis-text size", sz("9", "11", "13")),
      sel("legend_text_size", "Legend-text size", sz("9", "11", "13")),
      sel("strip_text_size", "Strip-text size", sz("10", "12", "14"))),
    .scroll_details("Legend", open = FALSE,
      sel("legend", "Position",
          c("Default" = "", "Right" = "right", "Left" = "left", "Top" = "top",
            "Bottom" = "bottom", "Hidden" = "none")),
      sel("legend_dir", "Direction",
          c("Default" = "", "Horizontal" = "horizontal", "Vertical" = "vertical")),
      sel("legend_title", "Legend title", showhide),
      sel("legend_key", "Key background",
          c("Default" = "", "White" = "white", "None" = "none"))),
    .scroll_details("Axes", open = FALSE,
      sel("axes", "Axis text", showhide),
      sel("axis_titles", "Axis titles", showhide),
      sel("axis_ticks", "Axis ticks", showhide),
      sel("axis_line", "Axis lines", showhide),
      col("axis_colour", "Axis colour"),                 # shared by lines + ticks
      sel("angle", "X label angle", c("Default" = "", "0" = "0", "45" = "45", "90" = "90")),
      sel("yangle", "Y label angle", c("Default" = "", "0" = "0", "90" = "90"))),
    .scroll_details("Panel", open = FALSE,
      sel("grid_major", "Major gridlines", onoff),
      sel("grid_minor", "Minor gridlines", onoff),
      col("grid_colour", "Gridline colour"),
      sel("line_size", "Line thickness",
          c("Default" = "", "Thin" = "thin", "Medium" = "medium", "Thick" = "thick")),
      col("line_colour", "Line colour"),
      sel("border", "Panel border", onoff),
      col("border_colour", "Border colour")),
    .scroll_details("Facets & spacing", open = FALSE,
      col("strip_bg", "Strip background"),
      sel("margin", "Plot margin",
          c("Default" = "", "Compact" = "compact", "Normal" = "normal", "Roomy" = "roomy"))))
}

# ---- labels family ----------------------------------------------------------

.scroll_label_keys <- function() c("title", "subtitle", "caption", "x", "y", "legend")
.SCROLL_LABEL_MAX <- 200L       # characters; a longer label is truncated

.scroll_label_inputs <- function(ns, values = list(), prefix = "st_labels_") {
  txt <- function(k, label, ph)
    div(class = "scroll-filter",
        textInput(ns(paste0(prefix, k)), label, value = values[[k]] %||% "", placeholder = ph))
  .scroll_details("Titles & labels", open = TRUE,
    txt("title", "Title", "none"),
    txt("subtitle", "Subtitle", "none"),
    txt("caption", "Caption", "none"),
    txt("legend", "Legend title", "automatic"),
    txt("x", "X axis label", "automatic"),
    txt("y", "Y axis label", "automatic"))
}

# ---- scales & axes family ---------------------------------------------------
# What a panel's scales allow is declared per panel as `style_caps` (a spec field):
#   list(limits = c("x", "y"),                  axes whose visible range can be set
#        trans  = list(x = c(...), y = c(...)), transforms offered per axis
#        breaks = c("x", "y"),                  axes with a break-count control
#        flip   = FALSE,                        offer "flip coordinates"
#        facet  = FALSE)                        offer fixed / free facet scales
# The default (a panel that declares nothing) is visible-range limits on both axes,
# which is safe for any plot; a panel opts into the rest. `list()` = no scales.
.SCROLL_DEFAULT_CAPS <- list(limits = c("x", "y"))
.SCROLL_TRANS_LABELS <- c("log10" = "Log10", "sqrt" = "Square root", "log1p" = "log(1 + x)",
                          "pseudo_log" = "Pseudo-log (signed)", "reverse" = "Reversed")

.scroll_style_caps <- function(sec) if (is.null(sec$style_caps)) .SCROLL_DEFAULT_CAPS else sec$style_caps

# The scale keys a panel's caps expose (the family's keys, in display order).
.scroll_scale_keys <- function(caps) {
  k <- character(0)
  for (a in c("x", "y")) {
    if (a %in% caps$limits) k <- c(k, paste0(a, c("min", "max")))
    if (length(caps$trans[[a]])) k <- c(k, paste0(a, "trans"))
    if (a %in% caps$breaks) k <- c(k, paste0(a, "breaks"))
  }
  if (isTRUE(caps$flip)) k <- c(k, "flip")
  if (isTRUE(caps$facet)) k <- c(k, "facet")
  k
}

.scroll_scale_inputs <- function(ns, values = list(), caps = .SCROLL_DEFAULT_CAPS,
                                 prefix = "st_scales_") {
  v <- function(k) values[[k]] %||% ""
  id <- function(k) ns(paste0(prefix, k))
  sel <- function(k, label, choices)
    div(class = "scroll-filter", selectInput(id(k), label, choices, selected = v(k), selectize = FALSE))
  num <- function(k, label)
    div(class = "scroll-filter", textInput(id(k), label, value = v(k), placeholder = "auto"))
  brk <- c("Default" = "", stats::setNames(c("3", "4", "5", "6", "8", "10"), c("3", "4", "5", "6", "8", "10")))
  axis <- function(a, name) {
    rows <- list(
      if (a %in% caps$limits) div(class = "scroll-sheet-pair",
                                  num(paste0(a, "min"), paste(name, "from")),
                                  num(paste0(a, "max"), paste(name, "to"))),
      if (length(caps$trans[[a]]))
        sel(paste0(a, "trans"), paste(name, "scale"),
            c("Linear" = "", stats::setNames(caps$trans[[a]], .SCROLL_TRANS_LABELS[caps$trans[[a]]]))),
      if (a %in% caps$breaks) sel(paste0(a, "breaks"), paste(name, "breaks"), brk))
    rows[!vapply(rows, is.null, logical(1))]
  }
  body <- c(axis("x", "X"), axis("y", "Y"),
    list(if (isTRUE(caps$flip))
           sel("flip", "Orientation", c("Default" = "", "Flipped (swap axes)" = "flip")),
         if (isTRUE(caps$facet))
           sel("facet", "Facet scales", c("Default" = "", "Fixed" = "fixed", "Free" = "free",
                                          "Free x" = "free_x", "Free y" = "free_y"))))
  body <- body[!vapply(body, is.null, logical(1))]
  if (!length(body)) return(NULL)
  do.call(.scroll_details, c(list("Scales & axes", open = FALSE), body))
}

# ---- family registry --------------------------------------------------------
# Each family: its input-id prefix, keys, the snapshot debounce (text is slower so
# typing a title doesn't redraw per keystroke), its sheet section, and how to
# push a value back into one of its controls.
.scroll_style_families <- function(caps = .SCROLL_DEFAULT_CAPS) list(
  labels = list(prefix = "st_labels_", keys = .scroll_label_keys(), delay = 500,
                ui = .scroll_label_inputs,
                update = function(session, id, key, value) updateTextInput(session, id, value = value)),
  theme  = list(prefix = .SCROLL_THEME_PREFIX, keys = .scroll_theme_keys(), delay = 150,
                ui = .scroll_theme_inputs,
                update = function(session, id, key, value) {
                  if (!key %in% .scroll_theme_colour_keys()) updateSelectInput(session, id, selected = value)
                  else if (requireNamespace("colourpicker", quietly = TRUE))
                    colourpicker::updateColourInput(session, id,
                      value = .scroll_nz(value) %||% .SCROLL_THEME_COLOUR_DEFAULTS[[key]])
                  else updateTextInput(session, id, value = value)
                }),
  scales = list(prefix = "st_scales_", keys = .scroll_scale_keys(caps), delay = 500,
                ui = function(ns, values) .scroll_scale_inputs(ns, values, caps),
                update = function(session, id, key, value) {
                  if (grepl("(min|max)$", key)) updateTextInput(session, id, value = value)
                  else updateSelectInput(session, id, selected = value)
                }))

# Canonical form of one family: only set keys, in key order, as trimmed-length
# strings. Two values that mean the same thing compare identical(), which is what
# stops the input <-> value sync from re-rendering a plot for a no-op.
.scroll_family_compact <- function(x, keys) {
  if (!length(x)) return(list())
  keys <- intersect(keys, names(x))
  out <- lapply(keys, function(k) {
    v <- x[[k]]
    if (is.null(v) || !length(v)) return(NULL)
    v <- as.character(v)[1]
    if (is.na(v) || !nzchar(v)) return(NULL)
    v <- if (k %in% names(.SCROLL_THEME_COLOUR_DEFAULTS)) .scroll_theme_colour(k, v)
         else .scroll_colour_value(v)
    if (is.null(v)) NULL else substr(v, 1L, .SCROLL_LABEL_MAX)
  })
  names(out) <- keys
  out[!vapply(out, is.null, logical(1))]
}

# ---- sheet UI ---------------------------------------------------------------

# The paintbrush button for a plot toolbar. `ns` is the panel's module namespace;
# the sheet and the first-open signal are addressed through it.
.scroll_style_button <- function(ns)
  tags$button(class = "scroll-style-btn", type = "button",
              `data-sheet` = ns("sheet"), `data-open` = ns("style_open"),
              `aria-controls` = ns("sheet"), `aria-expanded` = "false",
              `data-tip` = .SCROLL_TIPS$style, `aria-label` = .SCROLL_TIPS$style,
              onclick = "scrollToggleStyle(this)", shiny::icon("paintbrush"))

# One panel's side sheet. `ns` is the dataset namespace (identity for a flat app),
# so the panel's own namespace is NS(ns(sec$id)) -- the same one its server uses.
.scroll_style_sheet <- function(sec, data, ns = identity) {
  pns <- shiny::NS(ns(sec$id))
  geom <- if (is.function(sec$style_ui)) sec$style_ui(ns(sec$id), data)
  tags$aside(
    id = pns("sheet"), class = "scroll-sheet", role = "dialog", `aria-modal` = "false",
    `aria-label` = paste("Style:", sec$label %||% sec$id), inert = NA,
    div(class = "scroll-sheet-head",
        div(class = "scroll-sheet-titles",
            span(class = "scroll-sheet-kicker", "Style"),
            span(class = "scroll-sheet-title", sec$label %||% sec$id)),
        tags$button(class = "scroll-sheet-close", type = "button", `aria-label` = "Close",
                    onclick = "scrollCloseSheets()", shiny::icon("xmark"))),
    div(class = "scroll-sheet-body",
        if (!is.null(geom)) .scroll_details("Plot", open = TRUE, geom),
        # Theme (and later Labels / Scales) are inserted here on first open
        div(id = pns("style_body"), class = "scroll-sheet-sections")),
    div(class = "scroll-sheet-foot",
        if (!isFALSE(data$config$theme_controls))
          actionButton(pns("style_all"), "Apply theme to all plots",
                       class = "btn-sm btn-outline-primary"),
        actionLink(pns("style_reset"), "Reset")))
}

# All sheets for a dataset, in one container outside the layout grid.
.scroll_style_sheets <- function(data, panels, ns = identity)
  div(class = "scroll-sheets", lapply(panels, .scroll_style_sheet, data = data, ns = ns))

# ---- sheet server -----------------------------------------------------------

# The families a sheet edits: labels always; theme unless `theme_controls: false`;
# scales when the panel's caps expose any.
.scroll_sheet_families <- function(data, caps = .SCROLL_DEFAULT_CAPS) {
  fam <- .scroll_style_families(caps)
  if (isFALSE(data$config$theme_controls)) fam$theme <- NULL
  if (!length(fam$scales$keys)) fam$scales <- NULL
  fam
}

# Wire one panel's sheet. `rv` is this panel's style value; `all` is every panel's
# (for "apply to all"). Runs inside the panel's module namespace.
.scroll_style_server <- function(id, data, rv, all = list(rv), caps = .SCROLL_DEFAULT_CAPS) {
  force(rv); force(all); force(caps)   # called in a loop: bind this panel's values now
  shiny::moduleServer(id, function(input, output, session) {
    fams <- .scroll_sheet_families(data, caps)
    inserted <- shiny::reactiveVal(FALSE)

    # first open: insert the sections, seeded from the current style
    shiny::observeEvent(input$style_open, {
      if (isTRUE(shiny::isolate(inserted()))) return()
      st <- shiny::isolate(rv())
      ui <- lapply(names(fams), function(f) fams[[f]]$ui(session$ns, st[[f]]))
      shiny::insertUI(paste0("#", session$ns("style_body")), "beforeEnd",
                      ui = shiny::tagList(ui), immediate = TRUE)
      inserted(TRUE)
    })

    # inputs -> style, per family (debounced; and only once every control of the
    # family has bound on the client -- a bound control always reports a value, ""
    # for Default -- so a half-bound snapshot can't drop keys set by apply-to-all)
    for (f in names(fams)) local({
      fam <- fams[[f]]; fname <- f
      ids <- paste0(fam$prefix, fam$keys)
      snap <- shiny::debounce(shiny::reactive({
        if (!isTRUE(inserted())) return(NULL)
        vals <- lapply(ids, function(i) input[[i]])
        if (any(vapply(vals, is.null, logical(1)))) return(NULL)
        .scroll_family_compact(stats::setNames(vals, fam$keys), fam$keys)
      }), fam$delay)
      shiny::observeEvent(snap(), {
        cur <- shiny::isolate(rv())
        new <- snap()
        if (identical(.scroll_family_compact(cur[[fname]], fam$keys), new)) return()
        cur[[fname]] <- if (length(new)) new
        rv(cur)
      })
    })

    # style -> inputs (apply-to-all from another panel, Reset). Only keys whose
    # STORED value changed since the last sync are pushed: comparing against the
    # inputs instead would reset a field the user is still typing into (its value
    # not yet snapshotted) whenever another family's snapshot lands first.
    last <- shiny::isolate(rv())
    shiny::observeEvent(rv(), {
      st <- rv(); prev <- last; last <<- st
      if (!isTRUE(shiny::isolate(inserted()))) return()
      for (f in names(fams)) {
        fam <- fams[[f]]
        for (k in fam$keys) {
          want <- st[[f]][[k]] %||% ""
          if (identical(prev[[f]][[k]] %||% "", want)) next
          id <- paste0(fam$prefix, k)
          have <- shiny::isolate(input[[id]]) %||% ""
          if (k %in% names(.SCROLL_THEME_COLOUR_DEFAULTS))
            have <- .scroll_theme_colour(k, have) %||% ""
          if (!identical(have, want)) fam$update(session, id, k, want)
        }
      }
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$style_all, {
      th <- .scroll_theme_compact(shiny::isolate(rv())$theme)
      for (other in all) {
        cur <- shiny::isolate(other())
        if (identical(.scroll_theme_compact(cur$theme), th)) next
        cur$theme <- if (length(th)) th
        other(cur)
      }
      shiny::showNotification("Theme applied to every plot.", duration = 3)
    })

    shiny::observeEvent(input$style_reset, rv(list()))
  })
}
