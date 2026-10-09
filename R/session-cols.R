# --- New column from a level mapping (session-only) ------------------------------

# A "New column" section in the right-hand rail: pick a categorical column, type a
# mapping from its levels to new labels, and the result becomes a categorical column
# in every panel for this session (the same registry as Signature's "Add as
# column"; nothing is written to disk). Mapping syntax, entries separated by commas
# or new lines:
#
#   [D1, D2]: 2, [D3, D4]: 1      several levels -> one label
#   D5: 3                         one level -> a label
#   *: other                      every level not listed (default: missing)
#
# config `session_columns: false` hides the section.

.scroll_newcol_on <- function(data) !isFALSE(data$config$session_columns)

# Columns a mapping can start from: every categorical column, including each subset
# view's own columns (e.g. its cluster ids; cells outside the subset get no label)
# and the session's added columns.
.scroll_newcol_cols <- function(m)
  names(Filter(function(e) identical(e$type, "categorical"), .scroll_with_derived(m)$meta))

.scroll_newcol_ui <- function(data, ns = identity) {
  if (!.scroll_newcol_on(data)) return(NULL)
  cats <- .scroll_newcol_cols(data$manifest)
  if (!length(cats)) return(NULL)
  .scroll_details("New column", open = FALSE,
    div(class = "scroll-filter",
        selectInput(ns("scroll_newcol_from"), "From column", stats::setNames(cats, cats))),
    uiOutput(ns("scroll_newcol_levels")),
    div(class = "scroll-filter",
        textAreaInput(ns("scroll_newcol_map"), "Mapping", rows = 3,
                      placeholder = "[D1, D2]: 2, [D3, D4]: 1\n*: other")),
    div(class = "scroll-filter",
        textInput(ns("scroll_newcol_name"), "Column name", placeholder = "e.g. genotype")),
    div(class = "scroll-apply-row",
        actionButton(ns("scroll_newcol_add"), "Add", class = "btn-sm btn-primary")),
    uiOutput(ns("scroll_newcol_list")))
}

# Parse a mapping against a column's levels. Returns list(map = named character
# (level -> label), default = label for unlisted levels or NA). Stops with a
# readable message on anything it can't use.
.scroll_parse_level_map <- function(text, levels) {
  text <- paste(text %||% "", collapse = "\n")
  if (!nzchar(trimws(text))) stop("Type a mapping, e.g. [D1, D2]: 2, [D3, D4]: 1", call. = FALSE)
  unq <- function(x) gsub("^[\"'`]+|[\"'`]+$", "", trimws(x))
  # an entry is "[a, b]: label" or "a: label"; commas inside [] belong to the group
  rx <- "(\\[[^]]*\\]|[^,:\n\\[\\]]+)\\s*:\\s*([^,\n\\[\\]]+)"
  hits <- regmatches(text, gregexpr(rx, text, perl = TRUE))[[1]]
  rest <- gsub("[\\s,\\[\\]]", "", gsub(rx, "", text, perl = TRUE), perl = TRUE)
  if (!length(hits) || nzchar(rest))
    stop("Couldn't read the mapping", if (nzchar(rest)) paste0(" near '", substr(rest, 1, 20), "'"),
         ". Use [level, level]: label, separated by commas.", call. = FALSE)
  map <- character(0); default <- NA_character_
  for (h in hits) {
    lhs <- sub(rx, "\\1", h, perl = TRUE); lab <- unq(sub(rx, "\\2", h, perl = TRUE))
    if (!nzchar(lab)) stop("An entry has an empty label: ", trimws(h), call. = FALSE)
    vals <- unq(strsplit(gsub("^\\[|\\]$", "", trimws(lhs)), ",", fixed = TRUE)[[1]])
    vals <- vals[nzchar(vals)]
    if (identical(vals, "*")) { default <- lab; next }
    # exact level names first, then a case-insensitive match
    hit <- match(vals, levels)
    ci <- is.na(hit)
    hit[ci] <- match(tolower(vals[ci]), tolower(levels))
    if (anyNA(hit))
      stop("Not a level of this column: ", paste(vals[is.na(hit)], collapse = ", "), call. = FALSE)
    lv <- levels[hit]
    dup <- lv[lv %in% names(map) & map[lv] != lab]
    if (length(dup)) stop("Mapped twice: ", paste(unique(dup), collapse = ", "), call. = FALSE)
    map[lv] <- lab
  }
  if (!length(map) && is.na(default)) stop("The mapping names no levels.", call. = FALSE)
  list(map = map, default = default)
}

# The new column's values, aligned to `src` (the source column's values).
.scroll_apply_level_map <- function(src, parsed) {
  out <- unname(parsed$map[as.character(src)])
  out[is.na(out) & !is.na(src)] <- parsed$default
  out
}

# Server side: live level hint, Add (replaces a session column of the same name),
# and removal. `derived` is the session registry from .scroll_wire.
.scroll_bind_newcol <- function(input, output, session, data, derived) {
  if (!.scroll_newcol_on(data) || !length(.scroll_newcol_cols(data$manifest))) return(invisible())
  m <- data$manifest
  src_values <- function(col)
    data$cells[[col]] %||% .scroll_derived_cols(m)[[col]]$values

  # the column menu follows the session's added columns (map a mapped column again)
  observeEvent(derived()$version, {
    cats <- .scroll_newcol_cols(m)
    cur <- isolate(input$scroll_newcol_from)
    updateSelectInput(session, "scroll_newcol_from", choices = stats::setNames(cats, cats),
                      selected = if (!is.null(cur) && cur %in% cats) cur else cats[[1]])
  }, ignoreInit = TRUE)

  output$scroll_newcol_levels <- renderUI({
    col <- input$scroll_newcol_from
    req(col)
    lv <- .scroll_meta_levels(data, col)
    shown <- utils::head(lv, 30)
    helpText(class = "scroll-newcol-levels",
             paste0("Levels: ", paste(shown, collapse = ", "),
                    if (length(lv) > 30) sprintf(" (+%d more)", length(lv) - 30) else ""))
  })

  output$scroll_newcol_list <- renderUI({
    added <- names(derived()$cols)
    if (!length(added)) return(NULL)
    tagList(
      helpText(paste("Added this session:", paste(added, collapse = ", "))),
      selectizeInput(session$ns("scroll_newcol_drop"), NULL, added, multiple = TRUE,
                     options = list(placeholder = "Columns to remove...")),
      actionLink(session$ns("scroll_newcol_remove"), "Remove selected"))
  })

  observeEvent(input$scroll_newcol_add, {
    col <- input$scroll_newcol_from
    nm <- .scroll_sig_clean_name(input$scroll_newcol_name)
    fail <- function(msg) showNotification(msg, type = "error", duration = 6)
    if (is.null(col) || !nzchar(col)) return(fail("Pick a column to map from."))
    if (!nzchar(nm)) return(fail("Give the new column a name."))
    if (nm %in% names(data$cells))
      return(fail(paste0("'", nm, "' is already a dataset column; pick another name.")))
    src <- src_values(col)
    parsed <- tryCatch(.scroll_parse_level_map(input$scroll_newcol_map,
                                               .scroll_meta_levels(data, col)),
                       error = function(e) e)
    if (inherits(parsed, "error")) return(fail(conditionMessage(parsed)))
    vals <- .scroll_apply_level_map(src, parsed)
    cur <- derived(); cols <- cur$cols
    cols[[nm]] <- list(values = vals, meta = .scroll_meta_entry(vals))
    derived(list(cols = cols, version = cur$version + 1L))
    unlisted <- setdiff(.scroll_meta_levels(data, col), names(parsed$map))
    note <- if (length(unlisted) && is.na(parsed$default))
      sprintf(" %d level%s not listed (%s) %s missing.", length(unlisted),
              if (length(unlisted) == 1) "" else "s",
              paste(utils::head(unlisted, 5), collapse = ", "),
              if (length(unlisted) == 1) "is" else "are")
    else ""
    showNotification(paste0("Added '", nm, "': pick it in the panels' column menus.", note),
                     duration = 6)
  })

  observeEvent(input$scroll_newcol_remove, {
    drop <- input$scroll_newcol_drop
    req(length(drop))
    cur <- derived()
    derived(list(cols = cur$cols[setdiff(names(cur$cols), drop)], version = cur$version + 1L))
  })
  invisible()
}
