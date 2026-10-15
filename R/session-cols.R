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
#   ["B, memory", "T:CD4"]: x     quote a level with a comma, colon or bracket
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
# readable message on anything it can't use. A level (or label) that holds a comma,
# colon or bracket is written in quotes: ["B, memory", 'Treg [act]']: x.
.scroll_parse_level_map <- function(text, levels) {
  text <- paste(text %||% "", collapse = "\n")
  if (!nzchar(trimws(text))) stop("Type a mapping, e.g. [D1, D2]: 2, [D3, D4]: 1", call. = FALSE)
  map <- character(0); default <- NA_character_
  for (entry in .scroll_map_split(text, c(",", "\n"))) {
    entry <- trimws(entry)
    if (!nzchar(entry)) next
    at <- .scroll_map_colons(entry)
    if (!length(at))
      stop("Couldn't read '", substr(entry, 1, 30), "'. Use [level, level]: label, ",
           "separated by commas.", call. = FALSE)
    at <- at[[length(at)]]                      # the label follows the last colon
    lhs <- trimws(substr(entry, 1, at - 1L)); lab <- .scroll_map_unquote(substr(entry, at + 1L, nchar(entry)))
    if (!nzchar(lab)) stop("An entry has an empty label: ", entry, call. = FALSE)
    if (grepl("[][]", lab) && !.scroll_map_quoted(trimws(substr(entry, at + 1L, nchar(entry)))))
      stop("Couldn't read the label '", lab, "'. Put a label with brackets in quotes.", call. = FALSE)
    vals <- if (startsWith(lhs, "[") && endsWith(lhs, "]"))
              .scroll_map_split(substr(lhs, 2L, nchar(lhs) - 1L), ",")
            else lhs
    vals <- vapply(vals, .scroll_map_unquote, "", USE.NAMES = FALSE)
    vals <- vals[nzchar(vals)]
    if (!length(vals)) stop("An entry names no levels: ", entry, call. = FALSE)
    if (identical(vals, "*")) { default <- lab; next }
    # exact level names first, then a case-insensitive match
    hit <- match(vals, levels)
    ci <- is.na(hit)
    hit[ci] <- match(tolower(vals[ci]), tolower(levels))
    if (anyNA(hit)) {
      miss <- vals[is.na(hit)]
      stop("Not a level of this column: ", paste(miss, collapse = ", "),
           if (any(grepl(":", miss, fixed = TRUE))) " (a comma missing between entries?)"
           else ". Put a level name with a comma, colon or bracket in quotes",
           call. = FALSE)
    }
    lv <- levels[hit]
    dup <- lv[lv %in% names(map) & map[lv] != lab]
    if (length(dup)) stop("Mapped twice: ", paste(unique(dup), collapse = ", "), call. = FALSE)
    map[lv] <- lab
  }
  if (!length(map) && is.na(default)) stop("The mapping names no levels.", call. = FALSE)
  list(map = map, default = default)
}

# Walk `text` tracking quotes and [ ] depth: the positions of the characters in
# `seps` that sit outside both. A stray closing bracket at the top level is ignored.
.scroll_map_scan <- function(text, seps) {
  ch <- strsplit(text, "", fixed = TRUE)[[1]]
  q <- ""; depth <- 0L; at <- integer(0); drop <- integer(0)
  for (k in seq_along(ch)) {
    c <- ch[[k]]
    if (nzchar(q)) { if (c == q) q <- ""; next }
    if (c %in% c("\"", "'", "`")) q <- c
    else if (c == "[") depth <- depth + 1L
    else if (c == "]") { if (depth > 0L) depth <- depth - 1L else drop <- c(drop, k) }
    else if (depth == 0L && c %in% seps) at <- c(at, k)
  }
  list(chars = ch, at = at, drop = drop)
}

.scroll_map_split <- function(text, seps) {
  sc <- .scroll_map_scan(text, seps)
  ch <- sc$chars; ch[sc$drop] <- ""
  cut <- c(0L, sc$at, length(ch) + 1L)
  vapply(seq_len(length(cut) - 1L), function(i)
    paste(ch[seq_len(cut[i + 1L] - cut[i] - 1L) + cut[i]], collapse = ""), "")
}

.scroll_map_colons <- function(entry) .scroll_map_scan(entry, ":")$at

.scroll_map_quoted <- function(x)
  nchar(x) >= 2L && substr(x, 1, 1) %in% c("\"", "'", "`") &&
    substr(x, nchar(x), nchar(x)) == substr(x, 1, 1)

.scroll_map_unquote <- function(x) {
  x <- trimws(gsub("]+$", "", trimws(x)))   # a stray closing bracket
  if (.scroll_map_quoted(x)) substr(x, 2L, nchar(x) - 1L) else x
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
