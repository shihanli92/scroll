# The CDR3 logo panel: amino-acid sequence logos of the CDR3 (per group, at a chosen
# length) via ggseqlogo. Gated on manifest$vdj + ggseqlogo; needs the CDR3 aa columns
# (cdr3_<chain>) baked by a current scroll_build (older builds store only lengths).

# Table of CDR3 lengths present across the baked chains (named by length).
.scroll_cdr3_length_counts <- function(data) {
  rc <- .scroll_vdj_read(data, "rep_cells.parquet")
  cols <- paste0("cdr3_", unlist(data$manifest$vdj$cdr3_chains))
  cols <- cols[cols %in% names(rc)]
  if (is.null(rc) || !length(cols)) return(integer())
  lens <- unlist(lapply(cols, function(c) nchar(rc[[c]][!is.na(rc[[c]]) & nzchar(rc[[c]])])))
  if (!length(lens)) return(integer())
  table(lens)
}
# Lengths in numeric order (menu), and the modal length (the default).
.scroll_cdr3_lengths <- function(data) {
  t <- .scroll_cdr3_length_counts(data)
  if (!length(t)) character() else as.character(sort(as.integer(names(t))))
}
.scroll_cdr3_modal <- function(data) {
  t <- .scroll_cdr3_length_counts(data)
  if (!length(t)) NULL else names(t)[which.max(t)]
}

# Assembled as a built-in gated on manifest$vdj (placed after the Clone-map panel).
.scroll_cdr3_logo_panels <- function() {
  gate <- function(m) !is.null(m$vdj) && requireNamespace("ggseqlogo", quietly = TRUE)
  list(.scroll_gated_panel(
    "cdr3_logo", "CDR3 logo", "CDR3 sequence logo",
    "Amino-acid sequence logos of the CDR3 at a chosen length, one per group.",
    gate,
    controls = list(
      scroll_input_choice("chain", "Chain",
        choices = function(data) unlist(data$manifest$vdj$cdr3_chains)),
      scroll_input_choice("length", "CDR3 length",
        choices = function(data) .scroll_cdr3_lengths(data),   # numeric order
        selected = function(data) .scroll_cdr3_modal(data),    # default = modal length
        widget = "select"),
      # per-group logos; "None" pools all sequences into a single logo
      scroll_input_choice("group", "Group by",
        choices = function(data) c(.scroll_vdj_split_cols(data), "None" = ""), widget = "select"),
      scroll_input_levels("group_levels", "Groups", none = TRUE, watch = "group",
        choices = function(input, data) .scroll_vdj_col_levels(input, data, "group")),
      scroll_input_choice("count", "Count", c("Clones (dedup)", "Cells")),
      scroll_input_choice("units", "Logo units", c("Bits", "Probability")),
      scroll_input_slider("aspect", "Aspect ratio", 0.4, 3, 1, 0.1)),
    plot = function(cells, input, data) {
      if (!requireNamespace("ggseqlogo", quietly = TRUE))
        stop("Install the ggseqlogo package for CDR3 logos.")
      rc <- .scroll_vdj_read(data, "rep_cells.parquet")
      if (is.null(rc) || !nrow(rc)) stop("No repertoire store; rebuild with vdj =.")
      chain <- input$chain %||% unlist(data$manifest$vdj$cdr3_chains)[1]
      sc <- paste0("cdr3_", chain)
      if (!sc %in% names(rc))
        stop("CDR3 sequences not baked — rebuild the project with a current scroll_build().")
      rc <- .scroll_vdj_scope(rc, cells)          # honour the global filter / subset view
      grp <- if (!is.null(input$group) && nzchar(input$group) && input$group %in% names(rc) &&
                 any(!is.na(rc[[input$group]]))) input$group else NULL
      gcol <- grp %||% ".grp"
      if (is.null(grp)) rc$.grp <- "all"
      d <- rc[!is.na(rc[[sc]]) & nzchar(rc[[sc]]) & !is.na(rc[[gcol]]), , drop = FALSE]
      if (!is.null(grp) && length(input$group_levels))
        d <- d[as.character(d[[gcol]]) %in% input$group_levels, , drop = FALSE]
      if (identical(input$count, "Clones (dedup)")) {          # one CDR3 per clone per group
        d <- d[!is.na(d$clone_id), , drop = FALSE]
        d <- d[!duplicated(paste(d$clone_id, d[[gcol]])), , drop = FALSE]
      }
      if (!nrow(d)) stop("No CDR3 sequences for this selection.")
      d$seq <- d[[sc]]; d$len <- nchar(d$seq)
      # a logo needs equal-length sequences: use the chosen length, or the modal one
      len <- suppressWarnings(as.integer(input$length))
      if (is.na(len) || !any(d$len == len))
        len <- as.integer(names(sort(table(d$len), decreasing = TRUE))[1])
      d <- d[d$len == len, , drop = FALSE]
      if (!nrow(d)) stop("No length-", len, " CDR3s for this selection.")
      method <- if (identical(input$units, "Probability")) "prob" else "bits"
      # suppressWarnings: ggseqlogo builds its layer with the deprecated aes_string()
      p <- suppressWarnings(
        if (is.null(grp)) ggseqlogo::ggseqlogo(d$seq, method = method)
        else {
          gl <- split(d$seq, as.character(d[[gcol]])); gl <- gl[order(names(gl))]
          ggseqlogo::ggseqlogo(gl, method = method, ncol = min(length(gl), 3L))
        })
      unit <- if (identical(input$count, "Clones (dedup)")) "clones" else "cells"
      p <- p + ggplot2::labs(title = sprintf("CDR3 %s logo - length %d (%d %s)",
                                             chain, len, nrow(d), unit))
      attr(p, "scroll_source") <- data.frame(group = as.character(d[[gcol]]),
                                             cdr3 = d$seq, stringsAsFactors = FALSE)
      p
    }, csv = TRUE))
}
