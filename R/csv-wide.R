# --- Wide CSVs (for GraphPad Prism) -----------------------------------------------

# The categorical panels offer their CSV in two shapes: long (tidy, one row per cell
# or per group x category -- the original download) and wide, laid out the way Prism
# tables are entered:
#   Violin / Ridge / Signature   one column per group, each cell's value underneath
#                                (a Prism Column table); several genes -> "gene: group"
#   Proportions                  one row per group, one column per category
#   DotPlot (dots or tiles)      genes x groups: a block of mean expression, then a
#                                block of % expressing
# Columns and rows follow the plot's order (natural, or the manual group order).

# Per-cell values -> one column per group, padded with blanks to equal length.
# `series` is a named list of numeric vectors aligned to `group` (one per gene / value);
# with more than one, columns are named "<series>: <group>".
.scroll_wide_columns <- function(series, group, order = NULL) {
  group <- as.character(group)
  lv <- .scroll_level_order(group, order)
  cols <- list()
  for (s in names(series)) {
    v <- series[[s]]
    for (l in lv) {
      x <- v[!is.na(group) & group == l & !is.na(v)]
      cols[[if (length(series) > 1L) paste0(s, ": ", l) else l]] <- x
    }
  }
  n <- max(c(0L, lengths(cols)))
  out <- as.data.frame(lapply(cols, function(x) c(x, rep(NA_real_, n - length(x)))),
                       check.names = FALSE, stringsAsFactors = FALSE)
  names(out) <- names(cols)
  out
}

# The plotted values of a Violin / Ridge data reactive (`d`), one per gene (or the
# numeric column), aligned to `d$cells`.
.scroll_wide_series <- function(d) {
  if (!is.null(d$value_col))
    return(stats::setNames(list(suppressWarnings(as.numeric(d$cells[[d$value_col]]))), d$value_col))
  feats <- d$features %||% d$feature
  vl <- d$values
  stats::setNames(lapply(feats, function(g) {
    v <- if (!is.null(vl) && "feature" %in% names(vl))
           vl[vl$feature == g, c("cell", "value"), drop = FALSE] else vl
    .scroll_expr_vector(d$cells, v)
  }), feats)
}

# The group label per cell for a violin-style plot: the (possibly composite) group,
# or "group: split" when split.
.scroll_wide_groups <- function(cells, group_by, split_by = NULL) {
  g <- .scroll_combo_levels(cells, group_by)
  split_by <- split_by[nzchar(split_by %||% "")]
  if (length(split_by)) {
    s <- .scroll_combo_levels(cells, split_by)
    g <- ifelse(is.na(g) | is.na(s), NA_character_, paste0(g, ": ", s))
  }
  g
}

# A long table -> wide: one row per `row_col` value (in `row_order`), one column per
# `col_col` value (in `col_order`), cells from `value_col`.
.scroll_wide_pivot <- function(df, row_col, col_col, value_col, row_order, col_order,
                               row_label = row_col) {
  M <- matrix(NA_real_, length(row_order), length(col_order))
  i <- match(as.character(df[[row_col]]), row_order)
  j <- match(as.character(df[[col_col]]), col_order)
  ok <- !is.na(i) & !is.na(j)
  M[cbind(i[ok], j[ok])] <- df[[value_col]][ok]
  out <- data.frame(row_order, M, check.names = FALSE, stringsAsFactors = FALSE)
  names(out) <- c(row_label, col_order)
  out
}

# Proportions: one row per x group, one column per fill category; percent of the
# group when the bars are 100% stacked, else cell counts.
.scroll_proportions_wide <- function(cells, group_by, fill_by, position = "fill",
                                     group_order = NULL) {
  long <- .scroll_proportions_source(cells, group_by, fill_by)
  long$percent <- 100 * long$proportion
  .scroll_wide_pivot(long, "group", "category",
                     if (identical(position, "fill")) "percent" else "n_cells",
                     .scroll_level_order(long$group, group_order),
                     .scroll_natural_sort(long$category),
                     row_label = paste(group_by, collapse = " | "))
}

# DotPlot / tiles: genes x groups, as two blocks (mean expression, % expressing).
.scroll_dotplot_wide <- function(cells, features, group_by, expr_long,
                                 feature_order = features, group_order = NULL) {
  long <- .scroll_heatmap_source(cells, features, group_by, expr_long)
  long$pct_expressing <- 100 * long$pct_expressing
  g <- .scroll_level_order(long$group, group_order)
  list(
    "Mean expression (log-normalised)" =
      .scroll_wide_pivot(long, "feature", "group", "avg_expr", feature_order, g, "gene"),
    "% expressing" =
      .scroll_wide_pivot(long, "feature", "group", "pct_expressing", feature_order, g, "gene"))
}

# Write a wide table, or several titled blocks one after another, as CSV (blank cells
# for missing values, which Prism reads as empty).
.scroll_write_wide <- function(x, file) {
  if (is.data.frame(x)) return(utils::write.csv(x, file, row.names = FALSE, na = ""))
  con <- file(file, open = "w", encoding = "UTF-8")
  on.exit(close(con))
  blocks <- x
  for (k in seq_along(blocks)) {
    if (k > 1L) writeLines("", con)
    writeLines(sprintf("\"%s\"", names(blocks)[k]), con)
    utils::write.table(blocks[[k]], con, sep = ",", row.names = FALSE, na = "",
                       qmethod = "double")
  }
}

.scroll_wide_handler <- function(wide_r, name)
  downloadHandler(filename = function() name,
                  content = function(file) .scroll_write_wide(wide_r(), file))

# The toolbar's CSV control when a panel offers both shapes: a small menu.
.scroll_csv_menu <- function(ns)
  div(class = "dropdown scroll-dl-menu",
      tags$button(class = "btn btn-sm scroll-dl dropdown-toggle", type = "button",
                  `data-bs-toggle` = "dropdown", `aria-expanded` = "false",
                  shiny::icon("download"), "CSV"),
      tags$ul(class = "dropdown-menu dropdown-menu-end",
              tags$li(downloadLink(ns("csv"), "Long (tidy)", class = "dropdown-item")),
              tags$li(downloadLink(ns("csv_wide"), "Wide (for Prism)", class = "dropdown-item"))))
