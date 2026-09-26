# The VDJ (TCR/BCR) repertoire panels + their runtime accessors. The build-time
# spec + bake logic lives in R/vdj.R; these read the baked repertoire/ store and
# render the four gated panels (clone overview, gene usage, CDR3 length, diversity).

# ---- runtime: read the baked store ------------------------------------------

.scroll_vdj_read <- function(data, f)
  .scroll_read_asset(file.path(data$dir, "repertoire", f), "parquet")
# Narrow the baked per-cell repertoire table to the app's ACTIVE cells (global filter +
# subset view) by barcode, when the store carries a `cell` key. So the runtime-computed
# repertoire views (rank-abundance, expansion, gene frequency, CDR3 length, diversity
# metric) follow the filter/view. Stores baked before the `cell` key (or the baked
# chi-square / tissue-correlation tables) are whole-dataset and left untouched.
.scroll_vdj_scope <- function(rc, cells) {
  if (is.null(rc) || is.null(cells) || !"cell" %in% names(rc) || !"cell" %in% names(cells))
    return(rc)
  rc <- rc[rc$cell %in% cells$cell, , drop = FALSE]
  # Join any categorical cells.parquet column not already baked into the store onto the
  # per-cell repertoire table by barcode, so a repertoire view can group / colour / split
  # by ANY metadata column (e.g. cloneSize, a cluster resolution), not only the columns
  # carried at build time. Non-members of the active view have no barcode match -> NA.
  extra <- setdiff(names(cells)[vapply(cells, function(v) is.character(v) || is.factor(v),
                                       logical(1))], names(rc))
  if (length(extra) && nrow(rc)) {
    idx <- match(rc$cell, cells$cell)
    for (col in extra) rc[[col]] <- as.character(cells[[col]])[idx]
  }
  rc
}
# Distinct non-NA values of a vector.
.scroll_ndistinct <- function(v) length(unique(v[!is.na(v)]))

# Categorical rep_cells columns that are candidates for grouping or as a clone id
# (excludes the segment genes, all CDR3 columns (sequences + lengths), and the count).
.scroll_vdj_cat_cols <- function(rc, segments) {
  cand <- setdiff(names(rc), c("cell", "clone_id", "clone_count",
                               grep("^cdr3_", names(rc), value = TRUE), unlist(segments)))
  cand[vapply(cand, function(c) is.character(rc[[c]]) && any(!is.na(rc[[c]])), logical(1))]
}

# Split the categorical columns by cardinality, relative to the baked clone_id: a
# column near clone-scale cardinality is treated as a clone id (an alternate clonal
# definition), everything lower-cardinality is a grouping/colour axis. This keeps a
# 50k-level clone column out of the group-by menus and out of the colour scales.
.scroll_vdj_high_card <- function(rc, col, ref) .scroll_ndistinct(rc[[col]]) >= 0.5 * ref

# categorical rep_cells columns available to split/colour by (low-cardinality only —
# a high-cardinality clone column would mean tens of thousands of facets/colours).
.scroll_vdj_split_cols <- function(data) {
  rc <- .scroll_vdj_read(data, "rep_cells.parquet"); if (is.null(rc)) return(character())
  ref <- if ("clone_id" %in% names(rc)) .scroll_ndistinct(rc$clone_id) else nrow(rc)
  baked <- Filter(function(c) !.scroll_vdj_high_card(rc, c, ref),
                  .scroll_vdj_cat_cols(rc, data$manifest$vdj$segments))
  # Every other categorical metadata column is offered too (joined per-cell by barcode at
  # render, see .scroll_vdj_scope), capped to a sane level count so the group menus stay
  # usable and identifier-like columns (barcode, clone id) never appear.
  meta  <- data$manifest$meta
  extra <- Filter(function(c) { n <- length(unlist(meta[[c]]$levels)); n >= 1 && n <= 50 },
                  names(Filter(function(x) identical(x$type, "categorical"), meta)))
  unique(c(baked, setdiff(extra, c(baked, data$manifest$vdj$group_col))))
}

# Columns that can serve as a clone id: the baked `clone_id` plus any carried
# categorical column of comparable (high) cardinality — an alternate clonal
# definition, e.g. a nucleotide-vs-amino-acid CDR3 clonotype. Low-cardinality
# grouping columns (donor, antigen, cell type) are excluded.
.scroll_vdj_clone_cols <- function(data) {
  rc <- .scroll_vdj_read(data, "rep_cells.parquet")
  if (is.null(rc) || !"clone_id" %in% names(rc)) return("clone_id")
  ref <- .scroll_ndistinct(rc$clone_id)
  alt <- Filter(function(c) .scroll_vdj_high_card(rc, c, ref),
                .scroll_vdj_cat_cols(rc, data$manifest$vdj$segments))
  unique(c("clone_id", alt))
}

# The distinct levels of a group column named by the control `group_id` (for a group
# filter's choices). Falls back to the baked `group` column.
.scroll_vdj_col_levels <- function(input, data, group_id = "group") {
  g <- if (!is.null(input[[group_id]]) && nzchar(input[[group_id]])) input[[group_id]] else "group"
  rc <- .scroll_vdj_read(data, "rep_cells.parquet")
  if (!is.null(rc) && g %in% names(rc))                       # a baked repertoire column
    return(sort(unique(as.character(rc[[g]][!is.na(rc[[g]])]))))
  lv <- unlist(data$manifest$meta[[g]]$levels)                # a cells.parquet column, joined at render
  if (length(lv)) sort(as.character(lv)) else character()
}

# The active group levels = the group column's levels narrowed by the group filter
# (`levels_id`). Shared by the colour control (to render one picker per level) and the
# plot (to build the colour scale), so the two stay in lock-step. With no group selected
# ("None"), the plot is a single pooled series -- return one synthetic "all" level so the
# same palette / Manual colour control drives that one colour.
.scroll_vdj_level_set <- function(input, data, group_id = "group", levels_id = "group_levels") {
  gv <- input[[group_id]]
  if (is.null(gv) || !nzchar(gv)) return("all")
  lv <- .scroll_vdj_col_levels(input, data, group_id)
  f <- input[[levels_id]]
  if (length(f)) lv <- lv[lv %in% f]
  lv
}

# A single colour for an ungrouped (pooled) series: the first colour the palette / Manual
# control resolves for the active level set (which is the synthetic "all" when ungrouped).
.scroll_vdj_one_colour <- function(input, data, group_id, levels_id = "group_levels")
  unname(.scroll_vdj_colours(input, .scroll_vdj_level_set(input, data, group_id, levels_id)))[1]

# Zero lower-end expansion for continuous (value) position scales, so bars / curves sit
# flush on the axis; a small upper margin is kept for headroom.
.scroll_expand0 <- function() ggplot2::expansion(mult = c(0, 0.05))

# A palette + manual-per-group colour control for the VDJ panels, keyed to the active
# group levels. The Manual pickers are rendered dynamically (needs the builder to pass
# `output` to the bind), so they track the chosen group column + group filter.
.scroll_vdj_colour_control <- function(group_id = "group", levels_id = "group_levels") {
  scroll_input_custom("palette", "Colour",
    ui = function(ns, data)
      tagList(selectInput(ns("palette"), "Colour", .scroll_cat_palettes()),
              uiOutput(ns("palette_manual"))),
    bind = function(input, session, data, output) {
      ns <- session$ns
      output$palette_manual <- renderUI({
        if (!identical(input$palette, "Manual")) return(NULL)
        .scroll_manual_ui(ns, .scroll_vdj_level_set(input, data, group_id, levels_id))
      })
    })
}

# Named colour vector for the active group levels: the chosen palette, or the manual
# per-group pickers when palette == "Manual" (falling back until they populate).
.scroll_vdj_colours <- function(input, levels) {
  pal <- input$palette %||% "Tableau 10"
  if (identical(pal, "Manual")) {
    m <- .scroll_manual_colors(input, levels)
    if (!is.null(m) && all(levels %in% names(m))) return(m)
    pal <- "Tableau 10"
  }
  .scroll_discrete_colors(levels, pal)
}

# Expansion composition colours: the five ORDERED clone-size categories, a few
# sequential palettes suited to an ordered scale, plus Manual. Kept separate from the
# group colour control (distinct input ids "epalette"/"ecol_*") so both can live in the
# clone_overview panel without colliding; ordered, so it can't reuse the qualitative,
# alphabetically-sorted group colour machinery.
.SCROLL_EXPANSION_LEVELS <- c("Single", "Small", "Medium", "Large", "Hyperexpanded")
.scroll_expansion_palettes <- function() list(
  "YlOrRd"  = c("#FFFFB2", "#FECC5C", "#FD8D3C", "#F03B20", "#BD0026"),
  "Blues"   = c("#EFF3FF", "#BDD7E7", "#6BAED6", "#3182BD", "#08519C"),
  "Reds"    = c("#FEE5D9", "#FCAE91", "#FB6A4A", "#DE2D26", "#A50F15"),
  "Greens"  = c("#EDF8E9", "#BAE4B3", "#74C476", "#31A354", "#006D2C"),
  "Viridis" = c("#FDE725", "#7AD151", "#22A884", "#2A788E", "#414487"))

.scroll_expansion_colours <- function(input) {
  lv <- .SCROLL_EXPANSION_LEVELS; pals <- .scroll_expansion_palettes()
  pal <- input$epalette %||% "YlOrRd"
  if (identical(pal, "Manual")) {
    m <- .scroll_manual_colors(input, lv, prefix = "ecol")
    if (!is.null(m) && all(lv %in% names(m))) return(m)
    pal <- "YlOrRd"
  }
  stats::setNames(pals[[pal]] %||% pals[["YlOrRd"]], lv)
}

.scroll_expansion_colour_control <- function()
  scroll_input_custom("epalette", "Colour",
    ui = function(ns, data)
      tagList(selectInput(ns("epalette"), "Colour",
                          c(names(.scroll_expansion_palettes()), "Manual"), selected = "YlOrRd"),
              uiOutput(ns("epalette_manual"))),
    bind = function(input, session, data, output) {
      output$epalette_manual <- renderUI({
        if (!identical(input$epalette, "Manual")) return(NULL)
        .scroll_manual_ui(session$ns, .SCROLL_EXPANSION_LEVELS, prefix = "ecol",
          defaults = stats::setNames(.scroll_expansion_palettes()[["YlOrRd"]],
                                     .SCROLL_EXPANSION_LEVELS))
      })
    })

# The max-Small / max-Medium / max-Large cut points as QUANTILES of the expanded
# (size >= 2) clones. Single (size 1) stays a fixed category. Three plain sliders: a
# multi-handle noUiSlider throws "$x.noUiSlider is not a function" while hidden in a
# conditionalPanel (this control is gated to the Expansion-composition view), a client
# error that can flake the whole Shiny session.
.scroll_expansion_quantile_control <- function()
  scroll_input_custom("exp_q", "Category cut quantiles",
    ui = function(ns, data) tagList(
      sliderInput(ns("exp_q1"), "Small max (quantile)",  0, 1, 0.5,  step = 0.01),
      sliderInput(ns("exp_q2"), "Medium max (quantile)", 0, 1, 0.8,  step = 0.01),
      sliderInput(ns("exp_q3"), "Large max (quantile)",  0, 1, 0.95, step = 0.01)))

# The three cut quantiles from the control (multi-handle slider, or the three-slider
# fallback), sorted ascending.
.scroll_expansion_quantiles <- function(input) {
  q <- if (length(input$exp_q) == 3) as.numeric(input$exp_q)
       else c(input$exp_q1 %||% 0.5, input$exp_q2 %||% 0.8, input$exp_q3 %||% 0.95)
  sort(q)
}

# Assign each clone an expansion category by the quantile of its size among the EXPANDED
# (size >= 2) clones. Single (size 1) is fixed; clones sharing a size share a category
# (ties by lower edge). `q` = c(max_small, max_medium, max_large), strictly increasing.
.scroll_expansion_cut <- function(count, q) {
  cat <- rep("Single", length(count))
  ex <- which(count > 1)
  if (length(ex)) {
    pos <- (rank(count[ex], ties.method = "min") - 1) / length(ex)   # frac strictly smaller
    cat[ex] <- as.character(cut(pos, c(-Inf, q, Inf),
      labels = c("Small", "Medium", "Large", "Hyperexpanded")))
  }
  factor(cat, levels = .SCROLL_EXPANSION_LEVELS)
}

# Natural (numeric-aware) ordering of gene-segment names, approximating genomic /
# germline orientation: IMGT segment names are numbered along the locus, so TRBV2
# should precede TRBV10-1 -- a plain alphabetical sort puts TRBV10-1 first. Split each
# name into number / non-number chunks (in order), zero-pad the numbers, and sort on
# the reconstructed key. Non-numeric names fall back to their alphabetical position.
.scroll_gene_natural_levels <- function(genes) {
  u <- unique(as.character(genes))
  key <- vapply(u, function(s) {
    parts <- regmatches(s, gregexpr("[0-9]+|[^0-9]+", s))[[1]]
    num <- grepl("^[0-9]+$", parts)
    parts[num] <- formatC(as.integer(parts[num]), width = 8, flag = "0")
    paste0(parts, collapse = "")
  }, character(1))
  u[order(key)]
}

# ---- the four panels --------------------------------------------------------

# Assembled as built-in panels gated on manifest$vdj (see .scroll_assemble_panels).
.scroll_vdj_panels <- function() {
  gate <- function(m) !is.null(m$vdj)
  mk <- function(id, label, title, desc, controls, plot, csv = TRUE)
    .scroll_gated_panel(id, label, title, desc, gate, controls, plot, csv = csv)

  clone_overview <- mk("clone_overview", "Clone overview", "Clonal expansion overview",
    "Rank-abundance of clone sizes and expansion-category composition (clone-level).",
    list(scroll_input_choice("view", "View", c("Rank-abundance", "Expansion composition")),
         # which column defines a clone (the baked clone_id, or an alternate clonal
         # definition carried in the store, e.g. a nucleotide vs amino-acid clonotype)
         scroll_input_choice("clone_col", "Clone ID",
           choices = function(data) .scroll_vdj_clone_cols(data), widget = "select"),
         # "None" (default) pools all clones into a single ungrouped plot.
         scroll_input_choice("group", "Group by",
           choices = function(data) c("None" = "", .scroll_vdj_split_cols(data)),
           widget = "select"),
         # restrict to a subset of the Group-by column's levels (empty = all)
         scroll_input_levels("group_levels", "Groups", none = TRUE, watch = "group",
           choices = function(input, data) .scroll_vdj_col_levels(input, data, "group")),
         # rank-abundance clone-size axis: log10 (default) or raw counts
         scroll_show_when(scroll_input_choice("yscale", "Y axis", c("Log", "Raw")),
                          control = "view", equals = "Rank-abundance"),
         # palette + manual per-group colour picker (rank-abundance colours by group)
         scroll_style_input(scroll_show_when(.scroll_vdj_colour_control("group", "group_levels"),
                                             control = "view", equals = "Rank-abundance")),
         # expansion category cut points, as quantiles of the expanded clones (one slider,
         # three knobs). Single = 1 cell; Hyperexpanded = above the Large quantile.
         scroll_show_when(.scroll_expansion_quantile_control(),
                          control = "view", equals = "Expansion composition"),
         # ordered palette + manual pickers for the five expansion categories
         scroll_style_input(scroll_show_when(.scroll_expansion_colour_control(),
                                             control = "view", equals = "Expansion composition")),
         scroll_style_input(scroll_input_slider("aspect", "Aspect ratio", 0.4, 3, 1, 0.1))),
    function(cells, input, data) {
      rc <- .scroll_vdj_read(data, "rep_cells.parquet")
      if (is.null(rc) || !nrow(rc)) stop("No repertoire store; rebuild with vdj =.")
      rc <- .scroll_vdj_scope(rc, cells)          # honour the global filter / subset view
      # Clone id column drives what counts as a clone; Group-by drives the split/facet,
      # or NULL for a single ungrouped plot (the "None" default).
      cid <- if (!is.null(input$clone_col) && input$clone_col %in% names(rc)) input$clone_col else "clone_id"
      grp <- if (!is.null(input$group) && nzchar(input$group) && input$group %in% names(rc) &&
                 any(!is.na(rc[[input$group]]))) input$group else NULL
      gcol <- grp %||% ".grp"
      if (is.null(grp)) rc$.grp <- "all"
      rc <- rc[!is.na(rc[[cid]]) & !is.na(rc[[gcol]]), , drop = FALSE]
      if (!is.null(grp) && length(input$group_levels))         # optional group filter
        rc <- rc[as.character(rc[[gcol]]) %in% input$group_levels, , drop = FALSE]
      if (!nrow(rc)) stop("No clones for the selected groups.")
      if (identical(input$view, "Rank-abundance")) {
        d <- rc |> dplyr::group_by(.data[[cid]], .data[[gcol]]) |>
          dplyr::summarise(count = dplyr::n(), .groups = "drop") |>
          dplyr::arrange(.data[[gcol]], dplyr::desc(.data$count)) |>
          dplyr::group_by(.data[[gcol]]) |>
          dplyr::mutate(rank = dplyr::row_number()) |> dplyr::ungroup()
        raw <- identical(input$yscale, "Raw")
        ylab <- if (raw) "Clone size (cells)" else "Clone size (cells, log10)"
        # keep default axis expansion here: the log-scale rank-abundance scatter needs
        # the padding (zero expansion crowds points against the axes).
        yscale <- if (raw) ggplot2::scale_y_continuous() else ggplot2::scale_y_log10()
        if (is.null(grp)) {
          p <- ggplot2::ggplot(d, ggplot2::aes(.data$rank, .data$count)) +
            ggplot2::geom_point(size = 0.5,
              colour = .scroll_vdj_one_colour(input, data, "group")) + yscale +
            ggplot2::labs(x = "Clone rank", y = ylab, title = "Clone rank-abundance") +
            .scroll_base_theme(legend = FALSE)
        } else {
          lv <- .scroll_vdj_level_set(input, data, "group", "group_levels")
          p <- ggplot2::ggplot(d, ggplot2::aes(.data$rank, .data$count, colour = .data[[gcol]])) +
            ggplot2::geom_point(size = 0.5) + yscale +
            ggplot2::scale_colour_manual(values = .scroll_vdj_colours(input, lv), name = grp) +
            ggplot2::facet_wrap(stats::as.formula(paste0("~`", gcol, "`")), scales = "free_x") +
            ggplot2::labs(x = "Clone rank", y = ylab, colour = grp, title = "Clone rank-abundance") +
            .scroll_base_theme(legend = FALSE)
        }
        attr(p, "scroll_source") <- as.data.frame(d); p
      } else {
        # Clone size + each clone's dominant group. Vectorized: one grouped (clone x
        # group) count, then take the top-count group per clone via order + de-dup --
        # replaces a table() call per clone, which cost seconds over tens of thousands
        # of clones (52k here). Equivalent result, ~8x faster.
        cg  <- dplyr::count(rc, .data[[cid]], .data[[gcol]], name = "n")
        cg  <- cg[order(cg[[cid]], -cg$n), , drop = FALSE]
        dom <- cg[!duplicated(cg[[cid]]), c(cid, gcol), drop = FALSE]   # dominant group
        cl  <- dplyr::left_join(dplyr::count(rc, .data[[cid]], name = "count"),
                                dom, by = cid)
        names(cl)[names(cl) == gcol] <- "g"
        # cut into categories by quantile of clone size (Single = 1 cell is fixed)
        q <- .scroll_expansion_quantiles(input)
        if (!all(is.finite(q)) || !(0 < q[1] && q[1] < q[2] && q[2] < q[3] && q[3] < 1))
          stop("Cut quantiles must satisfy 0 < Small < Medium < Large < 1.")
        cl$expansion <- .scroll_expansion_cut(cl$count, q)
        p <- ggplot2::ggplot(cl, ggplot2::aes(.data$g, fill = .data$expansion)) +
          ggplot2::geom_bar(position = "fill", colour = "white", linewidth = 0.2) +
          ggplot2::scale_y_continuous(labels = scales::percent,
                                      expand = ggplot2::expansion(mult = c(0, 0))) +
          ggplot2::scale_fill_manual(values = .scroll_expansion_colours(input),
                                     name = "Expansion", drop = FALSE) +
          ggplot2::labs(x = grp, y = "Fraction of clones", fill = "Expansion",
                        title = "Expansion-category composition") +
          .scroll_base_theme()
        attr(p, "scroll_source") <- as.data.frame(cl); p
      }
    })

  gene_usage <- mk("gene_usage", "V/J gene usage", "TCR/BCR V/J gene usage",
    "Per-group V/J frequency, segment-pairing heatmap, or chi-square residuals.",
    list(scroll_input_choice("segment", "Segment",
           choices = function(data) unlist(data$manifest$vdj$segments)),
         scroll_input_choice("view", "View",
           c("Frequency", "Pairing", "Chi-square residuals")),
         # Pairing view: the second (y-axis) segment. Default to a DIFFERENT segment
         # than x (the choices fn lists the 2nd segment first, so it's pre-selected).
         scroll_show_when(
           scroll_input_choice("segment_y", "Segment (y)",
             choices = function(data) {
               s <- as.character(unlist(data$manifest$vdj$segments))
               if (length(s) > 1) c(s[2], s[-2]) else s
             }, widget = "select"),
           control = "view", equals = "Pairing"),
         # Count each cell, or deduplicate to one count per clone so a few large
         # (expanded) clones don't dominate the profile. Dedup needs a clone
         # definition, so the Clone-ID picker rides with it (Frequency + Pairing).
         scroll_show_when(
           scroll_input_choice("count", "Count", c("Cells", "Clones (dedup)")),
           control = "view", equals = c("Frequency", "Pairing")),
         scroll_show_when(
           scroll_input_choice("clone_col", "Clone ID",
             choices = function(data) .scroll_vdj_clone_cols(data), widget = "select"),
           control = "view", equals = c("Frequency", "Pairing")),
         # gene-axis order: Genomic (natural sort of IMGT gene numbers, approximating
         # germline/locus orientation) or plain Alphabetical.
         scroll_show_when(
           scroll_input_choice("gene_order", "Gene order", c("Genomic", "Alphabetical")),
           control = "view", equals = c("Frequency", "Pairing")),
         # Pairing heatmap fill: a continuous palette + min/max quantile cut-offs on the
         # colour scale (a two-knob range slider), so a few dominant pairings don't wash
         # out the rest -- values past the cut-offs saturate to the end colour.
         scroll_style_input(scroll_show_when(
           scroll_input_palette("cpalette", "Colour", type = "continuous",
                                selected = "viridis"),
           control = "view", equals = "Pairing")),
         scroll_style_input(scroll_show_when(
           scroll_input_slider("cquant", "Colour quantiles", 0, 1, c(0, 1), 0.01),
           control = "view", equals = "Pairing")),
         # Group-by applies to every view: re-groups Frequency and the (runtime-recomputed)
         # chi-square residuals, and facets the pairing heatmap. "None" (the default) pools
         # all cells into a single ungrouped plot (chi-square still requires a group).
         scroll_input_choice("group", "Group by",
           choices = function(data) c("None" = "", .scroll_vdj_split_cols(data)),
           widget = "select"),
         scroll_input_levels("group_levels", "Groups", none = TRUE, watch = "group",
           choices = function(input, data) .scroll_vdj_col_levels(input, data, "group")),
         # per-group point colours (Frequency only; the heatmap views use a fill scale)
         scroll_style_input(scroll_show_when(.scroll_vdj_colour_control("group", "group_levels"),
                                             control = "view", equals = "Frequency")),
         scroll_style_input(scroll_input_slider("aspect", "Aspect ratio", 0.4, 3, 1, 0.1))),
    function(cells, input, data) {
      rc <- .scroll_vdj_read(data, "rep_cells.parquet")
      if (is.null(rc) || !nrow(rc)) stop("No repertoire store; rebuild with vdj =.")
      rc <- .scroll_vdj_scope(rc, cells)          # honour the global filter / subset view
      seg <- input$segment
      if (!seg %in% names(rc)) stop("Segment '", seg, "' not baked.")
      view <- input$view %||% "Frequency"
      # group-by (all views): the chosen categorical, or NULL for a single ungrouped plot
      # (the "None" default). A synthetic constant column keeps the filter/dedup path uniform.
      grp <- if (!is.null(input$group) && nzchar(input$group) && input$group %in% names(rc) &&
                 any(!is.na(rc[[input$group]]))) input$group else NULL
      gcol <- grp %||% ".grp"
      if (is.null(grp)) rc$.grp <- "all"
      cid <- if (!is.null(input$clone_col) && input$clone_col %in% names(rc)) input$clone_col else "clone_id"

      if (identical(view, "Pairing")) {
        segs <- as.character(unlist(data$manifest$vdj$segments))
        segy <- if (!is.null(input$segment_y) && input$segment_y %in% names(rc))
                  input$segment_y else setdiff(segs, seg)[1]
        if (is.na(segy) || identical(segy, seg))
          stop("Pick two different segments for the pairing heatmap.")
        d <- rc[!is.na(rc[[seg]]) & !is.na(rc[[segy]]) & !is.na(rc[[gcol]]), , drop = FALSE]
        if (!is.null(grp) && length(input$group_levels))
          d <- d[as.character(d[[gcol]]) %in% input$group_levels, , drop = FALSE]
        if (identical(input$count, "Clones (dedup)")) {
          d <- d[!is.na(d[[cid]]), , drop = FALSE]
          d <- d[!duplicated(paste(d[[cid]], d[[gcol]])), , drop = FALSE]
        }
        if (!nrow(d)) stop("No paired segments for this selection.")
        # frequency of each (x, y) segment pairing WITHIN each group (facet)
        pf <- d |> dplyr::count(.data[[gcol]], .data[[seg]], .data[[segy]], name = "n") |>
          dplyr::group_by(.data[[gcol]]) |>
          dplyr::mutate(freq = .data$n / sum(.data$n)) |> dplyr::ungroup()
        names(pf)[names(pf) == seg]  <- "gx"
        names(pf)[names(pf) == segy] <- "gy"
        names(pf)[names(pf) == gcol] <- "grp"
        # complete the (x, y) grid within each group so absent pairings render as 0
        # (the low colour of the scale), not blank white tiles
        full <- expand.grid(gx = sort(unique(pf$gx)), gy = sort(unique(pf$gy)),
                            grp = unique(pf$grp), stringsAsFactors = FALSE)
        pf <- dplyr::left_join(full, pf, by = c("gx", "gy", "grp"))
        pf$n[is.na(pf$n)] <- 0; pf$freq[is.na(pf$freq)] <- 0
        if (!identical(input$gene_order, "Alphabetical")) {  # default: genomic (natural)
          pf$gx <- factor(pf$gx, levels = .scroll_gene_natural_levels(pf$gx))
          pf$gy <- factor(pf$gy, levels = .scroll_gene_natural_levels(pf$gy))
        }
        unit <- if (identical(input$count, "Clones (dedup)")) "clones" else "cells"
        # colour-scale quantile cut-offs (min/max range slider); c(0,1) = no clipping
        qr <- if (length(input$cquant) == 2) as.numeric(input$cquant) else c(0, 1)
        lims <- .scroll_expr_limits(pf$freq, sort(qr))
        p <- ggplot2::ggplot(pf, ggplot2::aes(.data$gx, .data$gy, fill = .data$freq)) +
          ggplot2::geom_tile() +
          .scroll_continuous_scale(input$cpalette %||% "viridis",
            name = paste0("Freq\n(", unit, ")"), limits = lims,
            aesthetic = "fill", labels = scales::percent) +
          (if (!is.null(grp)) ggplot2::facet_wrap(stats::as.formula("~grp"))) +
          ggplot2::labs(x = seg, y = segy, title = paste0(seg, "-", segy, " pairing")) +
          .scroll_base_theme() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, size = 7),
                         axis.text.y = ggplot2::element_text(size = 7))
        attr(p, "scroll_source") <- as.data.frame(pf); return(p)
      }

      if (identical(view, "Chi-square residuals")) {
        # residuals are of a group x gene table -- a grouping is intrinsic to this view
        if (is.null(grp))
          stop("Pick a 'Group by' column for the chi-square residuals view.")
        d <- rc[!is.na(rc[[grp]]), , drop = FALSE]
        if (length(input$group_levels))
          d <- d[as.character(d[[grp]]) %in% input$group_levels, , drop = FALSE]
        d <- d[!is.na(d[[cid]]), , drop = FALSE]
        d <- d[!duplicated(d[[cid]]), , drop = FALSE]     # one row per clone (usage, not size)
        # recomputed at runtime (not the baked chisq.parquet) so it re-groups by any
        # column AND honours the global cell filter / subset view.
        ch <- .scroll_vdj_chisq(d, grp, seg)
        if (is.null(ch) || !nrow(ch))
          stop("Chi-square needs >= 2 groups and >= 2 genes with data.")
        ch$residual <- pmax(pmin(ch$residual, 3), -3)
        ord <- names(sort(tapply(abs(ch$residual), ch$gene, max)))
        ch$gene <- factor(ch$gene, levels = ord)
        p <- ggplot2::ggplot(ch, ggplot2::aes(.data$gene, .data$group, fill = .data$residual)) +
          ggplot2::geom_tile(colour = "grey90") +
          ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                        midpoint = 0, limits = c(-3, 3)) +
          ggplot2::labs(x = seg, y = grp, fill = "std\nresidual",
                        title = paste(seg, "usage vs", grp)) +
          .scroll_base_theme() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, size = 8))
        attr(p, "scroll_source") <- as.data.frame(ch); return(p)
      }

      # Frequency (default)
      d <- rc[!is.na(rc[[seg]]) & !is.na(rc[[gcol]]), , drop = FALSE]
      d$gene <- as.character(d[[seg]])
      if (!is.null(grp) && length(input$group_levels))
        d <- d[as.character(d[[gcol]]) %in% input$group_levels, , drop = FALSE]
      if (!nrow(d)) stop("No cells for the selected groups.")
      # optional: collapse each clone to one row per group (one representative segment
      # call), so expanded clones count once toward the usage frequency
      if (identical(input$count, "Clones (dedup)")) {
        d <- d[!is.na(d[[cid]]), , drop = FALSE]
        d <- d[!duplicated(paste(d[[cid]], d[[gcol]])), , drop = FALSE]
        if (!nrow(d)) stop("No clones for the selected groups.")
      }
      f <- d |> dplyr::group_by(.data$gene, .data[[gcol]]) |>
        dplyr::summarise(count = dplyr::n(), .groups = "drop") |>
        dplyr::group_by(.data[[gcol]]) |>
        dplyr::mutate(freq = .data$count / sum(.data$count)) |> dplyr::ungroup()
      unit <- if (identical(input$count, "Clones (dedup)")) "clones" else "cells"
      if (!identical(input$gene_order, "Alphabetical"))      # default: genomic (natural)
        f$gene <- factor(f$gene, levels = .scroll_gene_natural_levels(f$gene))
      if (is.null(grp)) {
        # ungrouped: a single overall usage barplot (colour from the palette / Manual control)
        p <- ggplot2::ggplot(f, ggplot2::aes(.data$gene, .data$freq)) +
          ggplot2::geom_col(fill = .scroll_vdj_one_colour(input, data, "group")) +
          ggplot2::scale_y_continuous(expand = .scroll_expand0()) +
          ggplot2::labs(x = seg, y = paste0("Frequency (", unit, ", overall)"),
                        title = paste(seg, "usage")) +
          .scroll_base_theme(legend = FALSE) +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, size = 7))
      } else {
        lv <- .scroll_vdj_level_set(input, data, "group", "group_levels")
        p <- ggplot2::ggplot(f, ggplot2::aes(.data$gene, .data$freq)) +
          # a line joins one gene's groups: only for genes seen in 2+ groups (a lone point
          # has nothing to join, and ggplot would message on every render)
          ggplot2::geom_line(ggplot2::aes(group = .data$gene), colour = "grey70", linewidth = 0.3,
                             data = function(d) d[d$gene %in% d$gene[duplicated(d$gene)], , drop = FALSE]) +
          ggplot2::geom_point(ggplot2::aes(fill = .data[[gcol]]), shape = 21, size = 2.5) +
          ggplot2::scale_fill_manual(values = .scroll_vdj_colours(input, lv), name = grp) +
          ggplot2::scale_y_continuous(expand = .scroll_expand0()) +
          ggplot2::labs(x = seg, y = paste0("Frequency (", unit, ", within group)"), fill = grp,
                        title = paste(seg, "usage")) +
          .scroll_base_theme() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, size = 7))
      }
      attr(p, "scroll_source") <- as.data.frame(f); p
    })

  cdr3_length <- mk("cdr3_length", "CDR3 length", "CDR3 length distribution",
    "CDR3 amino-acid length per group (clone-deduplicated).",
    list(scroll_input_choice("chain", "Chain",
           choices = function(data) c(unlist(data$manifest$vdj$cdr3_chains), "Combined")),
         scroll_input_choice("style", "Style", c("Histogram", "Density")),
         scroll_input_choice("clone_col", "Clone ID",
           choices = function(data) .scroll_vdj_clone_cols(data), widget = "select"),
         # "None" (default) draws a single distribution over all clones.
         scroll_input_choice("colorby", "Colour by",
           choices = function(data) c("None" = "", .scroll_vdj_split_cols(data)),
           widget = "select"),
         scroll_input_levels("group_levels", "Groups", none = TRUE, watch = "colorby",
           choices = function(input, data) .scroll_vdj_col_levels(input, data, "colorby")),
         scroll_style_input(.scroll_vdj_colour_control("colorby", "group_levels")),
         scroll_style_input(scroll_input_slider("aspect", "Aspect ratio", 0.4, 3, 1, 0.1))),
    function(cells, input, data) {
      rc <- .scroll_vdj_read(data, "rep_cells.parquet")
      lc <- if (identical(input$chain, "Combined")) "cdr3_combined" else paste0("cdr3_", input$chain, "_len")
      if (is.null(rc) || !lc %in% names(rc)) stop("CDR3 length for '", input$chain, "' not baked.")
      rc <- .scroll_vdj_scope(rc, cells)          # honour the global filter / subset view
      # colour/split column, or NULL for a single pooled distribution (the "None" default)
      col <- if (!is.null(input$colorby) && nzchar(input$colorby) && input$colorby %in% names(rc))
               input$colorby else NULL
      ccol <- col %||% ".grp"
      if (is.null(col)) rc$.grp <- "all"
      cid <- if (!is.null(input$clone_col) && input$clone_col %in% names(rc)) input$clone_col else "clone_id"
      d <- rc[!is.na(rc[[lc]]) & !is.na(rc[[ccol]]), , drop = FALSE]
      d$len <- d[[lc]]; d$col <- as.character(d[[ccol]])
      if (!is.null(col) && length(input$group_levels))
        d <- d[d$col %in% input$group_levels, , drop = FALSE]
      d <- d[!duplicated(paste(d[[cid]], d$col)), , drop = FALSE]     # one length per clone
      if (!nrow(d)) stop("No clones with a length for this selection.")
      if (identical(input$style, "Density")) {
        if (is.null(col)) {
          p <- ggplot2::ggplot(d, ggplot2::aes(.data$len)) +
            ggplot2::geom_density(adjust = 2, linewidth = 1,
              colour = .scroll_vdj_one_colour(input, data, "colorby")) +
            ggplot2::scale_x_continuous(expand = .scroll_expand0()) +
            ggplot2::scale_y_continuous(expand = .scroll_expand0()) +
            ggplot2::labs(x = "CDR3 length", y = "Density",
                          title = paste(input$chain, "CDR3 length")) +
            .scroll_base_theme(legend = FALSE)
        } else {
          cols <- .scroll_vdj_colours(input, .scroll_vdj_level_set(input, data, "colorby", "group_levels"))
          p <- ggplot2::ggplot(d, ggplot2::aes(.data$len, colour = .data$col)) +
            ggplot2::geom_density(adjust = 2, linewidth = 1) +
            ggplot2::scale_colour_manual(values = cols, name = col) +
            ggplot2::scale_x_continuous(expand = .scroll_expand0()) +
            ggplot2::scale_y_continuous(expand = .scroll_expand0()) +
            ggplot2::labs(x = "CDR3 length", y = "Density", colour = col,
                          title = paste(input$chain, "CDR3 length")) +
            .scroll_base_theme()
        }
        attr(p, "scroll_source") <- data.frame(length = d$len, group = d$col,
                                               clone = d[[cid]]); p
      } else {
        h <- d |> dplyr::group_by(.data$len, .data$col) |>
          dplyr::summarise(count = dplyr::n(), .groups = "drop") |>
          dplyr::group_by(.data$col) |>
          dplyr::mutate(freq = .data$count / sum(.data$count)) |> dplyr::ungroup()
        if (is.null(col)) {
          p <- ggplot2::ggplot(h, ggplot2::aes(.data$len, .data$freq)) +
            ggplot2::geom_col(fill = .scroll_vdj_one_colour(input, data, "colorby"),
                              colour = "black", linewidth = 0.2) +
            ggplot2::scale_x_continuous(expand = .scroll_expand0()) +
            ggplot2::scale_y_continuous(expand = .scroll_expand0()) +
            ggplot2::labs(x = "CDR3 length", y = "Frequency (overall)",
                          title = paste(input$chain, "CDR3 length")) +
            .scroll_base_theme(legend = FALSE)
        } else {
          cols <- .scroll_vdj_colours(input, .scroll_vdj_level_set(input, data, "colorby", "group_levels"))
          p <- ggplot2::ggplot(h, ggplot2::aes(.data$len, .data$freq, fill = .data$col)) +
            ggplot2::geom_col(position = ggplot2::position_dodge2(preserve = "single"),
                              colour = "black", linewidth = 0.2) +
            ggplot2::scale_fill_manual(values = cols, name = col) +
            ggplot2::scale_x_continuous(expand = .scroll_expand0()) +
            ggplot2::scale_y_continuous(expand = .scroll_expand0()) +
            ggplot2::labs(x = "CDR3 length", y = "Frequency (within group)", fill = col,
                          title = paste(input$chain, "CDR3 length")) +
            .scroll_base_theme()
        }
        attr(p, "scroll_source") <- as.data.frame(h); p
      }
    })

  diversity <- mk("diversity", "Diversity", "Repertoire diversity",
    "Shannon / Simpson / clonality / Gini per group, for one or more clone or gene features.",
    list(# Diversity of one or more features: a clone definition (the baked clone_id or an
         # alternate clonotype) AND/OR a segment gene (TRBV/TRBJ/...). Recomputed at runtime
         # from the clone table. Several selections -> one facet per feature (own y-scale).
         scroll_input_choice("features", "Diversity of",
           choices = function(data) c(.scroll_vdj_clone_cols(data),
                                      as.character(unlist(data$manifest$vdj$segments))),
           selected = "clone_id", multiple = TRUE, widget = "select"),
         # With several features: facet them separately, or Combine into one composite
         # feature (e.g. TRBV+TRBJ = V-J pairing diversity).
         scroll_input_choice("combine", "Multiple features", c("Separate", "Combine")),
         # "None" (default) reports one overall value for the whole repertoire.
         scroll_input_choice("by", "Group by",
           choices = function(data) c("None" = "", .scroll_vdj_split_cols(data)),
           widget = "select"),
         scroll_input_levels("group_levels", "Groups", none = TRUE, watch = "by",
           choices = function(input, data) .scroll_vdj_col_levels(input, data, "by")),
         scroll_input_choice("metric", "Metric",
           c("shannon", "simpson", "clonality", "gini", "top_clone_prop")),
         scroll_style_input(.scroll_vdj_colour_control("by", "group_levels")),
         scroll_style_input(scroll_input_slider("aspect", "Aspect ratio", 0.4, 3, 1, 0.1))),
    function(cells, input, data) {
      rc <- .scroll_vdj_read(data, "rep_cells.parquet")
      if (is.null(rc) || !nrow(rc)) stop("No repertoire store; rebuild with vdj =.")
      rc <- .scroll_vdj_scope(rc, cells)          # honour the global filter / subset view
      feats <- input$features; if (!length(feats)) feats <- "clone_id"
      feats <- feats[feats %in% names(rc)]
      if (!length(feats)) stop("Pick at least one clone id or gene.")
      # group column, or NULL for one overall value (the "None" default)
      by_col <- if (!is.null(input$by) && nzchar(input$by) && input$by %in% names(rc)) input$by else NULL
      bcol <- by_col %||% ".grp"
      if (is.null(by_col)) rc$.grp <- "all"
      if (!is.null(by_col) && length(input$group_levels))           # optional group filter
        rc <- rc[!is.na(rc[[bcol]]) & as.character(rc[[bcol]]) %in% input$group_levels, , drop = FALSE]
      len_cols <- paste0("cdr3_", unlist(data$manifest$vdj$cdr3_chains), "_len")
      len_cols <- len_cols[len_cols %in% names(rc)]
      metric <- input$metric %||% "shannon"
      # each "feature" to score: a (column, label) pair. With several features and
      # Combine, paste them into one composite column (e.g. TRBV+TRBJ = V-J pairing);
      # otherwise score each separately (faceted). A composite is NA if any part is NA.
      feat_specs <- if (identical(input$combine, "Combine") && length(feats) > 1) {
        comps <- lapply(feats, function(f) as.character(rc[[f]]))
        key <- do.call(paste, c(comps, sep = "+"))
        key[Reduce(`|`, lapply(comps, is.na))] <- NA
        rc$.combo <- key
        list(list(col = ".combo", label = paste(feats, collapse = "+")))
      } else lapply(feats, function(f) list(col = f, label = f))
      # per-level diversity for each feature spec, tagged so the facets separate them
      # (a clone repertoire and a V-gene set live on different scales).
      parts <- lapply(feat_specs, function(s) {
        dd <- .scroll_vdj_diversity(rc, s$col, bcol, len_cols)
        if (is.null(dd) || !nrow(dd) || !metric %in% names(dd)) return(NULL)
        dd$feature <- s$label; dd
      })
      d <- do.call(rbind, Filter(Negate(is.null), parts))
      if (is.null(d) || !nrow(d)) stop("No '", metric, "' for this selection.")
      if (all(is.na(d[[metric]])))
        stop("'", metric, "' is undefined for this dataset (needs two CDR3 chains).")
      d$value <- d[[metric]]
      d$feature <- factor(d$feature, levels = vapply(feat_specs, function(s) s$label, character(1)))
      ord <- names(sort(tapply(d$value, d$level, function(v) mean(v, na.rm = TRUE))))
      d$level <- factor(d$level, levels = ord)                       # levels ordered by value
      lv <- .scroll_vdj_level_set(input, data, "by", "group_levels")
      p <- ggplot2::ggplot(d, ggplot2::aes(.data$level, .data$value, fill = .data$level)) +
        ggplot2::geom_col(colour = "black", linewidth = 0.3) +
        ggplot2::scale_fill_manual(values = .scroll_vdj_colours(input, lv)) +
        ggplot2::scale_y_continuous(expand = .scroll_expand0()) +
        ggplot2::labs(x = by_col, y = metric,
                      title = paste(metric, if (is.null(by_col)) "(overall)" else paste("by", by_col))) +
        .scroll_base_theme(legend = FALSE, x_angle = 30)
      if (length(feat_specs) > 1)                                    # one facet per feature
        p <- p + ggplot2::facet_wrap(stats::as.formula("~feature"), scales = "free_y")
      attr(p, "scroll_source") <- as.data.frame(d); p
    })

  list(clone_overview, gene_usage, cdr3_length, diversity)
}
