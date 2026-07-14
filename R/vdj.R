# VDJ (TCR/BCR immune repertoire) support. VDJ is per-cell metadata, so ingestion
# needs no new store type: at build time `scroll_build(..., vdj = vdj_spec(...))`
# bakes a small `repertoire/` Parquet store (clone table + pre-computed diversity /
# gene-usage residuals) alongside the expr store, and records a `vdj` block in the
# manifest. Four built-in panels (clone overview, V/J gene usage, CDR3 length,
# diversity) then auto-appear for any project built with a `vdj` spec. Generalized
# from the Zareie/Tablo repertoire panels; TCR vs BCR is a `chain_type` parameter.

# Default per-cell column maps by receptor. Users override any of these in vdj_spec().
.SCROLL_VDJ_DEFAULTS <- list(
  TCR = list(
    segments = c(TRBV = "v_gene_TRB", TRBJ = "j_gene_TRB",
                 TRAV = "v_gene_TRA", TRAJ = "j_gene_TRA"),
    cdr3     = c(beta = "cdr3_beta", alpha = "cdr3_alpha")),
  BCR = list(
    segments = c(IGHV = "v_gene_IGH", IGHJ = "j_gene_IGH",
                 IGKV = "v_gene_IGK", IGKJ = "j_gene_IGK",
                 IGLV = "v_gene_IGL", IGLJ = "j_gene_IGL"),
    cdr3     = c(heavy = "cdr3_heavy", light = "cdr3_light")))

#' Describe a dataset's TCR/BCR (VDJ) metadata for `scroll_build()`
#'
#' Maps a Seurat object's per-cell repertoire columns onto the schema scroll's
#' repertoire panels use. Pass the result as `scroll_build(..., vdj = vdj_spec(...))`
#' to bake the `repertoire/` store and enable the VDJ panels. Column defaults follow
#' common 10x naming and differ by `chain_type` (TCR vs BCR); override any that
#' differ in your object.
#'
#' @param chain_type `"TCR"` or `"BCR"` — sets default segment/CDR3 column names and
#'   panel labels (TRBV/TRAV… vs IGHV/IGKV…).
#' @param group_col A categorical per-cell column to group repertoire summaries by
#'   (e.g. cell type / cluster). Required.
#' @param clone_col Per-cell clonotype id. `NULL` (default) derives a clonotype from
#'   the pasted CDR3 columns.
#' @param count_col Optional per-cell clone-size column; `NULL` computes clone size
#'   from `clone_col` frequency.
#' @param segments Named character vector `label = column` of V/J gene columns
#'   (defaults by `chain_type`).
#' @param cdr3 Named character vector `chain = column` of CDR3 amino-acid columns
#'   (defaults by `chain_type`); drives CDR3-length views.
#' @param antigen_col,tissue_col,cluster_col,loc_col,broad_col Optional grouping /
#'   annotation columns (antigen split, tissue split, per-cluster diversity, and the
#'   location/broad pair used for the tissue-correlation view).
#' @param carry Extra categorical columns to keep in the clone table so panels can
#'   split / colour by them.
#' @param exclude Optional value of `antigen_col` to drop (e.g. a negative control).
#' @return A `scroll_vdj_spec` list.
#' @seealso [scroll_build()]
#' @export
vdj_spec <- function(chain_type = c("TCR", "BCR"), group_col, clone_col = NULL,
                     count_col = NULL, segments = NULL, cdr3 = NULL,
                     antigen_col = NULL, tissue_col = NULL, cluster_col = NULL,
                     loc_col = NULL, broad_col = NULL, carry = character(),
                     exclude = NULL) {
  chain_type <- match.arg(chain_type)
  if (missing(group_col) || is.null(group_col))
    stop("vdj_spec(): `group_col` is required.", call. = FALSE)
  d <- .SCROLL_VDJ_DEFAULTS[[chain_type]]
  structure(list(
    chain_type = chain_type, group_col = group_col, clone_col = clone_col,
    count_col = count_col,
    segments = if (is.null(segments)) d$segments else segments,
    cdr3 = if (is.null(cdr3)) d$cdr3 else cdr3,
    antigen_col = antigen_col, tissue_col = tissue_col, cluster_col = cluster_col,
    loc_col = loc_col, broad_col = broad_col, carry = carry, exclude = exclude
  ), class = "scroll_vdj_spec")
}

# ---- statistics (base R; no vegan) ------------------------------------------

# Standardized residuals of a group x gene contingency table (Monte-Carlo chi-sq,
# baked once at build). Returns long (group, gene, residual) + p_value, or NULL.
.scroll_vdj_chisq <- function(df, group_col, gene_col) {
  x <- df[!is.na(df[[gene_col]]) & !is.na(df[[group_col]]), , drop = FALSE]
  ct <- table(x[[group_col]], x[[gene_col]])
  if (nrow(ct) < 2 || ncol(ct) < 2) return(NULL)
  r <- suppressWarnings(tryCatch(
    stats::chisq.test(ct, simulate.p.value = TRUE, B = 2000), error = function(e) NULL))
  if (is.null(r)) return(NULL)
  out <- as.data.frame(as.table(r$stdres), stringsAsFactors = FALSE)
  names(out) <- c("group", "gene", "residual")
  out$p_value <- r$p.value
  out
}

# Per-level repertoire diversity metrics from a per-cell clone vector.
.scroll_vdj_diversity <- function(df, clone_col, by_col, cdr3_cols = character()) {
  d <- df[!is.na(df[[clone_col]]) & !is.na(df[[by_col]]), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  parts <- split(seq_len(nrow(d)), as.character(d[[by_col]]))
  rows <- lapply(names(parts), function(lv) {
    cl <- as.character(d[[clone_col]][parts[[lv]]])
    n <- length(cl); f <- as.numeric(table(cl)); p <- f / n; k <- length(f)
    shannon <- -sum(p * log(p))
    fs <- sort(f)
    gini <- 2 * sum(seq_len(k) * fs) / (k * sum(fs)) - (k + 1) / k
    paired <- if (length(cdr3_cols) == 2)
      mean(!is.na(d[[cdr3_cols[1]]][parts[[lv]]]) & !is.na(d[[cdr3_cols[2]]][parts[[lv]]]))
      else NA_real_
    data.frame(level = lv, n_cells = n, n_clones = k,
               shannon = shannon, simpson = 1 - sum(p^2),
               clonality = if (k > 1) 1 - shannon / log(k) else 0,
               gini = gini, top_clone_prop = max(f) / n, paired_rate = paired,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# Clone-frequency correlation across (broad x location), long form.
.scroll_vdj_tissue_corr <- function(df, clone_col, broad_col, loc_col, count_col) {
  d <- df[!is.na(df[[clone_col]]), , drop = FALSE]
  if (!is.null(count_col) && count_col %in% names(d))
    d <- d[!is.na(d[[count_col]]) & d[[count_col]] >= 2, , drop = FALSE]
  if (!nrow(d)) return(NULL)
  d$.grp <- paste0(as.character(d[[broad_col]]), "_", as.character(d[[loc_col]]))
  ct <- table(as.character(d[[clone_col]]), d$.grp)
  if (ncol(ct) < 2) return(NULL)
  m <- suppressWarnings(stats::cor(as.matrix(unclass(ct)), method = "pearson"))
  if (is.null(m) || any(!is.finite(m))) m[!is.finite(m)] <- NA_real_
  data.frame(g1 = rep(rownames(m), times = ncol(m)),
             g2 = rep(colnames(m), each = nrow(m)),
             cor = as.numeric(m), stringsAsFactors = FALSE)
}

# ---- bake -------------------------------------------------------------------

# Write the repertoire/ store (clone table + diversity/chisq/tissue_corr) from a
# per-cell metadata frame + a vdj_spec. Called by scroll_build(). Returns the
# `vdj` manifest block (or NULL if no clones survive).
.scroll_bake_vdj <- function(md, outdir, spec) {
  seg <- spec$segments[spec$segments %in% names(md)]
  cd3 <- spec$cdr3[spec$cdr3 %in% names(md)]
  gcol <- spec$group_col
  if (!gcol %in% names(md)) stop("vdj: group_col '", gcol, "' not in metadata.", call. = FALSE)

  # clone id: given, or derived from the pasted CDR3 columns
  clone <- if (!is.null(spec$clone_col) && spec$clone_col %in% names(md))
    as.character(md[[spec$clone_col]])
  else if (length(cd3))
    do.call(paste, c(lapply(cd3, function(c) as.character(md[[c]])), sep = "|"))
  else stop("vdj: no clone_col and no CDR3 columns to derive a clonotype from.", call. = FALSE)
  clone[clone %in% c("", "NA|NA", "NA")] <- NA

  keep <- !is.na(clone) & !is.na(md[[gcol]]) & nzchar(as.character(md[[gcol]]))
  if (!is.null(spec$exclude) && !is.null(spec$antigen_col) && spec$antigen_col %in% names(md))
    keep <- keep & (is.na(md[[spec$antigen_col]]) | md[[spec$antigen_col]] != spec$exclude)
  if (!any(keep)) return(NULL)
  m <- md[keep, , drop = FALSE]; clone <- clone[keep]

  gv <- function(col) if (!is.null(col) && col %in% names(m)) as.character(m[[col]]) else NA_character_
  rep_cells <- data.frame(clone_id = clone, group = as.character(m[[gcol]]),
                          antigen = gv(spec$antigen_col), tissue = gv(spec$tissue_col),
                          loc = gv(spec$loc_col), stringsAsFactors = FALSE)
  rep_cells$clone_count <- if (!is.null(spec$count_col) && spec$count_col %in% names(m))
    as.numeric(m[[spec$count_col]]) else as.numeric(stats::ave(clone, clone, FUN = length))
  for (lab in names(seg)) rep_cells[[lab]] <- as.character(m[[seg[[lab]]]])         # gene per segment
  for (ch in names(cd3)) rep_cells[[paste0("cdr3_", ch, "_len")]] <- nchar(as.character(m[[cd3[[ch]]]]))
  lens <- grep("^cdr3_.*_len$", names(rep_cells), value = TRUE)
  if (length(lens)) rep_cells$cdr3_combined <- rowSums(as.matrix(rep_cells[lens]), na.rm = TRUE)
  for (cc in setdiff(intersect(spec$carry, names(m)), names(rep_cells)))
    rep_cells[[cc]] <- as.character(m[[cc]])

  dir.create(file.path(outdir, "repertoire"), showWarnings = FALSE, recursive = TRUE)
  wp <- function(x, f) if (!is.null(x) && nrow(x))
    arrow::write_parquet(x, file.path(outdir, "repertoire", f), compression = "zstd")
  wp(rep_cells, "rep_cells.parquet")

  m_dedup <- m[!duplicated(clone), , drop = FALSE]
  chisq <- do.call(rbind, Filter(Negate(is.null), lapply(names(seg), function(lab) {
    r <- .scroll_vdj_chisq(m_dedup, gcol, seg[[lab]]); if (!is.null(r)) r$segment <- lab; r })))
  wp(chisq, "chisq.parquet")

  clone_dd <- clone[!duplicated(clone)]                       # for derived-clone diversity
  m_dd <- m; m_dd$.clone <- clone
  div <- do.call(rbind, Filter(Negate(is.null), list(
    { d <- .scroll_vdj_diversity(m_dd, ".clone", gcol, unname(cd3)); if (!is.null(d)) d$level_type <- "group"; d },
    if (!is.null(spec$cluster_col) && spec$cluster_col %in% names(m_dd)) {
      d <- .scroll_vdj_diversity(m_dd, ".clone", spec$cluster_col, unname(cd3))
      if (!is.null(d)) d$level_type <- "cluster"; d })))
  wp(div, "diversity.parquet")

  tcorr <- NULL
  if (!is.null(spec$loc_col) && !is.null(spec$broad_col) &&
      all(c(spec$loc_col, spec$broad_col) %in% names(m))) {
    m2 <- m; m2$.clone <- clone
    tcorr <- .scroll_vdj_tissue_corr(m2, ".clone", spec$broad_col, spec$loc_col, spec$count_col)
    wp(tcorr, "tissue_corr.parquet")
  }

  list(chain_type = spec$chain_type, segments = as.list(names(seg)),
       cdr3_chains = as.list(names(cd3)), group_col = gcol,
       has_antigen = !is.null(spec$antigen_col) && spec$antigen_col %in% names(md),
       has_tissue  = !is.null(spec$tissue_col)  && spec$tissue_col  %in% names(md),
       has_cluster = !is.null(div) && "cluster" %in% div$level_type,
       has_tissue_corr = !is.null(tcorr) && nrow(tcorr) > 0,
       n_clones = length(unique(clone)))
}

# ---- runtime: read the baked store ------------------------------------------

.scroll_vdj_read <- function(data, f) {
  p <- file.path(data$dir, "repertoire", f)
  if (!file.exists(p)) return(NULL)
  tryCatch(as.data.frame(arrow::read_parquet(p)), error = function(e) NULL)
}
# categorical rep_cells columns available to split/colour by
.scroll_vdj_split_cols <- function(data) {
  rc <- .scroll_vdj_read(data, "rep_cells.parquet"); if (is.null(rc)) return(character())
  cand <- setdiff(names(rc), c("clone_id", "clone_count", grep("_len$", names(rc), value = TRUE),
                               "cdr3_combined", names(data$manifest$vdj$segments)))
  cand <- setdiff(cand, unlist(data$manifest$vdj$segments))
  cand[vapply(cand, function(c) is.character(rc[[c]]) && any(!is.na(rc[[c]])), logical(1))]
}
.scroll_vdj_levels <- function(rc, col) if (col %in% names(rc))
  sort(unique(stats::na.omit(as.character(rc[[col]])))) else character()

# ---- the four panels --------------------------------------------------------

# Assembled as built-in panels gated on manifest$vdj (see .scroll_assemble_panels).
.scroll_vdj_panels <- function() {
  gate <- function(m) !is.null(m$vdj)
  mk <- function(id, label, title, desc, controls, plot)
    c(list(id = id, label = label, title = title, desc = desc, when = gate),
      .scroll_plot_panel_uiserver(plot, controls, compute = FALSE))

  clone_overview <- mk("clone_overview", "Clone overview", "Clonal expansion overview",
    "Rank-abundance of clone sizes and expansion-category composition (clone-level).",
    list(scroll_input_choice("view", "View", c("Rank-abundance", "Expansion composition")),
         scroll_input_choice("group", "Group by",
           choices = function(data) intersect(c("group", "tissue", "antigen"),
                                              names(.scroll_vdj_read(data, "rep_cells.parquet"))))),
    function(cells, input, data) {
      rc <- .scroll_vdj_read(data, "rep_cells.parquet")
      if (is.null(rc) || !nrow(rc)) stop("No repertoire store; rebuild with vdj =.")
      splitc <- if (!is.null(rc$antigen) && any(!is.na(rc$antigen))) "antigen" else "group"
      if (identical(input$view, "Rank-abundance")) {
        d <- rc |> dplyr::group_by(.data$clone_id, .data[[splitc]]) |>
          dplyr::summarise(count = dplyr::n(), .groups = "drop") |>
          dplyr::arrange(.data[[splitc]], dplyr::desc(.data$count)) |>
          dplyr::group_by(.data[[splitc]]) |>
          dplyr::mutate(rank = dplyr::row_number()) |> dplyr::ungroup()
        ggplot2::ggplot(d, ggplot2::aes(.data$rank, .data$count, colour = .data[[splitc]])) +
          ggplot2::geom_point(size = 0.5) + ggplot2::scale_y_log10() +
          ggplot2::facet_wrap(stats::as.formula(paste0("~`", splitc, "`")), scales = "free_x") +
          ggplot2::labs(x = "Clone rank", y = "Clone size (cells, log10)", colour = splitc,
                        title = "Clone rank-abundance") +
          ggplot2::theme_minimal() +
          ggplot2::theme(legend.position = "none", panel.grid.minor = ggplot2::element_blank())
      } else {
        grp <- if (input$group %in% names(rc) && any(!is.na(rc[[input$group]]))) input$group else "group"
        cl <- rc |> dplyr::group_by(.data$clone_id) |>
          dplyr::summarise(count = dplyr::n(),
            g = names(sort(table(.data[[grp]]), decreasing = TRUE))[1], .groups = "drop") |>
          dplyr::mutate(expansion = cut(.data$count, c(0, 1, 4, 19, 99, Inf),
            labels = c("Single", "Small", "Medium", "Large", "Hyperexpanded")))
        ggplot2::ggplot(cl, ggplot2::aes(.data$g, fill = .data$expansion)) +
          ggplot2::geom_bar(position = "fill", colour = "white", linewidth = 0.2) +
          ggplot2::scale_y_continuous(labels = scales::percent) +
          ggplot2::scale_fill_brewer(palette = "YlOrRd") +
          ggplot2::labs(x = grp, y = "Fraction of clones", fill = "Expansion",
                        title = "Expansion-category composition") +
          ggplot2::theme_minimal() +
          ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
      }
    })

  gene_usage <- mk("gene_usage", "V/J gene usage", "TCR/BCR V/J gene usage",
    "Per-group V/J segment frequency, or the chi-square standardized-residual heatmap.",
    list(scroll_input_choice("segment", "Segment",
           choices = function(data) unlist(data$manifest$vdj$segments)),
         scroll_input_choice("view", "View", c("Frequency", "Chi-square residuals"))),
    function(cells, input, data) {
      rc <- .scroll_vdj_read(data, "rep_cells.parquet"); seg <- input$segment
      if (is.null(rc) || !seg %in% names(rc)) stop("Segment '", seg, "' not baked.")
      d <- rc[!is.na(rc[[seg]]), , drop = FALSE]; d$gene <- as.character(d[[seg]])
      if (identical(input$view, "Chi-square residuals")) {
        ch <- .scroll_vdj_read(data, "chisq.parquet")
        ch <- if (!is.null(ch)) ch[ch$segment == seg, , drop = FALSE] else NULL
        if (is.null(ch) || !nrow(ch)) stop("No chi-square residuals for ", seg, ".")
        ch$residual <- pmax(pmin(ch$residual, 3), -3)
        ord <- names(sort(tapply(abs(ch$residual), ch$gene, max)))
        ch$gene <- factor(ch$gene, levels = ord)
        ggplot2::ggplot(ch, ggplot2::aes(.data$gene, .data$group, fill = .data$residual)) +
          ggplot2::geom_tile(colour = "grey90") +
          ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                        midpoint = 0, limits = c(-3, 3)) +
          ggplot2::labs(x = seg, y = NULL, fill = "std\nresidual",
                        title = paste(seg, "usage vs group")) +
          ggplot2::theme_minimal() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, size = 8),
                         panel.grid = ggplot2::element_blank())
      } else {
        f <- d |> dplyr::group_by(.data$gene, .data$group) |>
          dplyr::summarise(count = dplyr::n(), .groups = "drop") |>
          dplyr::group_by(.data$group) |>
          dplyr::mutate(freq = .data$count / sum(.data$count)) |> dplyr::ungroup()
        ggplot2::ggplot(f, ggplot2::aes(.data$gene, .data$freq)) +
          ggplot2::geom_line(ggplot2::aes(group = .data$gene), colour = "grey70", linewidth = 0.3) +
          ggplot2::geom_point(ggplot2::aes(fill = .data$group), shape = 21, size = 2.5) +
          ggplot2::labs(x = seg, y = "Frequency (within group)", fill = "group",
                        title = paste(seg, "usage")) +
          ggplot2::theme_minimal() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, size = 7),
                         panel.grid.major.x = ggplot2::element_blank())
      }
    })

  cdr3_length <- mk("cdr3_length", "CDR3 length", "CDR3 length distribution",
    "CDR3 amino-acid length per group (clone-deduplicated).",
    list(scroll_input_choice("chain", "Chain",
           choices = function(data) c(unlist(data$manifest$vdj$cdr3_chains), "Combined")),
         scroll_input_choice("style", "Style", c("Histogram", "Density")),
         scroll_input_choice("colorby", "Colour by",
           choices = function(data) .scroll_vdj_split_cols(data))),
    function(cells, input, data) {
      rc <- .scroll_vdj_read(data, "rep_cells.parquet")
      lc <- if (identical(input$chain, "Combined")) "cdr3_combined" else paste0("cdr3_", input$chain, "_len")
      if (is.null(rc) || !lc %in% names(rc)) stop("CDR3 length for '", input$chain, "' not baked.")
      col <- if (!is.null(input$colorby) && input$colorby %in% names(rc)) input$colorby else "group"
      d <- rc[!is.na(rc[[lc]]) & !is.na(rc[[col]]), , drop = FALSE]
      d$len <- d[[lc]]; d$col <- as.character(d[[col]])
      d <- d[!duplicated(paste(d$clone_id, d$col)), , drop = FALSE]
      if (!nrow(d)) stop("No clones with a length for this selection.")
      if (identical(input$style, "Density")) {
        ggplot2::ggplot(d, ggplot2::aes(.data$len, colour = .data$col)) +
          ggplot2::geom_density(adjust = 2, linewidth = 1) +
          ggplot2::labs(x = "CDR3 length", y = "Density", colour = col,
                        title = paste(input$chain, "CDR3 length")) +
          ggplot2::theme_minimal() + ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
      } else {
        h <- d |> dplyr::group_by(.data$len, .data$col) |>
          dplyr::summarise(count = dplyr::n(), .groups = "drop") |>
          dplyr::group_by(.data$col) |>
          dplyr::mutate(freq = .data$count / sum(.data$count)) |> dplyr::ungroup()
        ggplot2::ggplot(h, ggplot2::aes(.data$len, .data$freq, fill = .data$col)) +
          ggplot2::geom_col(position = ggplot2::position_dodge2(preserve = "single"),
                            colour = "black", linewidth = 0.2) +
          ggplot2::labs(x = "CDR3 length", y = "Frequency (within group)", fill = col,
                        title = paste(input$chain, "CDR3 length")) +
          ggplot2::theme_minimal() + ggplot2::theme(panel.grid = ggplot2::element_blank())
      }
    })

  diversity <- mk("diversity", "Diversity", "Repertoire diversity",
    "Shannon / Simpson / clonality / Gini per group (or cluster), and tissue correlation.",
    list(scroll_input_choice("view", "View", c("Diversity metric", "Tissue correlation")),
         scroll_input_choice("by", "Group by",
           choices = function(data) if (isTRUE(data$manifest$vdj$has_cluster)) c("group", "cluster") else "group"),
         scroll_input_choice("metric", "Metric",
           c("shannon", "simpson", "clonality", "gini", "top_clone_prop"))),
    function(cells, input, data) {
      if (identical(input$view, "Tissue correlation")) {
        tc <- .scroll_vdj_read(data, "tissue_corr.parquet")
        if (is.null(tc) || !nrow(tc)) stop("Tissue correlation not available (needs loc_col + broad_col).")
        lev <- unique(c(tc$g1, tc$g2))
        tc$g1 <- factor(tc$g1, levels = lev); tc$g2 <- factor(tc$g2, levels = rev(lev))
        return(ggplot2::ggplot(tc, ggplot2::aes(.data$g1, .data$g2, fill = .data$cor)) +
          ggplot2::geom_tile(colour = "white") +
          ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data$cor)), size = 3) +
          ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                        midpoint = 0, limits = c(-1, 1)) +
          ggplot2::labs(x = NULL, y = NULL, fill = "Pearson r", title = "Clone-frequency correlation") +
          ggplot2::coord_fixed() + ggplot2::theme_minimal() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                         panel.grid = ggplot2::element_blank()))
      }
      dv <- .scroll_vdj_read(data, "diversity.parquet")
      if (is.null(dv) || !nrow(dv)) stop("No diversity table baked.")
      d <- dv[dv$level_type == input$by, , drop = FALSE]
      if (!nrow(d) || !input$metric %in% names(d)) stop("No '", input$metric, "' for '", input$by, "'.")
      d$value <- d[[input$metric]]; d$level <- factor(d$level, levels = d$level[order(d$value)])
      ggplot2::ggplot(d, ggplot2::aes(.data$level, .data$value, fill = .data$value)) +
        ggplot2::geom_col(colour = "black", linewidth = 0.3) +
        ggplot2::scale_fill_viridis_c() +
        ggplot2::labs(x = input$by, y = input$metric, title = paste(input$metric, "by", input$by)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none",
                       axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
                       panel.grid.major.x = ggplot2::element_blank())
    })

  list(clone_overview, gene_usage, cdr3_length, diversity)
}
