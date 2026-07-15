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

