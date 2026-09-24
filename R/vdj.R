# VDJ (TCR/BCR immune repertoire) support. VDJ is per-cell metadata, so ingestion
# needs no new store type: at build time `scroll_build(..., vdj = vdj_spec(...))`
# bakes a small `repertoire/` Parquet store (clone table + pre-computed diversity /
# gene-usage residuals) alongside the expr store, and records a `vdj` block in the
# manifest. Four built-in panels (clone overview, V/J gene usage, CDR3 length,
# diversity) then auto-appear for any project built with a `vdj` spec. Generalized
# from the Zareie/Tablo repertoire panels; TCR vs BCR is a `chain_type` parameter.

# Per-cell column maps by upstream source and receptor. The runtime panels read a fixed
# normalized schema (rep_cells.parquet), so supporting a new package is purely a build-time
# mapping: `source=` (or auto-detection in .scroll_vdj_resolve_spec) picks one of these and
# fills any segment/CDR3/clone/count column the user didn't pass. `scroll` is the historical
# default; `airr` = dandelion per-cell obs (_VDJ/_VJ), `platypus` = VDJ_/VJ_ columns.
.SCROLL_VDJ_PRESETS <- list(
  scroll = list(
    TCR = list(
      segments = c(TRBV = "v_gene_TRB", TRBJ = "j_gene_TRB",
                   TRAV = "v_gene_TRA", TRAJ = "j_gene_TRA"),
      cdr3     = c(beta = "cdr3_beta", alpha = "cdr3_alpha")),
    BCR = list(
      segments = c(IGHV = "v_gene_IGH", IGHJ = "j_gene_IGH",
                   IGKV = "v_gene_IGK", IGKJ = "j_gene_IGK",
                   IGLV = "v_gene_IGL", IGLJ = "j_gene_IGL"),
      cdr3     = c(heavy = "cdr3_heavy", light = "cdr3_light"))),
  # dandelion / AIRR per-cell obs: _VDJ = heavy/beta chain, _VJ = light/alpha chain.
  airr = list(
    TCR = list(
      segments = c(TRBV = "v_call_VDJ", TRBJ = "j_call_VDJ",
                   TRAV = "v_call_VJ",  TRAJ = "j_call_VJ"),
      cdr3     = c(beta = "junction_aa_VDJ", alpha = "junction_aa_VJ"),
      clone_col = "clone_id", count_col = "duplicate_count"),
    BCR = list(
      segments = c(IGHV = "v_call_VDJ", IGHJ = "j_call_VDJ",
                   IGLV = "v_call_VJ",  IGLJ = "j_call_VJ"),
      cdr3     = c(heavy = "junction_aa_VDJ", light = "junction_aa_VJ"),
      clone_col = "clone_id", count_col = "duplicate_count")),
  # Platypus VDJ.GEX per-cell matrix: VDJ_ = heavy/beta chain, VJ_ = light/alpha chain.
  platypus = list(
    TCR = list(
      segments = c(TRBV = "VDJ_vgene", TRBJ = "VDJ_jgene",
                   TRAV = "VJ_vgene",  TRAJ = "VJ_jgene"),
      cdr3     = c(beta = "VDJ_cdr3s_aa", alpha = "VJ_cdr3s_aa"),
      clone_col = "clonotype_id_10x"),
    BCR = list(
      segments = c(IGHV = "VDJ_vgene", IGHJ = "VDJ_jgene",
                   IGLV = "VJ_vgene",  IGLJ = "VJ_jgene"),
      cdr3     = c(heavy = "VDJ_cdr3s_aa", light = "VJ_cdr3s_aa"),
      clone_col = "clonotype_id_10x")))

# Back-compat: eager defaulting in vdj_spec() reads the historical `scroll` map.
.SCROLL_VDJ_DEFAULTS <- .SCROLL_VDJ_PRESETS$scroll

#' Describe a dataset's TCR/BCR (VDJ) metadata for `scroll_build()`
#'
#' Maps a Seurat object's per-cell repertoire columns onto the schema scroll's
#' repertoire panels use. Pass the result as `scroll_build(..., vdj = vdj_spec(...))`
#' to bake the `repertoire/` store and enable the VDJ panels.
#'
#' `scroll` reads TCR/BCR data straight from the common upstream tools — set `source`
#' (or leave it `"auto"` to detect) and the segment/CDR3/clone columns are mapped for
#' you; override any of them explicitly when your object differs:
#'   * `"scRepertoire"` — Seurat objects from `combineTCR()`/`combineBCR()` +
#'     `combineExpression()`; the compound `CTgene`/`CTaa`/`CTstrict` columns are parsed
#'     into per-segment/per-chain columns.
#'   * `"airr"` — dandelion / AIRR per-cell obs (`v_call_VDJ`/`v_call_VJ`,
#'     `junction_aa_VDJ`/`_VJ`, `clone_id`).
#'   * `"platypus"` — Platypus VDJ.GEX per-cell columns (`VDJ_vgene`/`VJ_vgene`,
#'     `VDJ_cdr3s_aa`/`VJ_cdr3s_aa`).
#'   * `"scroll"` — scroll's own `v_gene_TRB`/`cdr3_beta`… convention (the historical
#'     default).
#' Long, one-row-per-contig tables (raw 10x `filtered_contig_annotations.csv` / long
#' AIRR `.tsv`) must be collapsed to per-cell columns upstream (e.g. with scRepertoire
#' or dandelion) before building.
#'
#' @param chain_type `"TCR"` or `"BCR"` — sets default segment/CDR3 column names and
#'   panel labels (TRBV/TRAV… vs IGHV/IGKV…).
#' @param source Upstream convention to read: `"auto"` (default; detect from the
#'   metadata columns), `"scroll"`, `"scRepertoire"`, `"airr"`, or `"platypus"`. Any
#'   column you pass explicitly (`segments`, `cdr3`, `clone_col`, …) overrides the preset.
#'   Column names are matched ignoring case when the exact name is absent (so
#'   `v_call_vdj` satisfies the AIRR preset's `v_call_VDJ`); each such match is
#'   reported with a message.
#' @param group_col A categorical per-cell column to group repertoire summaries by
#'   (e.g. cell type / cluster). Required.
#' @param clone_col Per-cell clonotype id. `NULL` (default) uses the preset's clone
#'   column, else (for a non-`"scroll"` source) the first present of `clone_id`,
#'   `strict_clone_id`, `clonotype_id`, `raw_clonotype_id`, `clonotype`, else derives
#'   a clonotype from the pasted CDR3 columns.
#' @param count_col Optional per-cell clone-size column. `NULL` uses
#'   `<clone_col>_count` / `<clone_col>_size` when present, else computes clone size
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
                     source = c("auto", "scroll", "scRepertoire", "airr", "platypus"),
                     antigen_col = NULL, tissue_col = NULL, cluster_col = NULL,
                     loc_col = NULL, broad_col = NULL, carry = character(),
                     exclude = NULL) {
  chain_type <- match.arg(chain_type)
  source <- match.arg(source)
  if (missing(group_col) || is.null(group_col))
    stop("vdj_spec(): `group_col` is required.", call. = FALSE)
  d <- .SCROLL_VDJ_DEFAULTS[[chain_type]]
  structure(list(
    chain_type = chain_type, source = source,
    group_col = group_col, clone_col = clone_col, count_col = count_col,
    # eager scroll defaults preserve back-compat; the resolver only augments/replaces
    # them when they don't match the object (see .scroll_vdj_resolve_spec).
    segments = if (is.null(segments)) d$segments else segments,
    cdr3 = if (is.null(cdr3)) d$cdr3 else cdr3,
    segments_explicit = !is.null(segments), cdr3_explicit = !is.null(cdr3),
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

# ---- source resolution (read scRepertoire / AIRR / Platypus / scroll) -------

# Tokens that mean "no value" in a clone id or a "|"-pasted clonotype field.
.SCROLL_VDJ_NA_TOKENS <- c("", "na", "none", "nan", ".", "unassigned")

# TRUE where a (possibly "|"-pasted) clone id is entirely missing — every field is an
# NA token. Broader than the old fixed c("", "NA", "NA|NA") so 3-field pastes and other
# conventions ("None", "NA_NA"-style) are caught.
.scroll_vdj_blank_clone <- function(clone) {
  parts <- strsplit(as.character(clone), "|", fixed = TRUE)
  vapply(parts, function(p) all(tolower(trimws(p)) %in% .SCROLL_VDJ_NA_TOKENS), logical(1))
}

# Expand scRepertoire's compound CT* columns into per-cell columns matching scroll's
# `scroll` map (v_gene_TRB / cdr3_beta …). scRepertoire stores, per cell, CTgene / CTaa
# with chains separated by "_" and segments within a chain by ".". Chain order is
# receptor-specific: TCR = alpha(TRA: V.J.C) _ beta(TRB: V.D.J.C); BCR = heavy(IGH) _ light.
.scroll_vdj_parse_screpertoire <- function(md, chain_type) {
  field <- function(col, chain_ix, seg_ix) {
    v <- as.character(md[[col]]); v[is.na(v)] <- ""
    vapply(strsplit(v, "_", fixed = TRUE), function(parts) {
      if (length(parts) < chain_ix) return(NA_character_)
      seg <- strsplit(parts[[chain_ix]], ".", fixed = TRUE)[[1]]
      x <- if (length(seg) >= seg_ix) seg[[seg_ix]] else NA_character_
      if (is.na(x) || !nzchar(x) || x %in% c("NA", "None")) NA_character_ else x
    }, character(1))
  }
  has <- function(col) col %in% names(md) && any(nzchar(as.character(md[[col]])), na.rm = TRUE)
  if (identical(chain_type, "BCR")) {
    if (has("CTgene")) {                                # heavy = chain 1, light = chain 2
      md$v_gene_IGH <- field("CTgene", 1, 1); md$j_gene_IGH <- field("CTgene", 1, 3)
      md$v_gene_IGL <- field("CTgene", 2, 1); md$j_gene_IGL <- field("CTgene", 2, 2)
    }
    if (has("CTaa")) { md$cdr3_heavy <- field("CTaa", 1, 1); md$cdr3_light <- field("CTaa", 2, 1) }
  } else {
    if (has("CTgene")) {                                # alpha = chain 1, beta = chain 2
      md$v_gene_TRA <- field("CTgene", 1, 1); md$j_gene_TRA <- field("CTgene", 1, 2)
      md$v_gene_TRB <- field("CTgene", 2, 1); md$j_gene_TRB <- field("CTgene", 2, 3)
    }
    if (has("CTaa")) { md$cdr3_alpha <- field("CTaa", 1, 1); md$cdr3_beta <- field("CTaa", 2, 1) }
  }
  clone_src <- Find(function(c) c %in% names(md), c("CTstrict", "CTaa", "CTnt"))
  md$.scroll_ct_clone <- if (!is.null(clone_src)) as.character(md[[clone_src]]) else NA_character_
  freq_src <- Find(function(c) c %in% names(md) && is.numeric(md[[c]]),
                   c("clonalFrequency", "Frequency", "clonalProportion"))
  if (!is.null(freq_src)) md$.scroll_ct_count <- as.numeric(md[[freq_src]])
  md
}

# Resolve a vdj_spec's column maps against the actual metadata `md` before baking.
# Explicit user maps always win; otherwise scRepertoire compounds are parsed and, for
# `source="auto"`, the preset whose columns are most present is chosen. Returns the
# (possibly augmented) `md` and the filled `spec`.
.scroll_vdj_resolve_spec <- function(spec, md) {
  ct <- spec$chain_type
  # case-insensitive, so a preset still scores when columns were exported as e.g.
  # `v_call_vdj` rather than dandelion's `v_call_VDJ`
  hits <- function(map) if (length(map)) sum(tolower(map) %in% tolower(names(md))) else 0L
  clone_explicit <- !is.null(spec$clone_col); count_explicit <- !is.null(spec$count_col)
  loosen <- function(spec) .scroll_vdj_loosen(spec, md, clone_explicit, count_explicit)

  # 1. scRepertoire: explicit, or auto-detected compound columns the current map misses
  if (identical(spec$source, "scRepertoire") ||
      (identical(spec$source, "auto") && "CTgene" %in% names(md) &&
       !spec$segments_explicit && hits(spec$segments) == 0)) {
    md <- .scroll_vdj_parse_screpertoire(md, ct)
    if (identical(ct, "BCR")) {
      seg <- c(IGHV = "v_gene_IGH", IGHJ = "j_gene_IGH", IGLV = "v_gene_IGL", IGLJ = "j_gene_IGL")
      cd  <- c(heavy = "cdr3_heavy", light = "cdr3_light")
    } else {
      seg <- c(TRBV = "v_gene_TRB", TRBJ = "j_gene_TRB", TRAV = "v_gene_TRA", TRAJ = "j_gene_TRA")
      cd  <- c(beta = "cdr3_beta", alpha = "cdr3_alpha")
    }
    if (!spec$segments_explicit) spec$segments <- seg
    if (!spec$cdr3_explicit)     spec$cdr3     <- cd
    if (is.null(spec$clone_col)) spec$clone_col <- ".scroll_ct_clone"
    if (is.null(spec$count_col) && ".scroll_ct_count" %in% names(md)) spec$count_col <- ".scroll_ct_count"
    spec$source <- "scRepertoire"
    return(list(spec = spec, md = md))
  }

  fill <- function(spec, p) {                            # fill only unset / zero-hit fields
    if (!spec$segments_explicit && !is.null(p$segments)) spec$segments <- p$segments
    if (!spec$cdr3_explicit && !is.null(p$cdr3))         spec$cdr3     <- p$cdr3
    if (is.null(spec$clone_col) && !is.null(p$clone_col)) spec$clone_col <- p$clone_col
    if (is.null(spec$count_col) && !is.null(p$count_col)) spec$count_col <- p$count_col
    spec
  }

  # 2. explicit name-map preset
  if (spec$source %in% c("scroll", "airr", "platypus"))
    return(list(spec = loosen(fill(spec, .SCROLL_VDJ_PRESETS[[spec$source]][[ct]])), md = md))

  # 3. auto: keep the (explicit or default) map when it matches; else sniff presets
  if (identical(spec$source, "auto")) {
    if (spec$segments_explicit || hits(spec$segments) > 0) {
      if (!spec$segments_explicit) spec$source <- "scroll"   # default map matched
      return(list(spec = loosen(spec), md = md))
    }
    scores <- vapply(names(.SCROLL_VDJ_PRESETS), function(nm)
      hits(.SCROLL_VDJ_PRESETS[[nm]][[ct]]$segments), integer(1))
    best <- names(scores)[which.max(scores)]
    if (scores[[best]] > 0) { spec <- fill(spec, .SCROLL_VDJ_PRESETS[[best]][[ct]]); spec$source <- best }
  }
  list(spec = loosen(spec), md = md)
}

# Common per-cell clonotype-id column names, tried (case-insensitively) when no
# clone_col was given and the preset's is absent.
.SCROLL_VDJ_CLONE_COLS <- c("clone_id", "strict_clone_id", "clonotype_id",
                            "raw_clonotype_id", "clonotype")

# Map wanted column names onto `have`: an exact name is kept; otherwise the one
# column equal ignoring case is used (ambiguous or no match -> left as is, so the
# usual missing-column warning still fires). Names of `want` are preserved.
.scroll_ci_cols <- function(want, have) {
  if (!length(want)) return(want)
  low <- tolower(have)
  out <- vapply(unname(want), function(w) {
    if (w %in% have) return(w)
    h <- have[low == tolower(w)]
    if (length(h) == 1L) h else w
  }, "")
  names(out) <- names(want)
  out
}

# Loosen a resolved spec against the actual columns: fix column-name case, find a
# clonotype column under a common name when none was given (preset/explicit maps
# only), and pick up a
# `<clone_col>_count` / `_size` clone-size column. An explicitly passed clone/count
# column is only case-corrected, never swapped for a different column. Reports
# what it loosened with a message.
.scroll_vdj_loosen <- function(spec, md, clone_explicit, count_explicit) {
  have <- names(md); notes <- character()
  fix <- function(x) {
    y <- .scroll_ci_cols(x, have)
    ch <- !is.na(x) & y != x
    if (any(ch)) notes <<- c(notes, paste0(x[ch], " -> ", y[ch]))
    y
  }
  spec$segments <- fix(spec$segments)
  spec$cdr3     <- fix(spec$cdr3)
  if (!is.null(spec$clone_col)) spec$clone_col <- fix(spec$clone_col)
  # (not for the scroll convention, whose documented NULL clone_col = derive from CDR3)
  if (!clone_explicit && !identical(spec$source, "scroll") &&
      (is.null(spec$clone_col) || !spec$clone_col %in% have)) {
    cand <- .scroll_ci_cols(.SCROLL_VDJ_CLONE_COLS, have)
    cand <- cand[cand %in% have]
    spec$clone_col <- if (length(cand)) cand[[1]] else NULL
    if (length(cand)) notes <- c(notes, paste0("clonotype = ", cand[[1]]))
  }
  if (!is.null(spec$count_col)) spec$count_col <- fix(spec$count_col)
  if (!count_explicit && !is.null(spec$clone_col) &&
      (is.null(spec$count_col) || !spec$count_col %in% have)) {
    cand <- .scroll_ci_cols(paste0(spec$clone_col, c("_count", "_size")), have)
    cand <- cand[cand %in% have]
    spec$count_col <- if (length(cand)) cand[[1]] else NULL
    if (length(cand)) notes <- c(notes, paste0("clone size = ", cand[[1]]))
  }
  if (length(notes))
    message("vdj (source '", spec$source, "'): matched ", paste(notes, collapse = "; "))
  spec
}

# ---- bake -------------------------------------------------------------------

# Write the repertoire/ store (clone table + diversity/chisq/tissue_corr) from a
# per-cell metadata frame + a vdj_spec. Called by scroll_build(). Returns the
# `vdj` manifest block (or NULL if no clones survive).
.scroll_bake_vdj <- function(md, outdir, spec) {
  r <- .scroll_vdj_resolve_spec(spec, md); spec <- r$spec; md <- r$md
  seg <- spec$segments[spec$segments %in% names(md)]
  cd3 <- spec$cdr3[spec$cdr3 %in% names(md)]
  gcol <- spec$group_col
  if (!gcol %in% names(md)) stop("vdj: group_col '", gcol, "' not in metadata.", call. = FALSE)
  miss <- setdiff(unname(spec$segments), names(md))      # mapped but absent -> warn, don't drop silently
  if (length(seg) && length(miss))
    warning("vdj: mapped segment/CDR3 columns not found in metadata: ",
            paste(miss, collapse = ", "), ". Set source= or pass segments=/cdr3= to fix.",
            call. = FALSE)
  if (!length(seg))
    warning("vdj: no V/J segment columns resolved (source='", spec$source,
            "'); gene-usage and diversity-by-gene views will be empty. ",
            "Set source= (\"scRepertoire\"/\"airr\"/\"platypus\") or pass segments=.",
            call. = FALSE)

  # clone id: given, or derived from the pasted CDR3 columns
  clone <- if (!is.null(spec$clone_col) && spec$clone_col %in% names(md))
    as.character(md[[spec$clone_col]])
  else if (length(cd3))
    do.call(paste, c(lapply(cd3, function(c) as.character(md[[c]])), sep = "|"))
  else stop("vdj: no clone_col and no CDR3 columns to derive a clonotype from.", call. = FALSE)
  clone[.scroll_vdj_blank_clone(clone)] <- NA

  keep <- !is.na(clone) & !is.na(md[[gcol]]) & nzchar(as.character(md[[gcol]]))
  if (!is.null(spec$exclude) && !is.null(spec$antigen_col) && spec$antigen_col %in% names(md))
    keep <- keep & (is.na(md[[spec$antigen_col]]) | md[[spec$antigen_col]] != spec$exclude)
  if (!any(keep)) return(NULL)
  m <- md[keep, , drop = FALSE]; clone <- clone[keep]

  gv <- function(col) if (!is.null(col) && col %in% names(m)) as.character(m[[col]]) else NA_character_
  # `cell` is the barcode (== cells.parquet$cell), so the runtime panels can narrow the
  # per-cell table to the app's active cells (global filter + subset view).
  rep_cells <- data.frame(cell = rownames(m), clone_id = clone, group = as.character(m[[gcol]]),
                          antigen = gv(spec$antigen_col), tissue = gv(spec$tissue_col),
                          loc = gv(spec$loc_col), stringsAsFactors = FALSE)
  rep_cells$clone_count <- if (!is.null(spec$count_col) && spec$count_col %in% names(m))
    as.numeric(m[[spec$count_col]]) else as.numeric(stats::ave(clone, clone, FUN = length))
  for (lab in names(seg)) rep_cells[[lab]] <- as.character(m[[seg[[lab]]]])         # gene per segment
  for (ch in names(cd3)) {                                       # CDR3 aa sequence + length
    aa <- as.character(m[[cd3[[ch]]]])
    rep_cells[[paste0("cdr3_", ch)]] <- aa
    rep_cells[[paste0("cdr3_", ch, "_len")]] <- nchar(aa)
  }
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

  list(chain_type = spec$chain_type, source = spec$source,
       segments = as.list(names(seg)), cdr3_chains = as.list(names(cd3)), group_col = gcol,
       has_antigen = !is.null(spec$antigen_col) && spec$antigen_col %in% names(md),
       has_tissue  = !is.null(spec$tissue_col)  && spec$tissue_col  %in% names(md),
       has_cluster = !is.null(div) && "cluster" %in% div$level_type,
       has_tissue_corr = !is.null(tcorr) && nrow(tcorr) > 0,
       n_clones = length(unique(clone)))
}

