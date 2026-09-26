# Pseudobulk differential expression. Cells are aggregated into sample-level
# count profiles (summed raw counts, via arrow) and compared with the
# edgeR / limma-voom pipeline, so the statistics are replicate-based rather than
# per-cell. Needs a counts store (scroll_build(counts = TRUE)).

# Combined-interaction level for each cell (e.g. "sgIcos | 0"); NA if any of the
# aggregate columns is NA for that cell.
.scroll_combo_levels <- function(cells, cols) {
  parts <- lapply(cols, function(c) as.character(cells[[c]]))
  combo <- do.call(paste, c(parts, sep = " | "))
  combo[Reduce(`|`, lapply(parts, is.na))] <- NA
  combo
}

# Group label ("group1" / "group2" / NA) for each cell given its combined level.
# NA = not in the contrast (its combo is NA, or it is outside both idents in a
# two-group contrast). Shared by scroll_pseudobulk_de() and the panel's live
# preview so the preview matches what the compute aggregates. Empty ident1 => NA.
.scroll_combo_group <- function(combo, ident1, ident2 = NULL) {
  ident1 <- as.character(ident1); ident1 <- ident1[nzchar(ident1)]
  if (!length(ident1)) return(rep(NA_character_, length(combo)))
  ident2 <- as.character(ident2); ident2 <- ident2[nzchar(ident2)]
  in1 <- !is.na(combo) & combo %in% ident1
  if (!length(ident2) || "rest" %in% ident2)
    ifelse(in1, "group1", ifelse(!is.na(combo), "group2", NA_character_))
  else
    ifelse(in1, "group1", ifelse(!is.na(combo) & combo %in% ident2, "group2", NA_character_))
}

.rbind_or_empty <- function(out)
  if (length(out)) do.call(rbind, out) else
    data.frame(cell = character(), psample = character(), group = character(),
               rep = character())

# Split a cell vector into n_pseudo DISJOINT pseudo-replicates (a shuffled
# partition, not overlapping draws), so no cell is shared between samples and the
# within-group variance is real. Returns list() when the pool is too small to
# give every bin >= min_cells (N < n_pseudo * min_cells); the caller turns that
# into a specific error. Each bin is capped at cells_per_pseudo (never below
# min_cells). Seeding is handled by the caller via .scroll_with_seed().
.scroll_partition_reps <- function(cellv, tag, group, min_cells, n_pseudo, cells_per_pseudo) {
  n <- length(cellv)
  bin <- n %/% n_pseudo                                  # smallest balanced bin (before cap)
  per <- min(bin, cells_per_pseudo)                      # cells emitted per bin (cap = cells_per_pseudo)
  if (n_pseudo < 1L || per < min_cells) return(list())   # can't meet the per-sample floor
  idx <- sample.int(n)                                   # shuffle once (seeded upstream)
  sizes <- rep(bin, n_pseudo)                            # balanced bin sizes (differ by <=1)
  rem <- n %% n_pseudo
  if (rem) sizes[seq_len(rem)] <- sizes[seq_len(rem)] + 1L
  ends <- cumsum(sizes); starts <- ends - sizes + 1L
  lapply(seq_len(n_pseudo), function(k) {
    # disjoint bin, capped from above at cells_per_pseudo (min_cells is only a floor,
    # it never inflates the sample size).
    take <- utils::head(idx[starts[k]:ends[k]], cells_per_pseudo)
    data.frame(cell = cellv[take], psample = paste0(tag, " :: pseudo", k),
               group = group, rep = NA_character_, stringsAsFactors = FALSE)
  })
}

# Build the cell -> pseudobulk-sample map, tagged group1/group2, and classify the
# design regime.
# - replicate_col NULL / "no_replicate": pseudo-replicate the two GROUPS directly
#   (disjoint partition of each group's cells into n_pseudo samples).
# - a real replicate column: ONE sample per replicate level per GROUP (a donor is
#   one observation per group, summing across combos), keeping levels with
#   >= min_cells; a group with < 2 real replicates falls back to a pseudo
#   partition of that group (a "mixed" run).
# Regime: "real-paired" (>= 2 replicate levels shared across both groups ->
# ~ replicate + group is fittable), "real-unpaired" (real reps but not shared),
# "mixed" (one group real, one pseudo), "pseudo" (both pseudo). NA-replicate cells
# are counted (dropped_na) and excluded, not silently dropped.
# Returns list(mapping = data.frame(cell, psample, group, rep), regime, pseudo,
# dropped_na).
.scroll_pseudobulk_map <- function(df, replicate_col, min_cells, n_pseudo, cells_per_pseudo) {
  no_rep <- is.null(replicate_col) || !nzchar(replicate_col) ||
    identical(replicate_col, "no_replicate")

  part <- function(g) {                       # disjoint pseudo-partition for a group
    ng <- sum(df$grp == g)
    made <- .scroll_partition_reps(df$cell[df$grp == g], g, g,
                                   min_cells, n_pseudo, cells_per_pseudo)
    if (!length(made))                        # can't meet the per-sample floor
      stop(sprintf(paste0("Group %s: cannot build %d pseudo-replicates (needs >= %d cells ",
                          "each; group has %d cells, Cells/pseudo-rep = %d). Lower Pseudo-reps ",
                          "or Min cells, raise Cells/pseudo-rep, widen the contrast, or use a ",
                          "replicate column."),
                   sub("group", "", g), n_pseudo, min_cells, ng, cells_per_pseudo),
           call. = FALSE)
    made
  }

  if (no_rep) {                               # pseudo-replicate the two groups
    out <- c(part("group1"), part("group2"))
    return(list(mapping = .rbind_or_empty(out), regime = "pseudo",
                pseudo = TRUE, dropped_na = 0L))
  }

  dropped_na <- sum(is.na(df$rep))            # real replicates: report NA, don't hide
  df <- df[!is.na(df$rep), , drop = FALSE]

  out <- list(); real_reps <- list()          # collapse to one sample per rep per GROUP
  for (g in c("group1", "group2")) {
    sub <- df[df$grp == g, , drop = FALSE]
    reps <- split(sub$cell, sub$rep)
    reps <- reps[vapply(reps, length, integer(1)) >= min_cells]
    if (length(reps) >= 2) {                   # real replicates for this group
      real_reps[[g]] <- names(reps)
      out <- c(out, lapply(names(reps), function(r)
        data.frame(cell = reps[[r]], psample = paste(g, r, sep = " :: "),
                   group = g, rep = r, stringsAsFactors = FALSE)))
    } else {                                    # too few replicates -> pseudo this group
      out <- c(out, part(g))
    }
  }
  g1_real <- !is.null(real_reps[["group1"]]); g2_real <- !is.null(real_reps[["group2"]])
  regime <- if (g1_real && g2_real) {
    shared <- intersect(real_reps[["group1"]], real_reps[["group2"]])
    if (length(shared) >= 2) "real-paired" else "real-unpaired"
  } else if (g1_real || g2_real) "mixed" else "pseudo"
  list(mapping = .rbind_or_empty(out), regime = regime,
       pseudo = regime %in% c("pseudo", "mixed"), dropped_na = dropped_na)
}

# Build the pseudobulk sample table (psample, group, rep) + regime from metadata
# ONLY -- no counts query, no model fit -- so the panel preview can show exactly
# the samples (and design) the compute will build. `cell_key` is the store join
# key (int32 global index / barcode) scroll_aggregate_counts needs; a preview may
# pass any per-cell id since it uses only the sample structure. Returns the map
# result plus `samp` (unique psample/group/rep).
.scroll_pseudobulk_samples <- function(cells, cell_key, aggregate_cols, ident1, ident2,
                                       replicate_col, min_cells, n_pseudo, cells_per_pseudo,
                                       seed = 1L) {
  aggregate_cols <- as.character(aggregate_cols)
  ident1 <- as.character(ident1); ident1 <- ident1[nzchar(ident1)]
  ident2 <- as.character(ident2); ident2 <- ident2[nzchar(ident2)]
  if (!length(aggregate_cols)) stop("Pick at least one aggregate-by column.", call. = FALSE)
  if (!length(ident1)) stop("Pick at least one Group 1 level.", call. = FALSE)
  real_rep <- !is.null(replicate_col) && nzchar(replicate_col) &&
    !identical(replicate_col, "no_replicate")
  if (real_rep && !(replicate_col %in% names(cells)))
    stop("Replicate column '", replicate_col, "' not found in the cell metadata.",
         call. = FALSE)
  combo <- .scroll_combo_levels(cells, aggregate_cols)
  grp <- .scroll_combo_group(combo, ident1, ident2)          # "group1"/"group2"/NA per cell
  rep_vals <- if (real_rep) as.character(cells[[replicate_col]]) else NA_character_
  df <- data.frame(cell = cell_key, combo = combo, grp = grp, rep = rep_vals,
                   stringsAsFactors = FALSE)
  df <- df[!is.na(df$combo) & !is.na(df$grp), , drop = FALSE]
  if (!nrow(df)) stop("No cells match the selected groups.", call. = FALSE)
  # seed the pseudo-replicate partition (reproducible) without perturbing the
  # session RNG; stability passes a per-run seed to vary the draws.
  pm <- .scroll_with_seed(seed,
    .scroll_pseudobulk_map(df, replicate_col, min_cells, n_pseudo, cells_per_pseudo))
  pm$samp <- unique(pm$mapping[, c("psample", "group", "rep")])
  pm
}

# Build the limma design for a sample table using the no-intercept (means)
# parameterization: `~ 0 + group` (each group is its own column), blocking on the
# replicate (`~ 0 + group + replicate`) when the regime is paired and the design
# stays estimable. The group1-vs-group2 effect is a `contrast` over the columns
# (+1 group1, -1 group2), tested via limma::contrasts.fit -- so a positive logFC
# is up in ident1. Rows follow samp row order.
.scroll_pseudobulk_design <- function(samp, regime, paired = "auto") {
  group <- factor(samp$group, levels = c("group2", "group1"))
  use_paired <- switch(paired, auto = , yes = identical(regime, "real-paired"), no = FALSE)
  warn <- NULL
  if (identical(paired, "yes") && !identical(regime, "real-paired"))
    warn <- "paired = 'yes' but no replicate is shared across both groups; using ~ 0 + group."
  design <- stats::model.matrix(~ 0 + group); formula <- "~ 0 + group"
  if (use_paired) {
    replicate <- factor(samp$rep)
    cand <- stats::model.matrix(~ 0 + group + replicate)
    if (qr(cand)$rank == ncol(cand) && (nrow(cand) - ncol(cand)) >= 1) {
      design <- cand; formula <- "~ 0 + group + replicate"
    } else warn <- "Paired design is not estimable here; using ~ 0 + group."
  }
  # contrast: group1 - group2 (other terms 0). Columns are "groupgroup1"/"groupgroup2".
  contrast <- stats::setNames(numeric(ncol(design)), colnames(design))
  contrast["groupgroup1"] <- 1; contrast["groupgroup2"] <- -1
  list(design = design, formula = formula, contrast = contrast,
       group = group, warn = warn)
}

#' Pseudobulk differential expression (edgeR / limma-voom)
#'
#' Aggregates cells into pseudobulk samples by summing raw counts, then compares
#' two groups with the limma-voom pipeline. Groups are built from
#' **combined-interaction levels**: the observed combinations of `aggregate_cols`
#' (e.g. `"sgIcos | 0"`); `ident1`/`ident2` select which combinations form each
#' group (either may name several). Pseudobulk samples come from `replicate_col`
#' when real replicates exist, otherwise from random pseudo-replicates.
#'
#' With a real `replicate_col`, each replicate contributes **one sample per group**
#' (summing across the combinations it spans). The model uses the no-intercept
#' (means) parameterization: `~ 0 + group`, or `~ 0 + group + replicate` when >= 2
#' replicate levels are shared across both groups (a **paired** design). The
#' group1-vs-group2 effect is tested as the `group1 - group2` contrast
#' (`limma::contrasts.fit`), so a positive `logFC` is up in `ident1`. A group with
#' fewer than two qualifying replicates falls back to pseudo-replicates for that
#' group (a "mixed" run).
#'
#' Pseudo-replicates are a **disjoint partition** of a group's cells into
#' `n_pseudo` non-overlapping samples (never sharing a cell), each requiring
#' `min_cells` cells; the group must therefore hold at least
#' `n_pseudo * min_cells` cells or the call errors. Pseudo-replication understates
#' biological variance and is anti-conservative — prefer a real replicate column.
#'
#' Requires a counts store — build with [scroll_build()] `counts = TRUE`.
#'
#' @param data A scroll data handle (`con`, `cells`, `manifest`).
#' @param assay Assay to test.
#' @param aggregate_cols Metadata columns whose combinations define the groups.
#' @param ident1 Combined levels forming the group of interest.
#' @param ident2 Combined levels for the comparison group, or `NULL`/`"rest"` for
#'   one-vs-rest (all cells not in `ident1`, restricted to observed combinations).
#' @param replicate_col Metadata column defining biological replicates. `NULL` or
#'   `"no_replicate"` pseudo-replicates the two groups directly (a disjoint
#'   partition of each group's cells) instead of using a real replicate column.
#'   Cells with a missing replicate value are excluded (and reported) rather than
#'   silently dropped; a named column that is absent is an error.
#' @param paired Whether to block on the replicate in the design when replicates
#'   are shared across both groups: `"auto"` (default) blocks when possible,
#'   `"yes"` forces it (falls back with a warning if not fittable), `"no"` always
#'   uses `~ 0 + group`.
#' @param min_cells Minimum cells for a pseudobulk sample to be kept (a per-sample
#'   floor, applied to real replicates and to each pseudo-partition bin).
#' @param n_pseudo,cells_per_pseudo Pseudo-replicate partition: number of
#'   pseudo-replicates and the maximum cells per pseudo-replicate.
#' @param cells Cells to test over (default `data$cells`).
#' @param seed Integer seed for the pseudo-replicate partition, so a compute is
#'   reproducible and the caller's RNG stream is left untouched. `NULL` uses the
#'   current RNG state (used internally to vary draws across stability runs).
#'   Real-replicate results are seed-independent (no random draw).
#' @param progress Optional `function(fraction, detail)` called at the aggregate
#'   and model-fit stages, for a UI progress bar. `NULL` (default) is a no-op.
#' @param runs Number of runs. `1` (default) is a single DE fit. `runs > 1`
#'   re-runs the fit with fresh random pseudo-replicates (seeds `1..runs`) and
#'   returns a **stability** table instead (see Value) -- which genes come up
#'   consistently rather than by a lucky draw. This checks robustness to the
#'   random sampling only; it does not fix the anti-conservative bias of
#'   pseudo-replication. Real-replicate fits are seed-independent, so stability is
#'   only informative with pseudo-replicates.
#' @param lfc,padj Stability only (`runs > 1`): the logFC and adjusted-p cutoffs
#'   defining a "hit" in each run.
#' @return A data.frame (`gene`, `logFC`, `avg_expr`, `p_val`, `p_val_adj`) ranked
#'   by adjusted p-value; positive `logFC` is up in `ident1`. `logFC`/`avg_expr`
#'   are log2 (limma's `logFC`/`AveExpr`). Attributes: `"pseudo"` (TRUE when any
#'   pseudo-replicates were used), `"regime"` (`"real-paired"`, `"real-unpaired"`,
#'   `"mixed"`, or `"pseudo"`), `"design"` (the model formula), `"n_samples"`
#'   (per group), `"sample_sizes"` (cells per sample), and `"dropped_na"` (cells
#'   excluded for a missing replicate value).
#'
#'   With `runs > 1`: a data.frame with `gene`, `sel_freq` (fraction of runs in
#'   which the gene is a hit), `median_logFC`, `sign_agree`, `median_padj`,
#'   `n_tested`, ranked by selection frequency (attribute `"runs"`: successful runs).
#' @export
scroll_pseudobulk_de <- function(data, assay, aggregate_cols, ident1, ident2 = NULL,
                                 replicate_col = NULL, paired = c("auto", "yes", "no"),
                                 min_cells = 10,
                                 n_pseudo = 3, cells_per_pseudo = 50, cells = NULL,
                                 seed = 1L, progress = NULL,
                                 runs = 1L, lfc = 1, padj = 0.05) {
  paired <- match.arg(paired)
  if (runs > 1) {
    runs_list <- .scroll_pseudobulk_runs(
      data, assay, aggregate_cols, ident1, ident2, replicate_col = replicate_col,
      paired = paired, min_cells = min_cells, n_pseudo = n_pseudo,
      cells_per_pseudo = cells_per_pseudo, cells = cells, runs = runs)
    return(.scroll_stability_aggregate(runs_list, lfc, padj))
  }
  pr <- function(f, d) if (is.function(progress)) progress(f, d)
  if (!requireNamespace("edgeR", quietly = TRUE) ||
      !requireNamespace("limma", quietly = TRUE))
    stop("Pseudobulk DE needs the 'edgeR' and 'limma' packages.", call. = FALSE)
  if (!isTRUE(data$manifest$has_counts))
    stop("This project has no counts store; rebuild with scroll_build(counts = TRUE).",
         call. = FALSE)
  cells <- cells %||% data$cells
  # `cell` here is the store's cell key that scroll_aggregate_counts joins on:
  # the int32 global index for a v2 store, the barcode string for v1.
  cell_key <- if (isTRUE(data$manifest$cell_index))
                (if (!is.null(cells$.gidx)) cells$.gidx else seq_len(nrow(cells))) else cells$cell
  sm <- .scroll_pseudobulk_samples(cells, cell_key, aggregate_cols, ident1, ident2,
                                   replicate_col, min_cells, n_pseudo, cells_per_pseudo, seed)
  mapping <- sm$mapping; samp <- sm$samp
  if (sm$dropped_na > 0)
    warning(sm$dropped_na, " cells with a missing '", replicate_col,
            "' value were excluded.", call. = FALSE)
  if (sum(samp$group == "group1") < 2 || sum(samp$group == "group2") < 2)
    stop("Each group needs >= 2 pseudobulk samples after the min-cell filter. ",
         "Lower Min cells, choose a replicate column, or add pseudo-replicates.",
         call. = FALSE)

  pr(0.2, "Aggregating pseudobulk counts...")
  agg <- scroll_aggregate_counts(data$con, assay, mapping[, c("cell", "psample")])
  # store's own feature names (encoding-robust; the manifest yaml can mangle a
  # non-ASCII feature name into an escaped form that would not match() the store).
  feats <- sort(unique(agg$feature))
  ps <- samp$psample
  M <- matrix(0L, nrow = length(feats), ncol = length(ps), dimnames = list(feats, ps))
  M[cbind(match(agg$feature, feats), match(agg$psample, ps))] <- as.integer(agg$count)

  pr(0.7, "Fitting model (edgeR / limma-voom)...")
  # no-intercept (means) design; the group1 - group2 effect is a contrast tested
  # via contrasts.fit, so positive logFC is up in ident1. design rows follow samp
  # order = the M columns (ps).
  ds <- .scroll_pseudobulk_design(samp, sm$regime, paired)
  if (!is.null(ds$warn)) warning(ds$warn, call. = FALSE)
  group <- ds$group; design <- ds$design; design_formula <- ds$formula

  dge <- edgeR::DGEList(counts = M, group = group)
  keep <- edgeR::filterByExpr(dge, group = group)
  dge <- edgeR::calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
  fit <- limma::lmFit(limma::voom(dge, design), design)
  fit <- limma::eBayes(limma::contrasts.fit(fit, ds$contrast))
  tt <- limma::topTable(fit, coef = 1, number = Inf, sort.by = "P")

  res <- data.frame(
    gene = rownames(tt),
    logFC = round(tt$logFC, 3),
    avg_expr = round(tt$AveExpr, 3),
    p_val = signif(tt$P.Value, 3),
    p_val_adj = signif(tt$adj.P.Val, 3),
    row.names = NULL, stringsAsFactors = FALSE)
  attr(res, "pseudo") <- sm$pseudo
  attr(res, "regime") <- sm$regime
  attr(res, "design") <- design_formula
  attr(res, "dropped_na") <- sm$dropped_na
  attr(res, "n_samples") <- c(group1 = sum(samp$group == "group1"),
                              group2 = sum(samp$group == "group2"))
  attr(res, "sample_sizes") <- tapply(mapping$cell, mapping$psample, length)
  res
}

# Run scroll_pseudobulk_de `runs` times (fresh pseudo-replicates each) and return
# the successful per-run result data.frames. A distinct per-run seed keeps the
# draws varied AND the whole stability computation reproducible, and (like a
# single compute) leaves the caller's RNG stream untouched.
.scroll_pseudobulk_runs <- function(..., runs) {
  out <- lapply(seq_len(runs),
                function(i) tryCatch(scroll_pseudobulk_de(..., seed = i),
                                     error = function(e) NULL))
  out <- Filter(Negate(is.null), out)
  if (length(out) < 2)
    stop("Stability needs >= 2 successful runs (", length(out), " succeeded); ",
         "check the contrast and min-cell settings.", call. = FALSE)
  out
}

# Summarise a list of per-run DE results into per-gene stability statistics: how
# often each gene clears the lfc + padj cutoffs (selection frequency), its median
# effect, and how consistently the effect points the same way.
.scroll_stability_aggregate <- function(runs, lfc = 1, padj = 0.05) {
  k <- length(runs)
  all <- do.call(rbind, runs)
  sig <- abs(all$logFC) >= lfc & all$p_val_adj <= padj
  ix <- split(seq_len(nrow(all)), all$gene)
  lf <- all$logFC; pj <- all$p_val_adj
  res <- data.frame(
    gene = names(ix),
    sel_freq = round(vapply(ix, function(j) sum(sig[j]) / k, numeric(1)), 3),
    median_logFC = round(vapply(ix, function(j) stats::median(lf[j]), numeric(1)), 3),
    sign_agree = round(vapply(ix, function(j) {
      s <- sign(lf[j]); max(mean(s > 0), mean(s < 0)) }, numeric(1)), 2),
    median_padj = signif(vapply(ix, function(j) stats::median(pj[j]), numeric(1)), 3),
    n_tested = vapply(ix, length, integer(1)),
    row.names = NULL, stringsAsFactors = FALSE)
  res <- res[order(-res$sel_freq, res$median_padj), , drop = FALSE]
  attr(res, "runs") <- k
  # aliases so view_volcano / the shared table still work
  res$logFC <- res$median_logFC; res$p_val_adj <- res$median_padj
  res
}

