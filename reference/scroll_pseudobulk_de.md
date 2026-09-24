# Pseudobulk differential expression (edgeR / limma-voom)

Aggregates cells into pseudobulk samples by summing raw counts, then
compares two groups with the limma-voom pipeline. Groups are built from
**combined-interaction levels**: the observed combinations of
`aggregate_cols` (e.g. `"sgIcos | 0"`); `ident1`/`ident2` select which
combinations form each group (either may name several). Pseudobulk
samples come from `replicate_col` when real replicates exist, otherwise
from random pseudo-replicates.

## Usage

``` r
scroll_pseudobulk_de(
  data,
  assay,
  aggregate_cols,
  ident1,
  ident2 = NULL,
  replicate_col = NULL,
  paired = c("auto", "yes", "no"),
  min_cells = 10,
  n_pseudo = 3,
  cells_per_pseudo = 50,
  cells = NULL,
  seed = 1L,
  progress = NULL,
  runs = 1L,
  lfc = 1,
  padj = 0.05
)
```

## Arguments

- data:

  A scroll data handle (`con`, `cells`, `manifest`).

- assay:

  Assay to test.

- aggregate_cols:

  Metadata columns whose combinations define the groups.

- ident1:

  Combined levels forming the group of interest.

- ident2:

  Combined levels for the comparison group, or `NULL`/`"rest"` for
  one-vs-rest (all cells not in `ident1`, restricted to observed
  combinations).

- replicate_col:

  Metadata column defining biological replicates. `NULL` or
  `"no_replicate"` pseudo-replicates the two groups directly (a disjoint
  partition of each group's cells) instead of using a real replicate
  column. Cells with a missing replicate value are excluded (and
  reported) rather than silently dropped; a named column that is absent
  is an error.

- paired:

  Whether to block on the replicate in the design when replicates are
  shared across both groups: `"auto"` (default) blocks when possible,
  `"yes"` forces it (falls back with a warning if not fittable), `"no"`
  always uses `~ 0 + group`.

- min_cells:

  Minimum cells for a pseudobulk sample to be kept (a per-sample floor,
  applied to real replicates and to each pseudo-partition bin).

- n_pseudo, cells_per_pseudo:

  Pseudo-replicate partition: number of pseudo-replicates and the
  maximum cells per pseudo-replicate.

- cells:

  Cells to test over (default `data$cells`).

- seed:

  Integer seed for the pseudo-replicate partition, so a compute is
  reproducible and the caller's RNG stream is left untouched. `NULL`
  uses the current RNG state (used internally to vary draws across
  stability runs). Real-replicate results are seed-independent (no
  random draw).

- progress:

  Optional `function(fraction, detail)` called at the aggregate and
  model-fit stages, for a UI progress bar. `NULL` (default) is a no-op.

- runs:

  Number of runs. `1` (default) is a single DE fit. `runs > 1` re-runs
  the fit with fresh random pseudo-replicates (seeds `1..runs`) and
  returns a **stability** table instead (see Value) – which genes come
  up consistently rather than by a lucky draw. This checks robustness to
  the random sampling only; it does not fix the anti-conservative bias
  of pseudo-replication. Real-replicate fits are seed-independent, so
  stability is only informative with pseudo-replicates.

- lfc, padj:

  Stability only (`runs > 1`): the logFC and adjusted-p cutoffs defining
  a "hit" in each run.

## Value

A data.frame (`gene`, `logFC`, `avg_expr`, `p_val`, `p_val_adj`) ranked
by adjusted p-value; positive `logFC` is up in `ident1`.
`logFC`/`avg_expr` are log2 (limma's `logFC`/`AveExpr`). Attributes:
`"pseudo"` (TRUE when any pseudo-replicates were used), `"regime"`
(`"real-paired"`, `"real-unpaired"`, `"mixed"`, or `"pseudo"`),
`"design"` (the model formula), `"n_samples"` (per group),
`"sample_sizes"` (cells per sample), and `"dropped_na"` (cells excluded
for a missing replicate value).

With `runs > 1`: a data.frame with `gene`, `sel_freq` (fraction of runs
in which the gene is a hit), `median_logFC`, `sign_agree`,
`median_padj`, `n_tested`, ranked by selection frequency (attribute
`"runs"`: successful runs).

## Details

With a real `replicate_col`, each replicate contributes **one sample per
group** (summing across the combinations it spans). The model uses the
no-intercept (means) parameterization: `~ 0 + group`, or
`~ 0 + group + replicate` when \>= 2 replicate levels are shared across
both groups (a **paired** design). The group1-vs-group2 effect is tested
as the `group1 - group2` contrast
([`limma::contrasts.fit`](https://rdrr.io/pkg/limma/man/contrasts.fit.html)),
so a positive `logFC` is up in `ident1`. A group with fewer than two
qualifying replicates falls back to pseudo-replicates for that group (a
"mixed" run).

Pseudo-replicates are a **disjoint partition** of a group's cells into
`n_pseudo` non-overlapping samples (never sharing a cell), each
requiring `min_cells` cells; the group must therefore hold at least
`n_pseudo * min_cells` cells or the call errors. Pseudo-replication
understates biological variance and is anti-conservative — prefer a real
replicate column.

Requires a counts store — build with
[`scroll_build()`](https://shihanli92.github.io/scroll/reference/scroll_build.md)
`counts = TRUE`.
