# Stability of pseudobulk DE across random pseudo-replicate draws

Deprecated: use `scroll_pseudobulk_de(..., runs = )`, which returns the
same stability table when `runs > 1`. This wrapper will be removed in
scroll 0.3.0.

## Usage

``` r
scroll_pseudobulk_stability(
  data,
  assay,
  aggregate_cols,
  ident1,
  ident2 = NULL,
  replicate_col = "no_replicate",
  min_cells = 10,
  n_pseudo = 3,
  cells_per_pseudo = 50,
  cells = NULL,
  runs = 25,
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

- min_cells:

  Minimum cells for a pseudobulk sample to be kept (a per-sample floor,
  applied to real replicates and to each pseudo-partition bin).

- n_pseudo, cells_per_pseudo:

  Pseudo-replicate partition: number of pseudo-replicates and the
  maximum cells per pseudo-replicate.

- cells:

  Cells to test over (default `data$cells`).

- runs:

  Number of random re-runs.

- lfc, padj:

  Stability only (`runs > 1`): the logFC and adjusted-p cutoffs defining
  a "hit" in each run.

## Value

See
[`scroll_pseudobulk_de()`](https://shihanli92.github.io/scroll/reference/scroll_pseudobulk_de.md)
(`runs > 1`).
