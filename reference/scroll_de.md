# Live differential expression for a contrast (presto / Wilcoxon)

Compares `ident1` against `ident2` or against the rest of the cells when
`ident2` is `NULL`/`"rest"`. Either side may name **several** levels of
`group_col` (they are pooled into one group). Reconstructs the
contrast's expression matrix from the store and runs
[`presto::wilcoxauc()`](https://rdrr.io/pkg/presto/man/wilcoxauc.html).

## Usage

``` r
scroll_de(
  data,
  assay,
  group_col,
  ident1,
  ident2 = NULL,
  min_pct = 0.1,
  cells = NULL,
  max_cells = NULL,
  progress = NULL
)
```

## Arguments

- data:

  A scroll data handle (`con`, `cells`, `manifest`).

- assay:

  Assay to test.

- group_col:

  Metadata column(s) defining the groups. Naming several columns
  compares their **interaction levels** (e.g.
  `c("genotype", "timepoint")` gives groups like `"KO | d7"`);
  `ident1`/`ident2` then select which combined levels form each side.

- ident1:

  One or more levels of `group_col` forming the group of interest.

- ident2:

  One or more comparison levels, or `NULL`/`"rest"` for one-vs-rest
  (every cell not in `ident1`).

- min_pct:

  Keep genes expressed in at least this fraction of either side.

- cells:

  Cells to test over (default `data$cells`); pass a subset to run DE
  within an active cell-subset filter.

- max_cells:

  Optional cap on cells **per group**. When set (a positive number),
  each side is randomly down-sampled to at most `max_cells` cells before
  testing, which bounds DE time and peak memory on large contrasts (a
  standard marker-detection shortcut; cf. Seurat's
  `max.cells.per.ident`). `NULL` (default) uses every cell. The
  subsample is deterministic (fixed seed) so repeated runs match, and
  the session RNG is left untouched.

- progress:

  Optional `function(fraction, detail)` called at the read and test
  stages, for a UI progress bar. `NULL` (default) is a no-op.

## Value

A data.frame of results, ranked by adjusted p-value, with columns
`gene`, `logFC`, `avg_log2FC`, `auc`, `pct.1`, `pct.2`, `p_val`,
`p_val_adj`. Two fold-change columns are reported because presto and
Seurat define it differently: **`logFC`** is presto's value – a
natural-log "mean of log" difference
(`mean(log-data in group1) - mean(...group2)`); **`avg_log2FC`** matches
[`Seurat::FindMarkers`](https://satijalab.org/seurat/reference/FindMarkers.html)
– the log2 of the mean of the *un-logged* normalized counts,
`log2((sum(expm1(x1))+1)/n1) - log2((sum(expm1(x2))+1)/n2)` (Seurat v5
`FoldChange`, pseudocount 1). For zero-inflated data `avg_log2FC` is
typically larger in magnitude than `logFC`.

## Details

With the default quantized build (`quantize = TRUE` in
[`scroll_build()`](https://shihanli92.github.io/scroll/reference/scroll_build.md)),
expression values below ~`max/510` round to zero. Such a floored value
is still stored (as an explicit zero), so `pct.1`/`pct.2` are
unaffected; instead those cells become tied at zero, which makes the
Wilcoxon p-values approximate, and `avg_log2FC` close to but not
bit-identical to Seurat's (max abs diff ~0.04, median ~0.005 in
practice; see `dev/validate_de_vs_findmarkers.R`). Build with
`quantize = FALSE` for exact statistics – `avg_log2FC`, `pct`, and the
p-value ranking then match
[`Seurat::FindMarkers`](https://satijalab.org/seurat/reference/FindMarkers.html)
to scroll's output rounding (Seurat's default `wilcox` test itself
dispatches to presto when it is installed, so the two run the same
Wilcoxon).
