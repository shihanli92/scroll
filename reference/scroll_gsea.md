# Preranked GSEA on a differential-expression result

Ranks genes from a
[`scroll_de()`](https://shihanli92.github.io/scroll/reference/scroll_de.md)
or
[`scroll_pseudobulk_de()`](https://shihanli92.github.io/scroll/reference/scroll_pseudobulk_de.md)
result and runs
[`fgsea::fgsea()`](https://rdrr.io/pkg/fgsea/man/fgsea.html)
(multilevel) against a list of gene sets. Uses the result's unrounded
metrics when present, and breaks rank ties deterministically.

## Usage

``` r
scroll_gsea(
  de,
  gene_sets,
  rank_by = NULL,
  min_size = 15,
  max_size = 500,
  seed = 1
)
```

## Arguments

- de:

  A result from
  [`scroll_de()`](https://shihanli92.github.io/scroll/reference/scroll_de.md)
  or
  [`scroll_pseudobulk_de()`](https://shihanli92.github.io/scroll/reference/scroll_pseudobulk_de.md).

- gene_sets:

  A named list of gene vectors (e.g. from a `.gmt` file).

- rank_by:

  Ranking metric: `"t"` (Pseudobulk's limma t), `"avg_log2FC"`,
  `"median_logFC"`, `"logFC"`, `"auc"` (AUC - 0.5), or `"signed_logp"`
  (sign x -log10 p). `NULL` picks the first available of that order.

- min_size, max_size:

  Gene-set size limits, after keeping ranked genes only.

- seed:

  Seed for fgsea's sampling (the session's RNG is restored).

## Value

A data.frame (`pathway`, `size`, `ES`, `NES`, `pval`, `padj`,
`leading_edge`) ordered by `padj`, with attributes `stats` (the ranked
statistic) and `sets` (the gene sets tested).
