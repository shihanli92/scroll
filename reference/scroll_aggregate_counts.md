# Aggregate raw counts into pseudobulk samples

Sums an assay's raw counts (the `counts/<assay>.parquet` store, written
by
[`scroll_build()`](https://shihanli92.github.io/scroll/reference/scroll_build.md)
with `counts = TRUE`) over a cell-to-sample mapping. The join and
summation run in arrow's engine, so the full genes x cells matrix is
never materialized in R. Used by
[`scroll_pseudobulk_de()`](https://shihanli92.github.io/scroll/reference/scroll_pseudobulk_de.md).

## Usage

``` r
scroll_aggregate_counts(con, assay, mapping)
```

## Arguments

- con:

  A handle from
  [`scroll_connect()`](https://shihanli92.github.io/scroll/reference/scroll_connect.md).

- assay:

  Assay name.

- mapping:

  A data.frame with columns `cell` and `psample` (a cell may map to
  several pseudobulk samples, e.g. overlapping pseudo-replicates).

## Value

A data.frame with columns `feature`, `psample`, `count` (summed).

## Deprecated

This is an internal step of
[`scroll_pseudobulk_de()`](https://shihanli92.github.io/scroll/reference/scroll_pseudobulk_de.md)
and will no longer be exported from scroll 0.3.0.
