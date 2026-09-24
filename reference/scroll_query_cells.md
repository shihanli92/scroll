# Query all features for a set of cells

Returns the long (feature, cell, value) records for the given cells
across every feature — used to reconstruct a contrast's expression
matrix for live DE. Unlike a feature lookup this scans all partitions
(the cell filter does not prune the feature-partitioned store), so it is
the one query that touches the whole assay; call it on demand, not per
interaction.

## Usage

``` r
scroll_query_cells(con, assay, cells, dict = FALSE)
```

## Arguments

- con:

  A handle from
  [`scroll_connect()`](https://shihanli92.github.io/scroll/reference/scroll_connect.md).

- assay:

  Assay name.

- cells:

  Cell keys to fetch: the int32 global row index into `cells.parquet`
  for a v2 store (i.e. `cells$.gidx`), or barcodes for a v1 store.

- dict:

  If `TRUE`, return the `feature` column dictionary-encoded (an R
  factor) instead of a character vector. arrow builds the dictionary
  during the scan, so downstream
  [`unique()`](https://rdrr.io/r/base/unique.html)/[`match()`](https://rdrr.io/r/base/match.html)
  over the tens of millions of repeated gene names become cheap integer
  ops — a large speedup when reconstructing a DE contrast's matrix. The
  factor's levels are the store's own feature names.

## Value

A data.frame with columns `feature`, `cell`, `value`.

## Deprecated

This is an internal step of
[`scroll_de()`](https://shihanli92.github.io/scroll/reference/scroll_de.md)
and will no longer be exported from scroll 0.3.0.
