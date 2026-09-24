# Query several features at once

Reads only the named features' partitions (partition pruning). Used by
aggregate views (e.g. dotplot) that draw a bounded marker panel; still
touches only those partitions, never the whole store.

## Usage

``` r
scroll_query_features(con, assay, features)
```

## Arguments

- con:

  A handle from
  [`scroll_connect()`](https://shihanli92.github.io/scroll/reference/scroll_connect.md).

- assay:

  Assay name.

- features:

  Character vector of feature names.

## Value

A data.frame with columns `feature`, `cell`, `value` (zero-valued cells
absent).
