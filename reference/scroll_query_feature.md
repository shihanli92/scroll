# Query one feature's expression from the Parquet store

Reads only the `feature=<feature>` partition of the assay subtree.
Values are returned as stored (quantized `uint8` when the project was
built with `quantize = TRUE`); use
[`scroll_dequantize()`](https://shihanli92.github.io/scroll/reference/scroll_dequantize.md)
to map back to normalized units.

## Usage

``` r
scroll_query_feature(con, assay, feature)
```

## Arguments

- con:

  A handle from
  [`scroll_connect()`](https://shihanli92.github.io/scroll/reference/scroll_connect.md).

- assay:

  Assay name (subtree under `expr/`).

- feature:

  Feature name to look up.

## Value

A data.frame with columns `cell` and `value`. Zero-valued cells are
absent (the matrix is sparse); callers treat missing cells as 0.
