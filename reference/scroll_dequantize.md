# Map stored (possibly quantized) values back to normalized expression

Map stored (possibly quantized) values back to normalized expression

## Usage

``` r
scroll_dequantize(values, manifest, assay)
```

## Arguments

- values:

  Numeric/integer vector of stored values.

- manifest:

  A manifest list (from
  [`scroll_manifest()`](https://shihanli92.github.io/scroll/reference/scroll_manifest.md)).

- assay:

  Assay the values came from.

## Value

Numeric vector in normalized units.
