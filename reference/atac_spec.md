# Describe a dataset's scATAC peaks assay for `scroll_build()`

Marks which exported assay holds chromatin-accessibility **peaks** and
bakes a peak-annotation table. Pass the result as
`scroll_build(..., atac = atac_spec())`. The peaks assay must still be
listed in `assays=` so its accessibility matrix is exported like any
other assay; `atac_spec()` additionally records genomic coordinates
(parsed from the peak names, e.g. `chr1-9776-10668`) and, when a Signac
`ChromatinAssay` annotation is present, each peak's nearest gene.

## Usage

``` r
atac_spec(assay = "peaks", nearest_gene = TRUE)
```

## Arguments

- assay:

  Name of the peaks assay (default `"peaks"`). Must be one of the
  exported `assays`.

- nearest_gene:

  Whether to compute each peak's nearest gene via
  [`Signac::ClosestFeature()`](https://stuartlab.org/signac/reference/ClosestFeature.html)
  when the assay carries an annotation (default `TRUE`). Falls back to
  coordinates-only if unavailable.

## Value

A `scroll_atac_spec` list.

## See also

[`scroll_build()`](https://shihanli92.github.io/scroll/reference/scroll_build.md)
