# Embedding colored by a metadata column

Embedding colored by a metadata column

## Usage

``` r
view_umap_colorby(cells, params, state = list())
```

## Arguments

- cells:

  The globally-loaded cells data.frame.

- params:

  List with `embedding` and `color_by`.

- state:

  Optional toggle state (`embedding`, `split_by`, `show_labels`).

## Value

A ggplot.

## Deprecated

This plot core is internal to the built-in panels and will no longer be
exported from scroll 0.3.0.
