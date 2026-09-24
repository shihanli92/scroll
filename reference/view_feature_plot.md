# Embedding colored by a feature's expression

Embedding colored by a feature's expression

## Usage

``` r
view_feature_plot(cells, params, values = NULL, state = list())
```

## Arguments

- cells:

  The globally-loaded cells data.frame.

- params:

  List with `embedding` (and `feature`, for the title).

- values:

  A data.frame(cell, value) in normalized units, or `NULL`.

- state:

  Optional toggle state.

## Value

A ggplot.

## Deprecated

This plot core is internal to the built-in panels and will no longer be
exported from scroll 0.3.0.
