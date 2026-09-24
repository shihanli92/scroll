# Violin: a feature's per-group distribution

Violin: a feature's per-group distribution

## Usage

``` r
view_violin(cells, params, values = NULL, state = list())
```

## Arguments

- cells:

  The cells data.frame.

- params:

  List with `group_by`, and a value source: either `feature` (+ optional
  `assay`, queried into `values`) or `value_col` (a numeric metadata
  column plotted directly).

- values:

  A data.frame(cell, value) for the feature, or `NULL`.

- state:

  Optional toggle state (`palette`, `jitter`, `legend`; `split_by`
  overrides `group_by`).

## Value

A ggplot.

## Deprecated

This plot core is internal to the built-in panels and will no longer be
exported from scroll 0.3.0.
