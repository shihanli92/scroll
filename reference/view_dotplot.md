# Dot plot: mean expression x fraction expressing, across groups

Dot plot: mean expression x fraction expressing, across groups

## Usage

``` r
view_dotplot(cells, params, expr_long, state = list(), assembly = NULL)
```

## Arguments

- cells:

  The cells data.frame (its group sizes are the denominators).

- params:

  List with `features`, `group_by`, and optionally `assay`.

- expr_long:

  A data.frame(feature, cell, value) in normalized units for the
  requested features (zero-valued cells absent).

- state:

  Optional toggle state (`split_by` overrides `group_by`).

- assembly:

  Optional precomputed result of the internal aggregation + clustering
  step; when supplied (as the Shiny app does), `cells`/`expr_long` are
  ignored for assembly and only cosmetics are applied.

## Value

A ggplot.

## Deprecated

This plot core is internal to the built-in panels and will no longer be
exported from scroll 0.3.0.
