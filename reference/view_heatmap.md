# Expression heatmap across groups (or subsampled cells)

A tile heatmap of mean (or % expressing) expression per gene per group,
with optional per-gene z-scoring, row/column clustering + dendrograms,
and a diverging palette. A `cells` mode instead draws genes x a random
subsample of cells (capped, grouped) via `geom_raster`. `group_by` may
name \>1 column (its `a | b` interaction).

## Usage

``` r
view_heatmap(cells, params, expr_long, state = list(), assembly = NULL)
```

## Arguments

- cells:

  The cells data.frame (its group sizes are the denominators).

- params:

  List with `features`, `group_by`, and optionally `assay`.

- expr_long:

  A data.frame(feature, cell, value) in normalized units.

- state:

  Optional toggle state (`mode`, `stat`, `scale`, `clip`, `cluster`,
  `palette`, `cell_cap`, `aspect`, `theme`).

- assembly:

  Optional precomputed assembly (as the Shiny app supplies).

## Value

A ggplot, or an aplot composite when dendrograms are attached.

## Deprecated

This plot core is internal to the built-in panels and will no longer be
exported from scroll 0.3.0.
