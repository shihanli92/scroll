# Reusable plotting helpers for custom views

The pieces the built-in panels use, exposed so a custom `plot` (or a raw
[`register_panel()`](https://shihanli92.github.io/scroll/reference/register_panel.md)
view) matches their look and performance:

## Usage

``` r
scroll_point_layer(mapping = NULL, size = 0.6, alpha = 1, raster = FALSE, ...)

scroll_discrete_colors(values, palette = "Tableau 10")

scroll_continuous_scale(palette = "viridis", name = NULL, limits = NULL)

scroll_group_labels(df, group, x, y)
```

## Arguments

- mapping, size, alpha, raster, ...:

  For `scroll_point_layer()`: a ggplot2 aesthetic mapping, point size,
  opacity, whether to rasterize on screen, and extra args passed to the
  geom.

- values:

  A vector of categorical values (its unique levels are coloured).

- palette:

  A palette name (discrete: e.g. `"Tableau 10"`; continuous: e.g.
  `"viridis"`, `"magma"`, a ColorBrewer name, or `"grey-red"`).

- name, limits:

  For `scroll_continuous_scale()`: legend title and optional
  `c(low, high)` clip (out-of-range values are squished, not dropped).

- df:

  A data.frame with the scatter coordinates and group column.

- group, x, y:

  For `scroll_group_labels()`: the grouping and coordinate column names
  in `df`.

## Value

A ggplot2 layer/scale (or, for `scroll_discrete_colors`, a named
vector).

## Details

- `scroll_point_layer()` — a scatter layer that rasterizes on screen for
  speed (via scattermore, when installed and `raster = TRUE`) but stays
  a true `geom_point` otherwise, so exports and small plots are
  unchanged.

- `scroll_discrete_colors()` — a named colour vector for a categorical
  column, keyed by level so a group keeps its colour across panels.

- `scroll_continuous_scale()` — a ggplot2 continuous colour scale
  (viridis, ColorBrewer, or grey-to-hue) matching FeaturePlot.

- `scroll_group_labels()` — per-group median-centroid text labels for a
  scatter.
