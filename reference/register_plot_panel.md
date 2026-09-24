# Register a plot panel declaratively (no Shiny boilerplate)

A higher-level companion to
[`register_panel()`](https://shihanli92.github.io/scroll/reference/register_panel.md):
describe the controls and a single plot function, and scroll generates
the whole module — UI, control population, error surfacing, a Compute
gate, downloads, spinner, and app-bar subset/view awareness. The
silent-failure footguns of a hand-written server (a stray `req()`, an
uncaught error) are handled for you.

## Usage

``` r
register_plot_panel(
  id,
  plot,
  controls = list(),
  label = id,
  title = label,
  desc = NULL,
  after = NULL,
  before = NULL,
  compute = TRUE,
  csv = FALSE,
  style_caps = NULL
)
```

## Arguments

- id:

  Panel id (a single name).

- plot:

  `function(cells, input, data)` returning a ggplot. `cells` is the
  active (subset/view-filtered) cell table; `input` exposes each control
  by its id; `data` is the handle (`data$query1`, `data$con`,
  `data$manifest`, …).

- controls:

  A list of `scroll_input_*()` specs (see
  [scroll_input](https://shihanli92.github.io/scroll/reference/scroll_input.md)).
  Wrap a look-only control in
  [`scroll_style_input()`](https://shihanli92.github.io/scroll/reference/scroll_style_input.md)
  to place it in the plot's Style sheet instead of the control column.

- label, title, desc, after, before:

  As in
  [`register_panel()`](https://shihanli92.github.io/scroll/reference/register_panel.md).

- compute:

  If `TRUE` (default), the plot recomputes on a Compute button; `FALSE`
  makes it live (recompute on any control change).

- csv:

  If `TRUE`, add a CSV button that exports the plot's source data. The
  `plot` function opts in by attaching the table to its result, e.g.
  `attr(p, "scroll_source") <- df; p`. Default `FALSE`.

- style_caps:

  Scales & axes options for the plot's Style sheet, as in
  [`register_panel()`](https://shihanli92.github.io/scroll/reference/register_panel.md).

## Value

Invisibly, `id`.

## Examples

``` r
if (FALSE) { # \dontrun{
register_plot_panel("counts",
  controls = list(scroll_input_column("grp", "Group", "categorical")),
  plot = function(cells, input, data) {
    ggplot2::ggplot(cells, ggplot2::aes(.data[[input$grp]])) + ggplot2::geom_bar()
  },
  label = "Counts", compute = FALSE)
} # }
```
