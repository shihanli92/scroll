# Wire a plot output with error surfacing and PNG/PDF downloads

One call replaces the usual `reactive`/`eventReactive` + `tryCatch` +
`renderPlot` + download-handler boilerplate. `fun` is evaluated inside a
`tryCatch`; any error (or a non-plot return) is shown **inline** in the
plot area rather than failing silently, and PNG/PDF handlers
(`output$png` / `output$pdf`, matching `scroll_render_plot()`'s
companion UI `.scroll_plot_area()`) are registered for the result.

## Usage

``` r
scroll_render_plot(
  output,
  id,
  fun,
  event = NULL,
  placeholder = NULL,
  lazy = FALSE,
  input = NULL
)
```

## Arguments

- output:

  The module's `output`.

- id:

  The panel id (used for download filenames).

- fun:

  A zero-argument function returning a ggplot (may
  [`stop()`](https://rdrr.io/r/base/stop.html) on bad input — the
  message is surfaced).

- event:

  `NULL` for a live plot, or a reactive gating recomputation (e.g.
  `shiny::reactive(input$go)` for a Compute button). When it is an
  action-button count, the plot waits for the first click.

- placeholder:

  Optional message shown before the first `event` fires.

- lazy:

  When `TRUE` (and `event` is `NULL`), gate a live plot on the panel's
  on-screen state, so a scrolled-away panel reuses its last render
  instead of redrawing on a global filter / View change. Requires
  `input`.

- input:

  The module's `input` (only needed when `lazy = TRUE`).

## Value

Invisibly, the plot reactive.
