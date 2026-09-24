# Preview a single panel while developing it

Mounts **one** panel in the full scroll harness — the app bar (with the
cell-subset filter and the subset-View selector) but none of the other
panels — so a custom panel can be driven for quick feedback without
building or scrolling through the whole app. Because the harness is the
real one, the panel is tested exactly as it will run: under the app-bar
filter and any subset views.

## Usage

``` r
scroll_preview_panel(panel, dir = ".")
```

## Arguments

- panel:

  The id of a panel registered with
  [`register_panel()`](https://shihanli92.github.io/scroll/reference/register_panel.md)
  /
  [`register_plot_panel()`](https://shihanli92.github.io/scroll/reference/register_plot_panel.md)
  (or a built-in id, e.g. `"dotplot"`).

- dir:

  A built scroll project directory to render the panel against.

## Value

A `shiny.appobj` — auto-prints/launches at the console, or pass to
[`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html).

## Details

The development loop is:

    source("panels/my_panel.R")               # (re-)registers the panel
    scroll_preview_panel("my_panel", "projects/demo")
    # edit the panel file -> re-source -> re-run

## See also

[`register_plot_panel()`](https://shihanli92.github.io/scroll/reference/register_plot_panel.md),
[`scroll_app()`](https://shihanli92.github.io/scroll/reference/scroll_app.md)
