# Put a control in the plot's Style sheet

Decorates a `scroll_input_*()` spec so the builder renders it in the
panel's **Style sheet** (opened from the paintbrush button in the plot
toolbar) instead of the control column. Use it for "how the plot looks"
options – palettes, point sizes, aspect ratio – and keep "what to plot"
in the column. The control keeps its id, so the `plot` function reads it
from `input` exactly as before, and it can be combined with
[`scroll_show_when()`](https://shihanli92.github.io/scroll/reference/scroll_show_when.md).

## Usage

``` r
scroll_style_input(control_spec)
```

## Arguments

- control_spec:

  A control spec from a `scroll_input_*()` constructor.

## Value

The decorated control spec.

## Examples

``` r
if (FALSE) { # \dontrun{
register_plot_panel("my_panel",
  controls = list(
    scroll_input_column("group", "Group by", type = "categorical"),
    scroll_style_input(scroll_input_slider("size", "Point size", 0.1, 3, 1, 0.1))),
  plot = function(cells, input, data) ...)
} # }
```
