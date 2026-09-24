# Show a control only when another control has a given value

Decorates a `scroll_input_*()` spec so the builder wraps it in a
[`shiny::conditionalPanel()`](https://rdrr.io/pkg/shiny/man/conditionalPanel.html)
— the control appears only when the control named `control` holds one of
`equals`. Enables declarative cascading (e.g. show a replicate picker
only when an engine choice is `"Pseudobulk"`).

## Usage

``` r
scroll_show_when(control_spec, control, equals)
```

## Arguments

- control_spec:

  A control spec from a `scroll_input_*()` constructor.

- control:

  Id of the control whose value gates visibility.

- equals:

  One or more values of `control` that reveal this control.

## Value

The decorated control spec.
