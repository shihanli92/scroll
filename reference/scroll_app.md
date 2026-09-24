# Build the scroll explorer app

Returns a
[`shiny::shinyApp`](https://rdrr.io/pkg/shiny/man/shinyApp.html) that
reads a built project directory. Deploy the scaffolded `app.R` (which
calls this) to a Shiny Server, or run it locally with
[`scroll_serve()`](https://shihanli92.github.io/scroll/reference/scroll_serve.md).

## Usage

``` r
scroll_app(dir = ".")
```

## Arguments

- dir:

  A built scroll project directory.

## Value

A `shiny.appobj`.
