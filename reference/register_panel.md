# Register a custom panel in the scroll explorer

Adds a user-defined analysis panel to every
[`scroll_app()`](https://shihanli92.github.io/scroll/reference/scroll_app.md)
built afterwards in this session. A panel follows the same module
contract as the built-ins — a pair of functions:

## Usage

``` r
register_panel(
  id,
  ui,
  server,
  label = id,
  title = label,
  desc = NULL,
  after = NULL,
  before = NULL,
  style_ui = NULL,
  style_caps = NULL
)
```

## Arguments

- id:

  Unique panel id: a letter followed by letters, digits, or underscores.
  Used as the panel anchor and the Shiny module namespace.

- ui, server:

  The panel's UI and server functions (see contract above).

- label:

  Short rail label (defaults to `id`).

- title, desc:

  Section heading and one-line description.

- after:

  Id of the panel to insert this one after; `NULL` (default) appends at
  the end. Ignored when replacing an existing id.

- before:

  Id of the panel to insert this one before (e.g. `"dimplot"` to put it
  at the top). Takes precedence over `after`. Ignored when replacing an
  existing id.

- style_ui:

  Optional `function(id, data)` returning look-only controls to show in
  the plot's **Style sheet** (the paintbrush button in a
  `.scroll_plot_area()` toolbar); namespace them with `NS(id)` like
  `ui`, and the server reads them as usual. The sheet also carries
  per-plot titles/labels and theme; a server receives them by declaring
  a `style_r` argument (the whole style) and/or `theme_r` (just the
  theme family).

- style_caps:

  Which "Scales & axes" options the Style sheet offers for this plot:
  `list(limits = c("x","y"), trans = list(x = , y = ), breaks = , flip = , facet = )`.
  `NULL` (default) offers visible-range limits on both axes;
  [`list()`](https://rdrr.io/r/base/list.html) offers none. Transforms
  are any of `"log10"`, `"sqrt"`, `"log1p"`, `"pseudo_log"`,
  `"reverse"`.

## Value

Invisibly, `id`.

## Details

- `ui(id, data)` returns a UI tag; namespace its inputs with
  [`shiny::NS()`](https://rdrr.io/pkg/shiny/man/NS.html)`(id)`.

- `server(id, data, cells_r)` wires the module, typically via
  [`shiny::moduleServer()`](https://rdrr.io/pkg/shiny/man/moduleServer.html).
  `cells_r()` is a reactive of the *active* cells (already narrowed by
  the app-bar subset filter); use it, not `data$cells`, so the panel
  honours the global filter.

`data` is the shared handle: `data$cells` (a data.frame of metadata +
embedding coordinates), `data$manifest`, `data$config`, and the query
helpers `data$query1(assay, feature)` / `data$queryN(assay, features)`,
which return dequantized long expression (`cell`, `value` / `feature`,
`cell`, `value`).

Registering an `id` that matches an existing panel — built-in or custom
— replaces it in place, so you can override a built-in. Clear all custom
panels with
[`scroll_reset_panels()`](https://shihanli92.github.io/scroll/reference/scroll_reset_panels.md).

## Examples

``` r
# A minimal custom panel: a cells-per-group bar chart.
count_ui <- function(id, data) {
  ns <- shiny::NS(id)
  cols <- names(Filter(function(x) identical(x$type, "categorical"),
                       data$manifest$meta))
  shiny::tagList(shiny::selectInput(ns("grp"), "Group", cols),
                 shiny::plotOutput(ns("plot")))
}
count_server <- function(id, data, cells_r = shiny::reactive(data$cells)) {
  shiny::moduleServer(id, function(input, output, session) {
    output$plot <- shiny::renderPlot({
      shiny::req(input$grp)
      barplot(table(cells_r()[[input$grp]]))
    })
  })
}
register_panel("counts", count_ui, count_server, label = "Counts",
               title = "Cells per group")
scroll_reset_panels()   # (undo, so the example leaves no state)
```
