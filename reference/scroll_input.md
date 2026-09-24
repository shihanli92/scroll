# Declarative controls for [`register_plot_panel()`](https://shihanli92.github.io/scroll/reference/register_plot_panel.md)

Each constructor returns a control spec. `scroll_input_column` picks a
metadata column of a type; `scroll_input_levels` is a level selector
bound to the column chosen by another control (`from`) or a fixed
column; `scroll_input_gene` is a server-side gene search over the
default assay; `scroll_input_assay` / `scroll_input_embedding` pick an
assay / reduction from the manifest (the embedding narrows to the active
subset view when `view_aware`); `scroll_input_custom` wraps your own
`ui`/`bind` as a control; the rest are thin wrappers over the matching
Shiny inputs. Wrap any of them in
[`scroll_show_when()`](https://shihanli92.github.io/scroll/reference/scroll_show_when.md)
for conditional visibility.

## Usage

``` r
scroll_input_column(
  id,
  label,
  type = "categorical",
  selected = 1,
  view_aware = FALSE,
  prefer = NULL
)

scroll_input_levels(
  id,
  label,
  from = NULL,
  choices = NULL,
  watch = NULL,
  multiple = TRUE,
  required = FALSE,
  none = FALSE,
  selected = 1
)

scroll_input_gene(id, label, multiple = FALSE, required = FALSE)

scroll_input_numeric(id, label, value, min = NA, max = NA, step = NA)

scroll_input_slider(id, label, min, max, value, step = NULL)

scroll_input_choice(
  id,
  label,
  choices,
  selected = NULL,
  inline = TRUE,
  multiple = FALSE,
  widget = c("auto", "radio", "select")
)

scroll_input_text(id, label, placeholder = NULL)

scroll_input_palette(
  id,
  label = "Palette",
  type = c("discrete", "continuous"),
  selected = NULL
)

scroll_input_assay(id, label = "Assay")

scroll_input_embedding(id, label = "Embedding", view_aware = TRUE)

scroll_input_custom(id, label, ui, bind = NULL, needs = NULL, required = FALSE)
```

## Arguments

- id:

  Control id (also the name under which its value reaches `plot`'s
  `input`).

- label:

  Control label.

- type:

  For `scroll_input_column`: `"categorical"`, `"numeric"`, or `"any"`;
  for `scroll_input_palette`: `"discrete"` or `"continuous"`.

- selected:

  Default selection: integer index/vector, or `"none"`. For
  `scroll_input_palette`, a palette name.

- view_aware:

  For `scroll_input_column`/`scroll_input_embedding`: repopulate choices
  from the active subset view (requires the panel to be view-aware).

- prefer:

  For `scroll_input_column`: one or more column names to default-select
  when present (first match wins), falling back to `selected` otherwise.

- from:

  For `scroll_input_levels`: a column-picking control's id, or a
  categorical column name. Supply either `from` or `choices`.

- choices:

  For most inputs a static vector; `scroll_input_choice` also accepts a
  `function(data)` and `scroll_input_levels` a `function(input, data)`,
  both evaluated to derive choices from the data handle (e.g. a baked
  asset under `data$dir`).

- watch:

  For `scroll_input_levels` with a `choices` function: control ids whose
  change repopulates the levels (empty = populate once at start-up).

- multiple:

  Allow multiple selections.

- required:

  If `TRUE`, an empty value shows a "Select ." message instead of
  computing.

- none:

  If `TRUE`, the level selector defaults to empty (e.g. a "vs. rest").

- value, min, max, step, placeholder, inline:

  Passed to the underlying Shiny input.

- widget:

  For `scroll_input_choice`: `"radio"`, `"select"`, or `"auto"` (radio
  for a short static list, a selectize dropdown otherwise).

- ui:

  For `scroll_input_custom`: a `function(ns, data)` returning the
  control's UI (namespace inputs with `ns()`).

- bind:

  For `scroll_input_custom`: an optional
  `function(input, session, data)` wiring server-side behaviour.

- needs:

  For `scroll_input_custom`: a required column type
  (`"categorical"`/`"numeric"`) that gates the panel's empty state, or
  `NULL`.

## Value

A control spec (a list) for `register_plot_panel(controls = )`.
