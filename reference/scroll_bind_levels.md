# Populate a selectize control from a metadata column's levels

Wires the observe-and-update pattern custom panels need: whenever the
column named by `from` changes, the `id` selectize is repopulated with
that column's levels. `from` may be another control's input id (dynamic)
or a literal categorical column name (static, populated once).

## Usage

``` r
scroll_bind_levels(
  input,
  session,
  id,
  from,
  data,
  selected = 1,
  control_ids = character()
)
```

## Arguments

- input, session:

  The module's `input` and `session`.

- id:

  Input id of the selectize to populate (unnamespaced).

- from:

  A control id whose value is a column name, or a column name.

- data:

  The scroll data handle.

- selected:

  Which levels to preselect: an integer index/vector (e.g. `1`,
  `c(1, 2)`), or `"none"` for an empty default.

- control_ids:

  Ids of the sibling controls, used only to catch a `from` that names
  neither a metadata column nor another control (a common typo); the
  builder passes this automatically.

## Value

Invisibly, `NULL`.

## Deprecated

Superseded by
[`scroll_input_levels()`](https://shihanli92.github.io/scroll/reference/scroll_input.md)
(`from = `), which wires this for you; it will no longer be exported
from scroll 0.3.0.
