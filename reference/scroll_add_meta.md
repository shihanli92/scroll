# Add (or replace) a metadata column in a built project – no Seurat object needed

Writes a new column into a built project's `cells.parquet` and registers
it in `manifest.yaml`, so it becomes selectable in the app on the next
restart (Color-by / group-by and DotPlot/Violin/Proportions for a
**categorical** column; FeaturePlot / Biaxial for a **numeric** one).
Unlike
[`scroll_update()`](https://shihanli92.github.io/scroll/reference/scroll_update.md)
this needs only the on-disk project, not the original Seurat object –
handy for editing a deployed app. The column must be something you can
compute from what the project already holds: other metadata columns, or
expression queried with
[`scroll_query_feature()`](https://shihanli92.github.io/scroll/reference/scroll_query_feature.md).
The expression store is untouched.

## Usage

``` r
scroll_add_meta(
  dir,
  name,
  values,
  scope = NULL,
  overwrite = FALSE,
  max_levels = .SCROLL_MAX_LEVELS
)
```

## Arguments

- dir:

  A built scroll project directory (with `cells.parquet` +
  `manifest.yaml`).

- name:

  Name of the column to add. Cannot be `"cell"` (the barcode column); an
  existing metadata column is replaced only when `overwrite = TRUE`.

- values:

  The column, one of: a vector of length `n_cells` in `cells.parquet`
  **row order**; a **named** vector keyed by cell barcode (unmatched
  cells become `NA`); a `data.frame` with a `cell` column and one value
  column; or a `function(cells)` returning one of those (it receives the
  loaded `cells` data.frame, so you can derive the column from existing
  columns). Character / factor / logical values register as categorical,
  numeric as numeric.

- scope:

  Optional name of a subset view (a key of the manifest `subsets`); when
  set, the column is scoped to that view – its levels are computed over
  the view's members only, so it surfaces only inside that View.

- overwrite:

  Replace an existing column of the same name (default `FALSE`).

- max_levels:

  Cap on cached categorical levels written to the manifest; above it the
  level list is omitted and recomputed at runtime (default 200).

## Value

`invisible(dir)`.

## See also

[`scroll_update()`](https://shihanli92.github.io/scroll/reference/scroll_update.md)
to add columns/embeddings from a Seurat object.
