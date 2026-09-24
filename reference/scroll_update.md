# Add or update reductions / metadata on a built project (no expr rebuild)

Rewrites only `cells.parquet` and `manifest.yaml`, aligning new columns
to the existing cells **by barcode** (cells absent from `object` get
`NA`, i.e. a partial-coverage subset embedding). The feature-partitioned
`expr/` and `counts/` stores are left byte-identical, so this is seconds
rather than a full
[`scroll_build()`](https://shihanli92.github.io/scroll/reference/scroll_build.md)
re-export.

## Usage

``` r
scroll_update(
  dir,
  object,
  embeddings = NULL,
  meta_cols = NULL,
  subsets = NULL,
  max_levels = .SCROLL_MAX_LEVELS,
  verbose = interactive()
)
```

## Arguments

- dir:

  A built scroll project directory (must contain `cells.parquet` and
  `manifest.yaml`).

- object:

  A Seurat object (or path to `.rds`) whose barcodes overlap the
  project. Only the requested reductions / metadata columns are read
  from it.

- embeddings:

  Character vector of reduction names in `object` to add as embedding
  coordinate columns (e.g. `c("umap_tcell")`). `NULL` adds none.

- meta_cols:

  Character vector of metadata columns in `object` to add / replace.
  `NULL` adds none.

- subsets:

  Optional named list of subset-view specs (as in
  [`scroll_build()`](https://shihanli92.github.io/scroll/reference/scroll_build.md));
  defaults to any spec attached by
  [`scroll_add_subset()`](https://shihanli92.github.io/scroll/reference/scroll_add_subset.md).
  Their embeddings / meta must be present after this update (existing or
  newly added).

- max_levels:

  Cap on cached categorical `levels:` per column, as in
  [`scroll_build()`](https://shihanli92.github.io/scroll/reference/scroll_build.md)
  (default 200). Columns above it record only `n_levels` and the app
  recomputes the level set at runtime.

- verbose:

  If `TRUE`, report what was written.

## Value

`dir`, invisibly.

## Details

The arguments mirror
[`scroll_build()`](https://shihanli92.github.io/scroll/reference/scroll_build.md),
so the *same* object you would rebuild from can instead be passed here
to append just the new pieces. Columns that already exist are
overwritten (e.g. to relabel a metadata column); new ones are added.
Subset views are registered exactly as in
[`scroll_build()`](https://shihanli92.github.io/scroll/reference/scroll_build.md)
/
[`scroll_add_subset()`](https://shihanli92.github.io/scroll/reference/scroll_add_subset.md).

## See also

[`scroll_build()`](https://shihanli92.github.io/scroll/reference/scroll_build.md),
[`scroll_add_subset()`](https://shihanli92.github.io/scroll/reference/scroll_add_subset.md)
