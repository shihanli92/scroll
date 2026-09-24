# Attach a reprocessed subset onto a parent object

Convenience assembly for
[`scroll_build()`](https://shihanli92.github.io/scroll/reference/scroll_build.md)'s
`subsets`: copies one or more dimensional reductions (and optional
metadata columns) from a separately reprocessed child object onto
`object`, aligned **by cell barcode**, and records the subset spec so a
later `scroll_build(object, ...)` exposes it as a linked **View**. This
is plumbing only — it performs no normalization, clustering, or
embedding; the child object must already be reprocessed.

## Usage

``` r
scroll_add_subset(
  object,
  sub_object,
  name,
  embeddings = "auto",
  label = name,
  meta = NULL,
  prefix = if (identical(meta, "auto")) paste0(name, "_") else ""
)
```

## Arguments

- object:

  The parent Seurat object (all cells).

- sub_object:

  A reprocessed child object (a subset of `object`'s cells, re-embedded
  in isolation). Its cell barcodes must match `object`'s.

- name:

  Short id for the subset (e.g. `"tcell"`).

- embeddings:

  Named character vector `new_name = source_reduction` naming
  reduction(s) in `sub_object` to copy (e.g. `c(umap_tcell = "umap")`).
  An unnamed value reuses the source name. Only the child's cells get
  coordinates; parent cells outside the subset are left uncovered
  (`NA`). The first is the view's primary embedding (membership = cells
  with coordinates in it). `"auto"` (default) copies every reduction
  that is new in `sub_object` or whose coordinates differ from
  `object`'s on the shared cells (i.e. it was re-run on the subset), as
  `<name>_<reduction>` (e.g. `tcell_umap`), ordered UMAP, then t-SNE,
  then the rest so a 2-D map is the primary. Reductions the child merely
  inherited unchanged are skipped. Every copied dim becomes a
  `cells.parquet` column, so a 50-dim PCA/Harmony adds 50 (mostly `NA`)
  columns – pass an explicit vector to keep only the maps you want.

- label:

  Human-readable view label (defaults to `name`).

- meta:

  Metadata columns in `sub_object` to copy onto `object` as
  subset-scoped columns (`NA` for non-members). `"auto"` detects them:
  every column that is new in `sub_object`, or whose values differ from
  `object`'s on the shared cells (e.g. re-run clusters or module
  scores). Unchanged inherited columns (`sample`, QC metrics, ...) are
  skipped, and the chosen columns are reported with a message – check
  them before building.

- prefix:

  Prefix for the copied column names. Defaults to `"<name>_"` with
  `meta = "auto"` (so a subset's `seurat_clusters` becomes
  `tcell_seurat_clusters` rather than overwriting the parent's) and to
  `""` (copy under the same name) for an explicit `meta`.

## Value

`object` with the reduction(s)/metadata added and a `scroll_subsets`
attribute recording the spec.
