# Build a scroll project from a Seurat object

Reads a processed Seurat object exactly once and emits a self-contained
project directory of lightweight on-disk artifacts: a small Parquet
table of cell metadata + embeddings (`cells.parquet`), a
feature-partitioned Parquet store of expression (`expr/<assay>/`), a
`manifest.yaml`, and app scaffolding (`config.yaml`, `app.R`). The
running app never touches the object again.

## Usage

``` r
scroll_build(
  object,
  outdir,
  assays = NULL,
  embeddings = NULL,
  meta_cols = NULL,
  quantize = TRUE,
  counts = FALSE,
  subsets = NULL,
  vdj = NULL,
  spatial = NULL,
  atac = NULL,
  panels = NULL,
  exclude_panels = NULL,
  max_levels = .SCROLL_MAX_LEVELS,
  overwrite = FALSE,
  verbose = interactive()
)
```

## Arguments

- object:

  A processed Seurat object, or a path to an `.rds` file.

- outdir:

  Output project directory (created if missing).

- assays:

  Assays to export. `NULL` exports the default assay only.

- embeddings:

  Reductions to export. `NULL` exports all reductions.

- meta_cols:

  Metadata columns to export. `NULL` infers a sensible set.

- quantize:

  If `TRUE`, expression is quantized to `uint8` (256 levels), with a
  per-assay `max` recorded in the manifest for dequantization. Values
  below ~`max/510` round to zero, so fraction-expressing statistics
  (dotplot dot size, DE `pct.1`/`pct.2`) slightly under-count
  low-expression cells; use `quantize = FALSE` for exact values.

- counts:

  If `TRUE`, also export each assay's raw `counts` layer to a
  `counts/<assay>.parquet` store (integer, unquantized). This roughly
  doubles expression storage but enables pseudobulk differential
  expression
  ([`scroll_pseudobulk_de()`](https://shihanli92.github.io/scroll/reference/scroll_pseudobulk_de.md),
  which needs counts for edgeR/limma-voom).

- subsets:

  Optional named list declaring reprocessed **subset views** — a slice
  of the data re-embedded in isolation (e.g. a T-cell-only UMAP) plus
  any subset-only metadata (subclusters). Each element is
  `list(label =, embeddings =, meta =)`: `embeddings` names one or more
  exported reductions that cover only the subset's cells (the first is
  the view's primary embedding; membership = cells with non-`NA`
  coordinates in it); `meta` names metadata columns that are meaningful
  only within the subset. In the app a subset becomes a **View** that
  restricts every panel to those cells and exposes its embedding +
  metadata. Defaults to any spec attached by
  [`scroll_add_subset()`](https://shihanli92.github.io/scroll/reference/scroll_add_subset.md).

- vdj:

  Optional
  [`vdj_spec()`](https://shihanli92.github.io/scroll/reference/vdj_spec.md)
  describing per-cell TCR/BCR columns. When supplied, a compact
  `repertoire/` store (clone table + diversity / gene-usage residuals /
  tissue correlation) is baked and the VDJ panels are enabled.

- spatial:

  Optional
  [`spatial_spec()`](https://shihanli92.github.io/scroll/reference/spatial_spec.md)
  (or `TRUE` for defaults) describing a tissue image, or a **list** of
  specs for multi-FOV / multi-slide objects. For each spec
  `GetTissueCoordinates()` is exported as a spatial embedding, the H&E
  image is baked as a raster asset, and the Spatial panel is enabled.
  Multiple specs must have distinct `name`s.

- atac:

  Optional
  [`atac_spec()`](https://shihanli92.github.io/scroll/reference/atac_spec.md)
  (or `TRUE` for defaults) naming a peaks assay. When supplied, the
  assay is marked `kind: peaks`, each peak's coordinates are parsed from
  its name, a peak-annotation table (with nearest gene when a Signac
  annotation is available) is baked, and the Peaks panel is enabled.

- panels:

  Optional character vector of panel ids to show **exactly**, in this
  order — an explicit override of the automatic data-driven gating (e.g.
  `c("dimplot", "featureplot", "spatial")`). Written to `config.yaml` as
  `panels:` and also hand-editable there. `NULL` (default) keeps the
  automatic behaviour.

- exclude_panels:

  Optional character vector of panel ids to hide (e.g.
  `c("de", "pseudobulk")`); applied after gating / the `panels`
  allowlist.

- max_levels:

  Cap on how many distinct values of a categorical column are cached as
  a `levels:` list in `manifest.yaml`. Columns above the cap (e.g. a
  clone id or barcode with thousands of values) still export normally
  and stay fully usable — the manifest just records their `n_levels`
  count and the app recomputes the level set from `cells.parquet` when a
  control needs it. Keeps the manifest small and hand-editable. Defaults
  to 200.

- overwrite:

  If `TRUE`, an existing `outdir` is removed first.

- verbose:

  If `TRUE`, report progress: a step message per phase and a progress
  bar over each assay's feature-partition export (the slow step).
  Defaults to
  [`interactive()`](https://rdrr.io/r/base/interactive.html), so it is
  quiet in scripts/CI; pass `TRUE` to force progress in a
  non-interactive session.

## Value

`outdir`, invisibly.
