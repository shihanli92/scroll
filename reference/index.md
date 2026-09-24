# Package index

## Build & assemble

Turn a processed Seurat object into an on-disk project, and attach
reprocessed sub-embeddings as linked views.

- [`scroll_build()`](https://shihanli92.github.io/scroll/reference/scroll_build.md)
  : Build a scroll project from a Seurat object

- [`scroll_build_stream()`](https://shihanli92.github.io/scroll/reference/scroll_build_stream.md)
  : Build a scroll project by streaming many sources (large-dataset
  build)

- [`scroll_update()`](https://shihanli92.github.io/scroll/reference/scroll_update.md)
  : Add or update reductions / metadata on a built project (no expr
  rebuild)

- [`scroll_add_meta()`](https://shihanli92.github.io/scroll/reference/scroll_add_meta.md)
  : Add (or replace) a metadata column in a built project – no Seurat
  object needed

- [`scroll_add_subset()`](https://shihanli92.github.io/scroll/reference/scroll_add_subset.md)
  : Attach a reprocessed subset onto a parent object

- [`vdj_spec()`](https://shihanli92.github.io/scroll/reference/vdj_spec.md)
  :

  Describe a dataset's TCR/BCR (VDJ) metadata for
  [`scroll_build()`](https://shihanli92.github.io/scroll/reference/scroll_build.md)

- [`spatial_spec()`](https://shihanli92.github.io/scroll/reference/spatial_spec.md)
  :

  Describe a dataset's spatial image + coordinates for
  [`scroll_build()`](https://shihanli92.github.io/scroll/reference/scroll_build.md)

- [`atac_spec()`](https://shihanli92.github.io/scroll/reference/atac_spec.md)
  :

  Describe a dataset's scATAC peaks assay for
  [`scroll_build()`](https://shihanli92.github.io/scroll/reference/scroll_build.md)

## Run the app

Launch the runtime explorer, or preview a single panel while developing
it.

- [`scroll_app()`](https://shihanli92.github.io/scroll/reference/scroll_app.md)
  : Build the scroll explorer app
- [`scroll_multi_app()`](https://shihanli92.github.io/scroll/reference/scroll_multi_app.md)
  : Build a multi-dataset scroll explorer
- [`scroll_serve()`](https://shihanli92.github.io/scroll/reference/scroll_serve.md)
  : Run the scroll explorer locally
- [`scroll_preview_panel()`](https://shihanli92.github.io/scroll/reference/scroll_preview_panel.md)
  : Preview a single panel while developing it

## Query layer

The flat-RAM arrow query interface, usable programmatically without the
app.

- [`scroll_connect()`](https://shihanli92.github.io/scroll/reference/scroll_connect.md)
  : Open a scroll project for querying
- [`scroll_disconnect()`](https://shihanli92.github.io/scroll/reference/scroll_disconnect.md)
  : Close a scroll project handle
- [`scroll_manifest()`](https://shihanli92.github.io/scroll/reference/scroll_manifest.md)
  : Read a scroll project manifest
- [`scroll_query_feature()`](https://shihanli92.github.io/scroll/reference/scroll_query_feature.md)
  : Query one feature's expression from the Parquet store
- [`scroll_query_features()`](https://shihanli92.github.io/scroll/reference/scroll_query_features.md)
  : Query several features at once
- [`scroll_dequantize()`](https://shihanli92.github.io/scroll/reference/scroll_dequantize.md)
  : Map stored (possibly quantized) values back to normalized expression

## Differential expression

Live Wilcoxon markers and replicate-aware pseudobulk testing.

- [`scroll_de()`](https://shihanli92.github.io/scroll/reference/scroll_de.md)
  : Live differential expression for a contrast (presto / Wilcoxon)
- [`scroll_pseudobulk_de()`](https://shihanli92.github.io/scroll/reference/scroll_pseudobulk_de.md)
  : Pseudobulk differential expression (edgeR / limma-voom)

## Extend the app

Register custom panels — declaratively from a plot function plus
declared controls, or with the low-level module contract the built-ins
use.

- [`register_plot_panel()`](https://shihanli92.github.io/scroll/reference/register_plot_panel.md)
  : Register a plot panel declaratively (no Shiny boilerplate)

- [`register_panel()`](https://shihanli92.github.io/scroll/reference/register_panel.md)
  : Register a custom panel in the scroll explorer

- [`scroll_input_column()`](https://shihanli92.github.io/scroll/reference/scroll_input.md)
  [`scroll_input_levels()`](https://shihanli92.github.io/scroll/reference/scroll_input.md)
  [`scroll_input_gene()`](https://shihanli92.github.io/scroll/reference/scroll_input.md)
  [`scroll_input_numeric()`](https://shihanli92.github.io/scroll/reference/scroll_input.md)
  [`scroll_input_slider()`](https://shihanli92.github.io/scroll/reference/scroll_input.md)
  [`scroll_input_choice()`](https://shihanli92.github.io/scroll/reference/scroll_input.md)
  [`scroll_input_text()`](https://shihanli92.github.io/scroll/reference/scroll_input.md)
  [`scroll_input_palette()`](https://shihanli92.github.io/scroll/reference/scroll_input.md)
  [`scroll_input_assay()`](https://shihanli92.github.io/scroll/reference/scroll_input.md)
  [`scroll_input_embedding()`](https://shihanli92.github.io/scroll/reference/scroll_input.md)
  [`scroll_input_custom()`](https://shihanli92.github.io/scroll/reference/scroll_input.md)
  :

  Declarative controls for
  [`register_plot_panel()`](https://shihanli92.github.io/scroll/reference/register_plot_panel.md)

- [`scroll_show_when()`](https://shihanli92.github.io/scroll/reference/scroll_show_when.md)
  : Show a control only when another control has a given value

- [`scroll_style_input()`](https://shihanli92.github.io/scroll/reference/scroll_style_input.md)
  : Put a control in the plot's Style sheet

- [`scroll_point_layer()`](https://shihanli92.github.io/scroll/reference/scroll_helpers.md)
  [`scroll_discrete_colors()`](https://shihanli92.github.io/scroll/reference/scroll_helpers.md)
  [`scroll_continuous_scale()`](https://shihanli92.github.io/scroll/reference/scroll_helpers.md)
  [`scroll_group_labels()`](https://shihanli92.github.io/scroll/reference/scroll_helpers.md)
  : Reusable plotting helpers for custom views

- [`scroll_render_plot()`](https://shihanli92.github.io/scroll/reference/scroll_render_plot.md)
  : Wire a plot output with error surfacing and PNG/PDF downloads

- [`scroll_columns()`](https://shihanli92.github.io/scroll/reference/scroll_columns.md)
  : Metadata columns of a given type

- [`scroll_embeddings()`](https://shihanli92.github.io/scroll/reference/scroll_embeddings.md)
  : Embeddings (reductions) available in a scroll project

- [`scroll_features()`](https://shihanli92.github.io/scroll/reference/scroll_features.md)
  : Features (e.g. genes) of an assay in a scroll project

- [`scroll_reset_panels()`](https://shihanli92.github.io/scroll/reference/scroll_reset_panels.md)
  :

  Clear all custom panels registered with
  [`register_panel()`](https://shihanli92.github.io/scroll/reference/register_panel.md)
