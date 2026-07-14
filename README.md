# scroll

*Interactive single-cell explorers, served from your own infrastructure.*

<!-- badges: start -->
[![R-CMD-check](https://github.com/shihanli1992/scroll/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/shihanli1992/scroll/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/shihanli1992/scroll/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/shihanli1992/scroll/actions/workflows/pkgdown.yaml)
<!-- badges: end -->

`scroll` turns a processed Seurat object into a polished, interactive web
explorer — a scrolling page of analysis panels (DimPlot, FeaturePlot, Biaxial,
DotPlot, Violin, Proportions, DE, Pseudobulk DE), each with its own fine-grained
controls. A heavy offline **build phase** extracts lightweight on-disk artifacts;
the **runtime** (a bslib Shiny app) reads only those, querying expression one
feature at a time via arrow — so runtime memory stays flat regardless of dataset
size. It deploys as a plain `app.R` on an open-source Shiny Server, with no render
step.

> **Status:** feature-complete for the MVP and validated at scale (a live 4.4M-cell
> atlas). Includes the two-phase build; eight analysis panels (incl. a live DE section
> and replicate-aware **Pseudobulk DE**); a compact **v2 storage format** (int32
> cell-index + float32 + zstd, ~6× smaller than v1 with no precision loss);
> **`scroll_build_stream()`** for streaming multi-million-cell builds on a laptop;
> **`scroll_multi_app()`** for several datasets behind one page; multimodal support
> (multi-assay, per-cell **VDJ / immune repertoire** via `vdj_spec()`, **spatial**
> tissue-image maps via `spatial_spec()`, and **scATAC** peaks via `atac_spec()`); **subset
> views** for reprocessed sub-embeddings; rich declarative controls (data-derived choices,
> cascading levels, preferred defaults, palette picker); PNG/PDF/CSV export; a
> rasterization/memoization performance pass; and the public `register_panel()`
> extension point. Documentation: [pkgdown site](https://shihanli1992.github.io/scroll/) +
> four vignettes (`browseVignettes("scroll")`).

## The two phases

```
Seurat .rds ──scroll_build()──▶  project/            ──scroll_serve()──▶  explorer
 (heavy, once)                    cells.parquet        (light, per session)
                                  expr/<assay>/feature=*/…
                                  manifest.yaml
                                  config.yaml   ← optional defaults
                                  app.R         ← deploy to a Shiny Server
```

The running app **never loads the Seurat object**. It loads `cells.parquet`
(metadata + embeddings, a few MB) once globally, and answers each feature lookup
with an arrow query that reads only that feature's Parquet partition.

## Tutorials

Four vignettes walk through the package (`browseVignettes("scroll")`):

- **Getting started** — from a Seurat object to a running app:
  `vignette("getting-started", package = "scroll")`.
- **Large datasets & streaming builds** — multi-million-cell projects one source at a
  time with `scroll_build_stream()`, the storage format, and `quantize`/`counts`
  tradeoffs: `vignette("large-datasets", package = "scroll")`.
- **Writing a custom panel** — extend the app declaratively with `register_plot_panel()`
  (richer controls + reusable helpers) or the low-level `register_panel()`:
  `vignette("custom-panels", package = "scroll")`.
- **Multimodal: VDJ, spatial & ATAC** — TCR/BCR repertoire with `vdj_spec()` (four
  repertoire panels), 10x Visium / imaging with `spatial_spec()` (a tissue-image
  Spatial panel with zoom), and scATAC peaks with `atac_spec()` (a gene/region peak
  search): `vignette("multimodal", package = "scroll")`.

## Install

```r
# runtime: install.packages(c("Matrix","arrow","dplyr","ggplot2",
#                             "scales","shiny","bslib","yaml","colourpicker"))
# build only: install.packages("SeuratObject")
R CMD INSTALL scroll        # from the repo root
```

## Build a project

```r
library(scroll)
library(SeuratData)                    # for the pbmc3k example object
data("pbmc3k.final", package = "pbmc3k.SeuratData")
obj <- SeuratObject::UpdateSeuratObject(pbmc3k.final)

scroll_build(obj, "pbmc3k")            # object -> project directory
```

`scroll_build()` writes `cells.parquet`, the feature-partitioned `expr/` store,
`manifest.yaml`, a starter `config.yaml`, and a deployable `app.R`. Expression is
quantized to `uint8` by default (`quantize = FALSE` to keep full precision).
Quantization floors values below ~`max/510` to zero, so fraction-expressing
stats (dotplot dot size, DE `pct.1`/`pct.2`) slightly under-count very low
expression — use `quantize = FALSE` when exact fractions matter.

**Storage & precision.** The expression store is compact and sparse, and answers
single-gene queries in tens of milliseconds even at several million cells, with no lossy
precision change (`quantize = FALSE`). Add `counts = TRUE` for a raw-counts store (needed
by Pseudobulk DE; roughly doubles build size). Projects built with an older `scroll` keep
working un-rebuilt. See `vignette("large-datasets")` for the `quantize`/`counts` tradeoffs.

**Build large datasets (streaming).** For many samples that won't fit in RAM as one
merged object, `scroll_build_stream()` reads **one source at a time** into a shared
store — peak memory ~per-source, not per-dataset — and is append-safe (re-run as data
arrives). It powers multi-million-cell atlases on a laptop:

```r
scroll_build_stream("atlas", sources = as.list(Sys.glob("data/*.h5ad")),
  reader = function(f) load_one_sample(f),   # -> a processed Seurat object
  assays = "RNA", embeddings = "umap", meta_cols = c("sample", "celltype"))
```

**Multimodal (CITE-seq / multiome).** By default only the object's default assay
is exported. Pass `assays =` to include more — each becomes its own `expr/<assay>/`
subtree with an independent dequantization scale:

```r
scroll_build(obj, "cite", assays = c("RNA", "ADT"))
```

Every expression panel (FeaturePlot, DotPlot, Violin, DE) then shows an **Assay**
selector; its gene search and DE run against the chosen assay. (Derived assays
like `SCT` or `integrated` are not exported unless you name them explicitly.)

## Explore

```r
scroll_serve("pbmc3k")                 # runs the app locally
```

Or deploy: copy the `pbmc3k/` directory to a Shiny Server — its `app.R`
(`scroll_app(".")`) launches the explorer. Every control auto-populates from the
manifest; `config.yaml` optionally sets defaults:

```yaml
title: "PBMC 3k"          # app-bar label
default_embedding: umap
default_assay: RNA
markers: [CD3D, CD8A, MS4A1, CD14, NKG7]   # DotPlot's starting panel
panels: [dimplot, featureplot, dotplot]    # optional: show exactly these, in order
exclude_panels: [pseudobulk]               # optional: hide specific sections
```

By default each section appears only when the data supports it (e.g. DE needs a
≥2-level grouping). To take explicit control, set `panels:` (an allowlist that shows
exactly those sections, in order, overriding the automatic gating) or
`exclude_panels:` (a denylist) — either in `config.yaml` (no rebuild needed) or via
`scroll_build(..., panels = , exclude_panels = )`.

**Multiple datasets.** `scroll_multi_app()` mounts several built projects behind one
page, each on its own tab (namespaced, independent):

```r
scroll_multi_app(c("PBMC 3k" = "pbmc3k", "CITE-seq" = "cite"))
```

### Panels (v1)

| Panel | Controls |
|-------|----------|
| **DimPlot** | reduction · color-by (metadata, categorical or numeric) · palette · point size · opacity · cluster labels · split-by |
| **FeaturePlot** | gene (server-side search over ~all genes) **or a numeric metadata column** (QC / hashtag-CLR / ADT) · assay · reduction · palette · point size · **colour-quantile clipping** · expressing-on-top · split-by |
| **Biaxial** | pairwise scatters of numeric metadata columns (e.g. hashtag / ADT / QC pairs), coloured by a categorical selection · point size · opacity · palette |
| **DotPlot** | marker genes (ordered multi-select) · group-by · assay · z-score scaling · palette · dot-size range · hclust rows/cols with dendrograms |
| **Violin** | gene · group-by · palette · jitter points |
| **Proportions** | group-by (x) · fill-by · palette · normalize-to-100% |
| **DE** | contrast (group vs group / vs rest) · live Wilcoxon via `presto`, computed once and shown as a **Table** (sortable, min-% / top-N filters) and a **Volcano** (logFC &amp; adj-p cutoffs · top-N labels via ggrepel) |
| **Pseudobulk DE** | aggregate cells into sample-level counts (combined-interaction groups × replicate) and test with **edgeR/limma-voom**; pseudo-replicate modes (`no_replicate` pools each group) + min-cell cutoffs; optional **stability** re-runs (report each gene's selection frequency across random draws). Needs a counts store (see below). Table + Volcano/stability output |

Every plot has an **aspect-ratio** control (reshapes within a fixed canvas) and a
clean black-box theme. Categorical colors are assigned **deterministically by
level name**, so a cell type keeps its color across every panel. Each plot
exports to **PNG** (raster) or **PDF** (vector, for figures/slides) from its
toolbar; the DE table also exports to **CSV**.

Live DE is the one memory-heavy operation: `presto::wilcoxauc` needs an in-memory
matrix, so a run reconstructs the contrast's expression matrix from the store
(on-demand, behind a Compute button). Precomputed `de/<contrast>` tables remain
supported via `view_de_table()` for full rigor / any test.

**Pseudobulk DE** needs raw counts (edgeR/limma-voom are count-based), which are
not in the default quantized store. Build with `counts = TRUE` to also export a
`counts/<assay>.parquet` store; the Pseudobulk DE panel then sums counts per
pseudobulk sample **in arrow's engine** (RAM stays flat) and runs limma-voom:

```r
scroll_build(obj, "proj", counts = TRUE)   # ~doubles expression storage
```

### Subset views (reprocessed sub-embeddings)

A slice of a dataset re-embedded in isolation (e.g. the T cells re-normalized,
re-PCA'd, and given their own UMAP, often with new subclusters) produces a
**partial embedding** — coordinates for the subset's cells only — and metadata
that is meaningful only within it. Declare these as **subset views**: assemble the
reprocessed child onto the parent by barcode with `scroll_add_subset()`, then
build.

```r
obj <- scroll_add_subset(obj, tcell_obj, name = "tcell", label = "T cells",
                         embeddings = c(umap_tcell = "umap"),  # child reduction -> new name
                         meta = "tcell_subcluster")            # subset-only metadata
scroll_build(obj, "proj")                                      # picks up the subset spec
```

The app bar then shows a **View** selector. Picking *T cells* switches to a single
coherent lens: every panel (DE, DotPlot, Proportions, FeaturePlot) restricts to
the subset's cells, the scatter panels default to the sub-UMAP, and the subset's
own columns (e.g. `tcell_subcluster`) become color-by / group-by options — they
stay hidden in the whole-dataset view. Membership is the set of cells with
coordinates in the sub-embedding; the manifest records each embedding's coverage
and each subset's cell count. (`scroll_add_subset()` is pure assembly — the
reprocessing itself happens in your own Seurat pipeline, outside `scroll`.)

## Add your own panel

The panel list is an extension point. The quickest way is `register_plot_panel()` —
declare the controls and one plot function, and scroll builds the whole Shiny
section (UI, control population, a Compute gate, inline error messages, PNG/PDF
export, and app-bar subset/view awareness):

```r
register_plot_panel("counts", label = "Counts", title = "Cells per group",
  after = "dimplot", compute = FALSE,
  controls = list(scroll_input_column("grp", "Group", "categorical")),
  plot = function(cells, input, data)
    ggplot2::ggplot(cells, ggplot2::aes(.data[[input$grp]])) + ggplot2::geom_bar())
scroll_serve("pbmc3k")           # the new section appears in scroll order
```

`plot = function(cells, input, data)` receives the **active** (subset/view-filtered)
cells, the control values by id, and the shared handle — `data$manifest`,
`data$config`, and `data$query1(assay, feature)` / `data$queryN(assay, features)`
for dequantized expression. Controls are declared with `scroll_input_column()`,
`scroll_input_levels()`, `scroll_input_gene()`, an `scroll_input_assay()` /
`scroll_input_embedding()` picker, `scroll_input_custom()` for your own widget, and
`scroll_input_numeric/slider/choice/text()`; wrap any in `scroll_show_when()` for
conditional visibility. To match the built-ins' look, reuse `scroll_point_layer()`
and the `scroll_discrete_colors()` / `scroll_continuous_scale()` palettes.

For full control (cross-output state, custom reactivity), `register_panel(id, ui,
server, …)` takes a raw `ui(id, data)` / `server(id, data, cells_r)` module pair —
the same contract the built-ins use — and the helpers `scroll_render_plot()` /
`scroll_bind_levels()` remove most of the boilerplate. Registering an existing `id`
overrides that panel in place; `scroll_reset_panels()` clears custom ones. See
`vignette("custom-panels")`.

## Query features directly (no app needed)

```r
con <- scroll_connect("pbmc3k")
hit <- scroll_query_feature(con, "RNA", "MS4A1")   # reads one partition
man <- scroll_manifest("pbmc3k")
hit$value <- scroll_dequantize(hit$value, man, "RNA")
scroll_disconnect(con)
```

## Layout

| Path | What |
|------|------|
| `R/build.R` | `scroll_build()` + Seurat v5 extraction |
| `R/manifest.R` / `R/scaffold.R` | manifest + `app.R`/`config.yaml` scaffolding |
| `R/query.R` | arrow query handle + `scroll_query_feature()` / `scroll_query_features()` |
| `R/views.R` | plot cores: `view_umap_colorby`, `view_feature_plot`, `view_dotplot`, … |
| `R/app.R` | `scroll_app()`, `scroll_serve()`, the panel modules |
| `R/styles.R` | the bslib app's design-system CSS + scroll-spy |
