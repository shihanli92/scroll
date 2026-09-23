# scroll

*Interactive single-cell explorers, served from your own infrastructure.*

<!-- badges: start -->
[![R-CMD-check](https://github.com/shihanli1992/scroll/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/shihanli1992/scroll/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/shihanli1992/scroll/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/shihanli1992/scroll/actions/workflows/pkgdown.yaml)
<!-- badges: end -->

`scroll` turns a processed Seurat object into a polished, interactive web
explorer — a scrolling page of analysis panels (DimPlot, FeaturePlot, Signature,
Biaxial, DotPlot, Heatmap, Violin, Proportions, DE, Pseudobulk DE), each with its own fine-grained
controls. A heavy offline **build phase** extracts lightweight on-disk artifacts;
the **runtime** (a bslib Shiny app) reads only those, querying expression one
feature at a time via arrow — so runtime memory stays flat regardless of dataset
size. It deploys as a plain `app.R` on an open-source Shiny Server, with no render
step.

> **Status:** feature-complete for the MVP and validated at scale (a live 4.4M-cell
> atlas). Includes the two-phase build; eight analysis panels (incl. a live DE panel
> and replicate-aware **Pseudobulk DE**); a compact **v2 storage format** (int32
> cell-index + zstd, ~5× smaller than v1; lossless with `quantize = FALSE`);
> **`scroll_build_stream()`** for streaming multi-million-cell builds on a laptop;
> **`scroll_multi_app()`** for several datasets behind one page; multimodal support
> (multi-assay, per-cell **VDJ / immune repertoire** via `vdj_spec()` — reads TCR/BCR
> from **scRepertoire / AIRR-dandelion / Platypus** out of the box, **spatial**
> tissue-image maps via `spatial_spec()`, and **scATAC** peaks via `atac_spec()`); **subset
> views** for reprocessed sub-embeddings; rich declarative controls (data-derived choices,
> cascading levels, preferred defaults, palette picker); PNG/PDF/CSV export; a
> rasterization/memoization performance pass; and the public `register_panel()`
> extension point. Documentation: [pkgdown site](https://shihanli1992.github.io/scroll/).

## The two phases

```
Seurat .rds ──scroll_build()──▶  project/            ──scroll_serve()──▶  explorer
 (heavy, once)                    cells.parquet        (light, per session)
                                  expr/<assay>/part-0.parquet
                                  manifest.yaml
                                  config.yaml   ← optional defaults
                                  app.R         ← deploy to a Shiny Server
```

The running app **never loads the Seurat object**. It loads `cells.parquet`
(metadata + embeddings, a few MB) once globally, and answers each feature lookup
with an arrow query that reads only that feature's Parquet partition.

## Install

Install from GitHub with **devtools** (or the lighter **remotes**) — it pulls the
runtime dependencies automatically:

```r
# install.packages("devtools")
devtools::install_github("shihanli92/scroll")

# to build projects from a Seurat object you also need SeuratObject (a Suggest,
# so it is not pulled automatically); the runtime app never needs it:
install.packages("SeuratObject")
```

`arrow` must be built with the **zstd** codec (scroll's store is zstd-compressed);
confirm with `arrow::arrow_info()$capabilities[["zstd"]]`. On HPC or other
environments where `install.packages("arrow")` yields a codec-less build, reinstall
with `Sys.setenv(LIBARROW_MINIMAL = "false"); install.packages("arrow")`, or use
`conda install -c conda-forge r-arrow`.

<details><summary>Install from a local checkout instead</summary>

```r
# runtime deps: install.packages(c("Matrix","arrow","dplyr","ggplot2",
#                                  "scales","shiny","bslib","yaml","colourpicker"))
devtools::install(".")      # from the repo root, or: R CMD INSTALL .
```
</details>

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
working un-rebuilt.

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
exclude_panels: [pseudobulk]               # optional: hide specific panels
```

By default each panel appears only when the data supports it (e.g. DE needs a
≥2-level grouping). To take explicit control, set `panels:` (an allowlist that shows
exactly those panels, in order, overriding the automatic gating) or
`exclude_panels:` (a denylist) — either in `config.yaml` (no rebuild needed) or via
`scroll_build(..., panels = , exclude_panels = )`.

**Multiple datasets.** `scroll_multi_app()` mounts several built projects behind one
page, each on its own tab (namespaced, independent):

```r
scroll_multi_app(c("PBMC 3k" = "pbmc3k", "CITE-seq" = "cite"))
```

### Panels

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
(on-demand, behind a Compute button).

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

Or let scroll work out what the subset changed. By default (`embeddings = "auto"`)
it copies every reduction that is new in the child or was re-run on it (e.g. a
subset UMAP) as `<name>_<reduction>`, UMAP first so it becomes the view's primary
map; `meta = "auto"` does the same for metadata: every column that is new or differs
from the parent on the shared cells (re-run clusters, recomputed scores), prefixed so
it can't overwrite the parent's (`seurat_clusters` becomes `tcells_seurat_clusters`).
Both report what they picked. Handy when attaching several lineages:

```r
for (n in c("tcells", "nk", "myeloid", "bcells"))
  obj <- scroll_add_subset(obj, lineages[[n]], name = n, meta = "auto")
scroll_build(obj, "proj")
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
panel (UI, control population, a Compute gate, inline error messages, PNG/PDF
export, and app-bar subset/view awareness):

```r
register_plot_panel("counts", label = "Counts", title = "Cells per group",
  after = "dimplot", compute = FALSE,
  controls = list(scroll_input_column("grp", "Group", "categorical")),
  plot = function(cells, input, data)
    ggplot2::ggplot(cells, ggplot2::aes(.data[[input$grp]])) + ggplot2::geom_bar())
scroll_serve("pbmc3k")           # the new panel appears in scroll order
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
the same contract the built-ins use — and `scroll_render_plot()` removes most of the
plot/download boilerplate. Registering an existing `id`
overrides that panel in place; `scroll_reset_panels()` clears custom ones.

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
| `R/views.R` | internal plot cores behind each built-in panel (`view_*`) |
| `R/app.R` | `scroll_app()` / `scroll_multi_app()` / `scroll_serve()`, the app shell + panel registry |
| `R/data-handle.R` | the runtime data handle: `.scroll_load()` + LRU-cached feature queries |
| `R/panel-helpers.R` | control-choice + plot-toolbar helpers shared across the built-in panels |
| `R/panel-<view>.R` | one file per built-in view (`panel-dimplot.R`, `panel-de.R`, …) |
| `R/panel-builder.R` | `register_plot_panel()` + `scroll_input_*` declarative controls |
| `R/vdj.R` / `R/spatial.R` / `R/atac.R` | multimodal `*_spec()` + build-time bake logic |
| `R/panel-vdj.R` / `R/panel-spatial.R` / `R/panel-atac.R` | the matching multimodal panels |
| `R/styles.R` | the bslib app's design-system CSS + scroll-spy |
