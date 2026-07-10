# scroll

*Interactive single-cell explorers, served from your own infrastructure.*

`scroll` turns a processed Seurat object into a polished, interactive web
explorer — a scrolling page of analysis panels (DimPlot, FeaturePlot, DotPlot,
Violin, Proportions, DE), each with its own fine-grained controls. A heavy
offline **build phase** extracts lightweight on-disk artifacts; the **runtime**
(a bslib Shiny app) reads only those, querying expression one feature at a time
via duckdb — so runtime memory stays flat regardless of dataset size. It deploys
as a plain `app.R` on an open-source Shiny Server, with no render step.

> **Status:** build phase + all six panels (incl. a live DE section — ranked
> table + volcano from one `presto` compute) + multimodal (multi-assay) support,
> wired end-to-end and validated live on pbmc3k. `register_panel()` and WebGL
> scatter for very large datasets are next.

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
with a duckdb query that reads only that feature's Parquet partition.

## Install

```r
# install.packages(c("SeuratObject","Matrix","arrow","duckdb","DBI",
#                     "ggplot2","scales","shiny","bslib","yaml"))
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
pseudobulk sample **in duckdb** (RAM stays flat) and runs limma-voom:

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

The panel list is an extension point. `register_panel()` adds a section to every
`scroll_app()` built afterwards — a panel is just a `ui(id, data)` /
`server(id, data, cells_r)` pair, the same contract the built-ins use:

```r
count_ui <- function(id, data) {
  ns <- shiny::NS(id)
  cats <- names(Filter(function(x) identical(x$type, "categorical"), data$manifest$meta))
  shiny::tagList(shiny::selectInput(ns("grp"), "Group", cats),
                 shiny::plotOutput(ns("plot")))
}
count_server <- function(id, data, cells_r = shiny::reactive(data$cells)) {
  shiny::moduleServer(id, function(input, output, session) {
    output$plot <- shiny::renderPlot(barplot(table(cells_r()[[input$grp]])))
  })
}

register_panel("counts", count_ui, count_server, label = "Counts",
               title = "Cells per group", after = "dimplot")
scroll_serve("pbmc3k")           # the new section appears in scroll order
```

`data` is the shared handle — `data$cells` (metadata + embeddings), `data$manifest`,
`data$config`, and `data$query1(assay, feature)` / `data$queryN(assay, features)`
for dequantized expression. Use `cells_r()` (not `data$cells`) so the panel honours
the app-bar subset filter. Registering an existing `id` overrides that panel in
place (including a built-in); `scroll_reset_panels()` clears custom ones.

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
| `R/query.R` | duckdb connection + `scroll_query_feature()` / `scroll_query_features()` |
| `R/views.R` | plot cores: `view_umap_colorby`, `view_feature_plot`, `view_dotplot`, … |
| `R/app.R` | `scroll_app()`, `scroll_serve()`, the panel modules |
| `R/styles.R` | the bslib app's design-system CSS + scroll-spy |
