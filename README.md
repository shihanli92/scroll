# scroll

*Interactive single-cell explorers, served from your own infrastructure.*

`scroll` turns a processed Seurat object into a polished, interactive web
explorer — a scrolling page of analysis panels (DimPlot, FeaturePlot, DotPlot),
each with its own fine-grained controls. A heavy offline **build phase** extracts
lightweight on-disk artifacts; the **runtime** (a bslib Shiny app) reads only
those, querying expression one feature at a time via duckdb — so runtime memory
stays flat regardless of dataset size. It deploys as a plain `app.R` on an
open-source Shiny Server, with no render step.

> **Status:** build phase + DimPlot / FeaturePlot / DotPlot panels, wired
> end-to-end and validated live on pbmc3k. Violin / Proportions / DE panels,
> `register_panel()`, multi-assay selectors, and WebGL scatter are next.

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
| **FeaturePlot** | gene (server-side search over ~all genes) · assay · reduction · palette · point size · expressing-on-top · split-by |
| **DotPlot** | marker genes (ordered multi-select) · group-by · assay · z-score scaling · palette · dot-size range |

Categorical colors are assigned **deterministically by level name**, so a cell
type keeps its color across every panel.

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
