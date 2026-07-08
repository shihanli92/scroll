# scroll

*Scroll-driven single-cell stories, served from your own infrastructure.*

`scroll` turns a processed Seurat object into an interactive, scroll-driven web
story for presenting single-cell analyses. A heavy **build phase** touches the
object once and emits lightweight on-disk artifacts; a memory-light **runtime
phase** (Quarto [closeread](https://closeread.dev) + Shiny) reads only those
artifacts, querying expression one feature at a time off disk — so runtime RAM
stays flat regardless of matrix size.

> **Status: Phase 2.** Single-assay build; the full view grammar
> (`umap_colorby`, `feature_plot`, `dotplot`, `violin`, `proportions`,
> `de_table`); toggles (embedding / split / subset / labels); duckdb feature
> lookups; and the scroll × reactive coupling. Multi-assay, `register_view()`,
> and WebGL scatter are planned for later phases.

## The two phases

```
Seurat .rds ──scroll_build()──▶  project/            ──scroll_serve()──▶  web story
 (heavy, once)                    cells.parquet        (light, per session)
                                  expr/<assay>/feature=*/…
                                  manifest.yaml
                                  config.yaml   ← authored
                                  story.qmd     ← authored
```

The running app **never loads the Seurat object**. It loads `cells.parquet`
(metadata + embeddings, a few MB) once globally, and answers each feature lookup
with a duckdb query that reads only that feature's Parquet partition.

## Install

```r
# install.packages(c("SeuratObject", "Matrix", "arrow", "duckdb", "DBI",
#                     "ggplot2", "yaml", "shiny"))
R CMD INSTALL scroll        # from the repo root
```

For the runtime app you also need the **Quarto CLI** and the closeread
extension (installed *inside a project directory*):

```sh
# Quarto: https://quarto.org
quarto add qmd-lab/closeread
```

## Build a project

```r
library(scroll)
library(SeuratData)                    # for the pbmc3k example object
data("pbmc3k.final", package = "pbmc3k.SeuratData")
obj <- SeuratObject::UpdateSeuratObject(pbmc3k.final)

scroll_build(obj, "pbmc3k-story")      # object -> project directory
```

`scroll_build()` writes `cells.parquet`, the `expr/RNA/` feature-partitioned
store, `manifest.yaml`, and authorable `config.yaml` + `story.qmd` scaffolding.
Expression is quantized to `uint8` by default (`quantize = FALSE` to keep full
precision).

## Query features directly (no app needed)

```r
con <- scroll_connect("pbmc3k-story")
hit <- scroll_query_feature(con, "RNA", "MS4A1")   # reads one partition
man <- scroll_manifest("pbmc3k-story")
hit$value <- scroll_dequantize(hit$value, man, "RNA")
scroll_disconnect(con)
```

## Author + serve the story

Edit `config.yaml` (ordered sections + view specs) and `story.qmd` (the
closeread narrative), then:

```r
scroll_serve("pbmc3k-story")     # quarto preview (localhost, for the meeting)
scroll_render("pbmc3k-story")    # rendered output to upload to a Shiny Server
```

## View grammar

Sections in `config.yaml` name a view type and its params:

| View | Params | Purpose |
|------|--------|---------|
| `umap_colorby` | `embedding`, `color_by` | embedding colored by a metadata column |
| `feature_plot` | `embedding`, `feature`, `assay` | expression on the embedding (live search target) |
| `dotplot` | `group_by`, `features`, `assay` | mean expression × fraction expressing across groups |
| `violin` | `group_by`, `feature`, `assay` | per-group distribution for a feature |
| `proportions` | `group_by`, `fill_by` | stacked composition of one categorical within another |
| `de_table` | `contrast`, `top` | a precomputed DE table (see below) |

**Toggles** (embedding, split-by, subset, labels) are Shiny inputs carried as
state on the active view; the feature search box retargets whatever view is in
focus (scatter recolors, violin switches feature, dotplot adds the gene).

### Precomputed DE tables (`de/`)

`de_table` reads a per-contrast file the author drops into the project's `de/`
directory — `de/<contrast>.parquet` (or `.csv` / `.tsv`) — and a section
references it by name:

```yaml
- id: bcell_de
  view: de_table
  params: { contrast: B_vs_rest, top: 50 }
```

`scroll` reads the table as-is (e.g. a `FindMarkers()` result written to
`de/B_vs_rest.parquet`); it does not compute DE.

## How the coupling works

The sticky visual answers to two drivers: **scroll position** (closeread) and
**the feature search box** (Shiny). `scroll` models it as one persistent Shiny
output redrawn by one reactive over `(active_section, feature)`. A small JS
bridge (`inst/app/cr-bridge.js`) reports the active closeread trigger to Shiny
as the global input `active_section`; the two drivers meet only inside
`.scroll_resolve()` and never contend for the DOM. The pure resolve logic and
the views are unit-tested without Quarto (`tests/testthat/`).

## Layout

| Path | What |
|------|------|
| `R/build.R` | `scroll_build()` + Seurat v5 extraction |
| `R/manifest.R` / `R/scaffold.R` | manifest + `config.yaml`/`story.qmd` scaffolding |
| `R/query.R` | duckdb connection + `scroll_query_feature()` |
| `R/views.R` | `umap_colorby`, `feature_plot`, `render_view()` dispatch |
| `R/coupling.R` | the scroll × reactive coupling + Shiny app UI/server |
| `R/serve.R` | `scroll_serve()` / `scroll_render()` (Quarto wrappers) |
| `inst/app/cr-bridge.js` | closeread → Shiny active-section bridge |
