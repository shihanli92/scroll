# scroll 0.0.0.9000

First development release. `scroll` turns a processed Seurat object into a
polished, flat-RAM interactive single-cell explorer.

## Architecture

* **Two phases.** A heavy offline `scroll_build()` extracts lightweight on-disk
  artifacts — cell metadata + embeddings as one small Parquet table
  (`cells.parquet`) and expression as a feature-partitioned Parquet store
  (`expr/<assay>/`). A light runtime (`scroll_app()` / `scroll_serve()`) reads
  only those, querying expression **one gene at a time** through arrow, so
  runtime memory stays flat regardless of matrix size.
* **Deploys as a plain `app.R`** on an open-source Shiny Server — no render step.
  The running app never loads the Seurat object.
* **Build progress.** `scroll_build(verbose = interactive())` reports each phase
  and shows a progress bar over the feature-partition export (the slow step);
  the export is written in feature batches, which also bounds peak memory.
* Expression is `uint8`-quantized by default (`quantize = FALSE` for exact
  precision); an opt-in raw `counts/<assay>.parquet` store (`counts = TRUE`)
  backs the pseudobulk panel.
* Multimodal: pass `assays =` to export more than the default assay; each panel
  gains an **Assay** selector.

## Panels

Eight built-in analysis sections, each with its own fine-grained controls:

* **DimPlot** — embedding coloured by any metadata column.
* **FeaturePlot** — coloured by a gene's expression **or** a numeric metadata
  column, with colour-quantile clipping.
* **Biaxial** — pairwise scatters of numeric metadata columns (hashtag / ADT / QC).
* **DotPlot** — marker panel with z-score scaling and hclust dendrograms.
* **Violin** and **Proportions** — distributions and composition.
* **DE** — live Wilcoxon markers via `presto`, one compute shown as a **Table**
  and a **Volcano**.
* **Pseudobulk DE** — replicate-aware `edgeR` / `limma-voom` on aggregated
  sample-level counts, with pseudo-replicate modes and optional stability re-runs.

## Explore & customize

* **Subset views** — expose a reprocessed slice (its own sub-embedding +
  subclusters) as a linked **View**; `scroll_add_subset()` assembles the child
  onto the parent by barcode, and selecting the view restricts every panel to it.
* **Global cell-subset filter** in the app bar threads through every panel.
* **Manual palettes** — per-level colour pickers on every categorical panel;
  colours are assigned deterministically by level name, so a cell type keeps its
  colour across panels.
* **Export** — every plot to PNG (raster) or PDF (vector); the DE table to CSV,
  with an "export all genes" option.
* **Performance** — large scatters rasterize on screen (`scattermore`) while
  exports stay vector; queries are memoized (LRU) and cosmetic controls debounced.

## Extend

* `register_panel(id, ui, server, ...)` adds a section to every app built
  afterwards, using the same module contract as the built-ins;
  `scroll_reset_panels()` clears custom registrations.

## Documentation

* Two vignettes: *Getting started* (Seurat object → running app) and *Writing a
  custom panel* (a worked centroid-map + minimum-spanning-tree example).
