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

## Multi-dataset

* **`scroll_multi_app(projects)`** mounts several built projects behind one app —
  a tab per dataset, each with its own app bar, manifest-driven panels, and query
  handle. The per-dataset page + server were refactored into a namespaced module
  (`shiny::NS()`), so independent datasets (even different organisms / feature
  spaces) coexist without id collisions; `scroll_app()` is the single-dataset case
  and is unchanged. Custom panels should read baked assets from the per-dataset
  handle (`data$dir`) so each tab loads its own files.

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
* **Stability** — killed scrollbar-driven resize loops on tall multi-panel pages.
  A scrollbar toggling as plots render (vertical bar → width change, horizontal bar
  → height change) fires window `resize`, which re-renders every fluid-width
  `plotOutput`, which nudges the size back — an endless loop where plots appeared to
  update on their own. The page now reserves the vertical gutter
  (`overflow-y:scroll; scrollbar-gutter:stable`) and hides page-level horizontal
  overflow (`overflow-x:hidden`; wide plots still scroll inside their card). The
  plot wrapper also pins `min-width:0` and snaps the plot output to a whole CSS
  pixel (`width: round(down, 100%, 1px)`). bslib's grid columns are fractional
  (e.g. 518.25px); on a HiDPI display (`devicePixelRatio` 2) that .25px is half a
  device pixel, so Shiny rendered the image at a rounded device size that displayed
  back a hair different, tripping the per-output `ResizeObserver` into re-rendering
  every plot forever (only at HiDPI + narrow widths). Snapping the width to an
  integer makes the device pixels land exactly, so nothing oscillates.

## Extend

* `register_panel(id, ui, server, ...)` adds a section to every app built
  afterwards, using the same module contract as the built-ins;
  `scroll_reset_panels()` clears custom registrations.
* `register_plot_panel(id, plot, controls, ...)` is a declarative wrapper: describe
  the controls (`scroll_input_column()`, `scroll_input_levels()`,
  `scroll_input_gene()`, …) and a single `plot(cells, input, data)` function, and
  scroll generates the whole module — UI, control population, a Compute gate,
  inline error messages, PNG/PDF export, and subset/view awareness. The lower-level
  helpers (`scroll_render_plot()`, `scroll_bind_levels()`, `scroll_columns()`) are
  exported for hand-written panels too.
* **Richer declarative controls** for `register_plot_panel()`:
  `scroll_input_assay()` (auto-hidden on single-assay data),
  `scroll_input_embedding()` and view-aware `scroll_input_column(view_aware=)`
  (narrow to the active subset view), `scroll_input_custom()` (wrap your own
  `ui`/`bind`), and `scroll_show_when()` for conditional visibility.
* **Reusable render + data helpers exported** so custom views match the built-ins:
  `scroll_point_layer()` (raster-aware scatter), `scroll_discrete_colors()` /
  `scroll_continuous_scale()` (shared palettes), `scroll_group_labels()`, and the
  accessors `scroll_columns(view=)`, `scroll_embeddings()`, `scroll_features()`.
* **Fewer silent failures.** A non-ggplot plot return now shows a clear message; a
  `scroll_input_levels(from=)` that names no column/control warns; contradictory
  `required` + `none` warns at construction.

## Documentation

* Two vignettes: *Getting started* (Seurat object → running app) and *Writing a
  custom panel* (a worked centroid-map + minimum-spanning-tree example).
