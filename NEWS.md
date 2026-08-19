# scroll 0.1.0

First tagged release. `scroll` turns a processed Seurat object into a
polished, flat-RAM interactive single-cell explorer.

## New: multi-gene and co-expression blend in the FeaturePlot panel

* The FeaturePlot gene box is now **multi-select**, and the numeric-metadata
  selector alongside it too. Pick any mix of genes and numeric columns and the
  panel lays out a **grid of feature plots**, one per feature, each with its own
  colour scale (as `Seurat::FeaturePlot(features = c(...))` does), up to 12.
  New exported view core `view_feature_multi()`.
* A **Blend** toggle reproduces `Seurat::FeaturePlot(blend = TRUE)` for the two
  selected genes: four views — each gene alone, their co-expression blend, and a
  2-D colour key — with an adjustable **Blend threshold** and pickers for the two
  gene colours. The blend colours are a faithful, Seurat-free port of Seurat's
  `BlendMatrix`/`BlendExpression` (byte-identical colours). New exported view core
  `view_feature_blend()`. The CSV export includes every plotted gene.
* Both features use the (Suggested) `patchwork` for their grid layout; without it
  the panel degrades to a single plot.
* FeaturePlot's **colour quantile** cutoffs are now taken over the expressing
  (non-zero) cells. Single-cell expression is zero-inflated, so a quantile of the
  full vector stayed pinned at 0 until the fraction passed the (often >90%) zero
  share — the lower cutoff had no visible effect. It now spans the real expression
  range.

## New: paste a gene list into the DotPlot marker box

* Paste a gene list **straight into the DotPlot marker box** — separated by
  spaces, tabs, commas or newlines — and the genes are added to the selection in
  the pasted order (typing one at a time still works in the same box). Symbols are
  matched case-insensitively and any that aren't in the assay are reported.

## New: reproducible CSV export for the built-in plot panels

* Every built-in plot panel (DimPlot, FeaturePlot, Biaxial, Violin, DotPlot,
  Proportions) now offers a **CSV** download beside its PNG/PDF buttons, containing
  the minimal source data behind the figure. For the per-cell panels the CSV leads
  with the **cell barcode** plus the plotted columns (embedding coords + colour-by /
  expression / group), so an exported point maps back to a tracked cell; the two
  aggregated panels export their group×feature (DotPlot: `avg_expr`, `pct_expressing`)
  and group×category (Proportions: `n_cells`, `proportion`) summary. The CSV reuses
  the exact data the plot draws, so it reproduces the figure. Custom panels built with
  `.scroll_plot_area(csv = TRUE)` + `.scroll_plot_downloads(csv_r = )` get the same.

## New: `scroll_preview_panel()` — fast custom-panel dev loop

* `scroll_preview_panel(panel, dir)` mounts **one** panel (by id) in the real app
  harness — the app bar with the cell-subset filter and subset-View selector, but
  none of the other panels — so a custom panel can be driven for quick feedback
  without building or scrolling through the whole app. Because it reuses the actual
  page + wiring, the panel is tested exactly as it will run (under the app-bar filter
  and subset views). The loop is `source()` the panel file → `scroll_preview_panel()`
  → edit → re-source → re-run.

## Correctness fixes

* **`scroll_update()` no longer re-escapes non-ASCII feature names.** It rewrites
  the manifest, but wrote it without `unicode = TRUE` (unlike `scroll_build()`), so
  an update on a project with a non-ASCII feature (e.g. an antibody `FcεRIa`)
  re-mangled the name to `<U+XXXX>` and broke its query/label match. Now consistent
  with the build path.
* **Pseudobulk DE is RNG-neutral and reproducible.** `scroll_pseudobulk_de()`'s
  random pseudo-replicate draw used the session's global RNG as a side effect and
  was non-reproducible. It now seeds the draw (new `seed =` argument, default `1L`)
  and restores the caller's `.Random.seed`, mirroring `scroll_de()`'s cap.
  `scroll_pseudobulk_stability()` still varies the draw per run (distinct per-run
  seeds), so it stays reproducible as a whole without perturbing the session.
* **`scroll_build_stream()` validates each source up front.** A source missing a
  declared assay / embedding / metadata column now errors with the offending source
  named, instead of failing deep in the final cells-frame merge.
* **`scroll_de()`** emits a heads-up when testing a very large contrast with no
  `max_cells` cap (the one runtime step whose memory scales with the matrix),
  pointing at `max_cells` to bound time and peak RAM.

## Fix: `scroll_add_subset()` preserves scoped-column types

* A subset's scoped metadata columns kept their **original type** instead of being
  coerced to character. Numeric scoped columns (e.g. `AddModuleScore` signature
  scores, `percent.*`) were being stringified, so the manifest classified them as
  high-cardinality **categorical** (one "level" per cell) rather than **numeric** —
  breaking numeric colour-by / violins / sliders on those columns. Rebuild any
  project with numeric scoped columns to pick up the correct types.
* The **Violin** (Metadata mode) and **Biaxial** numeric-column pickers are now
  view-aware, so scoped numeric columns appear when their subset view is active —
  matching the categorical Group-by, which was already view-aware. (Added
  `.scroll_bind_view_nums`, a numeric sibling of `.scroll_bind_view_cats`.)

## DE: group by multiple columns

* The live-DE panel's **Group by** is now multi-select: naming several categorical
  columns compares their **interaction levels** (e.g. `genotype` + `timepoint` gives
  groups like `"KO | d7"`), and the `ident1`/`ident2` selectors pick which combined
  levels form each side — matching how the Pseudobulk panel already groups.
  `scroll_de(group_col = c("genotype", "timepoint"))` accepts a column vector too.
  A single column behaves exactly as before.

## Live contrast preview (DE + Pseudobulk)

* Both differential panels now show a small **live mini-embedding** of the chosen
  contrast, updating as you pick controls — *before* you Compute. In the DE panel,
  `ident1` cells are drawn **red** and `ident2` (or "rest") **blue** over a faint grey
  outline of the whole UMAP; in Pseudobulk, one coloured point per pseudobulk **sample**
  is placed at that sample's medoid (showing how cells compact into their groups). It's
  pure metadata + coordinates (no store query), stratified-subsampled for speed, and
  driven by the same contrast-membership logic the compute uses, so the picture always
  matches what will be tested. New exported view cores `view_contrast_preview()` and
  `view_contrast_medoids()`.
* Fix: both differential panels' live previews now pick the embedding the way DimPlot
  does — the config `default_embedding` for the whole dataset, else the active subset
  view's primary embedding. Previously Pseudobulk ignored the active view (subset
  medoids landed on the whole-dataset UMAP), and both panels fell back to the *first*
  global reduction for the whole dataset (e.g. `pca` instead of the default `umap`).

## Faster, lighter live DE

* Live DE (`scroll_de()`) reconstructs the contrast's sparse matrix with a
  reverse-index gather instead of `match()` over tens of millions of cell ids, and
  frees the long (feature, cell, value) table **before** allocating the matrix so
  the two never coexist. On a 44k-cell one-vs-rest contrast this cut **peak memory
  ~2.3×** (≈5.9 GB → ≈2.6 GB) with byte-identical results — enough to stop the
  swapping that made DE crawl on a 16 GB server.
* New `scroll_de(max_cells = )` + a **"Max cells / group"** control in the DE panel
  (default off) down-sample each side of a contrast to bound DE time and memory on
  large or many-cell datasets — a standard marker-detection shortcut (cf. Seurat's
  `max.cells.per.ident`); e.g. capping at 5,000/group ran ~4× faster with an
  unchanged top-marker list. The subsample is deterministic and leaves the session
  RNG untouched.
* The DE and Pseudobulk **Compute** buttons are now `bslib::input_task_button`s:
  clicking one immediately turns it into a disabled spinner labelled "Computing…"
  until results return. Previously the only progress cue was an output-area spinner
  that needed the optional `shinycssloaders` package — absent on many deploy
  servers, users saw nothing happen during the multi-second compute and assumed the
  app had frozen. The task button needs only `bslib` (already required) and gives
  clear feedback right where the user clicked.

## Expression store: first-letter bucketing

* The expression store now partitions by the feature's **case-folded first
  character** (`expr/<assay>/bucket=<char>/`) instead of one directory per feature.
  On a 32k-gene store that is **~36 directories instead of ~25,600**, which removes
  the arrow Dataset-open crawl that dominated cold queries and could stall a server
  (fewer inodes / file descriptors), and — because rows are sorted by `feature`
  within each bucket — compresses **~1.8x smaller** with **no loss** (a single-gene
  lookup still reads little via Parquet row-group pruning; benchmarked *faster* than
  per-feature). Case is folded (so `Cd8a` and `ccdc198` share bucket `C`) because
  case-insensitive filesystems collide `C`/`c` directories. The runtime query layer
  is unchanged and reads old per-feature stores and new bucketed stores alike, so
  existing projects keep working un-rebuilt. Differential expression is unaffected
  in results and slightly faster (its whole-store scan opens far fewer files);
  pseudobulk uses the separate single-file counts store and is untouched.

## Manifest: cap cached levels for high-cardinality columns

* `manifest.yaml` no longer inlines **every** distinct value of a categorical column
  as a `levels:` list. Above `scroll_build(max_levels = 200)` (new argument, also on
  `scroll_update()`) a column records only its `n_levels` count; the app recomputes
  the level set from `cells.parquet` on demand. This keeps the manifest small and
  hand-editable when a build carries a clone id / barcode / sample id with thousands
  of values — those columns still export and stay fully usable everywhere (Color-by,
  Violin, DE gating). Manual per-level colour pickers now fall back to a note above
  30 levels instead of rendering a wall of widgets.

## Fixes

* **Non-ASCII feature names.** A feature whose name contains a non-ASCII character
  (e.g. an antibody like `FcεRIa`) broke DE and pseudobulk with `'i' and 'j' must
  not contain NA`, and silently returned no data in FeaturePlot/DotPlot: `yaml`
  escaped the name in the manifest (`Fc<U+03B5>RIa`) so it no longer matched the
  store's UTF-8 feature column. `scroll_de()` / `scroll_pseudobulk_de()` now index
  on the store's own feature names (encoding-robust; no rebuild needed), and the
  manifest is written with `unicode = TRUE` so names round-trip as UTF-8 (fixes the
  query/label path on rebuild).

## Incremental updates

* **`scroll_update(dir, object, embeddings =, meta_cols =, subsets =)`** adds or
  replaces reductions / metadata / subset views on an **already-built** project by
  rewriting only `cells.parquet` + `manifest.yaml`, aligned to the existing cells by
  barcode. The feature-partitioned `expr/` and `counts/` stores are left
  byte-identical, so appending a reduction is seconds rather than a full
  `scroll_build()` re-export (measured: ~2s vs ~8min to add a column to a 47k-cell,
  794MB project). Arguments mirror `scroll_build()`, so the same object flows through.

## Panels

* **Violin** and the ridge view can plot a **numeric metadata column** as the value
  axis (not just a queried feature), via a Gene/Metadata source toggle.
* **Biaxial** can put **genes on the axes** (query 2+ genes from any assay) in
  addition to numeric metadata columns; facets now default to **square** and gain
  **facet-columns / facet-rows** layout controls.
* Fixed a `conditionalPanel` namespacing bug that hid the source-dependent pickers
  (gene vs metadata) in the Violin, Biaxial, and Spatial panels.

## Storage format v2 + streaming builds

* **Leaner store (v2).** The expression store now keys `cell` on an **int32 global
  row-index** into `cells.parquet` (not the repeated barcode string), stores
  unquantized values as **float32** (not float64; ~lossless), compresses with
  **zstd**, and **compacts** to one part-file per feature. On a 4.4M-cell atlas this
  cut the store ~3× with no lossy precision change and faster cold queries. The
  runtime join became simpler and faster (positional integer indexing). Older v1
  (string-cell) stores keep working — a manifest flag (`cell_index`) selects the read
  path, so existing projects need no rebuild.
* **`scroll_build_stream()`** — build one project from many sources *one at a time*,
  so peak memory stays ~per-source instead of loading the whole dataset. Idempotent
  and append-safe (re-run as data arrives). Powers multi-million-cell atlases on a
  laptop.

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

Built-in panels **surface only when the project's data supports them** — a panel
that can't work on a given project simply doesn't appear (rather than showing an
empty-state message). DimPlot and FeaturePlot are always present; Biaxial needs ≥2
numeric columns; DotPlot/Violin need a categorical column; Proportions needs ≥2;
**DE / Pseudobulk require a categorical column with ≥2 levels** (a real contrast),
and Pseudobulk additionally needs a counts store. So, e.g., a single-region Visium
slide shows the Spatial, DimPlot, FeaturePlot, DotPlot, Violin and Biaxial panels
but not DE/Pseudobulk/Proportions, while a multi-region study gets them back.

For explicit control, `scroll_build(panels = c(...))` (or a `panels:` list in
`config.yaml`) shows **exactly** those panels in that order, overriding the automatic
gating; `exclude_panels` / `exclude_panels:` hides specific panels. `config.yaml` is
hand-editable, so the panel set can be changed without rebuilding.

Each built-in analysis panel carries its own fine-grained controls:

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

## Multimodal: scATAC

* **`atac_spec()` + `scroll_build(..., atac =)`** ingest a chromatin-accessibility
  **peaks** assay. Because peaks are just an assay, **FeaturePlot / DotPlot / Violin
  / DE work on accessibility unchanged**; `atac_spec()` additionally marks the assay
  `kind: peaks`, parses each peak's chr/start/end from its `chr-start-end` name, and
  bakes a `peaks.parquet` annotation table (with the nearest gene when a Signac
  `ChromatinAssay` annotation is present).
* **Peaks panel** — find a peak **by nearby gene** (when annotation was baked) or
  **by region** (a text filter over peak names), then colour the embedding by its
  accessibility. Auto-surfaces only for projects built with an `atac` spec.
  Validated on the 10x PBMC multiome scATAC dataset.

## Multimodal: spatial

* **`spatial_spec()` + `scroll_build(..., spatial =)`** ingest 10x Visium /
  imaging-based data. The build extracts `GetTissueCoordinates()` into a `spatial`
  embedding (image-pixel space, y-oriented for ggplot) and bakes the tissue image
  to a small raster asset; the manifest records an `images` block and marks the
  embedding `kind: spatial`. Because coordinates are just an embedding, **DimPlot
  and FeaturePlot work on spatial data unchanged**.
* **Spatial panel** — cells/spots in tissue space with a fixed aspect ratio, coloured
  by a gene or any metadata column, over the tissue image when one was baked.
  **Drag to zoom into a region, double-click to reset.** Auto-surfaces whenever the
  project has a spatial embedding — so it also works for **imaging platforms
  (Xenium / CosMx)** that carry cell centroids but no H&E image (the image toggle is
  simply hidden). Validated on a 10x Visium mouse-brain section (spots overlay the
  H&E exactly) and a Xenium-style FOV object.

## Multimodal: VDJ / immune repertoire

* **`vdj_spec()` + `scroll_build(..., vdj =)`** ingest per-cell TCR/BCR metadata.
  The build bakes a compact `repertoire/` Parquet store (clone table +
  pre-computed diversity, V/J gene-usage residuals, and tissue correlation) and
  records a `vdj` block in the manifest. `chain_type` parameterises TCR vs BCR
  segment names (TRBV/TRAV… vs IGHV/IGKV…); column defaults follow 10x naming.
* **Four repertoire panels** — *Clone overview* (rank-abundance + expansion
  composition), *V/J gene usage* (frequency + chi-square residual heatmap),
  *CDR3 length*, and *Diversity* (Shannon / Simpson / clonality / Gini +
  tissue correlation). They **auto-surface only for projects built with a
  `vdj` spec** and are absent otherwise, via a new `when(manifest)` gate on
  built-in panels — so RNA-only apps are unchanged.

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

* `register_panel(id, ui, server, ...)` adds a panel to every app built
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
* **Data-derived choices.** Controls can now compute options from the data handle
  instead of only the manifest: `scroll_input_choice(choices = function(data))` and
  `scroll_input_levels(choices = function(input, data), watch = c(...))` populate from
  e.g. a baked asset under `data$dir`, recomputing when a watched control changes.
  `scroll_input_column(prefer=)` sets a preferred default column, and
  `scroll_input_palette()` offers the built-in discrete/continuous palette names.

## Documentation

* Reference documentation for every export, published as a
  [pkgdown site](https://shihanli1992.github.io/scroll/).

## Internal

* Split the monolithic `R/app.R` into one file per built-in view (`R/panel-<view>.R`),
  plus `R/panel-helpers.R` (shared panel helpers) and `R/data-handle.R` (the runtime data
  handle). Pure reorganization — no behaviour change.
