# Changelog

## scroll 0.2.22

### Build

- **Pseudobulk on streamed builds.**
  [`scroll_build_stream()`](https://shihanli92.github.io/scroll/reference/scroll_build_stream.md)
  gains `counts = TRUE`, which exports each source’s raw `counts` layer
  under the running global cell index (one part per source in
  `counts/<assay>/`) and records `has_counts`, so the Pseudobulk DE
  panel now works on large multi-source projects. The counts are
  identical to a single `scroll_build(counts = TRUE)` of the same cells.
  An append run must match the project’s original `counts` setting (it
  errors otherwise, so pseudobulk never silently misses sources).

### Fixes

- **Windows:**
  [`scroll_update()`](https://shihanli92.github.io/scroll/reference/scroll_update.md),
  [`scroll_add_meta()`](https://shihanli92.github.io/scroll/reference/scroll_add_meta.md)
  and appending to a streamed build failed with “cannot be performed on
  a file with a user-mapped section open”, because `cells.parquet` was
  still memory-mapped when rewritten in place. Those reads no longer
  memory-map the file, and if something else in the session still holds
  it mapped, the rewrite frees unused mappings and retries once before
  failing with an error saying what to close.
- R CMD check is clean on all CI platforms: non-ASCII characters removed
  from R code, arrow’s `cast()` binding used unqualified, and `Seurat` /
  `cachem` declared in Suggests (used by tests).

### Docs

- Repository, badge and pkgdown links now point at `shihanli92/scroll`,
  and the default branch is `main` (CI runs on it).

## scroll 0.2.21

### Repertoire (VDJ)

- **AIRR / preset column matching is more forgiving.**
  [`vdj_spec()`](https://shihanli92.github.io/scroll/reference/vdj_spec.md)
  now matches a preset’s (or your explicit) column names ignoring case
  when the exact name is absent, so AIRR per-cell columns exported as
  `v_call_vdj` / `junction_aa_vj` satisfy the preset’s `v_call_VDJ` /
  `junction_aa_VJ` (auto-detection recognises them too). With no
  `clone_col` given and the preset’s clone column absent, it looks for a
  common clonotype column (`clone_id`, `strict_clone_id`,
  `clonotype_id`, `raw_clonotype_id`, `clonotype`), and picks up a
  `<clone_col>_count` / `_size` clone-size column. Each match is
  reported with a message. An explicit clone/count column is only
  case-corrected, never swapped; the `"scroll"` convention keeps
  deriving clonotypes from the CDR3s when `clone_col` is `NULL`.
  Previously such a project failed to build with “vdj: no clone_col and
  no CDR3 columns to derive a clonotype from”.

## scroll 0.2.20

### UI

- **App-bar buttons explain themselves on hover.** The dataset-info,
  two-up and controls buttons now show a styled label just below them
  (after a short delay, and on keyboard focus too) saying what they do
  in their current state, e.g. “Show two panels side by side” / “Back to
  one panel per row”, or “Open the filters & theme panel” when the
  controls live in a drawer. Replaces the native `title` tooltips, which
  appeared only after ~1s and were easy to miss. The label is computed
  from live state on hover, so it stays right after a resize changes the
  controls button’s role.

### Subset views

- **[`scroll_add_subset()`](https://shihanli92.github.io/scroll/reference/scroll_add_subset.md)
  works out what a subset changed.** Attaching several reprocessed
  lineages no longer needs per-subset `embeddings`/`meta` lists or
  manual renaming:
  - `embeddings = "auto"` (now the default; `embeddings` was previously
    required) copies every reduction that is new in the child or whose
    coordinates differ from the parent’s on the shared cells
    (i.e. re-run on the subset), as `<name>_<reduction>`, ordered UMAP
    \> t-SNE \> rest so a 2-D map is the view’s primary embedding.
    Inherited, unchanged reductions are skipped; if nothing was re-run
    it stops with a clear error.
  - `meta = "auto"` copies every metadata column that is new or differs
    from the parent on the shared cells (re-run clusters, recomputed
    scores). Numerics compare with a tolerance and factors by value, so
    float noise and re-ordered levels are not flagged. Columns get a
    `<name>_` prefix so a subset’s `seurat_clusters` can’t overwrite the
    parent’s; `prefix =` overrides it, and an explicit `meta = c(...)`
    still copies under the same names.
  - Both report what they picked with a message.
- [`scroll_build()`](https://shihanli92.github.io/scroll/reference/scroll_build.md)
  with inferred `meta_cols` now always exports the subsets’ declared
  columns, instead of failing validation when inference dropped one
  (e.g. a high-cardinality subset label).

## scroll 0.2.19

### API surface trimmed (67 -\> 43 documented exports)

An audit of the exported functions against the docs, tests, and real
downstream apps found many exports that were implementation details. The
public API is now the build / run / query / DE functions plus the
custom-panel extension API.

- **No longer exported** (internal; were only used by the built-in
  panels or tests): `render_view()`, `scroll_view_kind()`,
  `view_de_table()` (read a v1 `de/` folder that v2 stores do not have),
  `view_biaxial()`, `view_volcano()`, `view_stability()`,
  `view_contrast_preview()`, `view_contrast_medoids()`,
  `view_feature_multi()`, `view_feature_blend()`, `view_violin_multi()`,
  `view_violin_stacked()`, `scroll_scaffold_app()` (called by
  [`scroll_build()`](https://shihanli92.github.io/scroll/reference/scroll_build.md)),
  and `scroll_config()` (panels read `data$config`).
- **Deprecated, removed from the export list in 0.3.0.** Calling these
  from your own code now warns once per session; the built-in panels are
  unaffected:
  - [`view_umap_colorby()`](https://shihanli92.github.io/scroll/reference/view_umap_colorby.md),
    [`view_feature_plot()`](https://shihanli92.github.io/scroll/reference/view_feature_plot.md),
    [`view_dotplot()`](https://shihanli92.github.io/scroll/reference/view_dotplot.md),
    [`view_heatmap()`](https://shihanli92.github.io/scroll/reference/view_heatmap.md),
    [`view_violin()`](https://shihanli92.github.io/scroll/reference/view_violin.md),
    [`view_proportions()`](https://shihanli92.github.io/scroll/reference/view_proportions.md):
    use the built-in panels, or build the plot with ggplot2 plus
    [`scroll_point_layer()`](https://shihanli92.github.io/scroll/reference/scroll_helpers.md)
    /
    [`scroll_discrete_colors()`](https://shihanli92.github.io/scroll/reference/scroll_helpers.md)
    /
    [`scroll_continuous_scale()`](https://shihanli92.github.io/scroll/reference/scroll_helpers.md).
  - [`scroll_bind_levels()`](https://shihanli92.github.io/scroll/reference/scroll_bind_levels.md):
    use `scroll_input_levels(from = )`.
  - [`scroll_query_cells()`](https://shihanli92.github.io/scroll/reference/scroll_query_cells.md):
    use
    [`scroll_de()`](https://shihanli92.github.io/scroll/reference/scroll_de.md),
    or
    [`scroll_query_features()`](https://shihanli92.github.io/scroll/reference/scroll_query_features.md)
    for specific genes. (Its `cells` doc now correctly describes the v2
    int32 cell key.)
  - [`scroll_aggregate_counts()`](https://shihanli92.github.io/scroll/reference/scroll_aggregate_counts.md):
    use
    [`scroll_pseudobulk_de()`](https://shihanli92.github.io/scroll/reference/scroll_pseudobulk_de.md).
- **[`scroll_pseudobulk_stability()`](https://shihanli92.github.io/scroll/reference/scroll_pseudobulk_stability.md)
  merged into `scroll_pseudobulk_de(runs = )`.** With `runs > 1` it
  returns the per-gene stability table (`sel_freq`, `median_logFC`,
  `sign_agree`, …), matching the Pseudobulk panel’s stability mode; new
  `lfc` / `padj` arguments set the per-run hit cutoffs. The old name is
  a deprecated wrapper returning the identical table.
- Fixes the pkgdown reference build, which failed on four exported
  `view_*` topics missing from `_pkgdown.yml`.

## scroll 0.2.18

### Robustness

- **Build assays whose layer is a dense matrix.** `GetAssayData(layer=)`
  can return a dense base matrix (some v5 layers, small/ADT assays),
  which has no sparse `@x` slot and failed the build with
  `no applicable method for`@`applied to an object of class "matrix"`.
  The export now coerces any non-sparse `data`/`counts` layer to a
  `dgCMatrix` first.
- **Build large assays past 2 GB of feature strings.** The expr/counts
  store’s `feature` column (one gene name per nonzero) is now written as
  `large_utf8` (int64 offsets). On a large assay the column exceeds 2 GB
  of string bytes, overflowing plain `utf8`’s int32 offsets and failing
  the build with
  `Invalid: Failed casting from large_string to string: input array too large`.
  Reads are unchanged (the `open_dataset` + `feature ==` filter path
  handles either type).
- **Clear error when `arrow` lacks the zstd codec.** scroll reads and
  writes every store with `compression = "zstd"`; a “minimal” arrow
  build (common on HPC, where `install.packages("arrow")` can’t fetch
  the prebuilt libarrow and falls back to a codec-less build) would
  otherwise fail deep inside with a cryptic
  `NotImplemented: Support for codec 'zstd' not built`.
  [`scroll_build()`](https://shihanli92.github.io/scroll/reference/scroll_build.md),
  [`scroll_build_stream()`](https://shihanli92.github.io/scroll/reference/scroll_build_stream.md),
  [`scroll_update()`](https://shihanli92.github.io/scroll/reference/scroll_update.md),
  [`scroll_add_meta()`](https://shihanli92.github.io/scroll/reference/scroll_add_meta.md),
  and
  [`scroll_connect()`](https://shihanli92.github.io/scroll/reference/scroll_connect.md)
  (the read path for
  [`scroll_app()`](https://shihanli92.github.io/scroll/reference/scroll_app.md)/[`scroll_serve()`](https://shihanli92.github.io/scroll/reference/scroll_serve.md))
  now check up front and stop with a message naming the fix (reinstall a
  full arrow: `LIBARROW_MINIMAL=false`, or
  `conda install -c conda-forge r-arrow`).

### Header redesign

- **The app bar is now title-led and fits one row on a laptop.**
  Previously it packed the wordmark, dataset title, version badge, a
  four-metric stats strip and the toggles onto one line — which wrapped
  to two rows (~115px) on most laptops. Reworked:
  - The **dataset title is the dominant element** (bold, ink); the
    `scroll` wordmark is a small muted prefix, and the stat numbers no
    longer out-weigh the title.
  - The **stats strip is gone** — genes/assays/reductions move into an
    **ⓘ dataset popover**; only the live **cell count** stays, beside
    the View/Subset pills (the one value those controls change).
  - The **scroll version** moves from a bar badge to the section-rail
    foot.
  - The three actions (dataset info · two-up · controls) are a
    **right-hugging cluster** with real Font Awesome icons (info /
    table-columns / sliders) instead of font glyphs, a pressed state,
    and tooltips that update with state. Result: the bar is one row down
    to ~1000px (≈65px vs the old ~115px two rows), and no longer grows
    with the number of subset views.

### UI

- **Denser control column.** A tighter, slightly narrower control column
  (`--sc-ctl-w` now `clamp(232px, 22%, 320px)`, with smaller group/field
  gaps) so the inputs don’t sprawl on a wide card and the plot gets more
  room — without introducing any horizontal scroll in the controls.
- **`layout:` config block.** `config.yaml` can now tune the layout
  without editing CSS: `content_max` (wide-screen content cap),
  `rail_width`, `control_width` (a number → px, or a CSS string), and
  `two_up` (start with two panels per row by default; a per-browser
  toggle still overrides it).

## scroll 0.2.16

### UI (wide screens)

- **Full-bleed wide screens.** The layout no longer caps at a fixed
  centered width leaving big edge margins: the section rail (and filters
  column) hug the screen edges while the content fills the middle up to
  `--sc-content-max` (1600px), so a wide monitor gives the plot much
  more room (e.g. ~660px → ~1140px on a 2560px screen) — capped so a
  scatter/UMAP never grows to an unusable width.
- **Tighter chrome.** A narrower section rail (`--sc-rail-w` 160px, long
  labels ellipsize with a hover tooltip) and a capped text measure on
  panel titles/descriptions, so text no longer floats across a very wide
  card.
- **Two panels per row (wide screens).** A toggle in the app bar (shown
  when the screen is wide enough) packs two panel cards side by side; it
  self-degrades to one column when there isn’t room, remembers the
  choice across reloads, and raises the content cap so two panels fit.
  Each card is now its own layout container, so a half-width card still
  lays its controls out sensibly. When two-up is on, the global
  filters/theme rail becomes an off-canvas drawer (reachable from the
  controls button) so it doesn’t steal the width the second panel needs
  — two-up now works even with the filters column up.
- **Drag to reorder panels.** Each rail item shows a grip of raised dots
  to signal it’s movable; drag one and the list reflows live (siblings
  slide as the item is pulled out and re-inserted). The panel cards
  follow, everything renumbers, and the order persists across reloads (a
  “Reset order” link restores the default). Reordering never re-renders
  or rebinds a panel. (Disabled on the narrow horizontal rail.)

## scroll 0.2.15

### UI

- **Responsive panel layout.** Each panel now adapts to the width its
  card actually has (a container query on the panel area), so a
  minimized window or a portrait/sideways monitor no longer squeezes the
  plot to a sliver:
  - Wide enough (card ≥720px): controls sit beside the plot in a width
    that scales with the screen
    (`--sc-ctl-w = clamp(240px, 24%, 360px)`) — a usable minimum on
    laptops, capped on large monitors so the plot takes the extra room
    instead of an over-wide control column.
  - Narrower (card \<720px): controls **stack above** the plot in a
    compact, capped multi-column box, and the plot spans the full card
    width.
- **Bounded plot height.** The default plot height is now
  `clamp(320px, min(100vh−190px, 90cqw), 1100px)` — a tall/portrait
  monitor no longer stretches a plot to ~1700px; landscape is unchanged.
- **App bar wraps** and ellipsizes a long dataset title instead of
  overflowing; a **compact numeric section rail** (901–1199px) frees
  width for the plot, with hover labels.
- **Filters/theme rail becomes an off-canvas drawer** below 1400px (the
  app-bar toggle slides it in; Escape or a click outside closes it), so
  the plot keeps the full width and the controls stay reachable. The
  drawer has an explicit height so a long filter/theme list **scrolls**
  within it, is narrower (300px), and lays its theme controls out one
  per row. Sticky offsets (rail, drawer, scroll-to-panel) now track the
  app bar’s measured height (`--sc-appbar-h`), so they stay correct when
  it wraps or under a multi-app tab strip.
- **Per-panel Controls toggle.** When a panel is stacked (narrow), a
  small **Controls** button in the card header collapses that panel’s
  controls so the plot alone is visible; it is hidden when the panel is
  side-by-side. Defaults open, so inputs always initialise while
  visible.

## scroll 0.2.14

### Heatmap & DotPlot

- **Grouped/aggregated heatmap is now a DotPlot display.** The DotPlot
  panel’s **Display** toggle draws the same genes x groups data as dots
  or as heatmap **tiles**; the standalone **Heatmap** panel is
  single-cell (genes x a subsample of individual cells) only.
- **Mark a subset of genes with leader lines.** On dense heatmaps (row
  labels auto-hide past ~60 genes), a **Label genes** control draws
  ComplexHeatmap `anno_mark`-style elbow leaders on the right of both
  the DotPlot-tiles and single-cell heatmaps. Labels sit at their true
  gene row and are nudged apart only enough to avoid overlap, so the
  leaders stay short.
- **Group-by is optional** on the single-cell heatmap — leave it empty
  for one ungrouped block.
- **Cell ordering by a metadata column**, in addition to grouped and PC1
  orderings.
- **Compute-gated assembly.** The single-cell heatmap’s expensive step
  (query + subsample + PC1) runs only on a **Compute** button
  (auto-rendering once on first view); gene/scale/cluster/cap/ order
  tweaks stage without recomputing, while the cheap cosmetics
  (palette/clip/legend/label/ aspect) re-skin the memoized assembly
  live. An app-bar filter/View change refreshes it.

### Robustness (code review)

- Per-gene z-scoring (`.scroll_row_zscore`) is vectorised and
  shape-preserving — a single-level group column (e.g. a one-cluster
  filter) no longer collapses the matrix and breaks the tiles.
- The single-cell heatmap builds its genes x cells matrix in one O(nnz)
  scatter (was one full expr scan per gene), errors clearly on an empty
  cell selection, and never mis-shapes on 1 cell.
- The Heatmap **Legend** switch now works, and global right-rail
  **theme** overrides apply even when a dendrogram or gene-mark panel is
  attached (previously silently dropped).
- `view_dotplot` / `view_proportions` give a clear error (instead of an
  R length-coercion error) when `group_by`/`fill_by` is missing or
  multi-column.
- One shared `.scroll_cell_key()` for the v1/v2 join; violin jitter
  subsampling is now seeded (stable across redraws and export).
- Stacked violins now join expression through the global cell index, so
  they stay correct under an app-bar subset/filter (previously the v2
  join could misalign).
- Internal: the gene-box observers (assay repopulate + paste-list +
  long-list warning) and the view-aware multi-select binders are
  factored into `.scroll_bind_gene_box()` /
  `.scroll_bind_view_cats(multiple=)`, and the two violin views share
  `.scroll_violin_base()` / `.scroll_violin_fill()` — ~50 fewer lines
  across the panels, one place per behaviour.

## scroll 0.2.13

### New: Heatmap panel

- A built-in **Heatmap** panel (`view_heatmap`). Two modes:
  - **Groups (aggregated):** genes x groups tiles of mean or %
    expressing, optional per-gene z-scoring + clip, row/column
    clustering with dendrograms, diverging palette. Aggregation
    collapses the cell dimension, so it is flat-RAM at any cell count.
  - **Cells (subsampled):** genes x a random, per-group-proportional
    subsample of cells (capped for safety), rasterized and faceted by
    group; a **gene dendrogram** and a cheap **PC1 cell ordering**
    (“similar cells together”, no O(n^2) tree) are available.
- Multi-column group-by, paste-a-gene-list, assay selector, and CSV
  export (per-group means).
- Projects can set a `heatmap_markers:` config key (e.g. top variable
  genes) as the panel’s default gene set, independent of DotPlot’s
  `markers:`.
- The group-aggregation used by DotPlot and Heatmap is now a shared
  `rowsum` helper (`.scroll_group_expr_matrix`), ~3x faster than the old
  [`aggregate()`](https://rdrr.io/r/stats/aggregate.html) path.

### Layout

- **Plots fill the viewport.** All panels now use a viewport-relative
  plot height (`calc(100vh - 190px)`) by default – including custom
  panels built with
  [`register_plot_panel()`](https://shihanli92.github.io/scroll/reference/register_plot_panel.md)
  – so they use the available space and re-render on window resize.
- **Control columns scroll** within the viewport instead of stretching
  the card, and a multi-select with many chips (e.g. a 100-gene default)
  scrolls instead of growing unbounded.

## scroll 0.2.12

### Composition aesthetics

- **Order groups by a fill level’s share** – e.g. sort samples by their
  % of a chosen category, with a picker for which level.
- **Segment-label threshold** – hide labels on segments below a chosen
  percent (declutters small slices), plus a **label size** control.
- **Per-group totals** – the n cells for each group shown above its bar
  (placed correctly for stack / fill / dodge).

(The first aesthetics batch – ordering, count/percent labels, horizontal
bars, bar outline – landed under 0.2.11.)

## scroll 0.2.11

### Violin & composition panels

- **Violin: grouped, stacked, and cleaner points.**
  - **Split-by** draws side-by-side coloured violins within each group
    (Seurat `split.by`).
  - **Group-by** is multi-select – group by the `a | b` interaction of
    several columns.
  - **Stacked violin** (`view_violin_stacked`): one compact faceted row
    per gene sharing the group axis (points off in this mode).
  - **Point controls** when points are shown: size, opacity, and a
    subsample slider for a cleaner view.
  - **Violin width** slider.
- **Composition: bar shape and multi-fill layout.**
  - **Bar width** slider.
  - **Fill by** multiple columns can now be **combined** into one
    interaction plot or **faceted** – as a grid or **stacked rows**
    (like the stacked violin).
  - **Aesthetics:** bar outline width, group/fill ordering (alphabetical
    / total or abundance / reverse), segment labels (count or percent),
    and horizontal bars.
  - (Bar position stack/fill/dodge landed in 0.2.10.)

## scroll 0.2.10

### DE & composition panels

- **Volcano plots the Seurat-style fold change.** The DE (Wilcoxon)
  volcano now uses `avg_log2FC` (log2 fold change) on the x-axis and for
  the up/down threshold, instead of presto’s natural-log `logFC` – more
  interpretable and consistent with the DE table. The Pseudobulk volcano
  is unchanged (its `logFC` is already log2). `view_volcano()` gained
  `fc_col`/`fc_label` params (default `logFC`); non-finite fold changes
  are dropped.
- **Composition: fill by multiple columns.** The Proportions panel’s
  “Fill by” is now a multi-select that fills by the `a | b` interaction
  of the chosen categorical columns – the same pattern as DimPlot’s
  colour-by.
  [`view_proportions()`](https://shihanli92.github.io/scroll/reference/view_proportions.md)
  and its CSV export build the interaction; Manual colours use the
  observed interaction levels.

## scroll 0.2.9

### Pseudobulk DE: sounder replicate handling

- **Pseudo-replicates are now a disjoint partition.**
  [`scroll_pseudobulk_de()`](https://shihanli92.github.io/scroll/reference/scroll_pseudobulk_de.md)
  splits a group’s cells into `n_pseudo` *non-overlapping* samples
  instead of drawing overlapping random subsets. The previous draws
  could return identical samples whenever a group had
  `<= cells_per_pseudo` cells (within-group variance 0 -\> spurious
  significance); the partition makes the variance real. A group must now
  hold at least `n_pseudo * min_cells` cells or the call errors clearly;
  `min_cells` is a per-sample floor and `cells_per_pseudo` a per-sample
  cap.
- **Replicate-aware design.** With a real `replicate_col`, each
  replicate now contributes **one sample per group** (summing across the
  combinations it spans) rather than one per combination – removing the
  pseudo-replication of biological replicates. The model uses the
  no-intercept (means) parameterization `~ 0 + group`, or
  `~ 0 + group + replicate` when \>= 2 replicate levels are shared
  across both groups (a **paired** design); the group1-vs-group2 effect
  is tested as the `group1 - group2` contrast
  ([`limma::contrasts.fit`](https://rdrr.io/pkg/limma/man/contrasts.fit.html)).
  A new `paired = c("auto","yes","no")` argument controls the blocking.
- **No more silent fallbacks.** A named `replicate_col` that is absent
  is now an error; cells with a missing replicate value are counted and
  reported (a warning + a `dropped_na` attribute) instead of being
  dropped silently; a group that lacks real replicates falls back to
  pseudo-replicates as a reported “mixed” run.
- **Honest panel labelling.** The Pseudobulk panel’s note now reflects
  the actual regime (`real-paired` / `real-unpaired` / `mixed` /
  `pseudo`) with the sample counts, and the controls document the
  partition requirement.
- Result attributes gained `regime`, `design`, `sample_sizes`, and
  `dropped_na`.
- **Note:** these are deliberate statistical corrections, so pseudobulk
  results computed before 0.2.9 will change (previously inflated
  significance becomes more conservative).

### App

- The app bar shows the scroll package version (e.g. `v0.2.9`); hide it
  with `show_version: false` in a project’s `config.yaml`.

## scroll 0.2.8

### Startup warm-up & on-screen rendering

- **No more post-overlay cycling.** The startup view warm-up now
  advances as each view’s render *completes* (event-driven on the view
  echo) instead of on a fixed timer, so the DimPlot no longer visibly
  cycles through subset views after the “Preparing views” overlay
  clears. A per-step watchdog aborts and hides the overlay if a step
  stalls.
- **Panels gate correctly at first load.** The on-screen render gate is
  now a Shiny input binding, so each panel’s on-screen state is known
  *before* the first flush – only panels in view render up front,
  off-screen ones wait until scrolled to (previously all live panels
  rendered at once on load).
- **Overlay polish** – a labelled first stage (“Preparing the main
  view…”), a progress bar, page-scroll lock while warming (so scrolling
  can’t pollute the plot cache mid-cycle), and a fade-out. The client
  failsafe now fires on server *silence* (rearmed on each step) rather
  than a fixed 60s wall-clock, so a long but healthy warm-up is never
  cut short.
- **`prewarm_views` is an on/off flag.** `prewarm_views: true` (or any
  positive number) warms all subset views; `false`/`0`/absent is off.
  The number never was a per-view cap.
- **[`scroll_preview_panel()`](https://shihanli92.github.io/scroll/reference/scroll_preview_panel.md)**
  no longer shows a stranded overlay on a warm-up-configured project (it
  runs no warm-up), and
  **[`scroll_multi_app()`](https://shihanli92.github.io/scroll/reference/scroll_multi_app.md)**
  now caches rendered plots (per-dataset cache keys) like
  [`scroll_app()`](https://shihanli92.github.io/scroll/reference/scroll_app.md).

## scroll 0.2.7

### DE: Seurat-compatible `avg_log2FC` column

- [`scroll_de()`](https://shihanli92.github.io/scroll/reference/scroll_de.md)
  now also returns **`avg_log2FC`** alongside presto’s `logFC`. presto’s
  `logFC` is a natural-log “mean of log” difference; `avg_log2FC` is the
  `log2` of the mean of the un-logged normalized counts (Seurat v5
  `FoldChange`, pseudocount 1), so it matches
  [`Seurat::FindMarkers`](https://satijalab.org/seurat/reference/FindMarkers.html)
  — to floating-point tolerance on a `quantize = FALSE` build, and
  within ~0.01 on the default quantized store. (presto’s `logFC` is
  generally smaller in magnitude on zero-inflated data, which is why
  scroll’s fold changes looked smaller than Seurat’s.)

### Faster view switches (cached + downsampled scatter renders)

The DimPlot/FeaturePlot scatters no longer re-draw from scratch on every
view/control change — three layers, measured on a 237k-cell project:

- **Cache** — rendered images are memoized in Shiny’s shared app-level
  cache keyed on the active cells + controls, so returning to a (view +
  controls) state serves the PNG without re-drawing (~ms vs a ~300-700
  ms draw). On by default; `cache_plots: false` disables. The shared
  cache also self-warms across users.
- **On-screen downsampling** — above `options(scroll.onscreen_cap=)`
  (default 50k) the interactive rasterized draw uses a deterministic
  subsample (a 237k embedding renders ~2-3x faster); highlighted /
  expressing cells are always kept, and exports draw every point.
- **Startup warm-up** — set `prewarm_views: true` in `config.yaml` (any
  positive number also works; off by default) to cycle **all** subset
  views once at load so the first visit to each is instant too. It runs
  behind a progress overlay (shown from the initial page, so it covers
  the whole startup), steps as each view’s render completes, locks page
  scroll while it runs, and fades out when done.
- **Off-screen panels stay off** — the on-screen render gate is now
  wired reliably (the IntersectionObserver is set up on
  `shiny:connected`, not before Shiny is ready), and the live
  declarative panels (repertoire panels, clone-map) are lazy-gated too,
  so a global filter / View change no longer re-draws every
  scrolled-away panel and stalls the switch.

## scroll 0.2.6

### DimPlot: colour by multiple columns + configurable highlight background

- The DimPlot **Color by** control is now multi-select: pick more than
  one categorical column to colour cells by their `a | b` interaction
  (e.g. `donor | antigen`), mirroring the differential-expression “group
  by”. The **Highlight** box then lists the composite levels, so you can
  highlight a specific combination (e.g. `Donor 1 | pp65_CMV`).
- When a highlight is active, a **Background colour** picker sets the
  colour of the non-highlighted cells (previously a fixed grey), and the
  **Manual** palette now shows a colour picker only for the highlighted
  group(s) rather than every level — so you can recolour the few groups
  in focus even on a high-cardinality (or composite) column.

### Add a metadata column with no Seurat object

- New `scroll_add_meta(dir, name, values)` writes a column into a built
  project’s `cells.parquet` + `manifest.yaml` directly — no original
  Seurat object required, so you can add a column to an already-deployed
  app. `values` may be a length-`n_cells` vector, a barcode-named
  vector, a `data.frame(cell, value)`, or a `function(cells)` that
  derives the column from existing columns. Categorical vs numeric is
  inferred (as in `scroll_build`), and `scope =` scopes it to a subset
  view. Restart the app to pick it up. (Use
  [`scroll_update()`](https://shihanli92.github.io/scroll/reference/scroll_update.md)
  when you *do* have the object and want to pull
  columns/embeddings/subsets from it.)

### Expression store: one Parquet file per assay

- [`scroll_build()`](https://shihanli92.github.io/scroll/reference/scroll_build.md)
  now writes the expression store as a single feature-sorted file
  `expr/<assay>/part-0.parquet` instead of the former first-letter
  buckets (`expr/<assay>/bucket=<CHAR>/…`). A benchmark on a 237k-cell
  store (317M nonzero rows) showed the bucketing gave no retrieval
  advantage — a single-gene lookup is pruned by Parquet **row-group
  min/max statistics**, which work identically in one sorted file —
  while the single file plans slightly faster and is simpler on disk.
- The runtime read path is unchanged
  ([`arrow::open_dataset`](https://arrow.apache.org/docs/r/reference/open_dataset.html) +
  a `feature` filter is layout-agnostic), so **existing bucketed
  projects keep working with no rebuild**, and new single-file builds
  run on any 0.2.x runtime. Streaming builds
  ([`scroll_build_stream()`](https://shihanli92.github.io/scroll/reference/scroll_build_stream.md))
  still append one part per source and compact them into a single
  feature-sorted file per assay at the end of the run.

## scroll 0.2.5

### New Signature panel

- A new always-on built-in panel (**Signature**, after FeaturePlot)
  scores a gene **signature** per cell live — type or paste a gene list
  (case-insensitive, unknown symbols reported) and colour the embedding
  by the score, or draw a violin split by any categorical column. Two
  runtime scoring methods: **Mean** (mean of the genes’ log-norm
  expression, absent cell = 0), **Scaled** (z-score each gene across the
  shown cells, then average), and **AddModuleScore** — a faithful
  runtime port of Seurat’s `AddModuleScore` (signature mean minus an
  expression-matched control mean; validated at Pearson r ≈ 0.999
  against Seurat on real data). No rebuild needed: AddModuleScore’s
  per-gene background means are computed by a one-time full-store scan
  cached on the connection handle, and the other methods touch only the
  signature’s genes.

### Read VDJ from scRepertoire / AIRR / Platypus

- [`vdj_spec()`](https://shihanli92.github.io/scroll/reference/vdj_spec.md)
  gains a **`source`** argument (default `"auto"`) so the shipped
  repertoire panels read TCR/BCR data straight from the common upstream
  tools — no hand-written column map needed. `"auto"` detects the
  convention from the metadata columns; you can also name it:
  `"scRepertoire"` (parses the compound `CTgene`/`CTaa`/`CTstrict`
  columns into per-segment/per-chain columns), `"airr"` (dandelion
  per-cell obs — `v_call_VDJ`/ `v_call_VJ`, `junction_aa_VDJ`/`_VJ`,
  `clone_id`), `"platypus"` (`VDJ_`/`VJ_` columns), or `"scroll"` (the
  historical `v_gene_TRB`/`cdr3_beta` default). Any column you pass
  explicitly still overrides the preset. Long per-contig tables (raw 10x
  `filtered_contig_annotations.csv` / long AIRR `.tsv`) should be
  collapsed to per-cell upstream first. The resolved source is recorded
  in the manifest `vdj` block.
- Build-time ingestion is more robust: a mapped-but-absent segment/CDR3
  column now **warns** (naming it) instead of silently dropping, and the
  missing-clonotype sentinel recognises more conventions (`NA|NA|NA`,
  `None`, …), not just `NA`/`NA|NA`.

### V/J gene usage panel

- New **Pairing** view: a segment × segment co-occurrence heatmap
  (e.g. TRBV × TRBJ), faceted by the group-by, with a continuous colour
  palette + min/max quantile colour cut-offs. Absent pairings render as
  0 (the low colour), not blank tiles.
- The **Chi-square residuals** view is recomputed at runtime, so it
  re-groups by any categorical column and honours the global cell filter
  / subset view (previously locked to the build-time group).
- The **Frequency** view can **deduplicate expanded clones** (count
  clones, not cells) and orders genes in **genomic** (natural IMGT)
  order on the axis, not alphabetically.

### New CDR3 logo panel

- A new repertoire panel (**CDR3 logo**, gated on `vdj` + the
  `ggseqlogo` package): amino-acid sequence logos of the CDR3 at a
  chosen length, one per group (faceted). Controls for chain, CDR3
  length (defaults to the modal), group-by, cell vs clone-deduplicated
  counting, and bits vs probability height. Needs the CDR3 sequences,
  which
  [`scroll_build()`](https://shihanli92.github.io/scroll/reference/scroll_build.md)
  now bakes into `repertoire/` (`cdr3_<chain>` columns) alongside the
  lengths.

### New Clone-map panel

- A new repertoire panel (**Clone map**, gated on `vdj`): an interactive
  per-clone table (clone id, size, group, V/J genes, CDR3 lengths —
  searchable/sortable via `DT`) over a greyed embedding. Selecting
  clones in the table highlights their cells on the UMAP, one colour per
  clone. Honours the global cell filter / subset view; CSV + image
  export.

### Optional group-by across the VDJ panels

- Group-by / Colour-by now defaults to **None**: each VDJ panel (clone
  overview, gene usage, CDR3 length, diversity) opens as a single pooled
  plot until a column is chosen. Chi-square still requires a group.
  Ungrouped plots honour the palette / Manual **Colour** control (a
  single pooled colour).

### Plot toolbar

- Every panel’s plot toolbar gained an **export-scale slider** (1×–5×,
  with a live value readout) next to the download buttons. It feeds
  `ggsave(scale =)`, so a larger value exports a bigger figure without
  changing the on-screen plot. Custom `register_plot_panel` panels
  inherit it by default.

### Cosmetics

- Continuous value axes sit **flush on the axis** (zero lower-end
  expansion) across the VDJ value plots, Violin, and Proportions.

## scroll 0.2.4

### Collapsible control rail

- A toggle in the app bar collapses/expands the entire right-hand
  control rail (theme + filters), giving the plots full width — useful
  on narrow screens where the rail would otherwise crowd them out. Pure
  client-side (no redraw).

### Global theme controls

- The right control rail gained a comprehensive **Theme** section,
  organised into collapsible sub-sections: **Text & fonts**
  (base/title/axis/legend/strip sizes, font family, text colour, title
  style), **Legend** (position, direction, title show/hide, key
  background), **Axes** (text/titles/ticks/lines show-hide, a shared
  axis colour for lines + ticks, x & y label angle), **Panel**
  (major/minor gridlines, gridline colour, line thickness & colour,
  border + colour, panel & plot background), and **Facets & spacing**
  (strip background, plot margin). Colour options are **colour pickers**
  shown as small circular swatches (`colourpicker`, with a hex
  text-field fallback); an empty picker means “no override”. The
  controls lay out two-per-row in a widened rail, and are **applied on
  an Apply button** (with Reset) so the plots don’t redraw on every
  tweak. Every option is a ggplot
  [`theme()`](https://ggplot2.tidyverse.org/reference/theme.html)
  element applied to *every* plot (RNA views and the VDJ/spatial/ATAC
  panels alike) via a shared `.scroll_ggtheme()` override threaded
  through a `theme_r` reactive. Each control defaults to “Default” (a
  no-op), so rendering is unchanged until you pick something. Disable
  with `theme_controls: false` in `config.yaml`. The declarative panel
  builder now also threads `theme_r`, so custom
  [`register_plot_panel()`](https://shihanli92.github.io/scroll/reference/register_plot_panel.md)
  panels inherit the global theme.

### Global filter rail

- A right-hand **Filters** sidebar narrows the cells *every* panel sees.
  Controls are auto-generated from the manifest — a level multi-select
  for each categorical column and a range slider for each numeric column
  — and compose with AND. A high-cardinality column (above the level
  cap, e.g. a clone id) is skipped. Filters are **applied on an Apply
  button** (with Reset), so dragging sliders doesn’t re-narrow every
  panel until you commit. Includes a cell-count readout (“N of total”).
  Configure via `config.yaml`: `filters: false` hides the rail,
  `filters: [col, col]` curates (and orders) which columns appear.
  Composes with subset views and the existing app-bar subset filter.

### Multimodal panel overhaul (VDJ / spatial / ATAC)

The modality panels gained depth, most improvements working on
**existing built projects with no rebuild** (they read columns already
baked into the store):

- **CSV export on every modality panel.** The declarative panel builder
  ([`register_plot_panel()`](https://shihanli92.github.io/scroll/reference/register_plot_panel.md),
  and the built-in VDJ/ATAC panels) now takes `csv = TRUE`; a panel opts
  in by attaching its source table as `attr(p, "scroll_source")`, and a
  CSV button appears alongside PNG/PDF.
- **VDJ — runtime re-grouping.** Clone-overview rank-abundance and V/J
  gene-usage frequency now group by **any** baked categorical column
  (group / antigen / tissue / location / carried columns), not just the
  build-time `group_col`. The gene-usage group control is hidden in the
  chi-square view (which stays baked-group).
- **VDJ — clone-id column picker.** Clone overview, Diversity, and CDR3
  length gained a **Clone ID** control that chooses which column defines
  a clone — the baked `clone_id` or any carried high-cardinality
  alternate (e.g. a nucleotide- vs amino-acid-level CDR3 clonotype).
  Clone sizes, rank-abundance, expansion categories, diversity metrics,
  and the per-clone CDR3-length dedup all recompute under the chosen
  definition.
- **VDJ — group filter + group colour picker on all group panels.**
  Clone overview, V/J gene usage, CDR3 length, and Diversity each gained
  a “Groups” levels selector (restrict to a subset of the group column’s
  levels; empty = all) and a “Colour” control — a colourblind-safe
  palette dropdown plus a **Manual** option that renders one colour
  input per active group level (tracking the chosen group column +
  filter).
- The declarative panel builder now passes `output` to a control’s
  `bind` (when its formals declare it), so a control can render a
  dynamic `uiOutput` — used by the VDJ group colour picker’s per-level
  Manual swatches.
- **VDJ — group-by menus now exclude clone-scale columns.**
  Grouping/colour menus offer only lower-cardinality categorical
  columns; a high-cardinality clone column (tens of thousands of levels)
  is kept out of the group-by axes and offered as a Clone ID instead.
  Clone-id and grouping candidates are complementary (split by
  cardinality relative to the baked `clone_id`).
- **VDJ — Diversity is recomputed at runtime**, so it can group by any
  categorical column (not just the baked group / cluster) and exposes
  `paired_rate` as a metric; the metric / group-by controls are hidden
  in the tissue-correlation view.
- **Spatial — multi-FOV / multi-slide.** `scroll_build(spatial = )`
  accepts a **list** of
  [`spatial_spec()`](https://shihanli92.github.io/scroll/reference/spatial_spec.md)s
  (distinct `name`s), baking several tissue maps. The Spatial panel adds
  a tissue-map selector (when \>1 exists), continuous + categorical
  palette pickers, view-scoped metadata choices, and a CSV of the
  plotted coordinates.
- **ATAC — genomic-interval region search.** A `chr:start-end` query now
  selects peaks by coordinate **overlap** (substring match remains the
  fallback for other text); the peak list cap is raised and truncation
  is surfaced in the plot subtitle, not silent.

## scroll 0.2.3

### Fix: Pseudobulk “Aggregate by” / “Replicate” follow the active View

- The Pseudobulk DE panel built its **Aggregate by** and **Replicate**
  column menus once from the whole-dataset columns, so a subset View’s
  own scoped categorical columns (e.g. a re-clustered subset’s
  resolutions) were missing. Both menus now refresh with the active View
  — scoped columns appear in their view, global columns stay available
  everywhere — matching the DE panel.

## scroll 0.2.2

### Fix: DimPlot/FeaturePlot render once per View switch (was up to 4×)

- Switching the app-bar **View** used to redraw DimPlot/FeaturePlot
  several times: the active cells changed, then the view-driven
  reduction / colour-by selectors round-tripped through the client and
  each re-render fired again. The effective reduction and colour column
  are now deduped `reactiveVal`s resolved server-side (updated at high
  priority, before the plot renders), so a View switch redraws the panel
  exactly once. Measured 4 renders → 1 on the demo.

## scroll 0.2.1

### New: lazy on-screen rendering (snappier View/subset switching)

- The built-in plot panels (DimPlot, FeaturePlot, Biaxial, DotPlot,
  Violin, Proportions) now **recompute only while on screen**. An
  `IntersectionObserver` reports each panel’s visibility to its module,
  and the panel keeps its last render while off screen, refreshing when
  scrolled back into view. Switching the app-bar **View** (subset) or
  cell filter therefore redraws just the panels in view instead of every
  panel on the page — the on-screen panel updates immediately and the
  rest catch up as you scroll, rather than the whole page re-rendering
  serially at once. The build is memoized, so scrolling a panel out of
  view and back does not rebuild it unless its data or controls actually
  changed while it was away.

## scroll 0.2.0

Bigger multi-gene panels: the FeaturePlot, DotPlot and Violin panels now
take many genes at once (grids and co-expression blend), with
paste-a-list gene entry.

### New: multi-gene and co-expression blend in the FeaturePlot panel

- The FeaturePlot gene box is now **multi-select**, and the
  numeric-metadata selector alongside it too. Pick any mix of genes and
  numeric columns and the panel lays out a **grid of feature plots**,
  one per feature, each with its own colour scale (as
  `Seurat::FeaturePlot(features = c(...))` does), up to 12. New exported
  view core `view_feature_multi()`.
- A **Blend** toggle reproduces `Seurat::FeaturePlot(blend = TRUE)` for
  the two selected genes: four views — each gene alone, their
  co-expression blend, and a 2-D colour key — with an adjustable **Blend
  threshold** and pickers for the two gene colours. The blend colours
  are a faithful, Seurat-free port of Seurat’s
  `BlendMatrix`/`BlendExpression` (byte-identical colours). New exported
  view core `view_feature_blend()`. The CSV export includes every
  plotted gene.
- Both features use the (Suggested) `patchwork` for their grid layout;
  without it the panel degrades to a single plot.
- FeaturePlot’s **colour quantile** cutoffs are now taken over the
  expressing (non-zero) cells. Single-cell expression is zero-inflated,
  so a quantile of the full vector stayed pinned at 0 until the fraction
  passed the (often \>90%) zero share — the lower cutoff had no visible
  effect. It now spans the real expression range.

### New: multi-gene Violin, and paste a gene list into DotPlot/Violin

- The Violin gene box is now **multi-select**: pick or paste several
  genes and the panel draws a grid of violins, one per gene (each with
  its own expression axis), up to 12. New exported view core
  `view_violin_multi()`.

- Paste a gene list **straight into the DotPlot or Violin gene box** —
  separated by spaces, tabs, commas or newlines — and the genes are
  added to the selection in the pasted order (typing one at a time still
  works in the same box). Symbols are matched case-insensitively and any
  that aren’t in the assay are reported.

## scroll 0.1.0

First tagged release. `scroll` turns a processed Seurat object into a
polished, flat-RAM interactive single-cell explorer.

### New: reproducible CSV export for the built-in plot panels

- Every built-in plot panel (DimPlot, FeaturePlot, Biaxial, Violin,
  DotPlot, Proportions) now offers a **CSV** download beside its PNG/PDF
  buttons, containing the minimal source data behind the figure. For the
  per-cell panels the CSV leads with the **cell barcode** plus the
  plotted columns (embedding coords + colour-by / expression / group),
  so an exported point maps back to a tracked cell; the two aggregated
  panels export their group×feature (DotPlot: `avg_expr`,
  `pct_expressing`) and group×category (Proportions: `n_cells`,
  `proportion`) summary. The CSV reuses the exact data the plot draws,
  so it reproduces the figure. Custom panels built with
  `.scroll_plot_area(csv = TRUE)` + `.scroll_plot_downloads(csv_r = )`
  get the same.

### New: `scroll_preview_panel()` — fast custom-panel dev loop

- `scroll_preview_panel(panel, dir)` mounts **one** panel (by id) in the
  real app harness — the app bar with the cell-subset filter and
  subset-View selector, but none of the other panels — so a custom panel
  can be driven for quick feedback without building or scrolling through
  the whole app. Because it reuses the actual page + wiring, the panel
  is tested exactly as it will run (under the app-bar filter and subset
  views). The loop is [`source()`](https://rdrr.io/r/base/source.html)
  the panel file →
  [`scroll_preview_panel()`](https://shihanli92.github.io/scroll/reference/scroll_preview_panel.md)
  → edit → re-source → re-run.

### Correctness fixes

- **[`scroll_update()`](https://shihanli92.github.io/scroll/reference/scroll_update.md)
  no longer re-escapes non-ASCII feature names.** It rewrites the
  manifest, but wrote it without `unicode = TRUE` (unlike
  [`scroll_build()`](https://shihanli92.github.io/scroll/reference/scroll_build.md)),
  so an update on a project with a non-ASCII feature (e.g. an antibody
  `FcεRIa`) re-mangled the name to `<U+XXXX>` and broke its query/label
  match. Now consistent with the build path.
- **Pseudobulk DE is RNG-neutral and reproducible.**
  [`scroll_pseudobulk_de()`](https://shihanli92.github.io/scroll/reference/scroll_pseudobulk_de.md)’s
  random pseudo-replicate draw used the session’s global RNG as a side
  effect and was non-reproducible. It now seeds the draw (new `seed =`
  argument, default `1L`) and restores the caller’s `.Random.seed`,
  mirroring
  [`scroll_de()`](https://shihanli92.github.io/scroll/reference/scroll_de.md)’s
  cap.
  [`scroll_pseudobulk_stability()`](https://shihanli92.github.io/scroll/reference/scroll_pseudobulk_stability.md)
  still varies the draw per run (distinct per-run seeds), so it stays
  reproducible as a whole without perturbing the session.
- **[`scroll_build_stream()`](https://shihanli92.github.io/scroll/reference/scroll_build_stream.md)
  validates each source up front.** A source missing a declared assay /
  embedding / metadata column now errors with the offending source
  named, instead of failing deep in the final cells-frame merge.
- **[`scroll_de()`](https://shihanli92.github.io/scroll/reference/scroll_de.md)**
  emits a heads-up when testing a very large contrast with no
  `max_cells` cap (the one runtime step whose memory scales with the
  matrix), pointing at `max_cells` to bound time and peak RAM.

### Fix: `scroll_add_subset()` preserves scoped-column types

- A subset’s scoped metadata columns kept their **original type**
  instead of being coerced to character. Numeric scoped columns
  (e.g. `AddModuleScore` signature scores, `percent.*`) were being
  stringified, so the manifest classified them as high-cardinality
  **categorical** (one “level” per cell) rather than **numeric** —
  breaking numeric colour-by / violins / sliders on those columns.
  Rebuild any project with numeric scoped columns to pick up the correct
  types.
- The **Violin** (Metadata mode) and **Biaxial** numeric-column pickers
  are now view-aware, so scoped numeric columns appear when their subset
  view is active — matching the categorical Group-by, which was already
  view-aware. (Added `.scroll_bind_view_nums`, a numeric sibling of
  `.scroll_bind_view_cats`.)

### DE: group by multiple columns

- The live-DE panel’s **Group by** is now multi-select: naming several
  categorical columns compares their **interaction levels**
  (e.g. `genotype` + `timepoint` gives groups like `"KO | d7"`), and the
  `ident1`/`ident2` selectors pick which combined levels form each side
  — matching how the Pseudobulk panel already groups.
  `scroll_de(group_col = c("genotype", "timepoint"))` accepts a column
  vector too. A single column behaves exactly as before.

### Live contrast preview (DE + Pseudobulk)

- Both differential panels now show a small **live mini-embedding** of
  the chosen contrast, updating as you pick controls — *before* you
  Compute. In the DE panel, `ident1` cells are drawn **red** and
  `ident2` (or “rest”) **blue** over a faint grey outline of the whole
  UMAP; in Pseudobulk, one coloured point per pseudobulk **sample** is
  placed at that sample’s medoid (showing how cells compact into their
  groups). It’s pure metadata + coordinates (no store query),
  stratified-subsampled for speed, and driven by the same
  contrast-membership logic the compute uses, so the picture always
  matches what will be tested. New exported view cores
  `view_contrast_preview()` and `view_contrast_medoids()`.
- Fix: both differential panels’ live previews now pick the embedding
  the way DimPlot does — the config `default_embedding` for the whole
  dataset, else the active subset view’s primary embedding. Previously
  Pseudobulk ignored the active view (subset medoids landed on the
  whole-dataset UMAP), and both panels fell back to the *first* global
  reduction for the whole dataset (e.g. `pca` instead of the default
  `umap`).

### Faster, lighter live DE

- Live DE
  ([`scroll_de()`](https://shihanli92.github.io/scroll/reference/scroll_de.md))
  reconstructs the contrast’s sparse matrix with a reverse-index gather
  instead of [`match()`](https://rdrr.io/r/base/match.html) over tens of
  millions of cell ids, and frees the long (feature, cell, value) table
  **before** allocating the matrix so the two never coexist. On a
  44k-cell one-vs-rest contrast this cut **peak memory ~2.3×** (≈5.9 GB
  → ≈2.6 GB) with byte-identical results — enough to stop the swapping
  that made DE crawl on a 16 GB server.
- New `scroll_de(max_cells = )` + a **“Max cells / group”** control in
  the DE panel (default off) down-sample each side of a contrast to
  bound DE time and memory on large or many-cell datasets — a standard
  marker-detection shortcut (cf. Seurat’s `max.cells.per.ident`);
  e.g. capping at 5,000/group ran ~4× faster with an unchanged
  top-marker list. The subsample is deterministic and leaves the session
  RNG untouched.
- The DE and Pseudobulk **Compute** buttons are now
  [`bslib::input_task_button`](https://rstudio.github.io/bslib/reference/input_task_button.html)s:
  clicking one immediately turns it into a disabled spinner labelled
  “Computing…” until results return. Previously the only progress cue
  was an output-area spinner that needed the optional `shinycssloaders`
  package — absent on many deploy servers, users saw nothing happen
  during the multi-second compute and assumed the app had frozen. The
  task button needs only `bslib` (already required) and gives clear
  feedback right where the user clicked.

### Expression store: first-letter bucketing

- The expression store now partitions by the feature’s **case-folded
  first character** (`expr/<assay>/bucket=<char>/`) instead of one
  directory per feature. On a 32k-gene store that is **~36 directories
  instead of ~25,600**, which removes the arrow Dataset-open crawl that
  dominated cold queries and could stall a server (fewer inodes / file
  descriptors), and — because rows are sorted by `feature` within each
  bucket — compresses **~1.8x smaller** with **no loss** (a single-gene
  lookup still reads little via Parquet row-group pruning; benchmarked
  *faster* than per-feature). Case is folded (so `Cd8a` and `ccdc198`
  share bucket `C`) because case-insensitive filesystems collide `C`/`c`
  directories. The runtime query layer is unchanged and reads old
  per-feature stores and new bucketed stores alike, so existing projects
  keep working un-rebuilt. Differential expression is unaffected in
  results and slightly faster (its whole-store scan opens far fewer
  files); pseudobulk uses the separate single-file counts store and is
  untouched.

### Manifest: cap cached levels for high-cardinality columns

- `manifest.yaml` no longer inlines **every** distinct value of a
  categorical column as a `levels:` list. Above
  `scroll_build(max_levels = 200)` (new argument, also on
  [`scroll_update()`](https://shihanli92.github.io/scroll/reference/scroll_update.md))
  a column records only its `n_levels` count; the app recomputes the
  level set from `cells.parquet` on demand. This keeps the manifest
  small and hand-editable when a build carries a clone id / barcode /
  sample id with thousands of values — those columns still export and
  stay fully usable everywhere (Color-by, Violin, DE gating). Manual
  per-level colour pickers now fall back to a note above 30 levels
  instead of rendering a wall of widgets.

### Fixes

- **Non-ASCII feature names.** A feature whose name contains a non-ASCII
  character (e.g. an antibody like `FcεRIa`) broke DE and pseudobulk
  with `'i' and 'j' must not contain NA`, and silently returned no data
  in FeaturePlot/DotPlot: `yaml` escaped the name in the manifest
  (`Fc<U+03B5>RIa`) so it no longer matched the store’s UTF-8 feature
  column.
  [`scroll_de()`](https://shihanli92.github.io/scroll/reference/scroll_de.md)
  /
  [`scroll_pseudobulk_de()`](https://shihanli92.github.io/scroll/reference/scroll_pseudobulk_de.md)
  now index on the store’s own feature names (encoding-robust; no
  rebuild needed), and the manifest is written with `unicode = TRUE` so
  names round-trip as UTF-8 (fixes the query/label path on rebuild).

### Incremental updates

- **`scroll_update(dir, object, embeddings =, meta_cols =, subsets =)`**
  adds or replaces reductions / metadata / subset views on an
  **already-built** project by rewriting only `cells.parquet` +
  `manifest.yaml`, aligned to the existing cells by barcode. The
  feature-partitioned `expr/` and `counts/` stores are left
  byte-identical, so appending a reduction is seconds rather than a full
  [`scroll_build()`](https://shihanli92.github.io/scroll/reference/scroll_build.md)
  re-export (measured: ~2s vs ~8min to add a column to a 47k-cell, 794MB
  project). Arguments mirror
  [`scroll_build()`](https://shihanli92.github.io/scroll/reference/scroll_build.md),
  so the same object flows through.

### Panels

- **Violin** and the ridge view can plot a **numeric metadata column**
  as the value axis (not just a queried feature), via a Gene/Metadata
  source toggle.
- **Biaxial** can put **genes on the axes** (query 2+ genes from any
  assay) in addition to numeric metadata columns; facets now default to
  **square** and gain **facet-columns / facet-rows** layout controls.
- Fixed a `conditionalPanel` namespacing bug that hid the
  source-dependent pickers (gene vs metadata) in the Violin, Biaxial,
  and Spatial panels.

### Storage format v2 + streaming builds

- **Leaner store (v2).** The expression store now keys `cell` on an
  **int32 global row-index** into `cells.parquet` (not the repeated
  barcode string), stores unquantized values as **float32** (not
  float64; ~lossless), compresses with **zstd**, and **compacts** to one
  part-file per feature. On a 4.4M-cell atlas this cut the store ~3×
  with no lossy precision change and faster cold queries. The runtime
  join became simpler and faster (positional integer indexing). Older v1
  (string-cell) stores keep working — a manifest flag (`cell_index`)
  selects the read path, so existing projects need no rebuild.
- **[`scroll_build_stream()`](https://shihanli92.github.io/scroll/reference/scroll_build_stream.md)**
  — build one project from many sources *one at a time*, so peak memory
  stays ~per-source instead of loading the whole dataset. Idempotent and
  append-safe (re-run as data arrives). Powers multi-million-cell
  atlases on a laptop.

### Architecture

- **Two phases.** A heavy offline
  [`scroll_build()`](https://shihanli92.github.io/scroll/reference/scroll_build.md)
  extracts lightweight on-disk artifacts — cell metadata + embeddings as
  one small Parquet table (`cells.parquet`) and expression as a
  feature-partitioned Parquet store (`expr/<assay>/`). A light runtime
  ([`scroll_app()`](https://shihanli92.github.io/scroll/reference/scroll_app.md)
  /
  [`scroll_serve()`](https://shihanli92.github.io/scroll/reference/scroll_serve.md))
  reads only those, querying expression **one gene at a time** through
  arrow, so runtime memory stays flat regardless of matrix size.
- **Deploys as a plain `app.R`** on an open-source Shiny Server — no
  render step. The running app never loads the Seurat object.
- **Build progress.** `scroll_build(verbose = interactive())` reports
  each phase and shows a progress bar over the feature-partition export
  (the slow step); the export is written in feature batches, which also
  bounds peak memory.
- Expression is `uint8`-quantized by default (`quantize = FALSE` for
  exact precision); an opt-in raw `counts/<assay>.parquet` store
  (`counts = TRUE`) backs the pseudobulk panel.
- Multimodal: pass `assays =` to export more than the default assay;
  each panel gains an **Assay** selector.

### Panels

Built-in panels **surface only when the project’s data supports them** —
a panel that can’t work on a given project simply doesn’t appear (rather
than showing an empty-state message). DimPlot and FeaturePlot are always
present; Biaxial needs ≥2 numeric columns; DotPlot/Violin need a
categorical column; Proportions needs ≥2; **DE / Pseudobulk require a
categorical column with ≥2 levels** (a real contrast), and Pseudobulk
additionally needs a counts store. So, e.g., a single-region Visium
slide shows the Spatial, DimPlot, FeaturePlot, DotPlot, Violin and
Biaxial panels but not DE/Pseudobulk/Proportions, while a multi-region
study gets them back.

For explicit control, `scroll_build(panels = c(...))` (or a `panels:`
list in `config.yaml`) shows **exactly** those panels in that order,
overriding the automatic gating; `exclude_panels` / `exclude_panels:`
hides specific panels. `config.yaml` is hand-editable, so the panel set
can be changed without rebuilding.

Each built-in analysis panel carries its own fine-grained controls:

- **DimPlot** — embedding coloured by any metadata column.
- **FeaturePlot** — coloured by a gene’s expression **or** a numeric
  metadata column, with colour-quantile clipping.
- **Biaxial** — pairwise scatters of numeric metadata columns (hashtag /
  ADT / QC).
- **DotPlot** — marker panel with z-score scaling and hclust
  dendrograms.
- **Violin** and **Proportions** — distributions and composition.
- **DE** — live Wilcoxon markers via `presto`, one compute shown as a
  **Table** and a **Volcano**.
- **Pseudobulk DE** — replicate-aware `edgeR` / `limma-voom` on
  aggregated sample-level counts, with pseudo-replicate modes and
  optional stability re-runs.

### Multimodal: scATAC

- **[`atac_spec()`](https://shihanli92.github.io/scroll/reference/atac_spec.md) +
  `scroll_build(..., atac =)`** ingest a chromatin-accessibility
  **peaks** assay. Because peaks are just an assay, **FeaturePlot /
  DotPlot / Violin / DE work on accessibility unchanged**;
  [`atac_spec()`](https://shihanli92.github.io/scroll/reference/atac_spec.md)
  additionally marks the assay `kind: peaks`, parses each peak’s
  chr/start/end from its `chr-start-end` name, and bakes a
  `peaks.parquet` annotation table (with the nearest gene when a Signac
  `ChromatinAssay` annotation is present).
- **Peaks panel** — find a peak **by nearby gene** (when annotation was
  baked) or **by region** (a text filter over peak names), then colour
  the embedding by its accessibility. Auto-surfaces only for projects
  built with an `atac` spec. Validated on the 10x PBMC multiome scATAC
  dataset.

### Multimodal: spatial

- **[`spatial_spec()`](https://shihanli92.github.io/scroll/reference/spatial_spec.md) +
  `scroll_build(..., spatial =)`** ingest 10x Visium / imaging-based
  data. The build extracts
  [`GetTissueCoordinates()`](https://satijalab.github.io/seurat-object/reference/GetTissueCoordinates.html)
  into a `spatial` embedding (image-pixel space, y-oriented for ggplot)
  and bakes the tissue image to a small raster asset; the manifest
  records an `images` block and marks the embedding `kind: spatial`.
  Because coordinates are just an embedding, **DimPlot and FeaturePlot
  work on spatial data unchanged**.
- **Spatial panel** — cells/spots in tissue space with a fixed aspect
  ratio, coloured by a gene or any metadata column, over the tissue
  image when one was baked. **Drag to zoom into a region, double-click
  to reset.** Auto-surfaces whenever the project has a spatial embedding
  — so it also works for **imaging platforms (Xenium / CosMx)** that
  carry cell centroids but no H&E image (the image toggle is simply
  hidden). Validated on a 10x Visium mouse-brain section (spots overlay
  the H&E exactly) and a Xenium-style FOV object.

### Multimodal: VDJ / immune repertoire

- **[`vdj_spec()`](https://shihanli92.github.io/scroll/reference/vdj_spec.md) +
  `scroll_build(..., vdj =)`** ingest per-cell TCR/BCR metadata. The
  build bakes a compact `repertoire/` Parquet store (clone table +
  pre-computed diversity, V/J gene-usage residuals, and tissue
  correlation) and records a `vdj` block in the manifest. `chain_type`
  parameterises TCR vs BCR segment names (TRBV/TRAV… vs IGHV/IGKV…);
  column defaults follow 10x naming.
- **Four repertoire panels** — *Clone overview* (rank-abundance +
  expansion composition), *V/J gene usage* (frequency + chi-square
  residual heatmap), *CDR3 length*, and *Diversity* (Shannon / Simpson /
  clonality / Gini + tissue correlation). They **auto-surface only for
  projects built with a `vdj` spec** and are absent otherwise, via a new
  `when(manifest)` gate on built-in panels — so RNA-only apps are
  unchanged.

### Multi-dataset

- **`scroll_multi_app(projects)`** mounts several built projects behind
  one app — a tab per dataset, each with its own app bar,
  manifest-driven panels, and query handle. The per-dataset page +
  server were refactored into a namespaced module
  ([`shiny::NS()`](https://rdrr.io/pkg/shiny/man/NS.html)), so
  independent datasets (even different organisms / feature spaces)
  coexist without id collisions;
  [`scroll_app()`](https://shihanli92.github.io/scroll/reference/scroll_app.md)
  is the single-dataset case and is unchanged. Custom panels should read
  baked assets from the per-dataset handle (`data$dir`) so each tab
  loads its own files.

### Explore & customize

- **Subset views** — expose a reprocessed slice (its own sub-embedding +
  subclusters) as a linked **View**;
  [`scroll_add_subset()`](https://shihanli92.github.io/scroll/reference/scroll_add_subset.md)
  assembles the child onto the parent by barcode, and selecting the view
  restricts every panel to it.
- **Global cell-subset filter** in the app bar threads through every
  panel.
- **Manual palettes** — per-level colour pickers on every categorical
  panel; colours are assigned deterministically by level name, so a cell
  type keeps its colour across panels.
- **Export** — every plot to PNG (raster) or PDF (vector); the DE table
  to CSV, with an “export all genes” option.
- **Performance** — large scatters rasterize on screen (`scattermore`)
  while exports stay vector; queries are memoized (LRU) and cosmetic
  controls debounced.
- **Stability** — killed scrollbar-driven resize loops on tall
  multi-panel pages. A scrollbar toggling as plots render (vertical bar
  → width change, horizontal bar → height change) fires window `resize`,
  which re-renders every fluid-width `plotOutput`, which nudges the size
  back — an endless loop where plots appeared to update on their own.
  The page now reserves the vertical gutter
  (`overflow-y:scroll; scrollbar-gutter:stable`) and hides page-level
  horizontal overflow (`overflow-x:hidden`; wide plots still scroll
  inside their card). The plot wrapper also pins `min-width:0` and snaps
  the plot output to a whole CSS pixel
  (`width: round(down, 100%, 1px)`). bslib’s grid columns are fractional
  (e.g. 518.25px); on a HiDPI display (`devicePixelRatio` 2) that .25px
  is half a device pixel, so Shiny rendered the image at a rounded
  device size that displayed back a hair different, tripping the
  per-output `ResizeObserver` into re-rendering every plot forever (only
  at HiDPI + narrow widths). Snapping the width to an integer makes the
  device pixels land exactly, so nothing oscillates.

### Extend

- `register_panel(id, ui, server, ...)` adds a panel to every app built
  afterwards, using the same module contract as the built-ins;
  [`scroll_reset_panels()`](https://shihanli92.github.io/scroll/reference/scroll_reset_panels.md)
  clears custom registrations.
- `register_plot_panel(id, plot, controls, ...)` is a declarative
  wrapper: describe the controls
  ([`scroll_input_column()`](https://shihanli92.github.io/scroll/reference/scroll_input.md),
  [`scroll_input_levels()`](https://shihanli92.github.io/scroll/reference/scroll_input.md),
  [`scroll_input_gene()`](https://shihanli92.github.io/scroll/reference/scroll_input.md),
  …) and a single `plot(cells, input, data)` function, and scroll
  generates the whole module — UI, control population, a Compute gate,
  inline error messages, PNG/PDF export, and subset/view awareness. The
  lower-level helpers
  ([`scroll_render_plot()`](https://shihanli92.github.io/scroll/reference/scroll_render_plot.md),
  [`scroll_bind_levels()`](https://shihanli92.github.io/scroll/reference/scroll_bind_levels.md),
  [`scroll_columns()`](https://shihanli92.github.io/scroll/reference/scroll_columns.md))
  are exported for hand-written panels too.
- **Richer declarative controls** for
  [`register_plot_panel()`](https://shihanli92.github.io/scroll/reference/register_plot_panel.md):
  [`scroll_input_assay()`](https://shihanli92.github.io/scroll/reference/scroll_input.md)
  (auto-hidden on single-assay data),
  [`scroll_input_embedding()`](https://shihanli92.github.io/scroll/reference/scroll_input.md)
  and view-aware `scroll_input_column(view_aware=)` (narrow to the
  active subset view),
  [`scroll_input_custom()`](https://shihanli92.github.io/scroll/reference/scroll_input.md)
  (wrap your own `ui`/`bind`), and
  [`scroll_show_when()`](https://shihanli92.github.io/scroll/reference/scroll_show_when.md)
  for conditional visibility.
- **Reusable render + data helpers exported** so custom views match the
  built-ins:
  [`scroll_point_layer()`](https://shihanli92.github.io/scroll/reference/scroll_helpers.md)
  (raster-aware scatter),
  [`scroll_discrete_colors()`](https://shihanli92.github.io/scroll/reference/scroll_helpers.md)
  /
  [`scroll_continuous_scale()`](https://shihanli92.github.io/scroll/reference/scroll_helpers.md)
  (shared palettes),
  [`scroll_group_labels()`](https://shihanli92.github.io/scroll/reference/scroll_helpers.md),
  and the accessors `scroll_columns(view=)`,
  [`scroll_embeddings()`](https://shihanli92.github.io/scroll/reference/scroll_embeddings.md),
  [`scroll_features()`](https://shihanli92.github.io/scroll/reference/scroll_features.md).
- **Fewer silent failures.** A non-ggplot plot return now shows a clear
  message; a `scroll_input_levels(from=)` that names no column/control
  warns; contradictory `required` + `none` warns at construction.
- **Data-derived choices.** Controls can now compute options from the
  data handle instead of only the manifest:
  `scroll_input_choice(choices = function(data))` and
  `scroll_input_levels(choices = function(input, data), watch = c(...))`
  populate from e.g. a baked asset under `data$dir`, recomputing when a
  watched control changes. `scroll_input_column(prefer=)` sets a
  preferred default column, and
  [`scroll_input_palette()`](https://shihanli92.github.io/scroll/reference/scroll_input.md)
  offers the built-in discrete/continuous palette names.

### Documentation

- Reference documentation for every export, published as a [pkgdown
  site](https://shihanli92.github.io/scroll/).

### Internal

- Split the monolithic `R/app.R` into one file per built-in view
  (`R/panel-<view>.R`), plus `R/panel-helpers.R` (shared panel helpers)
  and `R/data-handle.R` (the runtime data handle). Pure reorganization —
  no behaviour change.
