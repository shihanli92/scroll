# scroll — dev notes & deferred items

Running record of known issues, deferred cleanups, and enhancement
ideas. Started 2026-07-09 after the first full code review (commit
`50e0452` fixed the top-3 blockers + documented the quantization
caveat). Items below were **consciously deferred**, not forgotten.

## Deferred issues (from code review)

### 1. Live DE is unbounded in the gene dimension

- **Where:**
  [`scroll_de()`](https://shihanli92.github.io/scroll/reference/scroll_de.md)
  (`R/de.R`) →
  [`scroll_query_cells()`](https://shihanli92.github.io/scroll/reference/scroll_query_cells.md)
  (`R/query.R`).
- **What:** A DE run reconstructs an *all-genes × contrast-cells* sparse
  matrix in memory for presto. Bounded by the contrast’s cell count and
  gated by the Compute button + a ≥3-cell floor, but there is **no upper
  cap**. A one-vs-rest on a large cluster in a 500k-cell dataset builds
  a big in-RAM matrix — the one intentional exception to the flat-RAM
  invariant, currently unguarded on size.
- **Impact:** memory spike / slowness on large contrasts. Fine for
  pbmc3k.
- **Options:** warn (or block) above N cells; offer downsampling per
  group; document the expected footprint.
- **Priority:** medium (matters only at scale).

### 2. `de_table` view / dispatch is app-unreachable (kept as public API)

- **Where:** `view_de_table()` + `render_view()`’s `de_table` branch.
  The DE panel uses live `scroll_de` instead, so these are reached only
  by unit tests.
- **Status:** the empty `de/` build dir was removed (nothing wrote/read
  it). The remaining `view_de_table` / `de_table` dispatch is part of
  the **exported** `render_view()` grammar, so it stays as public API
  rather than dead code.
- **Priority:** low (a future precomputed-DE source could reuse it).

### 3. `logFC` labeling is a mean-difference, not literal log2FC

- **Where:**
  [`scroll_de()`](https://shihanli92.github.io/scroll/reference/scroll_de.md)
  output column `logFC`; `view_volcano()` axis label and default cutoff
  `1`.
- **What:** presto’s `logFC` = difference of group means of the
  (lognorm) input — the standard presto/Seurat convention, **not a
  bug**, but the label + the cutoff-of-1 are scale-dependent and
  slightly imprecise.
- **Options:** relabel (e.g. “avg log-expr diff”), or convert to log2
  and label `avg_log2FC`; document the cutoff units.
- **Priority:** low.

### 4. Assay name interpolated into the glob path — RESOLVED

- `scroll_query_feature/features/cells` now build the glob via
  `.scroll_assay_glob()`, which rejects any `assay` that is not a plain
  path segment (no separators, no `..`, non-empty). Covered by
  test-multiassay.R.

### 5. Feature names with path separators are warned, not rejected

- **Where:** `.scroll_warn_unsafe_features()` (`build.R`) warns on
  `/`/`\` in feature names but does not reject/sanitize; such a name
  would break arrow’s hive partitioning at build.
- **Fix:** sanitize or stop on unsafe feature names.
- **Priority:** low (gene symbols are normally safe).

### 6. Unused Suggests — RESOLVED

- `shinycssloaders` is now wired via `.scroll_spin()` (guarded
  `requireNamespace`) around every panel’s plot/table output, so slow
  renders (DE presto, feature queries) show a spinner. No longer an
  unused Suggest.
- `SeuratData` is used only by tests/README (fine — keep).

## Test gaps to close later

- ~~**Multi-assay** build/query/UI~~ — DONE. The shared fixture is now
  RNA + a synthetic ADT; test-multiassay.R covers two `expr/` subtrees,
  per-assay `max` dequantize, the `length(assays) > 1` selector
  branches, an ADT FeaturePlot, and a two-group ADT DE.
- ~~**Two-group DE** end-to-end~~ — DONE (ADT `T` vs `B` in
  test-multiassay.R).
- **Quantization floor on `pct`**: a test pinning the near-zero drop
  behavior so the tradeoff is captured in tests.
- **`.scroll_warn_unsafe_features`** and the `overwrite = FALSE`
  non-empty-dir stop path.
- **App-level connection cleanup**: assert
  [`scroll_app()`](https://shihanli92.github.io/scroll/reference/scroll_app.md)
  registers the `onStop` (hard to introspect cleanly; currently only
  manually verified).
- **Browser/Shiny-UI layout**: `testServer` covers module logic well but
  not the rendered bslib layout / scroll-spy / selectize behavior —
  those are validated manually in a foreground browser, not in CI.

## Enhancement ideas (not bugs)

- ~~**Merge DE + Volcano**~~ — DONE. One DE section, `Table | Volcano`
  tabs from a single presto compute.
- ~~**Global cell-subset filter**~~ — DONE. App-bar pill threads
  `cells_r` through every panel; cells stat reads “N of total”.
- ~~**Per-plot export**~~ — DONE. PNG (all panels, via a png-device +
  [`print()`](https://rdrr.io/r/base/print.html) path that also handles
  the DotPlot aplot) + CSV (DE table).
- ~~**[`register_panel()`](https://shihanli92.github.io/scroll/reference/register_panel.md)**~~
  — DONE. Public `register_panel(id, ui, server, ...)`
  - [`scroll_reset_panels()`](https://shihanli92.github.io/scroll/reference/scroll_reset_panels.md);
    a session-global registry assembled with the built-ins
    (replace-by-id, `after` positioning, numbered by position).
- **WebGL / rasterized scatter** (`scattermore`/`ggrastr`) for DimPlot/
  FeaturePlot when cell counts climb past what canvas renders smoothly
  (~100k+).
- **Block-wise build** for very large objects (the export currently
  builds the triplet in one pass).
- **Pseudobulk `duplicateCorrelation`** (deferred, 0.2.9). The paired
  real-replicate path fits an explicit `~ replicate + group` block
  (deterministic, simple). For many unbalanced donors
  [`limma::duplicateCorrelation`](https://rdrr.io/pkg/limma/man/dupcor.html)
  (block = replicate) can be the better estimator; slot it as a
  `paired`/`method` option later. See
  [`scroll_pseudobulk_de()`](https://shihanli92.github.io/scroll/reference/scroll_pseudobulk_de.md).
- **Pseudobulk covariates.** Only `group` (+ optional replicate block)
  enter the design; a general `covariates =` (extra metadata terms) is a
  natural follow-up.
- ~~**Counts in streaming builds.**~~ — DONE (0.2.22).
  `scroll_build_stream(counts = TRUE)` writes one counts part per source
  under `counts/<assay>/` (no compaction needed: pseudobulk full-scans
  anyway, and `open_dataset()` reads the directory).

## Known behavior (documented, working as intended)

- **Quantization** (`quantize = TRUE` default) floors values below
  ~`max/510` to zero → fraction-expressing stats (dotplot dot size, DE
  `pct.1`/`pct.2`) slightly under-count low expression. Documented on
  `scroll_build`, `scroll_de`, README. Use `quantize = FALSE` for exact
  fractions.
- **Aspect ratio** is a `theme(aspect.ratio)` reshape within a fixed
  canvas (bounded, no page growth); on the DotPlot it applies only when
  clustering is off (with dendrograms the aplot composite fills the
  fixed canvas — aspect would detach the trees).
- **Flat-RAM invariant** holds on the runtime hot path (verified in
  review); live DE is the sole, intentional exception (see item 1).
