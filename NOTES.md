# scroll — dev notes & deferred items

Running record of known issues, deferred cleanups, and enhancement ideas.
Started 2026-07-09 after the first full code review (commit `50e0452` fixed the
top-3 blockers + documented the quantization caveat). Items below were
**consciously deferred**, not forgotten.

## Deferred issues (from code review)

### 1. Live DE is unbounded in the gene dimension
- **Where:** `scroll_de()` (`R/de.R`) → `scroll_query_cells()` (`R/query.R`).
- **What:** A DE run reconstructs an *all-genes × contrast-cells* sparse matrix
  in memory for presto. Bounded by the contrast's cell count and gated by the
  Compute button + a ≥3-cell floor, but there is **no upper cap**. A one-vs-rest
  on a large cluster in a 500k-cell dataset builds a big in-RAM matrix — the one
  intentional exception to the flat-RAM invariant, currently unguarded on size.
- **Impact:** memory spike / slowness on large contrasts. Fine for pbmc3k.
- **Options:** warn (or block) above N cells; offer downsampling per group;
  document the expected footprint.
- **Priority:** medium (matters only at scale).

### 2. Vestigial precomputed-DE path
- **Where:** `build.R` creates an empty `de/` dir; `view_de_table()` +
  `render_view()`'s `de_table` branch + `.scroll_read_de()` exist but are
  **unreachable from the app** (the DE panel uses live `scroll_de`). Exercised
  only by direct unit tests.
- **Options:** (a) wire it — let the DE panel offer a "Precomputed contrast"
  source that reads `de/<contrast>.parquet`; or (b) remove the dead path to
  reduce confusion.
- **Priority:** low (cleanliness).

### 3. `logFC` labeling is a mean-difference, not literal log2FC
- **Where:** `scroll_de()` output column `logFC`; `view_volcano()` axis label and
  default cutoff `1`.
- **What:** presto's `logFC` = difference of group means of the (lognorm) input
  — the standard presto/Seurat convention, **not a bug**, but the label + the
  cutoff-of-1 are scale-dependent and slightly imprecise.
- **Options:** relabel (e.g. "avg log-expr diff"), or convert to log2 and label
  `avg_log2FC`; document the cutoff units.
- **Priority:** low.

### 4. Assay name interpolated into the glob path (defense-in-depth)
- **Where:** `scroll_query_feature/features/cells` build
  `file.path(dir, "expr", assay, "**", "*.parquet")`. `feature`/`cell` are bound
  params / registered temp tables (injection-safe); `assay` comes from the
  manifest (not runtime user input), so practical risk is low.
- **Fix:** allowlist `assay` against `names(manifest$assays)` (or `basename()`)
  before building the path.
- **Priority:** low.

### 5. Feature names with path separators are warned, not rejected
- **Where:** `.scroll_warn_unsafe_features()` (`build.R`) warns on `/`/`\` in
  feature names but does not reject/sanitize; such a name would break arrow's
  hive partitioning at build.
- **Fix:** sanitize or stop on unsafe feature names.
- **Priority:** low (gene symbols are normally safe).

### 6. Unused Suggests
- `shinycssloaders` is in Suggests but never used in `R/`. Either wire a spinner
  on the DE/Volcano/feature renders, or drop it from DESCRIPTION.
- `SeuratData` is used only by tests/README (fine — keep).
- **Priority:** trivial.

## Test gaps to close later

- **Multi-assay** build/query/UI: every test uses a single `RNA` assay; the
  `if (length(assays) > 1)` UI branches and per-assay `max` dequantize are
  untested with ≥2 assays.
- **Two-group DE** end-to-end (only one-vs-rest is covered).
- **Quantization floor on `pct`**: a test pinning the near-zero drop behavior so
  the tradeoff is captured in tests.
- **`.scroll_warn_unsafe_features`** and the `overwrite = FALSE` non-empty-dir
  stop path.
- **App-level connection cleanup**: assert `scroll_app()` registers the `onStop`
  (hard to introspect cleanly; currently only manually verified).
- **Browser/Shiny-UI layout**: `testServer` covers module logic well but not the
  rendered bslib layout / scroll-spy / selectize behavior — those are validated
  manually in a foreground browser, not in CI.

## Enhancement ideas (not bugs)

- **Merge DE + Volcano** into one section with `Table | Volcano` tabs sharing a
  single presto compute (today they're separate panels, so a user computes DE
  twice for the same contrast).
- **Global cell-subset filter** in the app bar (restrict every panel to a cell
  type / subset at once). The design mockup called for this; not built.
- **Per-plot export**: PNG (all panels) + CSV (DE table).
- **`register_panel()`**: promote the internal panel registry to a public
  extension point so users can add custom panels.
- **WebGL / rasterized scatter** (`scattermore`/`ggrastr`) for DimPlot/
  FeaturePlot when cell counts climb past what canvas renders smoothly (~100k+).
- **Block-wise build** for very large objects (the export currently builds the
  triplet in one pass).

## Known behavior (documented, working as intended)

- **Quantization** (`quantize = TRUE` default) floors values below ~`max/510`
  to zero → fraction-expressing stats (dotplot dot size, DE `pct.1`/`pct.2`)
  slightly under-count low expression. Documented on `scroll_build`, `scroll_de`,
  README. Use `quantize = FALSE` for exact fractions.
- **Aspect ratio** is a `theme(aspect.ratio)` reshape within a fixed canvas
  (bounded, no page growth); on the DotPlot it applies only when clustering is
  off (with dendrograms the aplot composite fills the fixed canvas — aspect would
  detach the trees).
- **Flat-RAM invariant** holds on the runtime hot path (verified in review);
  live DE is the sole, intentional exception (see item 1).
