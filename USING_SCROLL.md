# Using `scroll` — a guide for an AI agent

`scroll` is an R package that turns a **processed Seurat object** into an interactive,
memory-light single-cell explorer (Shiny app). This document is written for an agent that has
been handed `scroll_0.2.5.tar.gz` and needs to install it, build a project from a Seurat object,
run the app, and extend it — without re-deriving the architecture.

---

## 1. The one invariant — two phases, runtime never loads Seurat

1. **Build (offline, heavy)** — `scroll_build()` reads a Seurat object and writes flat on-disk
   artifacts (Parquet + YAML).
2. **Runtime (light)** — `scroll_app()` reads only those artifacts. Runtime RAM scales with the
   cell-metadata table only; **expression stays on disk** and is queried one feature at a time.

`SeuratObject` is a **Suggests**, needed only at build time. If you load a Seurat object in app
code, you have broken the design.

---

## 2. Install

```r
install.packages("scroll_0.2.5.tar.gz", repos = NULL, type = "source")
```

Hard runtime deps (all CRAN, binary): `shiny`, `bslib`, `arrow`, `dplyr`, `tidyr`, `stringr`,
`ggplot2`, `yaml`, `scales`, `tibble`, `Matrix`, `colourpicker`. Build-time extra: `SeuratObject`.

Optional features degrade gracefully via Suggests — install only what you need:

| Feature | Package(s) |
|---|---|
| Live Wilcoxon DE | `presto` |
| Pseudobulk DE | `edgeR`, `limma` |
| Rasterized large scatters | `scattermore` |
| DotPlot dendrograms | `ggtree`, `aplot`, `ape` |
| ATAC | `Signac` |

---

## 3. Build → run in 30 seconds

```r
library(scroll)
# obj is a processed Seurat object: >=1 reduction + a non-empty `data` layer
scroll_build(obj, "myproject")   # writes ./myproject/
scroll_serve("myproject")        # launches the Shiny app
```

`scroll_build()` auto-infers assays, embeddings, and metadata columns when you don't name them.
Panels appear automatically based on what the data supports (see §6).

---

## 4. What a built project looks like

```
myproject/
├── cells.parquet          # one row per cell: cell, <meta…>, <reduction>_<dim>… (zstd)
├── expr/<assay>/part-0.parquet   # long (feature, cell, value); one file per assay
├── counts/<assay>.parquet  # only if counts = TRUE (integer, unquantized)
├── repertoire/             # only if vdj =
├── spatial/<name>.rds      # only if spatial =
├── atac/peaks.parquet      # only if atac =
├── manifest.yaml           # the contract the runtime reads
├── config.yaml             # title, defaults, markers, panel order (hand-editable)
└── app.R                   # library(scroll); scroll_app(".")
```

**Store format v2 facts that matter:**
- The `cell` column in `expr/` is an **int32, 1-based global row index into `cells.parquet`**, not
  a barcode.
- `value` is `uint8` when `quantize = TRUE` (default, ~4× smaller, floors values below ~max/510),
  else `float32` (exact). `quantize = FALSE` is required for `scroll_build_stream()`.
- Everything is **zstd**. One feature-sorted Parquet file per assay; row-group stats prune a
  single-gene query.
- Only **nonzero** entries are stored — an absent cell means 0.

---

## 5. `scroll_build()` — the arguments you'll actually use

```r
scroll_build(object, outdir,
             assays = NULL, embeddings = NULL, meta_cols = NULL,
             quantize = TRUE,      # FALSE = exact float32
             counts = FALSE,       # TRUE unlocks Pseudobulk DE
             subsets = NULL,       # named list of subset-view specs (see §7)
             vdj = NULL, spatial = NULL, atac = NULL,   # modality specs (§6)
             panels = NULL, exclude_panels = NULL,       # override panel gating
             overwrite = FALSE)    # ⚠ TRUE wipes the whole outdir
```

Related build/assemble helpers:
- `scroll_build_stream(outdir, sources, reader, ...)` — multi-source streaming build for datasets
  too big for RAM (one source at a time; `quantize = FALSE` only).
- `scroll_update(dir, object, ...)` — rewrite `cells.parquet` + manifest only (metadata /
  embeddings), **without** touching the expression store. Fast; doesn't wipe side-baked dirs.
  Needs the Seurat object.
- `scroll_add_meta(dir, name, values, scope=, overwrite=)` — add **one** metadata column to a
  **built** project (writes `cells.parquet` + `manifest.yaml`), **no Seurat object needed** — ideal
  for a deployed app. `values` can be a length-`n_cells` vector, a barcode-named vector, a
  `data.frame(cell, value)`, or a `function(cells)` deriving it from existing columns; type
  (categorical/numeric) is inferred. Restart the app to see it.
- `scroll_add_subset(object, sub_object, name, embeddings, label=, meta=)` — attach a reprocessed
  sub-embedding + scoped metadata to the parent by barcode.

⚠️ **`overwrite = TRUE` wipes the entire project dir**, including subdirectories a separate script
baked. Re-run side bakes after every rebuild, or use `scroll_update()` for metadata-only changes.

---

## 6. Panels are data-driven — you don't write panel code to get them

Built-in panels self-gate on the manifest:

| id | appears when |
|---|---|
| `dimplot`, `featureplot` | always |
| `biaxial` | ≥2 numeric columns |
| `dotplot`, `violin` | ≥1 categorical column |
| `proportions` | ≥2 categorical columns |
| `de` | a categorical with ≥2 levels (needs `presto`) |
| `pseudobulk` | `counts = TRUE` at build + a contrast (needs `edgeR`/`limma`) |
| `clone_overview`, `gene_usage`, `cdr3_length`, `diversity` | `vdj =` was supplied |
| `spatial` | a spatial embedding exists |
| `peaks` | `atac =` was supplied |

**To surface the four repertoire panels you pass a spec, not panel code:**

```r
scroll_build(obj, "proj", counts = TRUE,
  vdj = vdj_spec(chain_type = "TCR", group_col = "sample",
                 clone_col = "clone_id", count_col = "clone_count"))
```

`vdj_spec()` supports presets + auto-detection for scRepertoire / dandelion(AIRR) / Platypus
outputs (`source = "auto"` by default). `spatial_spec()` and `atac_spec()` are the other two.

To control which panels appear, edit `config.yaml` (`panels:` allowlist bypasses the gate;
`exclude_panels:` drops ids) — no rebuild needed — or pass `panels=`/`exclude_panels=` at build.

---

## 7. Subset views (reprocessed sub-embeddings)

A subset view restricts **every** panel to a slice of cells and defaults to that slice's
embedding. Membership is implicit: **cells with non-NA coords in the subset's primary embedding.**

```r
# obj already carries a re-embedded child's UMAP + subcluster column (NA elsewhere)
obj <- scroll_add_subset(obj, tcell_obj, name = "tcell", label = "T cells",
                         embeddings = c(umap_t = "umap"), meta = "tcell_subcluster")
scroll_build(obj, "proj")   # the subset is picked up automatically
```

Scoped meta columns get their levels computed over members only, so subcluster labels never leak
into whole-dataset selectors. The app bar composes subset View → ad-hoc filter → the cells every
panel sees.

---

## 8. Writing a custom panel — the 90% path

```r
register_plot_panel("my_panel",
  label = "My panel", title = "…", after = "dimplot", compute = TRUE,
  controls = list(
    scroll_input_column("group", "Group by", type = "categorical"),
    scroll_input_levels("lv", "Levels", from = "group", none = TRUE),
    scroll_input_gene("genes", "Genes", multiple = TRUE)),
  plot = function(cells, input, data) {
    # `cells` is ALREADY the active (view + filter narrowed) table — never use data$cells
    if (!length(input$lv)) stop("Pick at least one level.")   # errors render inline
    long <- data$queryN("RNA", input$genes)                   # sparse, already dequantized
    ggplot2::ggplot(...)                                       # must return a ggplot
  })
```

- **Throwing is the error channel** — the message renders in the panel. Never silently `req()`.
- `compute = TRUE` gates behind a Compute button (use for anything expensive); `FALSE` is live.
- Register **before** calling `scroll_app()`. Call `scroll_reset_panels()` first if re-sourcing.
- Controls: `scroll_input_column / _levels / _gene / _numeric / _slider / _choice / _text /
  _palette / _assay / _embedding / _custom`; `scroll_show_when(spec, control, equals)` for
  cascades. Dynamic choices: `scroll_input_choice(choices = function(data))` and
  `scroll_input_levels(choices = function(input, data), watch = c("ctl"))`.
- Need multiple outputs or cross-panel state? Use `register_panel(id, ui, server)` and add a 4th
  `view_r` formal to receive the active view.

---

## 9. The `data` handle + the `.gidx` join (the #1 gotcha)

Every panel receives a `data` handle:

| field | what |
|---|---|
| `data$dir` | project root — resolve baked assets from here |
| `data$cells` | `cells.parquet` + a `.gidx` column = `seq_len(nrow())` stamped on at load |
| `data$manifest` / `data$config` | parsed YAML |
| `data$query1(assay, feature)` / `data$queryN(assay, features)` | LRU-cached, **already dequantized**, sparse long (absent cell ⇒ 0) |

Do **not** call `scroll_dequantize()` on `query1`/`queryN` results — they're already dequantized.

**Joining expression to cells** (store `cell` is a global index; `cells` may be filtered — `.gidx`
rides along on every filtered copy):

```r
values <- data$queryN("RNA", genes)                 # df(feature, cell, value)
key  <- if (!is.null(cells$.gidx)) cells$.gidx else seq_len(nrow(cells))
expr <- values$value[match(key, values$cell)]
expr[is.na(expr)] <- 0                               # sparse ⇒ absent means zero
```

Reference implementation: `.scroll_expr_vector()` in the package source.

---

## 10. Headless querying (no app)

```r
con <- scroll_connect("myproject")
df  <- scroll_query_feature(con, "RNA", "CD8A")   # data.frame(cell, value) — raw, NOT dequantized
scroll_disconnect(con)
```

Raw `scroll_query_*` values are quantized; call `scroll_dequantize(values, manifest, assay)` if
you need real units. (The `data$query1/queryN` handle in panels dequantizes for you.)

---

## 11. Differential expression

```r
scroll_de(data, assay, group_col, ident1, ident2 = NULL, min_pct = 0.1, cells = NULL)
scroll_pseudobulk_de(data, assay, aggregate_cols, ident1, ident2 =,
                     replicate_col =, ...)   # needs counts = TRUE at build
```

---

## 12. Downstream app convention

Apps built on scroll live **outside** the package, one folder per study:

```
<study>/scroll_app/
├── build_<name>.R   # Seurat → scroll_build() → projects/<name>/
├── panels/*.R       # custom panels, each self-registering
├── assets/          # small baked tables the panels read
├── projects/        # BUILT output — git-ignore this, regenerate from build_*.R
└── app.R            # scroll_reset_panels() → source(panels) → scroll_app("projects/<name>")
```

Deploy = copy the whole folder (with a built `projects/`) to a Shiny Server, or zip it. Keep
dataset-specific code and built artifacts out of the package repo.

---

## 13. Verify your build worked

```r
m <- scroll_manifest("myproject")
m$n_cells                          # cell count
names(m$assays)                    # assays present
names(m$embeddings)                # reductions
m$has_counts                       # pseudobulk-ready?
scroll_app("myproject")            # returns a shiny.appobj without error
```

---

## 14. Quick gotcha list

- `overwrite = TRUE` wipes the whole dir including side bakes — re-bake after, or use
  `scroll_update()`.
- In a custom panel `plot()`, use the passed `cells`, never `data$cells` (that's unfiltered).
- `query1`/`queryN` are already dequantized; raw `scroll_query_*` are not.
- Throwing an error in a panel is the intended way to show a message; don't `req()` silently.
- Non-ASCII feature names round-trip (manifests written with `unicode = TRUE`).
