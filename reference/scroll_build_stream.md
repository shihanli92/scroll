# Build a scroll project by streaming many sources (large-dataset build)

Produces the same flat-RAM store as
[`scroll_build()`](https://shihanli92.github.io/scroll/reference/scroll_build.md),
but assembled from many sources one at a time, so peak memory stays
~per-source instead of loading the whole (possibly multi-million-cell)
dataset at once. Each source is read by `reader`, its expression
appended to a shared v2 store under a running global cell index, and the
manifest written once at the end. The store is **compacted** (one
part-file per feature) on completion.

## Usage

``` r
scroll_build_stream(
  outdir,
  sources,
  reader,
  assays = NULL,
  embeddings = NULL,
  meta_cols = NULL,
  quantize = FALSE,
  id_of = as.character,
  overwrite = FALSE,
  verbose = interactive()
)
```

## Arguments

- outdir:

  Project directory to create or append to.

- sources:

  A list/vector of sources (file paths, ids, …) passed one at a time to
  `reader`.

- reader:

  `function(source) -> Seurat object`, already processed (a
  log-normalized `data` layer, the reductions, and the metadata
  columns).

- assays, embeddings:

  Names to export; `NULL` takes the default assay / all reductions of
  the **first** source (then held constant across sources).

- meta_cols:

  Metadata columns to expose (required; declared explicitly rather than
  inferred). The declared assays, embeddings, and meta columns must be
  present in **every** source — a source missing any is an error naming
  it (extra columns a source happens to carry are simply ignored).

- quantize:

  Must be `FALSE` (streaming stores float32; see Details).

- id_of:

  `function(source) -> character` id used for append-tracking and
  part-file tags (default `as.character`).

- overwrite:

  If `TRUE`, wipe `outdir` first (a fresh build).

- verbose:

  Print per-source progress.

## Value

Invisibly, `outdir`.

## Details

Idempotent and **append-safe**: re-running with more `sources` adds only
the new ones (tracked in `outdir/.stream_sources.txt`), so a build can
grow as data arrives. Values are stored unquantized (`float32`); a
global quantization max is unknown while streaming, so `quantize = TRUE`
is not supported here.

## See also

[`scroll_build()`](https://shihanli92.github.io/scroll/reference/scroll_build.md)
