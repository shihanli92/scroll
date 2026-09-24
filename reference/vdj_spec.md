# Describe a dataset's TCR/BCR (VDJ) metadata for `scroll_build()`

Maps a Seurat object's per-cell repertoire columns onto the schema
scroll's repertoire panels use. Pass the result as
`scroll_build(..., vdj = vdj_spec(...))` to bake the `repertoire/` store
and enable the VDJ panels.

## Usage

``` r
vdj_spec(
  chain_type = c("TCR", "BCR"),
  group_col,
  clone_col = NULL,
  count_col = NULL,
  segments = NULL,
  cdr3 = NULL,
  source = c("auto", "scroll", "scRepertoire", "airr", "platypus"),
  antigen_col = NULL,
  tissue_col = NULL,
  cluster_col = NULL,
  loc_col = NULL,
  broad_col = NULL,
  carry = character(),
  exclude = NULL
)
```

## Arguments

- chain_type:

  `"TCR"` or `"BCR"` — sets default segment/CDR3 column names and panel
  labels (TRBV/TRAV… vs IGHV/IGKV…).

- group_col:

  A categorical per-cell column to group repertoire summaries by (e.g.
  cell type / cluster). Required.

- clone_col:

  Per-cell clonotype id. `NULL` (default) uses the preset's clone
  column, else (for a non-`"scroll"` source) the first present of
  `clone_id`, `strict_clone_id`, `clonotype_id`, `raw_clonotype_id`,
  `clonotype`, else derives a clonotype from the pasted CDR3 columns.

- count_col:

  Optional per-cell clone-size column. `NULL` uses `<clone_col>_count` /
  `<clone_col>_size` when present, else computes clone size from
  `clone_col` frequency.

- segments:

  Named character vector `label = column` of V/J gene columns (defaults
  by `chain_type`).

- cdr3:

  Named character vector `chain = column` of CDR3 amino-acid columns
  (defaults by `chain_type`); drives CDR3-length views.

- source:

  Upstream convention to read: `"auto"` (default; detect from the
  metadata columns), `"scroll"`, `"scRepertoire"`, `"airr"`, or
  `"platypus"`. Any column you pass explicitly (`segments`, `cdr3`,
  `clone_col`, …) overrides the preset. Column names are matched
  ignoring case when the exact name is absent (so `v_call_vdj` satisfies
  the AIRR preset's `v_call_VDJ`); each such match is reported with a
  message.

- antigen_col, tissue_col, cluster_col, loc_col, broad_col:

  Optional grouping / annotation columns (antigen split, tissue split,
  per-cluster diversity, and the location/broad pair used for the
  tissue-correlation view).

- carry:

  Extra categorical columns to keep in the clone table so panels can
  split / colour by them.

- exclude:

  Optional value of `antigen_col` to drop (e.g. a negative control).

## Value

A `scroll_vdj_spec` list.

## Details

`scroll` reads TCR/BCR data straight from the common upstream tools —
set `source` (or leave it `"auto"` to detect) and the segment/CDR3/clone
columns are mapped for you; override any of them explicitly when your
object differs:

- `"scRepertoire"` — Seurat objects from `combineTCR()`/`combineBCR()` +
  `combineExpression()`; the compound `CTgene`/`CTaa`/`CTstrict` columns
  are parsed into per-segment/per-chain columns.

- `"airr"` — dandelion / AIRR per-cell obs (`v_call_VDJ`/`v_call_VJ`,
  `junction_aa_VDJ`/`_VJ`, `clone_id`).

- `"platypus"` — Platypus VDJ.GEX per-cell columns
  (`VDJ_vgene`/`VJ_vgene`, `VDJ_cdr3s_aa`/`VJ_cdr3s_aa`).

- `"scroll"` — scroll's own `v_gene_TRB`/`cdr3_beta`… convention (the
  historical default). Long, one-row-per-contig tables (raw 10x
  `filtered_contig_annotations.csv` / long AIRR `.tsv`) must be
  collapsed to per-cell columns upstream (e.g. with scRepertoire or
  dandelion) before building.

## See also

[`scroll_build()`](https://shihanli92.github.io/scroll/reference/scroll_build.md)
