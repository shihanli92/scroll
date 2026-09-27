# Add gene sets to a built scroll project

Saves gene-set collections into a project (`genesets/sets.parquet`) so
the GSEA panel can use them without any extra package at runtime. Needs
no Seurat object, so it also works on a deployed project. Collections
already saved under the same name are replaced; others are kept.

## Usage

``` r
scroll_add_genesets(
  dir,
  msigdb = NULL,
  gmt = NULL,
  species = NULL,
  assay = NULL
)
```

## Arguments

- dir:

  A built scroll project directory.

- msigdb:

  Names of curated MSigDB collections to fetch with the msigdbr package:
  any of `"H"` (Hallmark), `"REACTOME"`, `"KEGG"`, `"WIKIPATHWAYS"`,
  `"GO:BP"`, `"GO:MF"`, `"GO:CC"`, `"IMMUNESIGDB"`, `"CELLTYPE"`.

- gmt:

  Paths to `.gmt` files. Each becomes a collection named after the file
  (or after the element name, when `gmt` is a named vector).

- species:

  Species for MSigDB (e.g. `"Homo sapiens"`, `"Mus musculus"`; mouse
  uses msigdbr's ortholog-mapped sets). `NULL` guesses from the gene
  symbols.

- assay:

  Assay whose genes the sets are trimmed to (default: the project's
  default assay). Genes absent from it are dropped, and so are sets left
  empty.

## Value

Invisibly, the names of the collections saved.

## See also

[`scroll_gsea()`](https://shihanli92.github.io/scroll/reference/scroll_gsea.md)
