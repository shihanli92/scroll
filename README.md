# scroll

*Interactive single-cell explorers, served from your own infrastructure.*

<!-- badges: start -->
[![R-CMD-check](https://github.com/shihanli92/scroll/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/shihanli92/scroll/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/shihanli92/scroll/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/shihanli92/scroll/actions/workflows/pkgdown.yaml)
<!-- badges: end -->

`scroll` turns a processed Seurat object into an interactive web explorer. It
builds the object once into small files on disk; the app reads only those, so it
never loads the Seurat object and its memory stays flat at any dataset size. The
result deploys as a plain `app.R` on a Shiny Server.

## Install

```r
remotes::install_github("shihanli92/scroll")
install.packages("SeuratObject")   # needed to build, not to run the app
```

## Build an app

```r
library(scroll)
scroll_build(seurat_obj, "myapp")
scroll_serve("myapp")
```

To deploy, copy the `myapp/` folder to a Shiny Server.

## Panels

DimPlot, FeaturePlot, Signature scores, Biaxial, DotPlot, Heatmap, Violin,
Proportions, DE, Pseudobulk DE and GSEA, plus repertoire (TCR/BCR), spatial and
ATAC panels when the data has them. Each plot has a Style sheet for its look and
exports to PNG or PDF.

## Documentation

- [Getting started](https://shihanli92.github.io/scroll/articles/getting-started.html):
  build options, `config.yaml`, subset views
- [Custom panel](https://shihanli92.github.io/scroll/articles/custom-panel-neighbours.html):
  add your own panel
- [Function reference](https://shihanli92.github.io/scroll/reference/)
