# scroll: Interactive Single-Cell Explorers from Seurat Objects

Turns a processed Seurat object into a polished, interactive single-cell
explorer: a Shiny app of eight analysis panels (DimPlot, FeaturePlot,
Biaxial, DotPlot, Violin, Proportions, and live differential expression
including replicate-aware pseudobulk), each with its own controls. A
heavy offline build phase extracts a compact, sparse on-disk store (cell
metadata and embeddings as one Parquet table, expression as a one
feature-sorted Parquet file per assay) that the app reads one feature at
a time via 'arrow', so runtime memory stays flat regardless of matrix
size. A streaming builder assembles multi-million-cell projects one
source at a time.

## See also

Useful links:

- <https://github.com/shihanli92/scroll>

- <https://shihanli92.github.io/scroll/>

- Report bugs at <https://github.com/shihanli92/scroll/issues>

## Author

**Maintainer**: Shihan Li <shihanli1992@gmail.com>
