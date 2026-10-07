# A tour of the app

This article walks through a scroll app built from the public 10x
Genomics dataset of CD8+ T cells from four donors (gene expression, TCR
and antigen binding; 237,883 cells). Everything shown comes from one
call to
[`scroll_build()`](https://shihanli92.github.io/scroll/reference/scroll_build.md);
the repertoire panels appear because the build was given a
[`vdj_spec()`](https://shihanli92.github.io/scroll/reference/vdj_spec.md).

## The page

![The scroll app: app bar with the View selector and cell count, the
panel rail on the left, the DimPlot panel in the centre and the filter
rail on the right](figures/tour-overview.png)

- **App bar** (top): the dataset title, the **View** selector, the
  number of cells currently shown, dataset details, a two-panels-per-row
  toggle and the control-rail toggle.
- **Panel rail** (left): one entry per panel. Click to jump; drag to
  reorder.
- **Panels** (centre): each card has *what to plot* on the left and the
  plot on the right, with a toolbar for the Style sheet, export size,
  and PNG / PDF / CSV downloads.
- **Control rail** (right): filters that narrow every panel at once, and
  the *New column* tool.

Panels only draw while they are on screen, and only show up when the
data supports them (DE needs a grouping column, the repertoire panels
need TCR data). The `panels:` and `exclude_panels:` keys in
`config.yaml` override this.

## Looking at cells

**DimPlot** colours any embedding by a metadata column, with optional
highlighting of chosen groups and a split into facets.

![DimPlot panel: a UMAP of 237,883 cells coloured by
donor](figures/tour-dimplot.png)

**FeaturePlot** takes several genes (or numeric columns) at once and
draws a grid, each with its own colour scale. With two genes, *Blend*
shows co-expression.

![FeaturePlot panel: a grid of four UMAPs coloured by CD8A, GZMB, CCR7
and GNLY expression](figures/tour-featureplot.png)

**Signature** scores a gene list per cell: a plain mean, a z-score,
Seurat’s `AddModuleScore` (matched exactly), UCell or an AUCell-style
score. The score can be shown on the embedding or as a violin, and *Add
as column* makes it (and optional high / low groups) a column in every
other panel.

![Signature panel: a UCell score for a five-gene cytotoxicity signature
on the UMAP, with the Add as column
controls](figures/tour-signature.png)

**Biaxial** plots pairs of numeric columns against each other, e.g. QC
metrics or antibody-derived tags.

## Comparing groups

**DotPlot** shows the fraction of cells expressing each gene and the
mean expression per group, with optional clustering and dendrograms.
**Heatmap**, **Violin** and **Ridge** show the same genes as a heatmap,
as violins or as one density ridge per group. The DotPlot and Violin
gene boxes accept a pasted list.

These panels can group by several columns at once (e.g. cluster and
donor), and list groups in natural order, so cluster 2 comes before
cluster 10. To set your own order, drag the group labels in the plot’s
Style sheet.

![DotPlot panel: eight marker genes by donor, dot size showing percent
expressing and colour the scaled mean
expression](figures/tour-dotplot.png)

**Proportions** shows how one grouping is made up of another, e.g. each
cluster’s donor composition.

![Proportions panel: stacked bars of donor composition within each
cluster](figures/tour-proportions.png)

## Differential expression

**DE** runs a Wilcoxon test (presto) for any contrast: one or more
groups against another set or against the rest. The result appears as a
table and a volcano plot. The CSV can hold the rows on screen, only the
genes passing the cutoffs, or every gene.

![DE panel: a volcano plot for cluster 3 against all other cells, with
up and down genes coloured](figures/tour-de.png)

**Pseudobulk DE** sums counts per sample and tests with edgeR /
limma-voom, using replicates (e.g. donors) when you name them and
pseudo-replicates when you don’t. It needs a project built with
`counts = TRUE`.

**GSEA** takes the latest DE or Pseudobulk result and runs preranked
GSEA (fgsea) against gene sets saved in the project
([`scroll_add_genesets()`](https://shihanli92.github.io/scroll/reference/scroll_add_genesets.md))
or fetched from MSigDB. It shows a table, the top pathways, and the
running enrichment score for any pathway.

![GSEA panel: top Hallmark pathways by normalized enrichment score for
the cluster 3 contrast](figures/tour-gsea.png)

## Immune repertoire, spatial and ATAC

Projects built with
[`vdj_spec()`](https://shihanli92.github.io/scroll/reference/vdj_spec.md)
get six repertoire panels: clone overview (rank-abundance and
expansion), V/J gene usage, CDR3 length, diversity, a clone map on the
embedding and a CDR3 sequence logo.
[`spatial_spec()`](https://shihanli92.github.io/scroll/reference/spatial_spec.md)
adds a tissue-image panel and
[`atac_spec()`](https://shihanli92.github.io/scroll/reference/atac_spec.md)
a peak panel.

![Clone overview panel: clone rank-abundance on a log
scale](figures/tour-repertoire.png)

## Styling a plot

The paintbrush in a plot’s toolbar opens its **Style sheet**. Changes
apply live: palettes and point sizes, titles and axis labels, the theme
(fonts, legend, gridlines), axis ranges and scales, and the fixed
settings of every layer the plot draws. *Apply theme to all plots*
copies one plot’s theme to the rest, and **Save** writes every plot’s
style to `style.yaml` in the project so the next session starts from it.

![The DimPlot Style sheet open beside the plot, showing palette, label
and theme controls with Save and Apply theme to all plots
buttons](figures/tour-style.png)

## Views, filters and new columns

A project built with subset views
([`scroll_add_subset()`](https://shihanli92.github.io/scroll/reference/scroll_add_subset.md))
gets a **View** selector in the app bar. Choosing a view restricts every
panel to that subset and switches the scatter plots to its own
embedding; here, one donor re-embedded on its own.

![The app switched to the Donor 1 view: the DimPlot shows only that
donor on its own UMAP and the cell count reads 55,206 of
237,883](figures/tour-view.png)

The control rail on the right holds **Filters**, which narrow every
panel to the chosen levels or ranges when you click *Apply*, and **New
column**, which maps the levels of a column to new labels for the
session. For example, `[Donor 1, Donor 2]: A, [Donor 3, Donor 4]: B`
adds a `batch` column that every panel can colour, group or split by.

![The control rail: Filters for each column, and New column with a donor
mapping to a batch column](figures/tour-rail.png)

## Next steps

- [Getting
  started](https://shihanli92.github.io/scroll/articles/getting-started.md):
  build and configure a project.
- [Custom
  panel](https://shihanli92.github.io/scroll/articles/custom-panel-neighbours.md):
  add a panel of your own.
