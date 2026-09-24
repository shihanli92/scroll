# Custom panel

## 1. Introduction

scroll allows for custom panels to be attached, extending usability
depending on the data context. This article details the process of
creating and attaching a custom panel.

For this panel we will leverage the k-nearest-neighbour graphs generated
by the Seurat `FindNeighbors` function. As this data is not exported in
the default build function, we will first create and save the graphs for
later accessibility (baking). The bake step is advised for new data
sources to keep the memory used during the app runtime as minimal as
possible. The alternative is to extract the principal components so we
can calculate the nearest neighbours on demand, but this would greatly
increase memory and load time when using the app.

## 2. Prerequisites

We will use the pbmc3k dataset from the first article.

- A built `pbmc3k` project — the `pbmc3k/` directory. The baked graph
  goes inside it.
- The `pbmc` Seurat object
- Seurat installed to run `FindNeighbors`

``` r

library(scroll)
library(Seurat)


pbmc <- SeuratData::LoadData("pbmc3k", type = "pbmc3k.final")
pbmc <- UpdateSeuratObject(pbmc)
```

## 3. Bake the neighbour graph

For this article we are interested in obtaining the Euclidean distances
between cells. To get these we need to rerun
[`FindNeighbors()`](https://satijalab.org/seurat/reference/FindNeighbors.html)
with `return.neighbor = TRUE`; otherwise the call only stores a binary
adjacency on `pbmc@graphs$RNA_nn`. The graph is flattened into a table
of cell_barcode \| neighbour_cell_barcode \| distance, then saved as a
Parquet table.

``` r

# 1. k-NN in expression/PCA space, keeping distances (creates pbmc@neighbors$RNA.nn)
pbmc <- FindNeighbors(pbmc, dims = 1:10, return.neighbor = TRUE)

# 2. flatten the neighbour matrices into a long edge list of barcodes
nn   <- pbmc@neighbors$RNA.nn              # a Neighbor object
idx  <- nn@nn.idx                          # cells x k : neighbour row indices
dst  <- nn@nn.dist                         # cells x k : matching distances
bc   <- nn@cell.names                      # barcodes, row-aligned to idx / dst
k    <- ncol(idx)

edges <- data.frame(
  cell      = rep(bc, times = k),          # the seed cell (column-major order)
  neighbour = bc[as.vector(idx)],          # one of its k neighbours (barcode)
  distance  = as.vector(dst)               # the k-NN distance between them
)
edges <- edges[edges$cell != edges$neighbour, ]   # drop self (distance 0)
nrow(edges)

dir.create("pbmc3k/neighbours", showWarnings = FALSE)
arrow::write_parquet(edges, "pbmc3k/neighbours/RNA_nn.parquet")
```

    #> [1] 50122

## 4. Write the neighbourhood panel

Each plot panel consists of two parts: a set of controls, which give
users control over plot adjustments and aesthetics, and a plotting
function that renders the ggplot. The
[`register_plot_panel()`](https://shihanli92.github.io/scroll/reference/register_plot_panel.md)
function then attaches these to the app. For this article, we will
slowly build up the nearest-neighbour panel with increasing complexity.

### 4.1 The simplest version

A range of Shiny widgets are already implemented in scroll for quick
use. To start, we add a column input, which lets a user select a
metadata column, and a level input linked to the column selector to
filter for specific factors within that column.

``` r

controls <- list(
  scroll_input_column("group", "Highlight column", type = "categorical"),
  scroll_input_levels("level", "Highlight group",  from = "group")
)
```

The plot function receives `cells`, the `input` values, and `data` — the
last used to find the baked graph via `data$dir` — and returns a ggplot.
Note that while we use a Parquet file to store the data, realistically
any readable data format works for custom panels, with the caveat that
memory constraints need to be considered.

``` r

neighbourhood_plot <- function(cells, input, data) {
  if (!length(input$level)) stop("Pick at least one highlight group.")

  # read the graph baked in step 3, resolved from the project root
  edges <- arrow::read_parquet(
    file.path(data$dir, "neighbours", "RNA_nn.parquet"))

  # seeds = active cells in the chosen group; highlight their neighbours
  seed <- cells$cell[as.character(cells[[input$group]]) %in% input$level]
  nb   <- unique(edges$neighbour[edges$cell %in% seed])

  df <- data.frame(
    x         = cells$umap_1,
    y         = cells$umap_2,
    neighbour = cells$cell %in% nb
  )

  ggplot(df, aes(x, y)) +
    geom_point(colour = "grey80", size = 0.4) +
    geom_point(data = df[df$neighbour, ], colour = "red", size = 0.7) +
    theme_classic() +
    labs(x = NULL, y = NULL,
         title = paste(input$group, paste(input$level, collapse = ", ")))
}
```

Call
[`register_plot_panel()`](https://shihanli92.github.io/scroll/reference/register_plot_panel.md)
to add the panel to the app. Adding `after = "featureplot"` places the
panel after FeaturePlot, and `compute = TRUE` adds a Compute button that
stops the graph from updating every time an option is changed.

``` r

library(scroll)
library(ggplot2)

register_plot_panel("neighbourhood",
  label    = "Neighbourhood",
  title    = "Neighbourhood explorer",
  desc     = "Highlight a group's k-NN graph neighbours on the UMAP.",
  after    = "featureplot",
  compute  = TRUE,
  controls = controls,
  plot     = neighbourhood_plot)

scroll_serve("pbmc3k")
```

### 4.2 Let the reader pick the embedding

The UMAP is hard-coded above. One more control lets the reader view the
same neighbourhood on any embedding in the project — PCA, or a subset
re-embedding. Add an embedding picker to `controls` (a third entry
alongside the two from 4.1):

``` r

controls <- list(
  scroll_input_column("group", "Highlight column", type = "categorical"),
  scroll_input_levels("level", "Highlight group",  from = "group"),
  scroll_input_embedding("emb", "Embedding")
)
```

scroll names an embedding’s coordinate columns `<name>_1` / `<name>_2`,
so read the chosen one instead of the fixed `umap_1` / `umap_2` — the
only change is the two coordinate lines in `df`:

``` r

  df <- data.frame(
    x         = cells[[paste0(input$emb, "_1")]],
    y         = cells[[paste0(input$emb, "_2")]],
    neighbour = cells$cell %in% nb
  )
```

### 4.3 Offer points or density

For a dense group, individual red points crowd together; 2-D density
contours read more clearly. This input control lets the user pick
whichever they prefer.

``` r

scroll_input_choice("style", "Style",
                    choices = c("Points", "Density"), selected = "Points")
```

Note that this control returns a value that we branch on in the plot
function.

``` r
  hits <- df[df$neighbour, ]
  if (!nrow(hits)) stop("This group has no neighbours among the shown cells.")

  p <- ggplot(df, aes(x, y)) +
    geom_point(colour = "grey80", size = 0.4) +
    theme_classic() +
    labs(x = NULL, y = NULL,
         title = paste(input$group, paste(input$level, collapse = ", ")))

  if (identical(input$style, "Density"))
    p + geom_density_2d(data = hits, colour = "red")
  else
    p + geom_point(data = hits, colour = "red", size = 0.7)
```

The Points branch still paints every neighbour the same flat red.

### 4.4 Colour neighbours by distance

Every edge carries the k-NN `distance` between a seed and its neighbour.
A neighbour reached from several seeds has several distances, so we will
use the smallest and carry it on `df` alongside the neighbour flag:

``` r

  # nearest-seed distance for each cell (NA when it isn't a neighbour)
  hit_edges <- edges[edges$cell %in% seed, ]
  dmin      <- tapply(hit_edges$distance, hit_edges$neighbour, min)

  df <- data.frame(
    x         = cells[[paste0(input$emb, "_1")]],
    y         = cells[[paste0(input$emb, "_2")]],
    neighbour = cells$cell %in% nb,
    dist      = dmin[cells$cell]
  )
```

Then, in the Points branch, map colour to `dist` on a continuous scale
instead of a flat red:

``` r

  p + geom_point(data = hits, aes(colour = dist), size = 0.7) +
      scale_colour_viridis_c(name = "distance")
```

### 4.5 Let the reader pick the gradient colours

Users may want their own colour ramps. A colour selector is not a
built-in control, so we use `scroll_input_custom`, which wraps any Shiny
input into a control. Here we add two
[`colourpicker::colourInput`](https://rdrr.io/pkg/colourpicker/man/colourInput.html)s,
one for each end of the gradient.

``` r
scroll_input_custom("lo", "Near colour",
  ui = function(ns, data)
    colourpicker::colourInput(ns("lo"), "Near colour", value = "#2C7FB8")),
scroll_input_custom("hi", "Far colour",
  ui = function(ns, data)
    colourpicker::colourInput(ns("hi"), "Far colour", value = "#EDF8B1"))
```

Each `ui` namespaces its input with `ns("lo")` / `ns("hi")`, so the
value reaches the plot function as `input$lo` / `input$hi`. Swap the
fixed viridis scale for a gradient between the two chosen colours:

``` r

      scale_colour_gradient(low = input$lo, high = input$hi, name = "distance")
```

> This registers in memory only.
> [`register_plot_panel()`](https://shihanli92.github.io/scroll/reference/register_plot_panel.md)
> just adds the panel to the current session’s registry. For a
> reproducible, deployable app, put it in its own file.

## 5. Assemble it as an app

For deployment, custom panels should be saved into a panels directory:

    your_app/
    ├── app.R                 # library(scroll); source the panel; scroll_app("pbmc3k")
    ├── panels/
    │   └── neighbourhood.R   # controls + plot fn + register_plot_panel()
    └── pbmc3k/               # the built project, with neighbours/ baked in

Create the folder — scroll won’t make it for you — and drop the panel
code into it:

``` r

dir.create("panels")
```

`panels/neighbourhood.R` holds the finished panel from section 4 — all
six controls, the plot function, and the registration:

``` r

# panels/neighbourhood.R — a neighbourhood explorer (SCneighbours port)

controls <- list(
  scroll_input_column("group", "Highlight column", type = "categorical"),
  scroll_input_levels("level", "Highlight group", from = "group"),
  scroll_input_embedding("emb",   "Embedding"),
  scroll_input_choice("style",  "Style", choices = c("Points", "Density"), selected = "Points"),
  scroll_input_custom("lo", "Near colour",
    ui = function(ns, data) colourpicker::colourInput(ns("lo"), "Near colour", value = "#2C7FB8")),
  scroll_input_custom("hi", "Far colour",
    ui = function(ns, data) colourpicker::colourInput(ns("hi"), "Far colour", value = "#EDF8B1"))
)

neighbourhood_plot <- function(cells, input, data) {
  # ... the finished body from section 4 ...
}

register_plot_panel("neighbourhood",
  label = "Neighbourhood", title = "Neighbourhood explorer",
  desc  = "Highlight a group's k-NN graph neighbours on the embedding.",
  after = "featureplot", compute = TRUE,
  controls = controls, plot = neighbourhood_plot)
```

`app.R` attaches the packages, sources the panel (which registers it),
then returns the app:

``` r

# app.R
library(scroll)
library(ggplot2)

scroll_reset_panels()                # clear the registry (safe to re-run)
source("panels/neighbourhood.R")     # registers the panel

scroll_app("pbmc3k")                 # the app object Shiny runs
```

Source order matters:
[`register_plot_panel()`](https://shihanli92.github.io/scroll/reference/register_plot_panel.md)
must run *before*
[`scroll_app()`](https://shihanli92.github.io/scroll/reference/scroll_app.md)
builds the UI, which is why the
[`source()`](https://rdrr.io/r/base/source.html) comes first. The
[`scroll_reset_panels()`](https://shihanli92.github.io/scroll/reference/scroll_reset_panels.md)
at the top clears any prior registration so re-running `app.R` doesn’t
add the panel twice. Launch it from the app folder:

``` r

shiny::runApp(".")   # or open app.R in RStudio and click "Run App"
```
