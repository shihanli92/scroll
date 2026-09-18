read_cells <- function(dir) as.data.frame(arrow::read_parquet(file.path(dir, "cells.parquet")))

test_that("scatter views render ggplots", {
  cells <- read_cells(test_project())
  expect_s3_class(view_umap_colorby(cells, list(embedding = "umap", color_by = "celltype")), "ggplot")
  vals <- data.frame(cell = cells$cell[1:10], value = runif(10))
  expect_s3_class(view_feature_plot(cells, list(embedding = "umap", feature = "CD3D"), vals), "ggplot")
})

test_that("view_feature_blend reproduces Seurat's blend colours and lays out 4 panels", {
  skip_if_not_installed("patchwork")
  cells <- read_cells(test_project())
  v1 <- data.frame(cell = cells$cell, value = pmax(0, stats::rnorm(nrow(cells), 1)))
  v2 <- data.frame(cell = cells$cell, value = pmax(0, stats::rnorm(nrow(cells), 1)))
  p <- view_feature_blend(cells, list(embedding = "umap", feature1 = "CD3D",
                                      feature2 = "MS4A1"), v1, v2)
  expect_s3_class(p, "patchwork")
  expect_length(p$patches$plots, 3)                     # + the 4th plot is the top-level object

  # the blend colour matrix is a faithful port of Seurat's internal BlendMatrix
  skip_if_not_installed("Seurat")
  ours <- scroll:::.scroll_blend_matrix(col.threshold = 0.5)
  seur <- getFromNamespace("BlendMatrix", "Seurat")(
    two.colors = c("#ff0000", "#00ff00"), col.threshold = 0.5, negative.color = "lightgrey")
  expect_identical(as.vector(ours), as.vector(seur))
})

test_that("view_feature_multi renders a grid, one panel per gene", {
  skip_if_not_installed("patchwork")
  cells <- read_cells(test_project())
  feats <- c("CD3D", "MS4A1", "CD8A")
  vl <- do.call(rbind, lapply(feats, function(g)
    data.frame(feature = g, cell = cells$cell, value = runif(nrow(cells)))))
  p <- view_feature_multi(cells, list(embedding = "umap", features = feats), vl)
  expect_s3_class(p, "patchwork")
  expect_length(p$patches$plots, length(feats) - 1)     # + top-level object = 3 panels
  # a single gene collapses to one plain feature plot, not a patchwork
  p1 <- view_feature_multi(cells, list(embedding = "umap", features = "CD3D"),
                           vl[vl$feature == "CD3D", ])
  expect_s3_class(p1, "ggplot")
  expect_false(inherits(p1, "patchwork"))
})

test_that(".scroll_parse_gene_list splits on any whitespace/comma, matches case-insensitively", {
  feats <- c("CD3D", "CD8A", "MS4A1", "NKG7")
  res <- scroll:::.scroll_parse_gene_list("CD3D, cd8a\tMS4A1\nNKG7 CD3D foo", feats)
  expect_equal(res$ok, c("CD3D", "CD8A", "MS4A1", "NKG7"))   # pasted order, deduped, canonical case
  expect_equal(res$missing, "foo")
  expect_equal(scroll:::.scroll_parse_gene_list("", feats)$ok, character(0))
})

test_that("dotplot aggregates fraction and mean per group", {
  cells <- read_cells(test_project())
  feats <- c("CD3D", "MS4A1")
  el <- data.frame(
    feature = c("CD3D", "CD3D", "MS4A1"),
    cell = cells$cell[1:3],
    value = c(1, 2, 3))
  p <- view_dotplot(cells, list(group_by = "celltype", features = feats), el)
  expect_s3_class(p, "ggplot")
  # grid is complete: every feature x group combination present
  expect_equal(nrow(p$data), length(feats) * length(unique(cells$celltype)))
})

test_that("dotplot accepts brewer palettes and hclust clustering", {
  cells <- read_cells(test_project())
  expect_s3_class(scroll:::.scroll_continuous_scale("RdBu"), "Scale")   # ColorBrewer
  feats <- c("CD3D", "CD8A", "MS4A1", "NKG7", "GNLY")
  el <- data.frame(feature = rep(feats, each = 4), cell = rep(cells$cell[1:4], 5),
                   value = runif(20))
  p <- view_dotplot(cells, list(group_by = "celltype", features = feats), el,
                    state = list(cluster = "both", palette = "RdBu", aspect = 1.2, scale = TRUE))
  expect_true(inherits(p, c("ggplot", "aplot", "patchwork")))   # trees -> aplot, else ggplot
})

test_that("view_violin_multi renders one violin panel per gene", {
  skip_if_not_installed("patchwork")
  cells <- read_cells(test_project())
  feats <- c("CD3D", "MS4A1", "NKG7")
  vl <- do.call(rbind, lapply(feats, function(g)
    data.frame(feature = g, cell = cells$cell, value = runif(nrow(cells)))))
  p <- view_violin_multi(cells, list(group_by = "celltype", features = feats), vl)
  expect_s3_class(p, "patchwork")
  expect_length(p$patches$plots, length(feats) - 1)      # + top-level = 3 panels
  p1 <- view_violin_multi(cells, list(group_by = "celltype", features = "CD3D"),
                          vl[vl$feature == "CD3D", ])
  expect_s3_class(p1, "ggplot")
  expect_false(inherits(p1, "patchwork"))
})

test_that("violin and proportions render ggplots", {
  cells <- read_cells(test_project())
  vals <- data.frame(cell = cells$cell[1:20], value = runif(20))
  expect_s3_class(view_violin(cells, list(group_by = "celltype", feature = "CD3D"), vals), "ggplot")
  expect_s3_class(view_proportions(cells, list(group_by = "condition", fill_by = "celltype")), "ggplot")
})

test_that("violin honours palette/jitter/legend; proportions honours normalize", {
  cells <- read_cells(test_project())
  vals <- data.frame(cell = cells$cell[1:20], value = runif(20))
  expect_s3_class(
    view_violin(cells, list(feature = "CD3D", group_by = "celltype"), vals,
                state = list(palette = "Set1", jitter = TRUE, legend = TRUE)),
    "ggplot")
  pc <- view_proportions(cells, list(group_by = "condition", fill_by = "celltype"),
                         state = list(normalize = FALSE))
  expect_s3_class(pc, "ggplot")
  expect_equal(pc$labels$y, "cells")     # raw counts, not "composition"
})

test_that("violin split-by draws grouped violins; group-by can be composite", {
  cells <- read_cells(test_project())
  vals <- data.frame(cell = cells$cell, value = runif(nrow(cells)))
  p <- view_violin(cells, list(feature = "CD3D", group_by = "celltype", split_by = "condition"), vals)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$fill, "condition")           # coloured by the split, not the group
  # composite group-by (interaction)
  p2 <- view_violin(cells, list(feature = "CD3D", group_by = c("celltype", "condition")), vals)
  expect_equal(p2$labels$x, "celltype | condition")
})

test_that("view_violin point controls subsample and apply size/opacity", {
  cells <- read_cells(test_project())
  vals <- data.frame(cell = cells$cell, value = runif(nrow(cells)))
  p <- view_violin(cells, list(feature = "X", group_by = "celltype"), vals,
                   state = list(jitter = TRUE, point_size = 1, point_alpha = 0.5, point_frac = 0.1))
  b <- ggplot2::ggplot_build(p)
  npts <- nrow(b$data[[2]])                 # the jitter layer
  expect_gt(npts, 0)
  expect_lt(npts, nrow(cells))              # subsampled to ~10%
})

test_that("view_violin_stacked facets one compact row per gene", {
  data <- scroll:::.scroll_load(test_project()); on.exit(scroll_disconnect(data$con))
  vl <- data$queryN("RNA", c("CD3D", "CD8A"))
  p <- view_violin_stacked(data$cells, list(group_by = "celltype",
                                            features = c("CD3D", "CD8A")), vl)
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  expect_gte(length(unique(b$layout$layout$feature)), 2)   # one facet row per gene
})

test_that("proportions position controls stack / fill / dodge", {
  cells <- read_cells(test_project())
  base <- list(group_by = "condition", fill_by = "celltype")
  expect_equal(view_proportions(cells, base, state = list(position = "stack"))$labels$y, "cells")
  expect_equal(view_proportions(cells, base, state = list(position = "fill"))$labels$y, "composition")
  expect_s3_class(view_proportions(cells, base, state = list(position = "dodge")), "ggplot")
  # back-compat: old `normalize` flag still maps to fill/stack
  expect_equal(view_proportions(cells, base, state = list(normalize = FALSE))$labels$y, "cells")
})

test_that("proportions facets one composition per fill column (grid/stacked)", {
  cells <- read_cells(test_project())
  fills <- c("celltype", "condition")
  is_multi <- function(p) if (requireNamespace("patchwork", quietly = TRUE))
    expect_s3_class(p, "patchwork") else expect_s3_class(p, "ggplot")
  is_multi(view_proportions(cells, list(group_by = "condition", fill_by = fills), state = list(fill_layout = "grid")))
  is_multi(view_proportions(cells, list(group_by = "condition", fill_by = fills), state = list(fill_layout = "stacked")))
  is_multi(view_proportions(cells, list(group_by = "condition", fill_by = fills), state = list(facet = TRUE)))  # back-compat
  # combine (default) stays a single combined-interaction plot
  pc <- view_proportions(cells, list(group_by = "condition", fill_by = fills), state = list(fill_layout = "combine"))
  expect_equal(pc$labels$fill, "celltype | condition")
})

test_that("proportions aesthetics: ordering, labels, orientation", {
  cells <- read_cells(test_project())
  base <- list(group_by = "condition", fill_by = "celltype")
  p <- view_proportions(cells, base, state = list(x_order = "total", fill_order = "abundance",
                                                  labels = "percent", horizontal = TRUE, outline = 0))
  expect_s3_class(p, "ggplot")
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_true("GeomText" %in% geoms)              # a label layer was added
  expect_s3_class(p$coordinates, "CoordFlip")     # horizontal
  # count labels also work in stack position
  pc <- view_proportions(cells, base, state = list(position = "stack", labels = "count"))
  expect_true("GeomText" %in% vapply(pc$layers, function(l) class(l$geom)[1], character(1)))
})

test_that("proportions honours bar width", {
  cells <- read_cells(test_project())
  p <- view_proportions(cells, list(group_by = "condition", fill_by = "celltype"), state = list(bar_width = 0.5))
  b <- ggplot2::ggplot_build(p)
  expect_s3_class(p, "ggplot")
  expect_true(any(abs(b$data[[1]]$xmax - b$data[[1]]$xmin - 0.5) < 1e-6))  # bars are 0.5 wide
})

test_that("proportions fills by the interaction of >1 column", {
  cells <- read_cells(test_project())
  p <- view_proportions(cells, list(group_by = "condition", fill_by = c("celltype", "condition")))
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$fill, "celltype | condition")          # composite fill
  # the fill categories are the observed interaction levels
  fills <- unique(as.character(ggplot2::ggplot_build(p)$plot$data$fill))
  expect_true(any(grepl(" | ", fills, fixed = TRUE)))
})

test_that("view_volcano plots the requested fold-change column", {
  de <- data.frame(gene = c("A", "B", "C"),
                   logFC = c(0.2, -0.2, 0.1), avg_log2FC = c(2, -2, 0),
                   p_val_adj = c(1e-5, 1e-5, 0.9))
  expect_equal(view_volcano(de)$labels$x, "logFC")                       # default
  p <- view_volcano(de, params = list(fc_col = "avg_log2FC", fc_label = "avg_log2FC", lfc = 1))
  expect_equal(p$labels$x, "avg_log2FC")
  # up/down classification now uses avg_log2FC (|2| >= 1), not logFC (|0.2| < 1)
  sig <- ggplot2::ggplot_build(p)$plot$data$sig
  expect_setequal(as.character(sig), c("up", "down", "ns"))
})

test_that("de_table returns a data.frame (and a message when absent)", {
  df <- data.frame(gene = c("A", "B"), avg_log2FC = c(1.2, -0.5))
  expect_s3_class(view_de_table(df, list(top = 1)), "data.frame")
  expect_equal(nrow(view_de_table(df, list(top = 1))), 1)
  expect_true("message" %in% names(view_de_table(NULL, list(contrast = "x"))))
})

test_that("render_view dispatches by ctx and honours feature override", {
  cells <- read_cells(test_project())
  ctx <- list(view = "umap_colorby", params = list(embedding = "umap", color_by = "condition"),
              cells = cells)
  expect_s3_class(render_view("umap_colorby", ctx), "ggplot")

  vals <- data.frame(cell = cells$cell[1:5], value = runif(5))
  ctx2 <- c(ctx, list(feature_values = vals, feature_name = "CD3D"))
  # feature override routes a colorby section through the feature plot
  expect_true(".expr" %in% names(render_view("umap_colorby", ctx2)$data))

  expect_equal(scroll_view_kind("de_table"), "table")
  expect_equal(scroll_view_kind("violin"), "plot")
  expect_error(render_view("nope", ctx), "Unknown view type")
})

test_that("view cores honour palette / size / order / scale controls", {
  cells <- read_cells(test_project())
  expect_s3_class(
    view_umap_colorby(cells, list(embedding = "umap", color_by = "celltype"),
                      state = list(palette = "Okabe-Ito", point_size = 1.2, alpha = 0.5)),
    "ggplot")
  vals <- data.frame(cell = cells$cell[1:10], value = runif(10))
  expect_s3_class(
    view_feature_plot(cells, list(embedding = "umap", feature = "CD3D"), vals,
                      state = list(palette = "magma", order = FALSE, point_size = 1)),
    "ggplot")
  el <- data.frame(feature = c("CD3D", "CD3D", "MS4A1"), cell = cells$cell[1:3],
                   value = c(1, 2, 3))
  expect_s3_class(
    view_dotplot(cells, list(group_by = "celltype", features = c("CD3D", "MS4A1")), el,
                 state = list(scale = TRUE, dot_size = c(2, 8), palette = "viridis")),
    "ggplot")
})

test_that("DimPlot legend toggle wins even when cluster labels are shown", {
  cells <- read_cells(test_project())
  base <- list(embedding = "umap", color_by = "celltype")
  # labels on + legend on -> legend visible (labels no longer suppress it)
  p_on <- view_umap_colorby(cells, base, list(show_labels = TRUE, legend = TRUE))
  expect_identical(p_on$theme$legend.position, "right")
  # labels on + legend off -> no legend
  p_off <- view_umap_colorby(cells, base, list(show_labels = TRUE, legend = FALSE))
  expect_identical(p_off$theme$legend.position, "none")
})

test_that("DimPlot highlight greys unselected groups; manual colors override", {
  cells <- read_cells(test_project())
  # highlight one cell type: plot keeps a grey background layer + the fg layer
  p <- view_umap_colorby(cells, list(embedding = "umap", color_by = "celltype"),
                         state = list(highlight = "B"))
  expect_s3_class(p, "ggplot")
  expect_gte(length(p$layers), 2)                 # background + foreground point layers
  # manual palette override changes a group's assigned color
  cols <- scroll:::.scroll_group_colors(c("B", "T", "NK"),
            list(palette = "Manual", manual_colors = c(B = "#123456")))
  expect_equal(unname(cols[["B"]]), "#123456")
})

test_that("FeaturePlot quantile caps clip the color scale", {
  cells <- read_cells(test_project())
  expect_null(scroll:::.scroll_expr_limits(1:100, c(0, 1)))          # no clip
  lim <- scroll:::.scroll_expr_limits(1:100, c(0.1, 0.9))
  expect_length(lim, 2)
  expect_lt(lim[1], lim[2])
  vals <- data.frame(cell = cells$cell, value = seq_len(nrow(cells)))
  p <- view_feature_plot(cells, list(embedding = "umap", feature = "CD3D"), vals,
                         state = list(clip = c(0.1, 0.9)))
  expect_s3_class(p, "ggplot")
})

test_that("quantile clip is taken over expressing cells, not the zero-inflated vector", {
  # a 95%-zero gene: the 0.5 quantile of the full vector is 0, so a naive clip
  # would pin the lower bound at 0; over non-zero cells it is a real value.
  ex <- c(rep(0, 950), seq_len(50))
  lim <- scroll:::.scroll_expr_limits(ex, c(0.5, 1.0))
  expect_false(is.null(lim))
  expect_gt(lim[1], 0)                               # lower bound moved off zero
  expect_equal(lim[1], unname(stats::quantile(1:50, 0.5)))
})

test_that("aspect ratio sets theme(aspect.ratio); default 1 leaves it unset", {
  cells <- read_cells(test_project())
  p <- view_umap_colorby(cells, list(embedding = "umap", color_by = "celltype"),
                         state = list(aspect = 2.5))
  expect_equal(p$theme$aspect.ratio, 2.5)
  p1 <- view_umap_colorby(cells, list(embedding = "umap", color_by = "celltype"),
                          state = list(aspect = 1))
  expect_equal(p1$theme$aspect.ratio, 1)          # embedding scatter is square at aspect 1
  # non-scatter panels leave aspect 1 unconstrained (fill the width), but honour != 1
  v1 <- view_violin(cells, list(feature = "CD3D", group_by = "celltype"), NULL,
                    state = list(aspect = 1))
  expect_null(v1$theme$aspect.ratio)
  v <- view_violin(cells, list(feature = "CD3D", group_by = "celltype"), NULL,
                   state = list(aspect = 0.5))
  expect_equal(v$theme$aspect.ratio, 0.5)
})

test_that("legend can be switched off", {
  cells <- read_cells(test_project())
  p <- view_umap_colorby(cells, list(embedding = "umap", color_by = "celltype"),
                         state = list(legend = FALSE))
  expect_equal(p$theme$legend.position, "none")
})

test_that("discrete colors are deterministic by level name", {
  a <- scroll:::.scroll_discrete_colors(c("B", "T", "NK"))
  b <- scroll:::.scroll_discrete_colors(c("NK", "B"))   # subset, different order
  expect_equal(a[["B"]], b[["B"]])                       # B keeps its color
  expect_equal(a[["NK"]], b[["NK"]])
})

test_that("toggles: split_by facets and overrides group_by", {
  cells <- read_cells(test_project())
  p <- view_umap_colorby(cells, list(embedding = "umap", color_by = "celltype"),
                         state = list(split_by = "condition"))
  expect_s3_class(p, "ggplot")
  # violin group column follows the split_by toggle
  v <- view_violin(cells, list(feature = "CD3D"), NULL, state = list(split_by = "condition"))
  expect_true(all(v$data$group %in% as.character(cells$condition)))
})

test_that("view_contrast_preview draws red/blue over a grey outline (placeholder on empty)", {
  set.seed(1); n <- 200
  cells <- data.frame(umap_1 = rnorm(n), umap_2 = rnorm(n),
                      ct = sample(c("A", "B", "C"), n, replace = TRUE),
                      stringsAsFactors = FALSE)
  lab <- scroll:::.scroll_contrast_labels(cells, "ct", "A", "B")
  expect_s3_class(view_contrast_preview(cells, "umap", lab), "ggplot")
  # all-NA labels -> outline only, still a ggplot
  expect_s3_class(view_contrast_preview(cells, "umap", rep(NA_character_, n)), "ggplot")
  # no coordinates -> friendly placeholder, no error
  expect_s3_class(view_contrast_preview(data.frame(ct = "A"), "umap", "group1"), "ggplot")
})

test_that("view_contrast_medoids emits exactly one point per sample", {
  set.seed(1); n <- 300
  cells <- data.frame(umap_1 = rnorm(n), umap_2 = rnorm(n),
                      samp = sample(paste0("s", 1:6), n, replace = TRUE),
                      stringsAsFactors = FALSE)
  grp <- ifelse(cells$samp %in% c("s1", "s2", "s3"), "group1", "group2")
  p <- view_contrast_medoids(cells, "umap", cells$samp, grp)
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  expect_equal(nrow(b$data[[length(b$data)]]), 6L)   # medoid layer: one point per sample
})

test_that("plot-source CSV builders carry barcodes + reproducer columns", {
  skip_if_not_installed("SeuratObject")
  dir <- test_project()
  cells <- as.data.frame(arrow::read_parquet(file.path(dir, "cells.parquet")))
  cells$.gidx <- seq_len(nrow(cells))
  emb <- "umap"

  # DimPlot source: cell barcode + the two embedding coords + the color-by column
  d <- .scroll_dimplot_source(cells, emb, "celltype")
  expect_true(all(c("cell", paste0(emb, c("_1", "_2")), "celltype") %in% names(d)))
  expect_true(all(d$cell %in% cells$cell))              # barcodes are real
  expect_equal(anyNA(d[[paste0(emb, "_1")]]), FALSE)    # NA-coord rows dropped

  # FeaturePlot source: barcode + coords + a feature-expression column, 0-filled
  con <- scroll_connect(dir); on.exit(scroll_disconnect(con))
  vals <- scroll_query_feature(con, "RNA", "CD3D")
  f <- .scroll_featureplot_source(cells, emb, "CD3D", vals)
  expect_true(all(c("cell", "CD3D") %in% names(f)))
  expect_false(anyNA(f$CD3D))

  # Violin source: barcode-first, group column, feature column
  v <- .scroll_violin_source(cells, "celltype", "CD3D", vals, value_col = NULL)
  expect_identical(names(v)[1], "cell")
  expect_true(all(c("celltype", "CD3D") %in% names(v)))

  # Aggregated panels: no per-cell barcode, the plotted summary instead
  p <- .scroll_proportions_source(cells, "celltype", "condition")
  expect_setequal(names(p), c("group", "category", "n_cells", "proportion"))
  expect_false("cell" %in% names(p))
})
