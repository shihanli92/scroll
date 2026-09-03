# On-screen scatter downsampling: caps the interactive (rasterized) draw while
# keeping the focus rows (highlighted / expressing) and drawing full-res on export.

test_that(".scroll_downsample caps on-screen points, keeps flagged rows, deterministic", {
  df <- data.frame(x = seq_len(300), y = rnorm(300),
                   grp = rep(c("a", "b", "c"), 100), stringsAsFactors = FALSE)
  st <- list(raster = TRUE)
  # off-screen / export (raster FALSE): untouched
  expect_equal(nrow(scroll:::.scroll_downsample(df, list(raster = FALSE), cap = 100)), 300)
  # below the cap: untouched
  expect_equal(nrow(scroll:::.scroll_downsample(df, st, cap = 500)), 300)
  # cap enforced + deterministic
  d1 <- scroll:::.scroll_downsample(df, st, cap = 100)
  d2 <- scroll:::.scroll_downsample(df, st, cap = 100)
  expect_equal(nrow(d1), 100)
  expect_identical(d1, d2)
  # flagged rows always kept, and fill the remaining budget with a subsample
  keep <- df$grp == "a"                              # 100 flagged rows
  d3 <- scroll:::.scroll_downsample(df, st, keep = keep, cap = 150)
  expect_equal(nrow(d3), 150)
  expect_true(all(df$x[keep] %in% d3$x))             # every flagged row retained
  # the global RNG is left undisturbed
  set.seed(42); r1 <- runif(1)
  set.seed(42); invisible(scroll:::.scroll_downsample(df, st, cap = 100)); r2 <- runif(1)
  expect_equal(r1, r2)
})

test_that("scatter views downsample on screen but keep highlighted / expressing cells", {
  old <- options(scroll.onscreen_cap = 150); on.exit(options(old))
  n <- 300
  cells <- data.frame(cell = paste0("c", seq_len(n)), umap_1 = rnorm(n), umap_2 = rnorm(n),
                      grp = rep(c("x", "y", "z"), 100), stringsAsFactors = FALSE)
  # DimPlot: highlighted group fully retained, total capped
  p <- view_umap_colorby(cells, list(embedding = "umap", color_by = "grp"),
                         list(raster = TRUE, highlight = "x"))
  expect_lte(nrow(p$data), 150)
  expect_true(all(cells$cell[cells$grp == "x"] %in% p$data$cell))   # all highlighted kept
  # export (raster FALSE) draws every point
  pf <- view_umap_colorby(cells, list(embedding = "umap", color_by = "grp"), list(raster = FALSE))
  expect_equal(nrow(pf$data), n)

  # FeaturePlot: expressing cells retained
  vals <- data.frame(cell = cells$cell[1:40], value = runif(40, 1, 5))   # 40 expressers
  pe <- view_feature_plot(cells, list(embedding = "umap", feature = "G"), vals,
                          list(raster = TRUE, order = TRUE))
  expect_lte(nrow(pe$data), 150)
  expect_true(all(vals$cell %in% pe$data$cell))                     # all expressers kept
})
