# Per-plot export: PNG (raster) + PDF (vector) for every panel, CSV for the DE
# table. These invoke the download handlers' content functions and validate the
# actual bytes written, incl. the DotPlot aplot composite (which ggsave rejects).

read_cells <- function(dir) as.data.frame(arrow::read_parquet(file.path(dir, "cells.parquet")))

# Reach a downloadHandler's content(file) function through its closure.
dl_content <- function(h) environment(environment(h)$renderFunc)$content

test_that(".scroll_img_handler writes real PNG and PDF files", {
  cells <- read_cells(test_project())
  p <- view_umap_colorby(cells, list(embedding = "umap", color_by = "celltype"))
  for (fmt in c("png", "pdf")) {
    h <- scroll:::.scroll_img_handler(function() p, paste0("x.", fmt), fmt)
    expect_s3_class(h, "shiny.render.function")
    f <- tempfile(fileext = paste0(".", fmt))
    dl_content(h)(f)
    expect_gt(file.info(f)$size, 1000)
    magic <- readBin(f, "raw", 4)                          # file-type signature
    if (fmt == "png") expect_equal(magic[2:4], as.raw(c(0x50, 0x4e, 0x47)))  # \x89PNG
    else expect_equal(rawToChar(magic), "%PDF")
  }
})

test_that("PNG and PDF export handle the DotPlot aplot composite", {
  cells <- read_cells(test_project())
  feats <- c("CD3D", "CD8A", "MS4A1", "NKG7", "GNLY")
  el <- data.frame(feature = rep(feats, each = 4), cell = rep(cells$cell[1:4], 5),
                   value = runif(20))
  p <- view_dotplot(cells, list(group_by = "celltype", features = feats), el,
                    state = list(cluster = "both", palette = "RdBu", aspect = 1.2, scale = TRUE))
  skip_if_not(inherits(p, "aplot"), "clustering did not produce an aplot")
  for (fmt in c("png", "pdf")) {
    h <- scroll:::.scroll_img_handler(function() p, paste0("dp.", fmt), fmt)
    f <- tempfile(fileext = paste0(".", fmt))
    dl_content(h)(f)
    expect_gt(file.info(f)$size, 1000)                     # aplot survives print()-to-device
  }
})

test_that(".scroll_csv_handler writes the DE table to CSV", {
  de <- data.frame(gene = c("A", "B"), logFC = c(1, -2), auc = c(.9, .1),
                   pct.1 = c(.8, .2), pct.2 = c(.1, .7),
                   p_val = c(1e-9, 1e-3), p_val_adj = c(1e-7, 1e-2))
  h <- scroll:::.scroll_csv_handler(function() de, "scroll_de.csv")
  expect_s3_class(h, "shiny.render.function")
  f <- tempfile(fileext = ".csv")
  dl_content(h)(f)
  back <- utils::read.csv(f)
  expect_equal(nrow(back), 2)
  expect_true(all(c("gene", "logFC", "p_val_adj") %in% names(back)))
})

test_that("panel servers register PNG and PDF download outputs", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  shiny::testServer(scroll:::dimplot_server, args = list(data = data), {
    session$setInputs(reduction = "umap", colorby = "celltype", palette = "Tableau 10",
                      size = 0.6, alpha = 0.85, labels = TRUE, legend = TRUE, split = "",
                      aspect = 1, highlight = character(0))
    expect_false(is.null(output$png))
    expect_false(is.null(output$pdf))
  })
})

test_that("de_server registers PNG, PDF and CSV downloads", {
  skip_if_not_installed("presto")
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  ct <- sort(unique(as.character(data$cells$celltype)))[[1]]
  shiny::testServer(scroll:::de_server, args = list(data = data), {
    session$setInputs(group = "celltype", ident1 = ct, ident2 = "rest",
                      minpct = 10, topn = 20, lfc = 1, padj = 0.05, labeln = 10,
                      aspect = 1, compute = 1)
    expect_false(is.null(output$png))
    expect_false(is.null(output$pdf))
    expect_false(is.null(output$csv))
  })
})
