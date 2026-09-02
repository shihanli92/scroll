# DimPlot: colour by multiple columns (composite) + configurable highlight background.

test_that("view_umap_colorby colours by a composite of >1 column", {
  cells <- data.frame(cell = paste0("c", 1:6),
                      umap_1 = rnorm(6), umap_2 = rnorm(6),
                      a = rep(c("x", "y"), 3), b = rep(c("p", "q", "r"), 2),
                      stringsAsFactors = FALSE)
  p <- view_umap_colorby(cells, list(embedding = "umap", color_by = c("a", "b")), list())
  expect_s3_class(p, "ggplot")
  expect_true(any(grepl(" \\| ", p$data$.col)))          # "x | p" interaction levels
  expect_match(p$labels$colour %||% p$labels$color, "a \\| b")
  # single-column path is unchanged (back-compat)
  p1 <- view_umap_colorby(cells, list(embedding = "umap", color_by = "a"), list())
  expect_false(any(grepl(" \\| ", p1$data$.col)))
})

test_that("highlight background colour is configurable (default grey85)", {
  cells <- data.frame(cell = paste0("c", 1:6),
                      umap_1 = rnorm(6), umap_2 = rnorm(6),
                      grp = rep(c("x", "y", "z"), 2), stringsAsFactors = FALSE)
  layer_cols <- function(p) unlist(lapply(p$layers, function(L) L$aes_params$colour))
  pd <- view_umap_colorby(cells, list(embedding = "umap", color_by = "grp"),
                          list(highlight = "x"))
  expect_true("grey85" %in% layer_cols(pd))              # default background
  pc <- view_umap_colorby(cells, list(embedding = "umap", color_by = "grp"),
                          list(highlight = "x", bg_color = "#123456"))
  expect_true("#123456" %in% layer_cols(pc))             # custom background
  expect_false("grey85" %in% layer_cols(pc))
})

test_that("dimplot_server: composite color-by, composite highlight levels, bg colour", {
  skip_if_not_installed("SeuratObject")
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  shiny::testServer(scroll:::dimplot_server, args = list(data = data), {
    session$setInputs(reduction = "umap", palette = "Tableau 10", size = 0.6,
                      alpha = 0.85, labels = TRUE, legend = TRUE, split = "", raster = TRUE)
    session$setInputs(colorby = c("celltype", "condition"))   # a distinct Color-by event
    expect_true(is_cat())                                # composite is categorical
    expect_equal(data_r()$color_by, c("celltype", "condition"))
    lv <- levels_of()
    expect_true(any(grepl(" \\| ", lv)))                 # composite highlight choices
    session$setInputs(highlight = lv[[1]], bg_color = "#EEDDDD")
    expect_false(is.null(output$plot))
    # Manual palette: pickers restrict to the highlighted level(s), not all 82
    session$setInputs(palette = "Manual")
    expect_equal(manual_levels(), lv[[1]])
    session$setInputs(highlight = lv[1:2])
    expect_setequal(manual_levels(), lv[1:2])
    # no highlight -> every level gets a picker
    session$setInputs(highlight = character(0))
    expect_setequal(manual_levels(), levels_of())
    # clearing the box keeps the last valid selection (color_by never empty)
    session$setInputs(colorby = character(0))
    expect_equal(cb_rv(), c("celltype", "condition"))
  })
})
