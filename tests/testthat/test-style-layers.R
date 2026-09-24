# Per-layer geom settings (R/style-layers.R): the catalog of a drawn plot's layers,
# applying overrides to clones of those layers, and the Style sheet's Layers section.

layers_of <- function(cat) vapply(cat, `[[`, "", "key")

test_that("the catalog lists each layer with its unmapped constants, keyed by role", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  p <- suppressWarnings(view_umap_colorby(data$cells, list(embedding = "umap", color_by = "celltype"),
                                          list(highlight = "T", raster = FALSE, show_labels = TRUE)))
  cat <- scroll:::.scroll_layer_catalog(p)
  expect_identical(layers_of(cat), c("Background cells", "Cells", "Cluster labels"))
  cells <- cat[[2]]$controls
  expect_false("colour" %in% names(cells))                 # mapped in aes(): never offered
  expect_true(all(c("size", "alpha", "shape", "stroke") %in% names(cells)))
  expect_identical(cells$size$current, "0.6")              # the hard-coded value as a hint
  expect_true("colour" %in% names(cat[[1]]$controls))       # background colour is a constant
  expect_identical(cat[[3]]$controls$size$current, "4")
  expect_identical(scroll:::.scroll_layer_catalog(p), cat) # stable (identical) across calls
})

test_that("overrides land in the drawn data; mapped aesthetics and the original are untouched", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  p <- suppressWarnings(view_umap_colorby(data$cells, list(embedding = "umap", color_by = "celltype"),
                                          list(highlight = "T", raster = FALSE, show_labels = TRUE)))
  q <- scroll:::.scroll_apply_style(p, list(layers = list(
    Cells = list(size = "2.5", shape = "17", colour = "#FF0000"),   # colour is mapped -> ignored
    `Background cells` = list(colour = "#00FF00"),
    `Cluster labels` = list(size = "6", fontface = "italic"))))
  b <- ggplot2::ggplot_build(q)
  expect_identical(unique(b$data[[1]]$colour), "#00FF00")
  expect_identical(unique(b$data[[2]]$size), 2.5)
  expect_identical(unique(b$data[[2]]$shape), 17L)
  expect_false("#FF0000" %in% b$data[[2]]$colour)           # the mapping still colours cells
  expect_identical(unique(b$data[[3]]$fontface), "italic")
  expect_identical(p$layers[[2]]$aes_params$size, 0.6)      # original plot unchanged
  # junk is dropped, not an error
  j <- scroll:::.scroll_apply_style(p, list(layers = list(Cells = list(size = "abc", alpha = "7",
                                                                       shape = "999"))))
  expect_no_error(ggplot2::ggplot_build(j))
  expect_identical(j$layers[[2]]$aes_params$size, 0.6)
})

test_that("geom / stat / position parameters apply (violin trim + width, jitter)", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  v <- suppressWarnings(view_violin(data$cells, list(feature = "x", group_by = "celltype",
                                                     value_col = "nCount_RNA"), NULL, list(jitter = TRUE)))
  expect_identical(layers_of(scroll:::.scroll_layer_catalog(v)), c("Violins", "Points"))
  q <- scroll:::.scroll_apply_style(v, list(layers = list(
    Violins = list(trim = "FALSE", width = "0.5", scale = "area"),
    Points = list(jitter_width = "0.05", size = "1"))))
  expect_false(q$layers[[1]]$stat_params$trim)
  expect_identical(q$layers[[1]]$stat_params$scale, "area")
  expect_identical(q$layers[[1]]$aes_params$width, 0.5)
  expect_identical(q$layers[[2]]$position$width, 0.05)
  expect_no_error(ggplot2::ggplot_build(q))
  expect_true(isTRUE(v$layers[[1]]$stat_params$trim))       # original unchanged
})

test_that("raster points take a pixel size; patchwork sub-plots are all styled", {
  skip_if_not_installed("scattermore")
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  s <- suppressWarnings(view_umap_colorby(data$cells, list(embedding = "umap", color_by = "celltype"),
                                          list(raster = TRUE)))
  q <- scroll:::.scroll_apply_style(s, list(layers = list(Cells = list(size = "3"))))
  expect_identical(q$layers[[1]]$geom_params$pointsize, 6)
  skip_if_not_installed("patchwork")
  mk <- function() ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point()
  pw <- patchwork::wrap_plots(mk(), mk())
  r <- scroll:::.scroll_apply_style(pw, list(layers = list(Points = list(size = "4"))))
  expect_identical(r$layers[[1]]$aes_params$size, 4)
  expect_identical(r$patches$plots[[1]]$layers[[1]]$aes_params$size, 4)
})

test_that("the style plot helper reports the catalog to the sheet, only when it changes", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  rv <- shiny::reactiveVal(list()); cat_rv <- shiny::reactiveVal(NULL)
  attr(rv, "scroll_catalog") <- cat_rv
  p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point(size = 2)
  shiny::isolate(scroll:::.scroll_style_plot(p, rv))
  expect_identical(shiny::isolate(layers_of(cat_rv())), "Points")
})

test_that("Layers inputs feed the style; undrawn layers keep their overrides; reset clears", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point(size = 2) +
    ggplot2::geom_hline(yintercept = 20)
  rv <- shiny::reactiveVal(list(layers = list(Gone = list(size = "9"))))
  cat_rv <- shiny::reactiveVal(scroll:::.scroll_layer_catalog(p))
  attr(rv, "scroll_catalog") <- cat_rv
  shiny::testServer(scroll:::.scroll_style_server,
                    args = list(data = data, rv = rv, all = list(rv)), {
    cat <- cat_rv(); vals <- list()
    for (e in cat) for (k in names(e$controls)) {
      id <- scroll:::.scroll_layer_input_id(e$key, k)
      if (identical(e$controls[[k]]$type, "colour")) vals[[paste0(id, "_mode")]] <- "" else vals[[id]] <- ""
    }
    vals[[scroll:::.scroll_layer_input_id("Points", "size")]] <- "3"
    vals[[scroll:::.scroll_layer_input_id("Reference lines", "linetype")]] <- "dashed"
    do.call(session$setInputs, vals)
    session$elapse(600)
    expect_identical(rv()$layers$Points, list(size = "3"))
    expect_identical(rv()$layers$`Reference lines`, list(linetype = "dashed"))
    expect_identical(rv()$layers$Gone, list(size = "9"))       # not drawn now: kept
    session$setInputs(style_reset = 1)
    expect_length(rv(), 0)
  })
})

test_that("duplicated look sliders are gone from the sheets; a Layers section is there", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  specs <- scroll:::.scroll_builtin_panels()
  get <- function(id) Filter(function(s) identical(s$id, id), specs)[[1]]
  gone <- list(dimplot = c("size", "alpha", "bg_color"), featureplot = "size",
               violin = c("vwidth", "psize", "palpha"),
               proportions = c("barwidth", "outline", "labelsize"))
  for (p in names(gone)) {
    html <- as.character(get(p)$style_ui(p, data))
    for (i in gone[[p]]) expect_false(grepl(sprintf('id="%s-%s"', p, i), html, fixed = TRUE),
                                      info = paste(p, i))
  }
  sheet <- as.character(scroll:::.scroll_style_sheet(get("dimplot"), data))
  expect_match(sheet, 'id="dimplot-layers_ui"', fixed = TRUE)
})
