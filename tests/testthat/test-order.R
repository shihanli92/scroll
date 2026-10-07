# Group label order (R/order.R) and the Ridge panel (R/panel-ridge.R).

test_that("natural sort orders digit runs by value", {
  ns <- scroll:::.scroll_natural_sort
  expect_equal(ns(c("10", "2", "1", NA, "2")), c("1", "2", "10"))
  expect_equal(ns(c("c10", "c2", "B cells", "b2", "T cells")),
               c("B cells", "b2", "c2", "c10", "T cells"))
  expect_equal(ns(c("D | 10", "D | 9", "A | 1")), c("A | 1", "D | 9", "D | 10"))
  expect_equal(ns(character(0)), character(0))
})

test_that("a manual order goes first, the rest follow naturally", {
  lo <- scroll:::.scroll_level_order
  expect_equal(lo(c("1", "2", "10", "3"), c("10", "gone")), c("10", "1", "2", "3"))
  expect_equal(lo(c("b", "a"), NULL), c("a", "b"))
  expect_equal(scroll:::.scroll_order_value(""), character(0))
  expect_equal(scroll:::.scroll_order_value(list("b", "a")), c("b", "a"))
})

test_that("views honour the group order", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  v <- data$query1("RNA", "CD3D")
  lv <- sort(unique(as.character(data$cells$celltype)))
  want <- rev(lv)
  p <- scroll:::view_violin(data$cells, list(feature = "CD3D", group_by = "celltype"), v,
                            list(group_order = want))
  expect_equal(levels(ggplot2::ggplot_build(p)$plot$data$group), want)

  el <- data$queryN("RNA", c("CD3D", "NKG7"))
  a <- scroll:::.scroll_dotplot_assemble(data$cells, c("CD3D", "NKG7"), "celltype", el)
  p <- scroll:::view_dotplot(data$cells, list(), NULL, list(group_order = want), assembly = a)
  expect_equal(levels(p$data$group), want)
  # several group columns -> their "a | b" interaction
  a2 <- scroll:::.scroll_dotplot_assemble(data$cells, c("CD3D"), c("celltype", "condition"), el)
  expect_true(all(grepl(" | ", levels(a2$agg$group), fixed = TRUE)))

  p <- scroll:::view_proportions(data$cells, list(group_by = "celltype", fill_by = "condition"),
                                 list(x_order = "manual", group_order = want))
  expect_equal(levels(p$data$x), want)
})

test_that("view_ridge draws one ridge per group, first group on top", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  v <- data$query1("RNA", "CD3D")
  lv <- sort(unique(as.character(data$cells$celltype)))
  p <- scroll:::view_ridge(data$cells, list(feature = "CD3D", group_by = "celltype"), v)
  b <- ggplot2::ggplot_build(p)
  expect_s3_class(p$layers[[1]]$geom, "GeomRibbon")
  expect_equal(levels(p$data$group), lv)
  top <- tapply(p$data$base, as.character(p$data$group), unique)
  expect_equal(names(sort(top, decreasing = TRUE)), lv)          # first level highest
  expect_true(max(p$data$top - p$data$base) <= 1.4 + 1e-9)        # tallest = ridge_scale
  # a manual order, expressing cells only, a numeric column, several genes
  p <- scroll:::view_ridge(data$cells, list(feature = "CD3D", group_by = "celltype", nonzero = TRUE),
                           v, list(group_order = rev(lv)))
  expect_equal(levels(p$data$group), rev(lv))
  num <- scroll:::.scroll_num_cols(data$manifest)[[1]]
  expect_s3_class(scroll:::view_ridge(data$cells, list(value_col = num, group_by = "celltype")), "ggplot")
  expect_s3_class(scroll:::view_ridge_multi(data$cells, list(features = c("CD3D", "NKG7"),
                    group_by = "celltype"), data$queryN("RNA", c("CD3D", "NKG7"))), "patchwork")
})

test_that("the Ridge panel renders and its Style sheet saves the order", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  shiny::testServer(scroll:::ridge_server, args = list(data = data), {
    session$setInputs(group = c("celltype", "condition"), feature = "CD3D", ridgescale = 1.4,
                      order = "", onscreen = TRUE)
    expect_s3_class(plot_r(), "ggplot")
    expect_true(all(grepl(" | ", levels(plot_r()$data$group), fixed = TRUE)))
    session$setInputs(feature = c("CD3D", "NKG7"))
    expect_s3_class(plot_r(), "patchwork")
  })
  sec <- Filter(function(p) p$id == "ridge", scroll:::.scroll_builtin_panels())[[1]]
  expect_true("order" %in% scroll:::.scroll_style_plot_ids(sec, data))
  for (id in c("violin", "dotplot", "heatmap", "proportions")) {
    sec <- Filter(function(p) p$id == id, scroll:::.scroll_builtin_panels())[[1]]
    expect_true("order" %in% scroll:::.scroll_style_plot_ids(sec, data), info = id)
  }
})

test_that("DotPlot and Proportions accept several group columns", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  shiny::testServer(scroll:::dotplot_server, args = list(data = data), {
    session$setInputs(group = c("celltype", "condition"), markers = c("CD3D", "NKG7"),
                      display = "dots", scale = TRUE, cluster = "off", onscreen = TRUE)
    expect_true(all(grepl(" | ", levels(data_r()$assembly$agg$group), fixed = TRUE)))
    expect_s3_class(plot_r(), "ggplot")
  })
  shiny::testServer(scroll:::proportions_server, args = list(data = data), {
    session$setInputs(group = c("celltype", "condition"), fill = "condition",
                      position = "fill", onscreen = TRUE)
    expect_s3_class(plot_r(), "ggplot")
  })
})

test_that("ridge layouts: overlapping, separate rows, or overlaid; line type by group", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  v <- data$query1("RNA", "CD3D")
  rp <- function(...) scroll:::view_ridge(data$cells, list(feature = "CD3D", group_by = "celltype"),
                                          v, list(...))
  p <- rp(ridge_mode = "separate")
  expect_lt(max(p$data$top - p$data$base), 1)                # never reaches the next row
  p <- rp(ridge_mode = "overlay")
  expect_equal(unique(p$data$base), 0)                       # one shared baseline
  expect_true(ggplot2::is_ggplot(p))
  lv <- sort(unique(as.character(data$cells$celltype)))
  p <- rp(linetypes = stats::setNames(c("dashed", "dotted"), lv[1:2]))
  lt <- unique(ggplot2::ggplot_build(p)$data[[1]][, c("group", "linetype")])$linetype
  expect_setequal(lt, c("dashed", "dotted", "solid"))        # unlisted groups stay solid
})

test_that("the Ridge panel offers per-group line types only for few groups", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  shiny::testServer(scroll:::ridge_server, args = list(data = data), {
    session$setInputs(group = "celltype", feature = "CD3D", ridgescale = 1.4, order = "",
                      ridgemode = "overlay", ltgroup = TRUE, onscreen = TRUE)
    lv <- lvl_r()
    expect_equal(names(linetypes_r()), lv)
    expect_equal(unname(linetypes_r()), scroll:::.SCROLL_LINETYPES[seq_along(lv)])
    session$setInputs(lt_1 = "longdash")
    expect_equal(unname(linetypes_r()[1]), "longdash")
    expect_s3_class(plot_r(), "ggplot")
    session$setInputs(ltgroup = FALSE)
    expect_null(linetypes_r())
  })
})
