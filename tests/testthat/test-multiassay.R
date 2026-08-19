# Multi-assay (multimodal) coverage. The shared test project is built with two
# assays (RNA + a synthetic ADT); these pin the per-assay build/query/dequant
# paths and the runtime assay-selector wiring that a single-assay fixture never
# exercises.

test_that("build emits one expr subtree per assay with distinct per-assay max", {
  dir <- test_project()
  expect_setequal(list.files(file.path(dir, "expr")), c("RNA", "ADT"))
  m <- scroll_manifest(dir)
  expect_setequal(names(m$assays), c("RNA", "ADT"))
  expect_equal(m$assays$ADT$n_features, 8)
  # the two assays normalize to different scales -> different dequant max
  expect_false(isTRUE(all.equal(m$assays$RNA$max, m$assays$ADT$max)))
})

test_that("queries and dequantization are per-assay", {
  dir <- test_project()
  con <- scroll_connect(dir); on.exit(scroll_disconnect(con))
  m <- scroll_manifest(dir)
  adt <- scroll_query_feature(con, "ADT", "CD3-P")
  expect_gt(nrow(adt), 0)                                   # ADT feature reachable
  dq <- scroll_dequantize(adt$value, m, "ADT")
  expect_lte(max(dq), m$assays$ADT$max + 1e-6)             # uses ADT max, not RNA's
  # same-named feature does not leak across assays: RNA has no "CD3-P"
  expect_equal(nrow(scroll_query_feature(con, "RNA", "CD3-P")), 0)
})

test_that(".scroll_assay_dir rejects unsafe assay names (defense-in-depth)", {
  con <- scroll_connect(test_project()); on.exit(scroll_disconnect(con))
  expect_error(scroll_query_feature(con, "../secret", "x"), "Invalid assay")
  expect_error(scroll_query_feature(con, "a/b", "x"), "Invalid assay")
  expect_error(scroll_query_feature(con, "", "x"), "Invalid assay")
  expect_s3_class(scroll_query_feature(con, "RNA", "CD3D"), "data.frame")   # valid still works
})

test_that("expression panels show an assay selector only with >=2 assays", {
  data <- scroll:::.scroll_load(test_project())               # 2 assays
  on.exit(scroll_disconnect(data$con))
  for (ui in list(scroll:::featureplot_ui, scroll:::dotplot_ui,
                  scroll:::violin_ui, scroll:::de_ui)) {
    expect_match(as.character(ui("x", data)), "x-assay")      # selector present
  }
  # a single-assay handle -> no assay selector
  fake1 <- list(config = list(), manifest = list(
    meta = list(celltype = list(type = "categorical", levels = list("A", "B"))),
    assays = list(RNA = list(features = list("A", "B"), max = 1, n_features = 2)),
    embeddings = list(umap = list(dims = 2)),
    default_assay = "RNA", default_embedding = "umap"))
  expect_false(grepl("x-assay", as.character(scroll:::featureplot_ui("x", fake1))))
})

test_that("featureplot_server renders a protein feature from the ADT assay", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  shiny::testServer(scroll:::featureplot_server, args = list(data = data), {
    session$setInputs(reduction = "umap", assay = "ADT", feature = "CD3-P",
                      palette = "grey-purple", size = 0.7, clip = c(0, 100),
                      order = TRUE, legend = TRUE, split = "", aspect = 1)
    expect_equal(assay(), "ADT")                             # assay reactive follows input
    expect_false(is.null(output$plot))
    expect_false(is.null(plot_r()))
  })
})

test_that("featureplot_server blends two genes into a co-expression layout", {
  skip_if_not_installed("patchwork")
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  shiny::testServer(scroll:::featureplot_server, args = list(data = data), {
    session$setInputs(reduction = "umap", assay = "RNA", feature = c("CD3D", "MS4A1"),
                      blend = TRUE, blend_threshold = 0.5, blend_c1 = "#FF0000",
                      blend_c2 = "#00FF00", palette = "grey-purple", size = 0.7,
                      clip = c(0, 100), order = TRUE, legend = TRUE, split = "", aspect = 1)
    d <- data_r()
    expect_true(isTRUE(d$blend))
    expect_false(is.null(d$values2))                       # second gene queried
    expect_s3_class(plot_r(), "patchwork")
    expect_true(all(c("CD3D", "MS4A1") %in% names(csv_r())))  # both genes exported
  })

  # multiple genes (no blend) -> a per-gene grid
  shiny::testServer(scroll:::featureplot_server, args = list(data = data), {
    session$setInputs(reduction = "umap", assay = "RNA", feature = c("CD3D", "MS4A1", "CD8A"),
                      blend = FALSE, metacol = "", palette = "grey-purple", size = 0.7,
                      clip = c(0, 100), order = TRUE, legend = TRUE, split = "", aspect = 1)
    d <- data_r()
    expect_true(isTRUE(d$multi))
    expect_length(d$features, 3)
    expect_s3_class(plot_r(), "patchwork")
    expect_true(all(c("CD3D", "MS4A1", "CD8A") %in% names(csv_r())))
  })

  # genes AND numeric metadata columns shown together in one grid
  shiny::testServer(scroll:::featureplot_server, args = list(data = data), {
    session$setInputs(reduction = "umap", assay = "RNA", feature = c("CD3D", "MS4A1"),
                      metacol = "nCount_RNA", blend = FALSE, palette = "grey-purple",
                      size = 0.7, clip = c(0, 100), order = TRUE, legend = TRUE,
                      split = "", aspect = 1)
    d <- data_r()
    expect_setequal(d$features, c("CD3D", "MS4A1", "nCount_RNA"))
    expect_setequal(unique(d$values$feature), c("CD3D", "MS4A1", "nCount_RNA"))
    expect_s3_class(plot_r(), "patchwork")
    expect_true(all(c("CD3D", "MS4A1", "nCount_RNA") %in% names(csv_r())))
  })
})

test_that("scroll_de runs a two-group contrast on the ADT assay", {
  skip_if_not_installed("presto")
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  res <- scroll_de(data, "ADT", "celltype", ident1 = "T", ident2 = "B")   # two-group
  expect_gt(nrow(res), 0)
  expect_true(all(res$gene %in% unlist(data$manifest$assays$ADT$features)))
  expect_false(is.unsorted(res$p_val_adj))
})
