# Pseudobulk DE: counts store, duckdb aggregation, edgeR/limma-voom, and the panel.
# The shared test_project() is built with counts = TRUE.

test_that("counts store is written and recorded when counts = TRUE", {
  dir <- test_project()
  expect_true(file.exists(file.path(dir, "counts", "RNA.parquet")))
  expect_true(isTRUE(scroll_manifest(dir)$has_counts))
})

test_that("scroll_aggregate_counts sums counts per pseudobulk sample", {
  dir <- test_project()
  con <- scroll_connect(dir); on.exit(scroll_disconnect(con))
  cells <- as.data.frame(arrow::read_parquet(file.path(dir, "cells.parquet")))$cell
  mapping <- data.frame(cell = cells[1:20], psample = rep(c("A", "B"), each = 10))
  agg <- scroll_aggregate_counts(con, "RNA", mapping)
  expect_setequal(names(agg), c("feature", "psample", "count"))
  expect_setequal(unique(agg$psample), c("A", "B"))
  # cross-check one gene's summed count against the raw store
  raw <- as.data.frame(arrow::read_parquet(file.path(dir, "counts", "RNA.parquet")))
  truth <- sum(raw$value[raw$feature == "CD3D" & raw$cell %in% cells[1:10]])
  got <- agg$count[agg$feature == "CD3D" & agg$psample == "A"]
  expect_equal(if (length(got)) got else 0, truth)
})

test_that("scroll_aggregate_counts errors without a counts store", {
  # a project built without counts
  dir <- file.path(tempdir(), "scroll-nocounts")
  if (!dir.exists(file.path(dir, "expr")))
    suppressMessages(scroll_build(make_test_object(), dir, assays = "RNA", overwrite = TRUE))
  con <- scroll_connect(dir); on.exit(scroll_disconnect(con))
  expect_error(scroll_aggregate_counts(con, "RNA", data.frame(cell = "x", psample = "A")),
               "counts store")
})

test_that("combined-interaction levels combine columns (and drop NAs)", {
  cells <- data.frame(g = c("A", "A", "B", NA), c = c("1", "2", "1", "1"),
                      stringsAsFactors = FALSE)
  expect_equal(scroll:::.scroll_combo_levels(cells, c("g", "c")),
               c("A | 1", "A | 2", "B | 1", NA))
  expect_setequal(scroll:::.scroll_combo_choices(cells, c("g", "c")),
                  c("A | 1", "A | 2", "B | 1"))
})

test_that("scroll_pseudobulk_de runs the real-replicate path", {
  skip_if_not_installed("edgeR"); skip_if_not_installed("limma")
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  res <- scroll_pseudobulk_de(data, "RNA", aggregate_cols = "celltype",
                              ident1 = "T", ident2 = "B",
                              replicate_col = "condition", min_cells = 5)
  expect_setequal(names(res), c("gene", "logFC", "avg_expr", "p_val", "p_val_adj"))
  expect_gt(nrow(res), 0)
  expect_false(is.unsorted(res$p_val_adj))
  expect_false(attr(res, "pseudo"))                 # real replicates used
})

test_that("scroll_pseudobulk_de falls back to pseudo-replicates", {
  skip_if_not_installed("edgeR"); skip_if_not_installed("limma")
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  set.seed(1)
  res <- scroll_pseudobulk_de(data, "RNA", aggregate_cols = "celltype",
                              ident1 = "T", ident2 = "B", replicate_col = NULL,
                              n_pseudo = 3, cells_per_pseudo = 12, min_cells = 5)
  expect_gt(nrow(res), 0)
  expect_true(attr(res, "pseudo"))                  # pseudo-reps used
  expect_equal(unname(attr(res, "n_samples")), c(3L, 3L))
})

test_that("no_replicate pseudo-replicates the two groups (pooled across combos)", {
  skip_if_not_installed("edgeR"); skip_if_not_installed("limma")
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  set.seed(2)
  combos <- scroll:::.scroll_combo_choices(data$cells, c("celltype", "condition"))
  g1 <- grep("^T ", combos, value = TRUE)   # T | ctrl, T | treat  (2 combos)
  g2 <- grep("^B ", combos, value = TRUE)
  res <- scroll_pseudobulk_de(data, "RNA", aggregate_cols = c("celltype", "condition"),
                              ident1 = g1, ident2 = g2, replicate_col = "no_replicate",
                              n_pseudo = 4, cells_per_pseudo = 20, min_cells = 5)
  expect_true(attr(res, "pseudo"))
  # group-based: 4 pseudo-reps per group regardless of the 2 combos each pools
  expect_equal(unname(attr(res, "n_samples")), c(4L, 4L))
})

test_that("scroll_pseudobulk_de errors clearly on too-few samples / no counts", {
  skip_if_not_installed("edgeR"); skip_if_not_installed("limma")
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  # min_cells too high -> not enough samples per group
  expect_error(
    scroll_pseudobulk_de(data, "RNA", aggregate_cols = "celltype",
                         ident1 = "T", ident2 = "B", replicate_col = "condition",
                         min_cells = 1000),
    "pseudobulk samples")
})

test_that(".scroll_stability_aggregate scores selection frequency + sign agreement", {
  # gene A: strong, always significant; gene B: null, never significant
  runs <- lapply(1:4, function(i) data.frame(
    gene = c("A", "B"), logFC = c(2 + i * 0.1, 0.05 * (-1)^i),
    avg_expr = 1, p_val = c(1e-5, 0.5), p_val_adj = c(1e-4, 0.6)))
  agg <- scroll:::.scroll_stability_aggregate(runs, lfc = 1, padj = 0.05)
  expect_equal(attr(agg, "runs"), 4)
  expect_equal(agg$sel_freq[agg$gene == "A"], 1)      # A passes every run
  expect_equal(agg$sel_freq[agg$gene == "B"], 0)      # B never passes
  expect_equal(agg$sign_agree[agg$gene == "A"], 1)    # A always positive
  expect_true(agg$median_logFC[agg$gene == "A"] > 1)
})

test_that("scroll_pseudobulk_stability returns per-gene stability over runs", {
  skip_if_not_installed("edgeR"); skip_if_not_installed("limma")
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  set.seed(1)
  st <- scroll_pseudobulk_stability(data, "RNA", aggregate_cols = "celltype",
                                    ident1 = "T", ident2 = "B", replicate_col = "no_replicate",
                                    runs = 6, n_pseudo = 3, cells_per_pseudo = 15, min_cells = 5)
  expect_true(all(c("gene", "sel_freq", "median_logFC", "sign_agree", "n_tested") %in% names(st)))
  expect_equal(attr(st, "runs"), 6)
  expect_true(all(st$sel_freq >= 0 & st$sel_freq <= 1))
})

test_that("pseudobulk panel gates off without a counts store", {
  fake <- list(manifest = list(
    has_counts = FALSE,
    meta = list(celltype = list(type = "categorical", levels = list("A", "B"))),
    assays = list(RNA = list(features = list("A", "B"), max = 1, n_features = 2)),
    embeddings = list(umap = list(dims = 2)),
    default_assay = "RNA", default_embedding = "umap"))
  tag <- scroll:::pseudobulk_de_ui("pb", fake)
  expect_match(as.character(tag), "counts store")
})

test_that("pseudobulk_de_server populates combined levels and computes", {
  skip_if_not_installed("edgeR"); skip_if_not_installed("limma")
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  shiny::testServer(scroll:::pseudobulk_de_server, args = list(data = data), {
    session$setInputs(aggregate_by = "celltype", ident1 = "T", ident2 = "B",
                      replicate = "condition", mincells = 5, npseudo = 3, cellsper = 12,
                      runs = 1, stabcut = 0.8,
                      topn = 20, lfc = 1, padj = 0.05, labeln = 10, aspect = 1, compute = 1)
    expect_false(is.null(result()$ok))
    expect_false(is.null(output$table))
    expect_false(is.null(output$plot))
  })
})

test_that("pseudobulk_de_server runs stability mode (runs > 1, no_replicate)", {
  skip_if_not_installed("edgeR"); skip_if_not_installed("limma")
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  shiny::testServer(scroll:::pseudobulk_de_server, args = list(data = data), {
    session$setInputs(aggregate_by = "celltype", ident1 = "T", ident2 = "B",
                      replicate = "no_replicate", mincells = 5, npseudo = 3, cellsper = 15,
                      runs = 5, stabcut = 0.8,
                      topn = 20, lfc = 1, padj = 0.05, labeln = 10, aspect = 1, compute = 1)
    expect_length(result()$runs, 5)               # 5 per-run results stored
    expect_true("sel_freq" %in% names(de_df()))   # aggregated to stability stats
    expect_false(is.null(output$plot))
  })
})
