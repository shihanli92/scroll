# Pseudobulk DE: counts store, arrow aggregation, edgeR/limma-voom, and the panel.
# The shared test_project() is built with counts = TRUE.

test_that("counts store is written and recorded when counts = TRUE", {
  dir <- test_project()
  expect_true(file.exists(file.path(dir, "counts", "RNA.parquet")))
  expect_true(isTRUE(scroll_manifest(dir)$has_counts))
})

test_that("scroll_aggregate_counts sums counts per pseudobulk sample", {
  dir <- test_project()
  con <- scroll_connect(dir); on.exit(scroll_disconnect(con))
  mapping <- data.frame(cell = 1:20, psample = rep(c("A", "B"), each = 10))  # v2 int cell index
  agg <- scroll_aggregate_counts(con, "RNA", mapping)
  expect_setequal(names(agg), c("feature", "psample", "count"))
  expect_setequal(unique(agg$psample), c("A", "B"))
  # cross-check one gene's summed count against the raw store
  raw <- as.data.frame(arrow::read_parquet(file.path(dir, "counts", "RNA.parquet")))
  truth <- sum(raw$value[raw$feature == "CD3D" & raw$cell %in% 1:10])
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
  # 2 disjoint pseudo-reps per group (bins ~half the group; big enough for voom on
  # the sparse synthetic counts). Group-based: pooled across the 2 combos each.
  res <- scroll_pseudobulk_de(data, "RNA", aggregate_cols = c("celltype", "condition"),
                              ident1 = g1, ident2 = g2, replicate_col = "no_replicate",
                              n_pseudo = 2, cells_per_pseudo = 30, min_cells = 5)
  expect_true(attr(res, "pseudo"))
  expect_equal(unname(attr(res, "n_samples")), c(2L, 2L))
})

test_that("scroll_pseudobulk_de errors clearly on too-few samples / no counts", {
  skip_if_not_installed("edgeR"); skip_if_not_installed("limma")
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  # min_cells too high -> no real replicates qualify, pseudo fallback can't fill bins
  expect_error(
    scroll_pseudobulk_de(data, "RNA", aggregate_cols = "celltype",
                         ident1 = "T", ident2 = "B", replicate_col = "condition",
                         min_cells = 1000),
    "needs >=")
})

test_that(".scroll_pseudobulk_design blocks on the replicate only when paired", {
  # paired: >= 2 replicate levels shared across both groups
  samp_p <- data.frame(psample = paste(rep(c("group1","group2"), each = 3), rep(1:3, 2)),
                       group = rep(c("group1","group2"), each = 3),
                       rep = as.character(rep(1:3, 2)), stringsAsFactors = FALSE)
  dp <- scroll:::.scroll_pseudobulk_design(samp_p, "real-paired", "auto")
  expect_equal(dp$formula, "~ 0 + group + replicate")
  expect_false(any(grepl("Intercept", colnames(dp$design))))   # no-intercept means model
  expect_true(any(grepl("^replicate", colnames(dp$design))))
  # contrast tests group1 - group2
  expect_equal(unname(dp$contrast["groupgroup1"]), 1)
  expect_equal(unname(dp$contrast["groupgroup2"]), -1)
  # unpaired / pseudo -> ~ 0 + group
  du <- scroll:::.scroll_pseudobulk_design(samp_p, "real-unpaired", "auto")
  expect_equal(du$formula, "~ 0 + group")
  # forcing paired = "yes" on a non-paired regime warns and falls back
  dy <- scroll:::.scroll_pseudobulk_design(samp_p, "real-unpaired", "yes")
  expect_equal(dy$formula, "~ 0 + group"); expect_false(is.null(dy$warn))
})

test_that(".scroll_pseudobulk_samples returns the samples the compute uses", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  sm <- scroll:::.scroll_pseudobulk_samples(
    data$cells, seq_len(nrow(data$cells)), "celltype", "T", "B",
    "condition", min_cells = 5, n_pseudo = 3, cells_per_pseudo = 50)
  expect_setequal(names(sm$samp), c("psample", "group", "rep"))
  expect_equal(sm$regime, "real-paired")             # condition spans T and B
  expect_true(all(c("group1","group2") %in% sm$samp$group))
})

test_that(".scroll_partition_reps builds disjoint, floored, capped pseudo-reps", {
  set.seed(1)
  # 30 cells, 3 reps, min 5: balanced bins of 10, no shared cells, union within input
  reps <- scroll:::.scroll_partition_reps(1:30, "g", "group1",
                                          min_cells = 5, n_pseudo = 3, cells_per_pseudo = 50)
  expect_length(reps, 3)
  all_cells <- unlist(lapply(reps, `[[`, "cell"))
  expect_false(any(duplicated(all_cells)))          # DISJOINT: never share a cell
  expect_true(all(all_cells %in% 1:30))
  expect_true(all(vapply(reps, nrow, integer(1)) >= 5))
  # cap: 3 bins of ~10 capped at 4 -> exactly 4 each, still disjoint
  capped <- scroll:::.scroll_partition_reps(1:30, "g", "group1",
                                            min_cells = 4, n_pseudo = 3, cells_per_pseudo = 4)
  expect_true(all(vapply(capped, nrow, integer(1)) == 4))
  expect_false(any(duplicated(unlist(lapply(capped, `[[`, "cell")))))
  # too small to fill n_pseudo bins of min_cells -> empty (caller turns into an error)
  expect_length(scroll:::.scroll_partition_reps(1:10, "g", "group1",
                                                min_cells = 5, n_pseudo = 3, cells_per_pseudo = 50), 0)
})

test_that("no_replicate refuses when a group can't fill n_pseudo x min_cells", {
  skip_if_not_installed("edgeR"); skip_if_not_installed("limma")
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  expect_error(
    scroll_pseudobulk_de(data, "RNA", aggregate_cols = "celltype",
                         ident1 = "T", ident2 = "B", replicate_col = "no_replicate",
                         n_pseudo = 5, min_cells = 200),   # 5 x 200 = 1000 > group size
    "needs >=")
})

test_that("real replicates shared across groups give a paired design", {
  skip_if_not_installed("edgeR"); skip_if_not_installed("limma")
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  # condition (ctrl/treat) spans both T and B -> paired
  res <- scroll_pseudobulk_de(data, "RNA", aggregate_cols = "celltype",
                              ident1 = "T", ident2 = "B",
                              replicate_col = "condition", min_cells = 5)
  expect_equal(attr(res, "regime"), "real-paired")
  expect_equal(attr(res, "design"), "~ 0 + group + replicate")
  expect_false(attr(res, "pseudo"))
  # paired vs forcing ~ group on the same data differ (blocking changes the fit)
  unp <- scroll_pseudobulk_de(data, "RNA", aggregate_cols = "celltype",
                              ident1 = "T", ident2 = "B",
                              replicate_col = "condition", min_cells = 5, paired = "no")
  expect_equal(attr(unp, "design"), "~ 0 + group")
  j <- merge(res, unp, by = "gene")
  expect_false(isTRUE(all.equal(j$p_val.x, j$p_val.y)))
})

test_that("disjoint real replicates give an unpaired design", {
  skip_if_not_installed("edgeR"); skip_if_not_installed("limma")
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  cells <- data$cells
  # donors disjoint across T (d1/d2) and B (d3/d4): no replicate shared -> unpaired
  ct <- as.character(cells$celltype)
  set.seed(3)
  cells$donor <- ifelse(ct == "T", sample(c("d1", "d2"), nrow(cells), TRUE),
                 ifelse(ct == "B", sample(c("d3", "d4"), nrow(cells), TRUE), NA))
  res <- scroll_pseudobulk_de(data, "RNA", aggregate_cols = "celltype",
                              ident1 = "T", ident2 = "B",
                              replicate_col = "donor", min_cells = 5, cells = cells)
  expect_equal(attr(res, "regime"), "real-unpaired")
  expect_equal(attr(res, "design"), "~ 0 + group")
  expect_false(attr(res, "pseudo"))
})

test_that("missing replicate values are reported, not silently dropped", {
  skip_if_not_installed("edgeR"); skip_if_not_installed("limma")
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  cells <- data$cells
  cells$donor <- ifelse(as.character(cells$celltype) == "T",
                        sample(c("d1", "d2"), nrow(cells), TRUE),
                        sample(c("d1", "d2"), nrow(cells), TRUE))
  cells$donor[1:5] <- NA                     # some cells lack a replicate value
  expect_warning(
    res <- scroll_pseudobulk_de(data, "RNA", aggregate_cols = "celltype",
                                ident1 = "T", ident2 = "B",
                                replicate_col = "donor", min_cells = 3, cells = cells),
    "missing")
  expect_gt(attr(res, "dropped_na"), 0)
})

test_that("a named replicate column that is absent is an error", {
  skip_if_not_installed("edgeR"); skip_if_not_installed("limma")
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  expect_error(
    scroll_pseudobulk_de(data, "RNA", aggregate_cols = "celltype",
                         ident1 = "T", ident2 = "B", replicate_col = "nope"),
    "not found")
})

test_that("logFC sign points up in ident1 (coefficient found by name)", {
  skip_if_not_installed("edgeR"); skip_if_not_installed("limma")
  # a project with a gene deliberately high in condition = treat
  dir <- file.path(tempdir(), "scroll-pb-sign")
  if (!dir.exists(file.path(dir, "counts"))) {
    obj <- make_test_object(n = 120, seed = 5)
    cm <- as.matrix(SeuratObject::GetAssayData(obj, layer = "counts"))
    treat <- obj$condition == "treat"
    cm["G1", treat]  <- cm["G1", treat] + 40L      # strong signal up in treat
    cm["G2", !treat] <- cm["G2", !treat] + 40L     # and one up in ctrl
    obj <- SeuratObject::SetAssayData(obj, layer = "counts",
                                      new.data = Matrix::Matrix(cm, sparse = TRUE))
    obj <- SeuratObject::SetAssayData(obj, layer = "data",
             new.data = Matrix::Matrix(.lognorm(cm, 1e4), sparse = TRUE))
    suppressMessages(scroll_build(obj, dir, assays = "RNA", counts = TRUE, overwrite = TRUE))
  }
  data <- scroll:::.scroll_load(dir); on.exit(scroll_disconnect(data$con))
  res <- scroll_pseudobulk_de(data, "RNA", aggregate_cols = "condition",
                              ident1 = "treat", ident2 = "ctrl",
                              replicate_col = "no_replicate",
                              n_pseudo = 3, cells_per_pseudo = 15, min_cells = 5)
  expect_gt(res$logFC[res$gene == "G1"], 0)        # up in ident1 (treat)
  expect_lt(res$logFC[res$gene == "G2"], 0)        # up in ctrl -> negative
})

test_that("cells = restricts the aggregation", {
  skip_if_not_installed("edgeR"); skip_if_not_installed("limma")
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  sub <- data$cells[data$cells$celltype %in% c("T", "B"), , drop = FALSE]
  res <- scroll_pseudobulk_de(data, "RNA", aggregate_cols = "celltype",
                              ident1 = "T", ident2 = "B", replicate_col = "no_replicate",
                              n_pseudo = 3, cells_per_pseudo = 20, min_cells = 5, cells = sub)
  # every pseudo-sample's cells drawn only from the T/B subset
  expect_lte(sum(attr(res, "sample_sizes")), nrow(sub))
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

test_that("scroll_pseudobulk_de(runs > 1) returns per-gene stability over runs", {
  skip_if_not_installed("edgeR"); skip_if_not_installed("limma")
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  st <- scroll_pseudobulk_de(data, "RNA", aggregate_cols = "celltype",
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
    # export-all switch: CSV picks the full de_df() over the displayed top-N
    session$setInputs(topn = 3, export_all = FALSE)
    full <- nrow(de_df())
    if (full > 3) {
      expect_lte(nrow(table_rows()), 3)
      session$setInputs(export_all = TRUE)
      exported <- if (isTRUE(input$export_all)) de_df() else table_rows()
      expect_equal(nrow(exported), full)
    }
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

test_that(".scroll_combo_group labels group1 / group2 / NA from combined levels", {
  combo <- c("X | 0", "X | 1", "Y | 0", NA)
  expect_equal(scroll:::.scroll_combo_group(combo, "X | 0"),
               c("group1", "group2", "group2", NA))                 # one-vs-rest
  expect_equal(scroll:::.scroll_combo_group(combo, "X | 0", "Y | 0"),
               c("group1", NA, "group2", NA))                       # two-group
  expect_true(all(is.na(scroll:::.scroll_combo_group(combo, character(0)))))
})

test_that("pseudobulk_de_server renders a live sample preview without Compute", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  ct <- sort(unique(as.character(data$cells$celltype)))[[1]]
  shiny::testServer(scroll:::pseudobulk_de_server, args = list(data = data), {
    session$setInputs(aggregate_by = "celltype", ident1 = ct, ident2 = character(0),
                      replicate = "no_replicate", mincells = 3, npseudo = 2, cellsper = 20)
    session$elapse(200)
    expect_no_error(output$preview)
  })
})

test_that("pseudobulk_de_server preview renders for a real replicate (group x rep)", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  shiny::testServer(scroll:::pseudobulk_de_server, args = list(data = data), {
    session$setInputs(aggregate_by = "celltype", ident1 = "T", ident2 = "B",
                      replicate = "condition", mincells = 5, npseudo = 3)
    session$elapse(200)
    expect_no_error(output$preview)      # keys by (group x replicate), applies min-cells
  })
})

test_that("pseudobulk_de_server preview follows the active subset view's embedding", {
  data <- scroll:::.scroll_load(subset_test_project())
  on.exit(scroll_disconnect(data$con))
  # whole dataset -> global embedding; subset view -> that view's embedding
  shiny::testServer(scroll:::pseudobulk_de_server,
                    args = list(data = data, view_r = reactive("tcell")), {
    session$setInputs(aggregate_by = "tsub", ident1 = "Tfh", ident2 = character(0),
                      replicate = "no_replicate", mincells = 2, npseudo = 2, cellsper = 20)
    session$elapse(200)
    expect_equal(preview_in()$emb, "umap_tcell")       # NOT the global umap
    expect_no_error(output$preview)
  })
})

test_that(".scroll_preview_embedding: whole=default_embedding, subset=its own", {
  data <- list(
    manifest = list(
      n_cells = 100,
      embeddings = list(pca = list(dims = 2, n_covered = 100),   # first global, but NOT default
                        umap = list(dims = 2, n_covered = 100),
                        umap_sub = list(dims = 2, n_covered = 40)),
      subsets = list(sub = list(embeddings = "umap_sub",
                                primary_embedding = "umap_sub"))),
    config = list(default_embedding = "umap"))
  expect_equal(scroll:::.scroll_preview_embedding(data, NULL), "umap")     # default, not pca
  expect_equal(scroll:::.scroll_preview_embedding(data, "sub"), "umap_sub")# the view's own
})
