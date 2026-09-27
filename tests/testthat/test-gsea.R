# Gene sets + preranked GSEA (R/genesets.R) and the GSEA panel (R/panel-gsea.R).

test_that("species is guessed from gene-symbol case", {
  expect_identical(scroll:::.scroll_guess_species(c("CD8A", "GZMB", "MS4A1", "RPL3")), "Homo sapiens")
  expect_identical(scroll:::.scroll_guess_species(c("Cd8a", "Gzmb", "Ms4a1", "mt-Co1")), "Mus musculus")
})

test_that(".gmt files are read and saved into a project, trimmed to its genes", {
  proj <- tempfile("gs"); dir.create(proj)
  file.copy(list.files(test_project(), full.names = TRUE), proj, recursive = TRUE)
  gmt <- tempfile(fileext = ".gmt")
  writeLines(c("T_CELL\tdesc\tCD3D\tCD8A\tNOTAGENE", "NK\t\tNKG7\tGNLY", "EMPTY\tx\tNOPE1"), gmt)
  sets <- scroll:::.scroll_read_gmt(gmt)
  expect_identical(names(sets), c("T_CELL", "NK", "EMPTY"))
  expect_identical(sets$T_CELL, c("CD3D", "CD8A", "NOTAGENE"))
  expect_message(scroll_add_genesets(proj, gmt = c(immune = gmt)), "saved immune")
  tab <- as.data.frame(arrow::read_parquet(file.path(proj, "genesets", "sets.parquet")))
  expect_setequal(unique(tab$set), c("T_CELL", "NK"))            # EMPTY dropped
  expect_false("NOTAGENE" %in% tab$gene)                          # trimmed to the assay
  m <- scroll_manifest(proj)
  expect_identical(m$genesets$collections$immune, 2L)
  expect_identical(m$genesets$species, "Homo sapiens")
  # the panel reads them back
  data <- scroll:::.scroll_load(proj); on.exit(scroll_disconnect(data$con), add = TRUE)
  expect_true("saved:immune" %in% scroll:::.scroll_geneset_choices(data))
  expect_identical(sort(names(scroll:::.scroll_geneset_get(data, "saved:immune"))), c("NK", "T_CELL"))
})

test_that("MSigDB collections are fetched through msigdbr", {
  skip_if_not_installed("msigdbr")
  h <- scroll:::.scroll_msigdb_sets("H", "Homo sapiens")
  expect_length(h, 50)
  expect_true("HALLMARK_INTERFERON_GAMMA_RESPONSE" %in% names(h))
  expect_error(scroll:::.scroll_msigdb_sets("NOPE", "Homo sapiens"), "Unknown MSigDB")
})

# A synthetic DE result: 2000 genes, set `top` placed at the head of the ranking.
fake_de <- function() {
  set.seed(1); g <- paste0("G", 1:2000)
  de <- data.frame(gene = g, avg_log2FC = round(rnorm(2000), 3), auc = round(runif(2000), 3),
                   p_val = runif(2000), stringsAsFactors = FALSE)
  attr(de, "ranking") <- data.frame(gene = g, avg_log2FC = c(seq(3, 2, length.out = 40), rnorm(1960)),
                                    auc = runif(2000), logFC = rnorm(2000), p_val = runif(2000))
  list(de = de, sets = list(top = g[1:40], rnd1 = sample(g, 50), rnd2 = sample(g, 30), tiny = g[1:5]))
}

test_that("scroll_gsea finds a planted set, deterministically, and matches fgsea", {
  skip_if_not_installed("fgsea")
  x <- fake_de()
  r <- expect_no_warning(scroll_gsea(x$de, x$sets))
  expect_identical(r$pathway[1], "top")
  expect_gt(r$NES[1], 2)
  expect_false("tiny" %in% r$pathway)                             # below min_size
  expect_identical(attr(r, "rank_by"), "avg_log2FC")
  expect_identical(scroll_gsea(x$de, x$sets)[, 1:6], r[, 1:6])    # seeded
  # same numbers as a direct fgsea call on the same statistic
  stats <- attr(r, "stats")
  set.seed(1)
  ref <- NULL
  utils::capture.output(ref <- fgsea::fgsea(attr(r, "sets"), stats, minSize = 15, maxSize = 500,
                                            eps = 0, nproc = 1))
  expect_equal(r$NES[match(ref$pathway, r$pathway)], ref$NES)
  # the ranking has no ties, even from rounded columns
  de2 <- x$de; attr(de2, "ranking") <- NULL
  expect_no_warning(r2 <- scroll_gsea(de2, x$sets, rank_by = "auc"))
  expect_false(anyDuplicated(attr(r2, "stats")) > 0)
  # a signed -log10 p ranking, and a clear error for a missing metric
  expect_s3_class(scroll_gsea(x$de, x$sets, rank_by = "signed_logp"), "data.frame")
  expect_error(scroll_gsea(x$de, x$sets, rank_by = "t"), "no 't' column")
})

test_that("the rank metrics offered follow the result type", {
  expect_identical(unname(scroll:::.scroll_rank_choices(
    data.frame(gene = "a", t = 1, logFC = 1, p_val = 0.1)))[1], "t")
  ch <- scroll:::.scroll_rank_choices(data.frame(gene = "a", avg_log2FC = 1, auc = 0.5, logFC = 1, p_val = 0.1))
  expect_identical(unname(ch[1]), "avg_log2FC")
  expect_false("logFC" %in% ch)                                   # presto's mean-of-logs is not offered
})

test_that("DE results reach the GSEA panel, which runs and plots", {
  skip_if_not_installed("fgsea"); skip_if_not_installed("presto")
  proj <- tempfile("gsp"); dir.create(proj)
  file.copy(list.files(test_project(), full.names = TRUE), proj, recursive = TRUE)
  genes <- scroll:::.scroll_features_of(scroll_manifest(proj), "RNA")
  gmt <- tempfile(fileext = ".gmt")
  writeLines(c(paste(c("SET_A", "", genes[1:10]), collapse = "\t"),
               paste(c("SET_B", "", genes[11:20]), collapse = "\t")), gmt)
  suppressMessages(scroll_add_genesets(proj, gmt = c(test = gmt)))
  data <- scroll:::.scroll_load(proj); on.exit(scroll_disconnect(data$con), add = TRUE)
  panels <- Filter(function(p) p$id %in% c("de", "gsea"), scroll:::.scroll_builtin_panels())
  shiny::testServer(function(input, output, session) {
    reg <- scroll:::.scroll_wire(input, output, session, data, panels)
  }, {
    expect_null(reg$results()$de)
    session$setInputs(`de-group` = "celltype", `de-ident1` = "T", `de-compute` = 1)
    s <- reg$results()$de
    expect_false(is.null(s))
    expect_match(s$label, "celltype: T vs rest")
    expect_true(all(c("avg_log2FC", "auc", "p_val") %in% names(s$ranking)))
    session$setInputs(`gsea-source` = "de", `gsea-collection` = "saved:test",
                      `gsea-rank_by` = "avg_log2FC", `gsea-min_size` = 5, `gsea-max_size` = 500,
                      `gsea-padj` = 1, `gsea-direction` = "both", `gsea-topn` = 10,
                      `gsea-show` = "top", `gsea-run` = 1)
  })
  # the panel server itself (inside its module) with a published result
  res <- shiny::reactiveVal(list())
  d2 <- data; d2$results <- res
  x <- fake_de()
  gmt2 <- tempfile(fileext = ".gmt")
  writeLines(vapply(names(x$sets), function(n) paste(c(n, "", x$sets[[n]]), collapse = "\t"), ""), gmt2)
  sets_id <- "saved:test"
  local_get <- function(data, id) x$sets
  testthat::local_mocked_bindings(.scroll_geneset_get = local_get, .package = "scroll")
  panel <- Filter(function(p) identical(p$id, "gsea"), scroll:::.scroll_builtin_panels())[[1]]
  shiny::testServer(panel$server, args = list(data = d2), {
    session$setInputs(source = "de", collection = sets_id, rank_by = "avg_log2FC",
                      min_size = 15, max_size = 500, padj = 1, direction = "both", topn = 10,
                      show = "top", palette = "viridis", run = 1)
    expect_error(res_r(), "DE or Pseudobulk")                     # nothing published yet
    res(list(de = list(df = x$de, ranking = attr(x$de, "ranking"), label = "fake",
                       n_genes = 2000, time = Sys.time())))
    session$setInputs(run = 2)
    expect_identical(res_r()$pathway[1], "top")
    expect_s3_class(plot_r(), "ggplot")                            # top pathways
    session$setInputs(show = "enrich", pathway = "top")
    expect_s3_class(plot_r(), "ggplot")                            # enrichment plot
    expect_gt(nrow(table_r()), 0)
  })
})

test_that("Pseudobulk results carry limma's t for ranking", {
  skip_if_not_installed("edgeR"); skip_if_not_installed("limma")
  data <- scroll:::.scroll_load(test_project()); on.exit(scroll_disconnect(data$con), add = TRUE)
  res <- scroll_pseudobulk_de(data, "RNA", aggregate_cols = "celltype", ident1 = "T",
                              ident2 = "B", replicate_col = "condition", min_cells = 5)
  rk <- attr(res, "ranking")
  expect_true(all(c("gene", "t", "logFC", "p_val") %in% names(rk)))
  expect_identical(unname(scroll:::.scroll_rank_choices(rk)[1]), "t")
})
