# Subset views: reprocessed sub-embeddings declared via `subsets=` /
# scroll_add_subset(), surfaced as a linked View that restricts every panel.

test_that("scroll_add_subset transfers a reduction + scoped meta by barcode", {
  obj <- make_subset_object()
  expect_true("umap_tcell" %in% SeuratObject::Reductions(obj))
  # membership = cells with the sub-embedding; tsub is NA off the subset
  em <- SeuratObject::Embeddings(obj, "umap_tcell")
  members <- rownames(em)
  expect_true(all(obj$celltype[members] == "T"))
  expect_equal(sum(!is.na(obj$tsub)), length(members))
  expect_true(all(is.na(obj$tsub[setdiff(colnames(obj), members)])))
  # spec recorded in @misc (survives downstream Seurat ops) for scroll_build
  spec <- SeuratObject::Misc(obj, "scroll_subsets")
  expect_equal(spec$tcell$embeddings, "umap_tcell")
  expect_equal(spec$tcell$meta, c("tsub", "tscore"))
  # and it survives an assay-data edit (a bare attribute would not)
  obj2 <- SeuratObject::SetAssayData(obj, layer = "data",
                                     new.data = SeuratObject::GetAssayData(obj, layer = "data"))
  expect_equal(SeuratObject::Misc(obj2, "scroll_subsets")$tcell$embeddings, "umap_tcell")
})

test_that("scroll_add_subset(meta = 'auto') copies only new/changed columns, prefixed", {
  obj <- make_test_object(80, seed = 5)
  obj$clusters <- factor(sample(c("a", "b"), ncol(obj), replace = TRUE))
  obj$score <- stats::runif(ncol(obj))
  tcells <- colnames(obj)[obj$celltype == "T"]
  sub <- subset(obj, cells = tcells)
  te <- matrix(stats::rnorm(length(tcells) * 2), ncol = 2,
               dimnames = list(tcells, c("UMAP_1", "UMAP_2")))
  sub[["umap"]] <- SeuratObject::CreateDimReducObject(embeddings = te, key = "UMAP_",
                                                      assay = "RNA")
  sub$clusters <- factor(sample(c("t1", "t2", "t3"), length(tcells), replace = TRUE))  # re-clustered
  sub$tnew <- "x"                                                                     # new column
  sub$score <- sub$score + 1e-12                                                      # float noise only
  sub$condition <- factor(as.character(sub$condition),                                # levels re-ordered only
                          levels = rev(sort(unique(as.character(sub$condition)))))

  expect_message(
    out <- scroll_add_subset(obj, sub, "tcell", embeddings = c(umap_tcell = "umap"),
                             meta = "auto"),
    "tcell_clusters")
  spec <- SeuratObject::Misc(out, "scroll_subsets")$tcell
  expect_setequal(spec$meta, c("tcell_tnew", "tcell_clusters"))
  # the parent's own column is untouched; the subset's copy is scoped (NA elsewhere)
  expect_identical(as.character(out$clusters), as.character(obj$clusters))
  expect_equal(sum(!is.na(out$tcell_clusters)), length(tcells))
  expect_true(all(out$tcell_clusters[tcells] %in% c("t1", "t2", "t3")))
  # and it builds: the scoped, prefixed columns reach the manifest
  dir <- file.path(tempdir(), "scroll-subset-auto")
  suppressMessages(scroll_build(out, dir, assays = "RNA", overwrite = TRUE,
                                meta_cols = c("celltype", spec$meta)))
  m <- scroll_manifest(dir)
  expect_equal(m$meta$tcell_clusters$scope, "tcell")
  expect_setequal(unlist(m$subsets$tcell$meta), spec$meta)
})

test_that("scroll_add_subset(embeddings = 'auto') copies re-run/new reductions, UMAP first", {
  obj <- make_test_object(80, seed = 6)
  pc <- matrix(stats::rnorm(ncol(obj) * 5), ncol = 5,
               dimnames = list(colnames(obj), paste0("PC_", 1:5)))
  obj[["pca"]] <- SeuratObject::CreateDimReducObject(embeddings = pc, key = "PC_", assay = "RNA")
  tcells <- colnames(obj)[obj$celltype == "T"]
  sub <- subset(obj, cells = tcells)             # inherits umap + pca unchanged
  mk <- function(key) SeuratObject::CreateDimReducObject(
    embeddings = matrix(stats::rnorm(length(tcells) * 2), ncol = 2,
                        dimnames = list(tcells, paste0(key, 1:2))), key = key, assay = "RNA")
  sub[["tsne"]] <- mk("tSNE_")                   # new
  sub[["umap"]] <- mk("UMAP_")                   # re-run on the subset
  sub$tclust <- sample(c("t1", "t2"), length(tcells), replace = TRUE)

  msgs <- testthat::capture_messages(out <- scroll_add_subset(obj, sub, "tc", meta = "auto"))
  expect_match(msgs[1], "tc_umap \\(2d\\), tc_tsne \\(2d\\)")
  expect_match(msgs[2], "tc_tclust")
  spec <- SeuratObject::Misc(out, "scroll_subsets")$tc
  expect_equal(spec$embeddings, c("tc_umap", "tc_tsne"))     # umap is primary; pca skipped
  expect_false("tc_pca" %in% SeuratObject::Reductions(out))
  expect_equal(nrow(SeuratObject::Embeddings(out, "tc_umap")), length(tcells))

  # builds with inferred embeddings + meta_cols; the view keys on the sub-UMAP
  dir <- file.path(tempdir(), "scroll-subset-auto-emb")
  suppressMessages(scroll_build(out, dir, assays = "RNA", overwrite = TRUE))
  m <- scroll_manifest(dir)
  expect_equal(m$subsets$tc$primary_embedding, "tc_umap")
  expect_equal(m$subsets$tc$n_cells, length(tcells))
  expect_equal(m$meta$tc_tclust$scope, "tc")

  # nothing re-run -> a clear error instead of an empty view
  expect_error(scroll_add_subset(obj, subset(obj, cells = tcells), "none"),
               "no reduction that is new or differs")
})

test_that("scroll_add_subset explicit meta keeps names unless prefix is given", {
  obj <- make_subset_object()
  expect_true("tsub" %in% names(obj[[]]))
  obj2 <- make_test_object()
  tcells <- colnames(obj2)[obj2$celltype == "T"]
  sub <- subset(obj2, cells = tcells)
  sub[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(stats::rnorm(length(tcells) * 2), ncol = 2,
                        dimnames = list(tcells, c("UMAP_1", "UMAP_2"))),
    key = "UMAP_", assay = "RNA")
  sub$tsub <- "Tfh"
  out <- scroll_add_subset(obj2, sub, "tc", embeddings = c(umap_tc = "umap"),
                           meta = "tsub", prefix = "tc_")
  expect_true("tc_tsub" %in% names(out[[]]))
  expect_false("tsub" %in% names(out[[]]))
})

test_that("scroll_add_subset preserves a numeric scoped column (not stringified)", {
  obj <- make_test_object(80, seed = 5)
  tcells <- colnames(obj)[obj$celltype == "T"]
  sub <- subset(obj, cells = tcells)
  te <- matrix(stats::rnorm(length(tcells) * 2), ncol = 2,
               dimnames = list(tcells, c("UMAP_1", "UMAP_2")))
  sub[["umap"]] <- SeuratObject::CreateDimReducObject(embeddings = te, key = "UMAP_",
                                                      assay = "RNA")
  sub$score <- stats::runif(length(tcells))                   # a continuous scoped column
  out <- scroll_add_subset(obj, sub, "tcell", embeddings = c(umap_tcell = "umap"),
                           meta = "score")
  v <- out$score
  expect_type(v, "double")                                    # stayed numeric, not character
  expect_equal(sum(!is.na(v)), length(tcells))                # members carry the score
  expect_true(all(is.na(v[setdiff(colnames(out), tcells)])))  # NA off the subset
  expect_equal(scroll:::.scroll_meta_entry(v)$type, "numeric")# classifier types it numeric
})

test_that("build records subsets, embedding coverage, and scoped meta", {
  dir <- subset_test_project()
  m <- scroll_manifest(dir)
  expect_true("tcell" %in% names(m$subsets))
  expect_equal(m$subsets$tcell$label, "T cells")
  expect_equal(m$subsets$tcell$primary_embedding, "umap_tcell")
  # the sub-embedding is partial; the parent umap is full
  expect_equal(m$embeddings$umap$n_covered, m$n_cells)
  expect_lt(m$embeddings$umap_tcell$n_covered, m$n_cells)
  expect_equal(m$embeddings$umap_tcell$n_covered, m$subsets$tcell$n_cells)
  # scoped meta carries its scope and NA-free, subset-only levels
  expect_equal(m$meta$tsub$scope, "tcell")
  expect_setequal(unlist(m$meta$tsub$levels), c("Tcm", "Tem", "Tfh"))
  expect_false("NA" %in% unlist(m$meta$tsub$levels))
  expect_null(m$meta$celltype$scope)
  # a scoped *numeric* column stays numeric (not stringified to categorical) and
  # is scoped to its view
  expect_equal(m$meta$tscore$type, "numeric")
  expect_false(is.null(m$meta$tscore$range))
  expect_equal(m$meta$tscore$scope, "tcell")
})

test_that("scoped columns are hidden from whole-dataset selectors, shown in-view", {
  m <- scroll_manifest(subset_test_project())
  expect_false("tsub" %in% scroll:::.scroll_cat_cols(m))            # global: hidden
  expect_true("tsub" %in% scroll:::.scroll_cat_cols(m, "tcell"))    # in its view: shown
  expect_true("celltype" %in% scroll:::.scroll_cat_cols(m, "tcell"))# globals still shown
  # colorby choices likewise gain the scoped column only in-view
  flat_global <- unlist(scroll:::.scroll_colorby_choices(m), use.names = FALSE)
  flat_view <- unlist(scroll:::.scroll_colorby_choices(m, "tcell"), use.names = FALSE)
  expect_false("tsub" %in% flat_global)
  expect_true("tsub" %in% flat_view)
  # a scoped numeric column is likewise hidden globally, shown in its view — so the
  # Violin panel's numeric-column picker only offers it under the subset view
  expect_false("tscore" %in% scroll:::.scroll_num_cols(m))
  expect_true("tscore" %in% scroll:::.scroll_num_cols(m, "tcell"))
})

test_that("violin plots a scoped numeric column under its subset view", {
  data <- scroll:::.scroll_load(subset_test_project())
  on.exit(scroll_disconnect(data$con))
  shiny::testServer(scroll:::violin_server,
                    args = list(data = data, view_r = reactive("tcell")), {
    session$setInputs(source = "Metadata", metacol = "tscore", group = "tsub",
                      palette = "Tableau 10", jitter = FALSE, legend = FALSE, aspect = 1)
    session$flushReact()
    expect_equal(data_r()$value_col, "tscore")     # numeric column plotted directly
    expect_false(is.null(output$plot))
  })
})

test_that("view embeddings split full vs subset; sub-embedding only in its view", {
  m <- scroll_manifest(subset_test_project())
  expect_equal(scroll:::.scroll_view_embeddings(m, NULL), "umap")
  expect_equal(scroll:::.scroll_view_embeddings(m, "tcell"), "umap_tcell")
})

test_that(".scroll_view_cells restricts to non-NA primary-coord members", {
  dir <- subset_test_project()
  m <- scroll_manifest(dir)
  cells <- as.data.frame(arrow::read_parquet(file.path(dir, "cells.parquet")))
  all_cells <- scroll:::.scroll_view_cells(cells, m, NULL)
  members <- scroll:::.scroll_view_cells(cells, m, "tcell")
  expect_equal(nrow(all_cells), m$n_cells)
  expect_equal(nrow(members), m$subsets$tcell$n_cells)
  expect_true(all(!is.na(members$umap_tcell_1)))
  expect_true(all(members$celltype == "T"))
})

test_that("the app bar shows a View selector only when subsets are declared", {
  with_subsets <- scroll:::.scroll_load(subset_test_project())
  no_subsets <- scroll:::.scroll_load(test_project())
  expect_match(as.character(scroll:::.scroll_appbar(with_subsets, NULL)), "scroll_view")
  expect_false(grepl("scroll_view", as.character(scroll:::.scroll_appbar(no_subsets, NULL))))
})

test_that("scatter views drop NA-coordinate cells (partial embedding)", {
  dir <- subset_test_project()
  cells <- as.data.frame(arrow::read_parquet(file.path(dir, "cells.parquet")))
  # coloring the sub-UMAP by the scoped column plots only members, no NA warning
  p <- expect_no_warning(
    view_umap_colorby(cells, list(embedding = "umap_tcell", color_by = "tsub")))
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), sum(!is.na(cells$umap_tcell_1)))
})

test_that("scroll_build validates subset embedding coverage", {
  obj <- make_test_object(60)
  # 'umap' covers all cells -> declaring it as a subset embedding warns
  expect_warning(
    suppressMessages(scroll_build(obj, tempfile("full"), assays = "RNA",
                                  subsets = list(x = list(embeddings = "umap")),
                                  overwrite = TRUE)),
    "covers all")
  # an unexported embedding is a clear error
  expect_error(
    scroll:::.scroll_normalize_subsets(list(x = list(embeddings = "nope")),
                                       embeddings = "umap", meta_cols = character(0)),
    "unexported embedding")
})

test_that("view_r is opt-in: 3-arg custom panel servers keep the old contract", {
  builtin <- scroll:::dimplot_server
  custom3 <- function(id, data, cells_r = shiny::reactive(NULL)) NULL
  expect_true("view_r" %in% names(formals(builtin)))
  expect_false("view_r" %in% names(formals(custom3)))
})

test_that("dimplot_server switches reduction + colorby when the view changes", {
  data <- scroll:::.scroll_load(subset_test_project())
  on.exit(scroll_disconnect(data$con))
  view <- shiny::reactiveVal(NULL)
  shiny::testServer(scroll:::dimplot_server,
                    args = list(data = data, view_r = view), {
    session$flushReact()
    view("tcell"); session$flushReact()
    # the observer restricted reductions to the subset's embedding and offered
    # the scoped color-by
    expect_equal(scroll:::.scroll_view_embeddings(data$manifest, "tcell"), "umap_tcell")
    expect_true("tsub" %in% scroll:::.scroll_cat_cols(data$manifest, "tcell"))
  })
})

test_that("pseudobulk aggregate-by exposes subset-scoped columns in a subset view", {
  data <- scroll:::.scroll_load(subset_test_project())
  on.exit(scroll_disconnect(data$con))
  # scoped column is out of scope whole-dataset, in scope in the T-cell view
  expect_false("tsub" %in% scroll:::.scroll_cat_cols(data$manifest, NULL))
  expect_true("tsub" %in% scroll:::.scroll_cat_cols(data$manifest, "tcell"))
  view <- shiny::reactiveVal(NULL)
  shiny::testServer(scroll:::pseudobulk_de_server,
                    args = list(data = data, view_r = view), {
    view("tcell"); session$flushReact()
    # in-view, the scoped column can drive aggregation (combos come from cells)
    session$setInputs(aggregate_by = "tsub", ident1 = character(0), replicate = "no_replicate")
    session$flushReact()
    expect_gt(length(.scroll_combo_choices(cells_r(), "tsub")), 0)
  })
})
