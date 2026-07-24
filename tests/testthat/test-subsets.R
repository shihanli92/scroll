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
  expect_equal(spec$tcell$meta, "tsub")
  # and it survives an assay-data edit (a bare attribute would not)
  obj2 <- SeuratObject::SetAssayData(obj, layer = "data",
                                     new.data = SeuratObject::GetAssayData(obj, layer = "data"))
  expect_equal(SeuratObject::Misc(obj2, "scroll_subsets")$tcell$embeddings, "umap_tcell")
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
