# scroll_add_meta: add a column to a built project with no Seurat object.
# Builds fresh projects in tempdir (it writes cells.parquet in place -- never touch
# the shared test_project() fixture).

fresh_project <- function(tag) {
  skip_if_not_installed("SeuratObject")
  dir <- file.path(tempdir(), paste0("scroll-addmeta-", tag))
  suppressMessages(scroll_build(make_test_object(), dir, assays = "RNA", overwrite = TRUE))
  dir
}

test_that("adds a categorical column (function form) selectable as categorical", {
  dir <- fresh_project("cat")
  scroll_add_meta(dir, "grp2",
                  function(cells) ifelse(cells$celltype == "T", "Tcell", "other"))
  cells <- as.data.frame(arrow::read_parquet(file.path(dir, "cells.parquet")))
  expect_true("grp2" %in% names(cells))
  m <- scroll_manifest(dir)
  expect_equal(m$meta$grp2$type, "categorical")
  expect_setequal(unlist(m$meta$grp2$levels), c("Tcell", "other"))
  expect_true("grp2" %in% scroll:::.scroll_cat_cols(m))
  # query layer / expr store untouched
  con <- scroll_connect(dir); on.exit(scroll_disconnect(con))
  expect_gt(nrow(scroll_query_feature(con, "RNA", "MS4A1")), 0)
})

test_that("adds a numeric column (row-order vector) with a manifest range", {
  dir <- fresh_project("num")
  n <- nrow(as.data.frame(arrow::read_parquet(file.path(dir, "cells.parquet"))))
  scroll_add_meta(dir, "score", seq_len(n) / n)
  m <- scroll_manifest(dir)
  expect_equal(m$meta$score$type, "numeric")
  # manifest ranges round-trip through yaml at ~7 sig figs (they seed slider bounds)
  expect_equal(m$meta$score$range$min, 1 / n, tolerance = 1e-5)
  expect_equal(m$meta$score$range$max, 1, tolerance = 1e-5)
  expect_true("score" %in% scroll:::.scroll_num_cols(m))
})

test_that("named vector aligns by barcode; unmatched cells become NA", {
  dir <- fresh_project("named")
  cells <- as.data.frame(arrow::read_parquet(file.path(dir, "cells.parquet")))
  half <- cells$cell[seq_len(floor(nrow(cells) / 2))]
  vals <- setNames(rep("hit", length(half)), half)
  scroll_add_meta(dir, "flag", vals)
  got <- as.data.frame(arrow::read_parquet(file.path(dir, "cells.parquet")))$flag
  expect_equal(sum(!is.na(got)), length(half))
  expect_setequal(unique(stats::na.omit(got)), "hit")
})

test_that("overwrite guard, reserved name, and replace-in-place", {
  dir <- fresh_project("ow")
  scroll_add_meta(dir, "lab", rep("a", nrow(
    as.data.frame(arrow::read_parquet(file.path(dir, "cells.parquet"))))))
  expect_error(scroll_add_meta(dir, "lab", "x"), "already exists")
  expect_error(scroll_add_meta(dir, "cell", "x"), "cannot be 'cell'")
  n <- nrow(as.data.frame(arrow::read_parquet(file.path(dir, "cells.parquet"))))
  scroll_add_meta(dir, "lab", rep("b", n), overwrite = TRUE)
  m <- scroll_manifest(dir)
  expect_setequal(unlist(m$meta$lab$levels), "b")
})

test_that("scoped column carries scope and member-only levels", {
  skip_if_not_installed("SeuratObject")
  dir <- file.path(tempdir(), "scroll-addmeta-scope")
  obj <- make_test_object()
  suppressMessages(scroll_build(obj, dir, assays = "RNA", overwrite = TRUE))
  tcells <- colnames(obj)[obj$celltype == "T"]
  sub <- subset(obj, cells = tcells)
  emb <- matrix(rnorm(length(tcells) * 2), ncol = 2,
                dimnames = list(tcells, c("UMAP_1", "UMAP_2")))
  sub[["umap_t"]] <- SeuratObject::CreateDimReducObject(embeddings = emb, key = "umapt_", assay = "RNA")
  obj2 <- scroll_add_subset(obj, sub, "tcell", embeddings = c(t_umap = "umap_t"), label = "T cells")
  suppressMessages(scroll_update(dir, obj2, embeddings = "t_umap"))

  # column defined for members (T cells), NA elsewhere; scope it to the view
  cells <- as.data.frame(arrow::read_parquet(file.path(dir, "cells.parquet")))
  v <- ifelse(cells$cell %in% tcells,
              sample(c("x", "y"), nrow(cells), replace = TRUE), NA_character_)
  scroll_add_meta(dir, "tsub", v, scope = "tcell")
  m <- scroll_manifest(dir)
  expect_equal(m$meta$tsub$scope, "tcell")
  expect_setequal(unlist(m$meta$tsub$levels), c("x", "y"))     # member-only, no NA leak
})
