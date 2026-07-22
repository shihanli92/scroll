test_that("scroll_build emits the expected artifacts", {
  dir <- test_project()
  expect_true(file.exists(file.path(dir, "cells.parquet")))
  expect_true(file.exists(file.path(dir, "manifest.yaml")))
  expect_true(file.exists(file.path(dir, "config.yaml")))
  expect_true(file.exists(file.path(dir, "app.R")))
  expect_true(dir.exists(file.path(dir, "expr", "RNA")))

  cells <- as.data.frame(arrow::read_parquet(file.path(dir, "cells.parquet")))
  expect_true(all(c("cell", "condition", "celltype", "umap_1", "umap_2") %in% names(cells)))
  expect_equal(nrow(cells), 120)
})

test_that("manifest records assays, embeddings, and defaults", {
  man <- scroll_manifest(test_project())
  expect_equal(man$default_assay, "RNA")
  expect_equal(man$default_embedding, "umap")
  expect_true(man$quantize)
  expect_true("CD3D" %in% unlist(man$assays$RNA$features))
  expect_equal(man$embeddings$umap$dims, 2)
})

test_that("expr store is first-letter bucket-partitioned and features query", {
  dir <- test_project()
  parts <- list.files(file.path(dir, "expr", "RNA"))
  expect_true(all(grepl("^bucket=", parts)))     # v3 layout: partition by 1st char
  expect_true("bucket=M" %in% parts)             # MS4A1 -> bucket M
  con <- scroll_connect(dir); on.exit(scroll_disconnect(con))
  expect_gt(nrow(scroll_query_feature(con, "RNA", "MS4A1")), 0)
})

# A minimal object where one gene's only nonzero data value sits below the
# per-assay quantization floor (~max/510), so it survives full precision but
# rounds to zero under uint8 quantization. Pins the documented tradeoff.
.lowexpr_object <- function(seed = 1) {
  set.seed(seed)
  genes <- c("HIGH", "LOW"); n <- 6
  cm <- matrix(1L, nrow = 2, ncol = n, dimnames = list(genes, paste0("c", seq_len(n))))
  obj <- SeuratObject::CreateSeuratObject(counts = Matrix::Matrix(cm, sparse = TRUE),
                                          min.cells = 0, min.features = 0)
  dat <- matrix(0, nrow = 2, ncol = n, dimnames = dimnames(cm))
  dat["HIGH", 1] <- 100     # sets the per-assay max -> floor is 100/510 ~ 0.196
  dat["LOW", 1]  <- 0.05    # below the floor
  obj <- SeuratObject::SetAssayData(obj, layer = "data",
                                    new.data = Matrix::Matrix(dat, sparse = TRUE))
  emb <- matrix(rnorm(n * 2), ncol = 2, dimnames = list(colnames(obj), c("UMAP_1", "UMAP_2")))
  obj[["umap"]] <- SeuratObject::CreateDimReducObject(embeddings = emb, key = "UMAP_", assay = "RNA")
  obj$grp <- factor(rep(c("a", "b"), length.out = n))
  obj
}

test_that("quantization floors sub-max/510 values to zero; quantize=FALSE keeps them", {
  skip_if_not_installed("SeuratObject")
  obj <- .lowexpr_object()

  qdir <- file.path(tempdir(), "scroll-quant-on")
  suppressMessages(scroll_build(obj, qdir, assays = "RNA", meta_cols = "grp",
                                quantize = TRUE, overwrite = TRUE))
  con <- scroll_connect(qdir); on.exit(scroll_disconnect(con), add = TRUE)
  low_q <- scroll_query_feature(con, "RNA", "LOW")
  expect_true(all(low_q$value == 0))                 # floored away by quantization
  expect_true(any(scroll_query_feature(con, "RNA", "HIGH")$value > 0))
  scroll_disconnect(con); on.exit(NULL)

  fdir <- file.path(tempdir(), "scroll-quant-off")
  suppressMessages(scroll_build(obj, fdir, assays = "RNA", meta_cols = "grp",
                                quantize = FALSE, overwrite = TRUE))
  con2 <- scroll_connect(fdir); on.exit(scroll_disconnect(con2), add = TRUE)
  low_f <- scroll_query_feature(con2, "RNA", "LOW")
  expect_true(any(low_f$value > 0))                  # preserved at full precision
})

test_that("unsafe feature names (path separators) warn", {
  expect_warning(scroll:::.scroll_warn_unsafe_features(c("OK", "BAD/NAME")),
                 "path separators")
  expect_silent(scroll:::.scroll_warn_unsafe_features(c("CD3D", "MS4A1")))
})

test_that("building into a non-empty dir without overwrite errors", {
  skip_if_not_installed("SeuratObject")
  dir <- file.path(tempdir(), "scroll-nonempty")
  dir.create(dir, showWarnings = FALSE)
  writeLines("x", file.path(dir, "sentinel.txt"))    # make it non-empty
  expect_error(
    suppressMessages(scroll_build(.lowexpr_object(), dir, assays = "RNA",
                                  meta_cols = "grp", overwrite = FALSE)),
    "use overwrite = TRUE")
})

# A tiny object carrying one high-cardinality categorical column (`clone`, 300
# distinct values) and one low-cardinality one (`grp`, 3), passed via explicit
# meta_cols so the inference cap is bypassed.
.highcard_object <- function(n = 300, seed = 1) {
  set.seed(seed)
  genes <- c("AAA", "BBB")
  cm <- matrix(rpois(length(genes) * n, 1), nrow = length(genes),
               dimnames = list(genes, paste0("c", seq_len(n))))
  obj <- SeuratObject::CreateSeuratObject(counts = Matrix::Matrix(cm, sparse = TRUE),
                                          min.cells = 0, min.features = 0)
  obj <- SeuratObject::SetAssayData(obj, layer = "data",
                                    new.data = Matrix::Matrix(.lognorm(cm, 1e4), sparse = TRUE))
  emb <- matrix(rnorm(n * 2), ncol = 2, dimnames = list(colnames(obj), c("UMAP_1", "UMAP_2")))
  obj[["umap"]] <- SeuratObject::CreateDimReducObject(embeddings = emb, key = "UMAP_", assay = "RNA")
  obj$clone <- paste0("cl", seq_len(n))                    # 300 distinct
  obj$grp   <- factor(rep(c("a", "b", "c"), length.out = n))
  obj
}

test_that("manifest caps cached levels for high-cardinality columns", {
  skip_if_not_installed("SeuratObject")
  obj <- .highcard_object()
  dir <- file.path(tempdir(), "scroll-levels-cap")
  expect_message(
    scroll_build(obj, dir, assays = "RNA", meta_cols = c("clone", "grp"),
                 overwrite = TRUE),
    "high-cardinality")
  man <- scroll_manifest(dir)
  # high-cardinality column: type + n_levels recorded, cached levels omitted
  expect_equal(man$meta$clone$type, "categorical")
  expect_equal(man$meta$clone$n_levels, 300L)
  expect_null(man$meta$clone$levels)
  # low-cardinality column keeps its cached levels
  expect_setequal(unlist(man$meta$grp$levels), c("a", "b", "c"))
  expect_equal(man$meta$grp$n_levels, 3L)
})

test_that("max_levels raised caches the full level list", {
  skip_if_not_installed("SeuratObject")
  obj <- .highcard_object()
  dir <- file.path(tempdir(), "scroll-levels-cap-hi")
  suppressMessages(scroll_build(obj, dir, assays = "RNA",
                                meta_cols = c("clone", "grp"),
                                max_levels = 1000, overwrite = TRUE))
  man <- scroll_manifest(dir)
  expect_equal(length(man$meta$clone$levels), 300L)
})
