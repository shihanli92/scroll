# Build a small scroll project once per test run from a synthetic Seurat object,
# so tests need no external data download. Cached in a temp dir.

# Hand log-normalize a counts matrix (avoids a Seurat dependency in tests).
.lognorm <- function(cm, scale) log1p(sweep(cm, 2, pmax(colSums(cm), 1), "/") * scale)

make_test_object <- function(n = 120, seed = 1) {
  set.seed(seed)
  genes <- c("CD3D", "CD8A", "MS4A1", "CD14", "NKG7", "GNLY",
             paste0("G", 1:20))
  cm <- matrix(rpois(length(genes) * n, lambda = 0.8),
               nrow = length(genes), dimnames = list(genes, paste0("cell", 1:n)))
  norm <- .lognorm(cm, 1e4)
  counts <- Matrix::Matrix(cm, sparse = TRUE)
  obj <- SeuratObject::CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0)
  obj <- SeuratObject::SetAssayData(obj, layer = "data",
                                    new.data = Matrix::Matrix(norm, sparse = TRUE))

  # A second (protein / ADT) assay so multi-assay paths are exercised. Distinct
  # feature namespace and a different normalization scale, so its per-assay `max`
  # differs from RNA's (pins the per-assay dequantization).
  adt <- c("CD3-P", "CD4-P", "CD8-P", "CD19-P", "CD14-P", "CD56-P", "CD16-P", "PTPRC-P")
  am <- matrix(rpois(length(adt) * n, lambda = 3),
               nrow = length(adt), dimnames = list(adt, colnames(obj)))
  obj[["ADT"]] <- SeuratObject::CreateAssay5Object(counts = Matrix::Matrix(am, sparse = TRUE))
  obj <- SeuratObject::SetAssayData(obj, assay = "ADT", layer = "data",
                                    new.data = Matrix::Matrix(.lognorm(am, 100), sparse = TRUE))

  emb <- matrix(rnorm(n * 2), ncol = 2, dimnames = list(colnames(obj), c("UMAP_1", "UMAP_2")))
  obj[["umap"]] <- SeuratObject::CreateDimReducObject(embeddings = emb, key = "UMAP_", assay = "RNA")
  obj$condition <- factor(rep(c("ctrl", "treat"), length.out = n))
  obj$celltype <- factor(sample(c("T", "B", "NK"), n, replace = TRUE))
  obj
}

test_project <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      dir <- file.path(tempdir(), "scroll-test-proj")
      if (!dir.exists(file.path(dir, "expr"))) {
        suppressMessages(scroll_build(make_test_object(), dir,
                                      assays = c("RNA", "ADT"), overwrite = TRUE))
      }
      cached <<- dir
    }
    cached
  }
})
