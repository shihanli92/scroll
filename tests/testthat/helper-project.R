# Build a small scroll project once per test run from a synthetic Seurat object,
# so tests need no external data download. Cached in a temp dir.

make_test_object <- function(n = 120, seed = 1) {
  set.seed(seed)
  genes <- c("CD3D", "CD8A", "MS4A1", "CD14", "NKG7", "GNLY",
             paste0("G", 1:20))
  cm <- matrix(rpois(length(genes) * n, lambda = 0.8),
               nrow = length(genes), dimnames = list(genes, paste0("cell", 1:n)))
  # Log-normalize by hand (avoids a Seurat dependency in tests).
  libsize <- pmax(colSums(cm), 1)
  norm <- log1p(sweep(cm, 2, libsize, "/") * 1e4)
  counts <- Matrix::Matrix(cm, sparse = TRUE)
  obj <- SeuratObject::CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0)
  obj <- SeuratObject::SetAssayData(obj, layer = "data",
                                    new.data = Matrix::Matrix(norm, sparse = TRUE))
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
        suppressMessages(scroll_build(make_test_object(), dir, overwrite = TRUE))
      }
      cached <<- dir
    }
    cached
  }
})
