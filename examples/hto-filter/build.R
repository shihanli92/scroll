#!/usr/bin/env Rscript
# Build a scroll project from two mouse 10x samples (RNA + hashtag/HTO), for the
# QC + hashtag-gating example app. Run once, offline:
#
#   Rscript examples/hto-filter/build.R
#
# Produces examples/hto-filter/project/ (git-ignored). The app (app.R) reads it.

suppressMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})
devtools::load_all("/Users/shihanl1/scroll", quiet = TRUE)   # the scroll package

root   <- "/Users/shihanl1/scroll"
inputs <- c(s1 = file.path(root, "data/s1_sample_filtered_feature_bc_matrix.h5"),
            s2 = file.path(root, "data/s2_sample_filtered_feature_bc_matrix.h5"))
outdir <- file.path(root, "examples/hto-filter/project")
ckpt   <- file.path(root, "examples/hto-filter/.processed.rds")  # cache (git-ignored)

# One Seurat object per sample: RNA (Gene Expression) + an HTO assay (the 3
# Antibody-Capture hashtags). Cells are labelled by sample.
make_obj <- function(h5, sample) {
  x   <- Read10X_h5(h5)
  rna <- x[["Gene Expression"]]
  hto <- x[["Antibody Capture"]]                     # Hashtag_1 / _2 / _3
  obj <- CreateSeuratObject(counts = rna, project = sample)
  obj[["HTO"]] <- CreateAssay5Object(counts = hto)
  obj$sample <- sample
  obj
}

# The read + normalize + PCA/UMAP pipeline is the slow part; cache the processed
# object so re-runs (e.g. iterating on the build settings) skip it. Delete
# .processed.rds to force a full recompute.
if (file.exists(ckpt)) {
  message("Loading cached processed object: ", ckpt)
  obj <- readRDS(ckpt)
} else {

message("Reading ", inputs[["s1"]])
o1 <- make_obj(inputs[["s1"]], "s1")
message("Reading ", inputs[["s2"]])
o2 <- make_obj(inputs[["s2"]], "s2")

# Merge; add.cell.ids prefixes barcodes (s1_/s2_) so ids stay unique across
# samples. JoinLayers consolidates the per-sample layers into one matrix so the
# downstream pipeline (and scroll_build's GetAssayData(layer="data")) see a
# single joined layer.
obj <- merge(o1, o2, add.cell.ids = c("s1", "s2"))
obj <- JoinLayers(obj)                               # RNA
obj <- JoinLayers(obj, assay = "HTO")
message("Merged: ", ncol(obj), " cells x ", nrow(obj), " genes")

# --- RNA QC ------------------------------------------------------------------
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")   # mouse mito

# --- HTO: CLR-normalize, then expose per-cell values as metadata for gating ---
# margin = 2 normalizes each hashtag across cells, so an independent threshold
# per hashtag cleanly separates positive from background.
obj <- NormalizeData(obj, assay = "HTO", normalization.method = "CLR", margin = 2,
                     verbose = FALSE)
hto_clr <- as.matrix(GetAssayData(obj, assay = "HTO", layer = "data"))   # 3 x cells
for (h in rownames(hto_clr)) obj[[paste0("hto_", h)]] <- hto_clr[h, ]

# scroll_build's exporter expects a sparse (dgCMatrix) `data` layer; the small
# CLR'd HTO layer comes back dense, so store it back as sparse.
obj <- SetAssayData(obj, assay = "HTO", layer = "data",
                    new.data = Matrix::Matrix(hto_clr, sparse = TRUE))

# --- RNA processing for the built-in viz panels (+ satisfies scroll_build) ----
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)

  saveRDS(obj, ckpt)
}

# --- Build the scroll project ------------------------------------------------
hto_cols <- grep("^hto_", colnames(obj[[]]), value = TRUE)
scroll_build(
  obj, outdir,
  assays    = c("RNA", "HTO"),
  meta_cols = c("sample", "nCount_RNA", "nFeature_RNA", "percent.mt", hto_cols),
  overwrite = TRUE)

message("DONE: ", outdir)
