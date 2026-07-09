#!/usr/bin/env Rscript
# Reprocess the mouse RNA+HTO dataset down to the final selection: HTO singlets
# (from the gating export) that also pass the RNA QC filters (percent.mt < 10,
# nFeature_RNA < 6000). The kept cells are re-normalized and re-embedded from
# scratch, so the PCA/UMAP reflect only the clean population. Run after gating:
#
#   Rscript examples/hto-filter/reprocess.R
#
# Outputs (git-ignored): filtered.rds (the reprocessed Seurat object) and a
# project-filtered/ scroll project you can serve to explore the clean data.

suppressMessages({
  library(Seurat); library(SeuratObject); library(Matrix)
})
devtools::load_all("/Users/shihanl1/scroll", quiet = TRUE)

root   <- "/Users/shihanl1/scroll"
ckpt   <- file.path(root, "examples/hto-filter/.processed.rds")   # full processed object
csv    <- file.path(root, "data/hto_keep_barcodes.csv")           # gating export
rds    <- file.path(root, "examples/hto-filter/filtered.rds")
outdir <- file.path(root, "examples/hto-filter/project-filtered")

MT_MAX     <- 10       # percent.mt filter
NFEAT_MAX  <- 6000     # nFeature_RNA filter

stopifnot(file.exists(ckpt), file.exists(csv))
obj <- readRDS(ckpt)

# HTO singlets from the saved gating (assignment is one hashtag).
keep <- utils::read.csv(csv, stringsAsFactors = FALSE)
singlet_labels <- grep("^Hashtag", unique(keep$assignment), value = TRUE)
sing <- keep[keep$assignment %in% singlet_labels, , drop = FALSE]
obj$hto <- sing$assignment[match(colnames(obj), sing$cell)]        # NA for non-singlets

# Final selection: singlet AND passing both QC filters.
keep_cells <- colnames(obj)[
  !is.na(obj$hto) & obj$percent.mt < MT_MAX & obj$nFeature_RNA < NFEAT_MAX]
message(sprintf("Keeping %d of %d cells (singlets, mt < %g%%, nFeature < %d)",
                length(keep_cells), ncol(obj), MT_MAX, NFEAT_MAX))
sub <- subset(obj, cells = keep_cells)

# Re-process RNA from counts on the filtered population.
DefaultAssay(sub) <- "RNA"
sub <- NormalizeData(sub, verbose = FALSE)
sub <- FindVariableFeatures(sub, verbose = FALSE)
sub <- ScaleData(sub, verbose = FALSE)
sub <- RunPCA(sub, npcs = 30, verbose = FALSE)
sub <- RunUMAP(sub, dims = 1:30, verbose = FALSE)
# keep the HTO data layer sparse for scroll_build's exporter
sub <- SetAssayData(sub, assay = "HTO", layer = "data",
                    new.data = as(GetAssayData(sub, assay = "HTO", layer = "data"),
                                  "CsparseMatrix"))

saveRDS(sub, rds)
message("Saved reprocessed object: ", rds)

# Build a scroll project of the clean data.
hto_cols <- grep("^hto_", colnames(sub[[]]), value = TRUE)
scroll_build(
  sub, outdir, assays = c("RNA", "HTO"),
  meta_cols = c("sample", "hto", "nCount_RNA", "nFeature_RNA", "percent.mt", hto_cols),
  overwrite = TRUE)
message("DONE: ", outdir)
