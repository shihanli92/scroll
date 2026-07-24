#!/usr/bin/env Rscript
# Unified Cao et al. scroll project: the RAW all-cells object (46,606 cells) is the
# base "whole dataset"; every downstream population is a SUBSET VIEW attached by
# barcode. Reuses the already-processed objects (no re-normalize/UMAP):
#   .processed.rds       -> raw 46,606-cell base (RNA + HTO, umap/pca, QC + HTO CLR)
#   filtered-subsets.rds -> the clean singlets (32,611) carrying umap + the six
#                           sub-population UMAPs and all their subcluster columns
#
# Clustering (resolutions 0.1..1.0) exists for EVERY object:
#   res_<r>          whole dataset (46,606)  -- computed here on the raw base
#   clean_res_<r>    Clean singlets (32,611) -- from filtered-subsets.rds
#   <id>_res_<r>     each sub-population      -- from filtered-subsets.rds
# percent.ribo is computed on the raw base (global, all 46,606 cells).
#
# Subset views: clean / cd101 / tfh / cx3cr1 / cycling / sgctrl / tpex_sgctrl
#
#   Rscript examples/hto-filter/build_unified.R

suppressMessages({ library(Seurat); library(SeuratObject) })
library(scroll)                      # installed package (NOT load_all: avoids data/ hang)

root   <- "/Users/shihanl1/scroll/examples/hto-filter"
raw    <- readRDS(file.path(root, ".processed.rds"))         # 46,606 base
sub    <- readRDS(file.path(root, "filtered-subsets.rds"))   # 32,611 clean + sub-UMAPs
outdir <- file.path(root, "project-unified")

RES <- sprintf("%.1f", seq(0.1, 1.0, 0.1))

# --- base: percent.ribo + whole-dataset clustering at every resolution ----------
DefaultAssay(raw) <- "RNA"
raw[["percent.ribo"]] <- PercentageFeatureSet(raw, pattern = "^Rp[sl]")   # mouse ribo
raw <- FindNeighbors(raw, reduction = "pca", dims = 1:30, verbose = FALSE)
for (r in RES) {
  raw <- FindClusters(raw, resolution = as.numeric(r), verbose = FALSE)
  raw[[sprintf("res_%s", r)]] <- paste0("c", as.character(raw$seurat_clusters))
}

# --- clean singlets: the whole reprocessed population as one subset view --------
# umap renamed to umap_clean so it doesn't collide with the raw base umap. Its
# clustering (renamed clean_res_*) + biology columns exist only for the 32,611
# singlets, so they are scoped to this view (absent for the raw-only cells).
for (r in RES) sub[[sprintf("clean_res_%s", r)]] <- sub[[sprintf("res_%s", r)]][[1]]
clean_meta <- c("genotype", "treatment", "celltype", "hto",
                "sig_tfh2", "sig_memory", "percent.hsp",
                sprintf("clean_res_%s", RES))
raw <- scroll_add_subset(raw, sub, name = "clean", label = "Clean singlets",
                         embeddings = c(umap_clean = "umap"), meta = clean_meta)

# --- six sub-populations (each a subset of the clean singlets) -------------------
ids <- c(cd101 = "Cd101+", tfh = "Tpex", cx3cr1 = "Cx3cr1+",
         cycling = "Cycling", sgctrl = "sgCtrl", tpex_sgctrl = "Tpex (sgCtrl)")
sub_meta <- list()
for (id in names(ids)) {
  cols <- sprintf("%s_res_%s", id, RES)
  emb  <- paste0("umap_", id)
  raw  <- scroll_add_subset(raw, sub, name = id, label = unname(ids[id]),
                            embeddings = setNames(emb, emb), meta = cols)
  sub_meta[[id]] <- cols
}

# --- genotype-pair views: control + one knockout, re-embedded together ----------
# Each pair (sgCtrl + a KO guide) is re-normalized / re-PCA / re-UMAP so the two
# genotypes share one embedding, with a scoped `<name>_genotype` column (colour
# control vs KO) and a subcluster sweep. Subset of the clean singlets.
pairs <- list(
  sgctrl_cxcr5 = list(label = "sgCtrl + sgCxcr5", genos = c("sgCtrl", "sgCxcr5")),
  sgctrl_icos  = list(label = "sgCtrl + sgIcos",  genos = c("sgCtrl", "sgIcos")))
pair_meta <- list()
for (nm in names(pairs)) {
  pr <- pairs[[nm]]
  p  <- subset(sub, cells = colnames(sub)[sub$genotype %in% pr$genos])
  DefaultAssay(p) <- "RNA"
  p <- p |>
    NormalizeData(verbose = FALSE) |> FindVariableFeatures(verbose = FALSE) |>
    ScaleData(verbose = FALSE) |> RunPCA(npcs = 30, verbose = FALSE) |>
    RunUMAP(dims = 1:30, verbose = FALSE) |> FindNeighbors(dims = 1:30, verbose = FALSE)
  gcol    <- paste0(nm, "_genotype"); p[[gcol]] <- as.character(p$genotype)
  rescols <- sprintf("%s_res_%s", nm, RES)
  for (i in seq_along(RES)) {
    p <- FindClusters(p, resolution = as.numeric(RES[i]), verbose = FALSE)
    p[[rescols[i]]] <- paste0("c", as.character(p$seurat_clusters))
  }
  emb <- paste0("umap_", nm)
  raw <- scroll_add_subset(raw, p, name = nm, label = pr$label,
                           embeddings = setNames("umap", emb), meta = c(gcol, rescols))
  pair_meta[[nm]] <- c(gcol, rescols)
  message("  pair ", nm, ": ", ncol(p), " cells")
}

# --- build ----------------------------------------------------------------------
base_meta <- c("sample", "nCount_RNA", "nFeature_RNA", "nCount_HTO", "nFeature_HTO",
               "percent.mt", "percent.ribo", sprintf("res_%s", RES),
               grep("^hto_Hashtag", colnames(raw[[]]), value = TRUE))
meta_all  <- c(base_meta, clean_meta, unlist(sub_meta, use.names = FALSE),
               unlist(pair_meta, use.names = FALSE))

message("Building unified store: base ", ncol(raw), " cells, ",
        length(SeuratObject::Misc(raw, "scroll_subsets")), " subset views")
scroll_build(
  raw, outdir, assays = c("RNA", "HTO"), counts = TRUE,
  meta_cols = meta_all, overwrite = TRUE)         # subsets read from @misc

# --- config: brand + curated markers --------------------------------------------
cfg <- yaml::read_yaml(file.path(outdir, "config.yaml"))
cfg$brand <- "Cao et al."
cfg$title <- "Cao et al. — all cells + subset views"
cfg$default_embedding <- "umap"
cfg$markers <- c("Cd3d","Cd8a","Il7r","Ccr7","Cd101","Cx3cr1","Tcf7","Mki67")
yaml::write_yaml(cfg, file.path(outdir, "config.yaml"))

m <- scroll_manifest(outdir)
message("DONE: ", outdir)
message("  n_cells=", m$n_cells, "  has_counts=", m$has_counts)
message("  base clustering: ", paste(grep("^res_", names(m$meta), value = TRUE), collapse = ", "))
message("  percent.ribo in meta: ", "percent.ribo" %in% names(m$meta))
message("  subsets: ", paste(names(m$subsets), collapse = ", "))
