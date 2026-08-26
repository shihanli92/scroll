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
                                      assays = c("RNA", "ADT"), counts = TRUE,
                                      overwrite = TRUE))
      }
      cached <<- dir
    }
    cached
  }
})

# A parent object carrying a reprocessed "T cell" subset: a partial embedding
# `umap_tcell` (covers only T cells) + a subset-scoped metadata column `tsub`,
# attached via scroll_add_subset(). Exercises the subset-view build path.
make_subset_object <- function(n = 120, seed = 1) {
  obj <- make_test_object(n, seed)
  tcells <- colnames(obj)[obj$celltype == "T"]
  sub <- subset(obj, cells = tcells)
  te <- matrix(rnorm(length(tcells) * 2), ncol = 2,
               dimnames = list(tcells, c("UMAP_1", "UMAP_2")))
  sub[["umap"]] <- SeuratObject::CreateDimReducObject(embeddings = te, key = "UMAP_", assay = "RNA")
  set.seed(seed + 7)
  sub$tsub <- factor(sample(c("Tfh", "Tcm", "Tem"), length(tcells), replace = TRUE))
  sub$tscore <- stats::runif(length(tcells))          # a scoped *numeric* column (module-score-like)
  scroll_add_subset(obj, sub, "tcell", embeddings = c(umap_tcell = "umap"),
                    label = "T cells", meta = c("tsub", "tscore"))
}

subset_test_project <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      dir <- file.path(tempdir(), "scroll-subset-proj")
      if (!dir.exists(file.path(dir, "expr")))
        suppressMessages(scroll_build(make_subset_object(), dir, assays = "RNA",
                                      meta_cols = c("condition", "celltype", "tsub", "tscore"),
                                      overwrite = TRUE))
      cached <<- dir
    }
    cached
  }
})

# A single synthetic TCR draw (skewed clone sizes + V/J genes + CDR3 aa + annotations),
# shared by make_vdj_object() and the alternate-source fixtures so they all express the
# SAME data in different upstream conventions — the parity tests then assert that
# auto-detection normalizes each to an identical rep_cells table.
.vdj_draw <- function(n = 120, seed = 3) {
  set.seed(seed)
  nclone <- 25L
  clone <- sample(seq_len(nclone), n, replace = TRUE,
                  prob = (nclone:1) / sum(nclone:1))          # a skewed clone-size dist
  vb <- paste0("TRBV", sample(1:8, nclone, replace = TRUE))[clone]
  jb <- paste0("TRBJ", sample(1:4, nclone, replace = TRUE))[clone]
  va <- paste0("TRAV", sample(1:8, nclone, replace = TRUE))[clone]
  ja <- paste0("TRAJ", sample(1:6, nclone, replace = TRUE))[clone]
  aa <- function(len) vapply(len, function(k)
    paste(sample(strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]], k, replace = TRUE), collapse = ""),
    character(1))
  list(clone = clone, vb = vb, jb = jb, va = va, ja = ja,
       cdr3b = aa(sample(10:18, nclone, replace = TRUE))[clone],
       cdr3a = aa(sample(8:15, nclone, replace = TRUE))[clone],
       clonotype = paste0("clone", clone),
       clonotype_nt = paste0("nt", clone, "_", sample(1:2, n, replace = TRUE)),
       antigen = factor(sample(c("gB", "B8R"), n, replace = TRUE)),
       tissue  = factor(sample(c("skin", "spleen"), n, replace = TRUE)))
}

# A synthetic TCR object: the base object plus per-cell repertoire columns
# (V/J segments, CDR3 aa strings, clone id, antigen/tissue) named per scroll's own
# conventions so a minimal vdj_spec() picks most of them up by default.
make_vdj_object <- function(n = 120, seed = 3) {
  obj <- make_test_object(n, seed); d <- .vdj_draw(n, seed)
  obj$v_gene_TRB <- d$vb; obj$j_gene_TRB <- d$jb
  obj$v_gene_TRA <- d$va; obj$j_gene_TRA <- d$ja
  obj$cdr3_beta  <- d$cdr3b; obj$cdr3_alpha <- d$cdr3a
  obj$clonotype  <- d$clonotype
  # an alternate, finer clonal definition (each clone splits into up to 2) — stands in
  # for a nucleotide-vs-amino-acid clonotype, so the clone-id picker has >1 candidate
  obj$clonotype_nt <- d$clonotype_nt
  obj$antigen    <- d$antigen
  obj$tissue     <- d$tissue
  obj
}

# The same TCR draw expressed in scRepertoire's compound-column convention (CTgene/CTaa/
# CTstrict + a numeric clonalFrequency), Platypus (VDJ_/VJ_ columns), and dandelion/AIRR
# per-cell obs (_VDJ/_VJ, clone_id). Auto-detection must normalize all three identically.
make_vdj_screpertoire_object <- function(n = 120, seed = 3) {
  obj <- make_test_object(n, seed); d <- .vdj_draw(n, seed)
  obj$CTgene   <- paste0(d$va, ".", d$ja, ".TRAC", "_", d$vb, ".TRBD1.", d$jb, ".TRBC2")
  obj$CTaa     <- paste0(d$cdr3a, "_", d$cdr3b)
  obj$CTstrict <- d$clonotype
  obj$clonalFrequency <- as.numeric(stats::ave(d$clone, d$clone, FUN = length))
  obj$antigen  <- d$antigen; obj$tissue <- d$tissue
  obj
}

make_vdj_platypus_object <- function(n = 120, seed = 3) {
  obj <- make_test_object(n, seed); d <- .vdj_draw(n, seed)
  obj$VDJ_vgene <- d$vb; obj$VDJ_jgene <- d$jb
  obj$VJ_vgene  <- d$va; obj$VJ_jgene  <- d$ja
  obj$VDJ_cdr3s_aa <- d$cdr3b; obj$VJ_cdr3s_aa <- d$cdr3a
  obj$clonotype_id_10x <- d$clonotype
  obj$antigen <- d$antigen; obj$tissue <- d$tissue
  obj
}

make_vdj_airr_object <- function(n = 120, seed = 3) {
  obj <- make_test_object(n, seed); d <- .vdj_draw(n, seed)
  obj$v_call_VDJ <- d$vb; obj$j_call_VDJ <- d$jb
  obj$v_call_VJ  <- d$va; obj$j_call_VJ  <- d$ja
  obj$junction_aa_VDJ <- d$cdr3b; obj$junction_aa_VJ <- d$cdr3a
  obj$clone_id <- d$clonotype
  obj$duplicate_count <- as.numeric(stats::ave(d$clone, d$clone, FUN = length))
  obj$antigen <- d$antigen; obj$tissue <- d$tissue
  obj
}

# A synthetic BCR object: the TCR draw relabelled as heavy/light IG segments, in scroll's
# BCR convention, so the chain_type = "BCR" path (IGHV/IGLV, heavy/light) is exercised.
make_vdj_bcr_object <- function(n = 120, seed = 3) {
  obj <- make_test_object(n, seed); d <- .vdj_draw(n, seed)
  obj$v_gene_IGH <- sub("^TRB", "IGH", d$vb); obj$j_gene_IGH <- sub("^TRB", "IGH", d$jb)
  obj$v_gene_IGL <- sub("^TRA", "IGL", d$va); obj$j_gene_IGL <- sub("^TRA", "IGL", d$ja)
  obj$cdr3_heavy <- d$cdr3b; obj$cdr3_light <- d$cdr3a
  obj$clonotype  <- d$clonotype
  obj$antigen <- d$antigen; obj$tissue <- d$tissue
  obj
}

vdj_test_project <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      dir <- file.path(tempdir(), "scroll-vdj-proj")
      if (!dir.exists(file.path(dir, "repertoire")))
        suppressMessages(scroll_build(
          make_vdj_object(), dir, assays = "RNA",
          meta_cols = c("condition", "celltype", "antigen", "tissue"),
          vdj = vdj_spec("TCR", group_col = "celltype", clone_col = "clonotype",
                         antigen_col = "antigen", tissue_col = "tissue",
                         cluster_col = "celltype", carry = "clonotype_nt"),
          overwrite = TRUE))
      cached <<- dir
    }
    cached
  }
})

# Build a VDJ project from an alternate-source object via source-auto detection, cached
# per source. `spec_args` overrides go to vdj_spec (e.g. chain_type = "BCR").
.vdj_alt_project <- local({
  cache <- new.env(parent = emptyenv())
  function(tag, obj, spec_args = list()) {
    if (is.null(cache[[tag]])) {
      dir <- file.path(tempdir(), paste0("scroll-vdj-", tag))
      spec <- do.call(vdj_spec, c(list(group_col = "celltype",
        antigen_col = "antigen", tissue_col = "tissue", cluster_col = "celltype"), spec_args))
      if (!dir.exists(file.path(dir, "repertoire")))
        suppressWarnings(suppressMessages(scroll_build(
          obj, dir, assays = "RNA",
          meta_cols = c("condition", "celltype", "antigen", "tissue"),
          vdj = spec, overwrite = TRUE)))
      cache[[tag]] <- dir
    }
    cache[[tag]]
  }
})
vdj_screpertoire_test_project <- function() .vdj_alt_project("screp", make_vdj_screpertoire_object())
vdj_platypus_test_project     <- function() .vdj_alt_project("platypus", make_vdj_platypus_object())
vdj_airr_test_project         <- function() .vdj_alt_project("airr", make_vdj_airr_object())
vdj_bcr_test_project <- function()
  .vdj_alt_project("bcr", make_vdj_bcr_object(), list(chain_type = "BCR", clone_col = "clonotype"))

# A synthetic spatial project. Constructing a real Seurat image class in a test is
# impractical, so we build a normal project with a `spatial` 2-D reduction, then
# bake a small tissue raster via the real .scroll_write_spatial() and patch the
# manifest exactly as scroll_build(spatial=) would (images block + kind:spatial).
make_spatial_object <- function(n = 120, seed = 5) {
  obj <- make_test_object(n, seed)
  set.seed(seed)
  coord <- cbind(runif(n, 5, 55), runif(n, 5, 45))            # inside a 60x50 image
  colnames(coord) <- c("spatial_1", "spatial_2"); rownames(coord) <- colnames(obj)
  obj[["spatial"]] <- suppressWarnings(
    SeuratObject::CreateDimReducObject(embeddings = coord, key = "sp_", assay = "RNA"))
  obj
}

# A synthetic scATAC object: an RNA assay plus a `peaks` assay whose features are
# genomic peak names ("chr-start-end"). SeuratObject-only (data layers set directly).
make_atac_object <- function(n = 120, seed = 9) {
  obj <- make_test_object(n, seed)
  set.seed(seed)
  npk <- 30
  peaks <- sprintf("chr1-%d-%d", (1:npk) * 10000, (1:npk) * 10000 + 800)
  pm <- matrix(rpois(npk * n, 0.5), nrow = npk, dimnames = list(peaks, colnames(obj)))
  obj[["peaks"]] <- SeuratObject::CreateAssay5Object(counts = Matrix::Matrix(pm, sparse = TRUE))
  obj <- SeuratObject::SetAssayData(obj, assay = "peaks", layer = "data",
                                    new.data = Matrix::Matrix(.lognorm(pm, 100), sparse = TRUE))
  obj
}

atac_test_project <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      dir <- file.path(tempdir(), "scroll-atac-proj")
      if (!dir.exists(file.path(dir, "atac"))) {
        suppressMessages(scroll_build(
          make_atac_object(), dir, assays = c("RNA", "peaks"),
          meta_cols = c("condition", "celltype"),
          atac = atac_spec("peaks"), overwrite = TRUE))
        # a real Signac ChromatinAssay annotation would fill nearest_gene; here we
        # patch the baked table to exercise the gene -> peak lookup (3 RNA markers).
        pf <- file.path(dir, "atac", "peaks.parquet")
        pk <- as.data.frame(arrow::read_parquet(pf))
        pk$nearest_gene <- rep(c("CD3D", "CD8A", "MS4A1"), length.out = nrow(pk))
        pk$distance <- rep(c(0, 1200, 5000), length.out = nrow(pk))
        arrow::write_parquet(pk, pf)
      }
      cached <<- dir
    }
    cached
  }
})

# A synthetic imaging-based (Xenium/CosMx-style) object: cell centroids in a FOV,
# no H&E image. Exercises the coordinates-only spatial build path. SeuratObject-only
# (no Seurat / NormalizeData) so the data layer is set directly.
make_xenium_object <- function(n = 200, seed = 7) {
  set.seed(seed)
  genes <- paste0("Gene", 1:40)
  cm <- matrix(rpois(length(genes) * n, 1.2), nrow = length(genes),
               dimnames = list(genes, paste0("cell", 1:n)))
  obj <- SeuratObject::CreateSeuratObject(counts = Matrix::Matrix(cm, sparse = TRUE),
                                          assay = "Xenium")
  obj <- SeuratObject::SetAssayData(obj, assay = "Xenium", layer = "data",
                                    new.data = Matrix::Matrix(.lognorm(cm, 1e4), sparse = TRUE))
  obj$celltype <- factor(sample(c("tumor", "stroma", "immune"), n, replace = TRUE))
  coords <- data.frame(x = runif(n, 0, 5000), y = runif(n, 0, 4000), cell = colnames(obj))
  fov <- SeuratObject::CreateFOV(
    list(centroids = SeuratObject::CreateCentroids(coords)),
    type = "centroids", assay = "Xenium")
  obj[["fov"]] <- fov
  obj
}

spatial_test_project <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      dir <- file.path(tempdir(), "scroll-spatial-proj")
      if (!dir.exists(file.path(dir, "spatial"))) {
        suppressMessages(scroll_build(
          make_spatial_object(), dir, assays = "RNA",
          embeddings = c("umap", "spatial"),
          meta_cols = c("condition", "celltype"), overwrite = TRUE))
        W <- 60L; H <- 50L
        arr <- array(stats::runif(H * W * 3), dim = c(H, W, 3))
        blk <- scroll:::.scroll_write_spatial(
          dir, list(name = "spatial", image = "slice1",
                    raster = grDevices::as.raster(arr), width = W, height = H))
        man <- scroll_manifest(dir)
        man$images <- blk
        man$embeddings$spatial$kind <- "spatial"
        yaml::write_yaml(man, file.path(dir, "manifest.yaml"))
      }
      cached <<- dir
    }
    cached
  }
})

# A two-tissue-map spatial object (spatial + spatial2), to exercise the multi-FOV /
# multi-slide embedding selector. Constructing real Seurat images in a test is
# impractical, so — like spatial_test_project — we build two spatial reductions and
# patch two `images` blocks + kind:spatial marks into the manifest.
make_multifov_object <- function(n = 120, seed = 5) {
  obj <- make_spatial_object(n, seed)                        # carries `spatial`
  set.seed(seed + 1)
  coord2 <- cbind(runif(n, 5, 35), runif(n, 5, 25))
  colnames(coord2) <- c("spatial2_1", "spatial2_2"); rownames(coord2) <- colnames(obj)
  obj[["spatial2"]] <- suppressWarnings(
    SeuratObject::CreateDimReducObject(embeddings = coord2, key = "sp2_", assay = "RNA"))
  obj
}

multifov_test_project <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      dir <- file.path(tempdir(), "scroll-multifov-proj")
      if (!dir.exists(file.path(dir, "spatial"))) {
        suppressMessages(scroll_build(
          make_multifov_object(), dir, assays = "RNA",
          embeddings = c("umap", "spatial", "spatial2"),
          meta_cols = c("condition", "celltype"), overwrite = TRUE))
        W <- 60L; H <- 50L; blk <- NULL
        for (nm in c("spatial", "spatial2")) {
          arr <- array(stats::runif(H * W * 3), dim = c(H, W, 3))
          blk <- c(blk, scroll:::.scroll_write_spatial(
            dir, list(name = nm, image = nm,
                      raster = grDevices::as.raster(arr), width = W, height = H)))
        }
        man <- scroll_manifest(dir)
        man$images <- blk
        man$embeddings$spatial$kind <- "spatial"
        man$embeddings$spatial2$kind <- "spatial"
        yaml::write_yaml(man, file.path(dir, "manifest.yaml"))
      }
      cached <<- dir
    }
    cached
  }
})
