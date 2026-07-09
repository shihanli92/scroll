# HTO / QC gating example

A self-contained `scroll` app that live-gates two mouse 10x samples (RNA + 3
cell-hashing antibodies) on RNA QC + **manual** hashtag thresholds, and exports
the barcodes to keep. Built entirely on the public `register_panel()` API — the
gating panel in `panel.R` is a worked example of the extension point.

## Run

1. **Build the project** once (reads `../../data/s{1,2}_*.h5`, ~a few minutes for
   normalize/PCA/UMAP on ~46k cells):

   ```sh
   Rscript examples/hto-filter/build.R
   ```

   This writes `examples/hto-filter/project/` (git-ignored): a two-assay
   (`RNA` + `HTO`) scroll project with per-cell QC (`nCount_RNA`, `nFeature_RNA`,
   `percent.mt`) and CLR hashtag values (`hto_Hashtag_1/2/3`) as metadata, plus a
   `sample` column.

2. **Serve the app:**

   ```sh
   R -e 'shiny::runApp("examples/hto-filter")'
   ```

## What it does

Two sections are inserted after DimPlot:

**RNA QC filter** — `nCount_RNA` / `nFeature_RNA` range sliders and a `percent.mt`
cap, with live histograms (thresholds drawn) and an "N of M pass" readout. It
publishes the QC-passing barcodes for the gating section to use.

**Hashtag gating** — one interactive **lasso** biaxial per hashtag pair (with 3
hashtags: `H1×H2`, `H1×H3`, `H2×H3`). Draw a lasso around a hashtag's
**single-positive** cloud on whichever plot shows it best, pick that hashtag in
*Assign lasso to*, and click *Set gate from lasso*. Repeat for each hashtag.
Classification is then automatic: a cell in **exactly one** singlet gate is a
**Singlet** (assigned that hashtag), in **two or more** a **Doublet**, and in
**none** a **Negative**. You can also gate doublets **by hand** — lasso a doublet
cloud and assign it to **Doublet** (manual doublets are forced to Doublet on top
of the double-gating rule). Choose which classes to **Keep** (singlets by default)
and **Download kept barcodes (CSV)** → `cell, sample, assignment` for kept cells
that also pass the QC filter.

The two sections are independent Shiny modules that cooperate through
`session$userData` (the QC section publishes its pass-set; the gating section
reads it for the export) — no changes to scroll's core.

Hand-drawn gating is deliberate: with only 3 hashtags and a noisy background,
neither automated demultiplexing (HTODemux) nor simple per-axis cutoffs separate
the populations cleanly — a polygon does.

The app-bar **Subset** pill restricts every panel — including this one — to one
`sample`. The built-in DimPlot / FeaturePlot / Violin panels also work (e.g. UMAP
by `sample`, or FeaturePlot of `Hashtag_1` on the `HTO` assay).
