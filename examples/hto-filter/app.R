# QC + hashtag-gating example app, built entirely on the public register_panel()
# extension point. Build the project first:
#
#   Rscript examples/hto-filter/build.R
#
# then serve this directory:
#
#   R -e 'shiny::runApp("examples/hto-filter")'

library(scroll)
source("panel.R")                 # qc_filter_* / hto_gate_*

scroll_reset_panels()             # start from a clean registry
register_panel(
  "qc", qc_filter_ui, qc_filter_server,
  label = "QC filter", title = "RNA QC filter",
  desc  = "Threshold nCount / nFeature / percent.mt; publishes the QC-passing cells.",
  after = "dimplot")
register_panel(
  "hto", hto_gate_ui, hto_gate_server,
  label = "HTO gate", title = "Hashtag gating",
  desc  = paste("Draw a lasso around a population on the biaxial hashtag plot,",
                "label it, and export the kept barcodes (QC-passing) as CSV."),
  after = "qc")

scroll_app("project")
