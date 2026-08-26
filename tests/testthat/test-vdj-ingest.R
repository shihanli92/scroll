# Generalised VDJ ingestion: presets + auto-detection + the scRepertoire parser let the
# shipped repertoire panels read TCR/BCR data from different upstream tools. Because the
# bake normalizes every source into one rep_cells schema, the key guarantee is that
# scRepertoire / AIRR / Platypus objects produce the SAME normalized table as scroll's own.

test_that("scRepertoire compound columns parse into per-segment / per-chain columns", {
  md <- data.frame(
    CTgene = c("TRAV1.TRAJ2.TRAC_TRBV3.TRBD1.TRBJ4.TRBC2",   # paired
               "TRAV5.TRAJ6.TRAC_NA",                         # alpha only
               "NA_TRBV7.TRBD1.TRBJ8.TRBC1"),                 # beta only
    CTaa     = c("CAVR_CASSL", "CAVX_NA", "NA_CASSZ"),
    CTstrict = c("cl1", "cl2", "cl3"),
    clonalFrequency = c(2, 1, 1),
    stringsAsFactors = FALSE)
  out <- scroll:::.scroll_vdj_parse_screpertoire(md, "TCR")
  expect_equal(out$v_gene_TRB, c("TRBV3", NA, "TRBV7"))
  expect_equal(out$j_gene_TRB, c("TRBJ4", NA, "TRBJ8"))
  expect_equal(out$v_gene_TRA, c("TRAV1", "TRAV5", NA))
  expect_equal(out$j_gene_TRA, c("TRAJ2", "TRAJ6", NA))
  expect_equal(out$cdr3_beta,  c("CASSL", NA, "CASSZ"))
  expect_equal(out$cdr3_alpha, c("CAVR", "CAVX", NA))
  expect_equal(out$.scroll_ct_clone, c("cl1", "cl2", "cl3"))
  expect_equal(out$.scroll_ct_count, c(2, 1, 1))
})

test_that("resolver auto-detects AIRR / Platypus by column presence; explicit maps win", {
  airr <- make_vdj_airr_object(40)
  r <- scroll:::.scroll_vdj_resolve_spec(vdj_spec("TCR", group_col = "celltype"), airr[[]])
  expect_equal(r$spec$source, "airr")
  expect_equal(unname(r$spec$segments[["TRBV"]]), "v_call_VDJ")
  expect_equal(r$spec$clone_col, "clone_id")
  expect_equal(r$spec$count_col, "duplicate_count")

  platy <- make_vdj_platypus_object(40)
  r2 <- scroll:::.scroll_vdj_resolve_spec(vdj_spec("TCR", group_col = "celltype"), platy[[]])
  expect_equal(r2$spec$source, "platypus")
  expect_equal(unname(r2$spec$segments[["TRAV"]]), "VJ_vgene")

  # scroll-convention object under auto keeps the default map, labelled "scroll"
  r3 <- scroll:::.scroll_vdj_resolve_spec(vdj_spec("TCR", group_col = "celltype"),
                                          make_vdj_object(40)[[]])
  expect_equal(r3$spec$source, "scroll")

  # an explicit segments= is never overwritten by a preset
  s <- vdj_spec("TCR", group_col = "celltype", segments = c(TRBV = "v_call_VDJ"))
  r4 <- scroll:::.scroll_vdj_resolve_spec(s, airr[[]])
  expect_true(r4$spec$segments_explicit)
  expect_equal(unname(r4$spec$segments), "v_call_VDJ")
})

test_that("broadened NA sentinel blanks all-missing clone ids", {
  expect_equal(
    scroll:::.scroll_vdj_blank_clone(c("NA|NA|NA", "None", "clone1", "TRBV|NA", "", "cl_2")),
    c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE))
})

test_that("mapped-but-absent segment columns warn instead of dropping silently", {
  obj <- make_vdj_object(40)
  expect_warning(
    suppressMessages(scroll_build(
      obj, file.path(tempdir(), "scroll-vdj-warn"), assays = "RNA", meta_cols = "celltype",
      vdj = vdj_spec("TCR", group_col = "celltype", clone_col = "clonotype",
                     segments = c(TRBV = "v_gene_TRB", TRBJ = "does_not_exist")),
      overwrite = TRUE)),
    "not found in metadata")
})

test_that("scRepertoire / AIRR / Platypus normalize to the same rep_cells as scroll", {
  keep <- c("cell", "clone_id", "group", "clone_count",
            "TRBV", "TRBJ", "TRAV", "TRAJ", "cdr3_beta", "cdr3_alpha")
  norm <- function(dir) {
    x <- as.data.frame(arrow::read_parquet(file.path(dir, "repertoire", "rep_cells.parquet")))
    x <- x[order(x$cell), keep]; rownames(x) <- NULL; x
  }
  r0 <- norm(vdj_test_project())
  expect_equal(norm(vdj_airr_test_project()), r0)
  expect_equal(norm(vdj_platypus_test_project()), r0)
  expect_equal(norm(vdj_screpertoire_test_project()), r0)
})

test_that("the resolved source is recorded in the manifest vdj block", {
  expect_equal(scroll_manifest(vdj_airr_test_project())$vdj$source, "airr")
  expect_equal(scroll_manifest(vdj_platypus_test_project())$vdj$source, "platypus")
  expect_equal(scroll_manifest(vdj_screpertoire_test_project())$vdj$source, "scRepertoire")
})

test_that("BCR projects build and expose heavy/light labels + a computable paired_rate", {
  dir <- vdj_bcr_test_project()
  m <- scroll_manifest(dir)$vdj
  expect_equal(m$chain_type, "BCR")
  expect_true("IGHV" %in% unlist(m$segments))
  expect_setequal(unlist(m$cdr3_chains), c("heavy", "light"))
  div <- as.data.frame(arrow::read_parquet(file.path(dir, "repertoire", "diversity.parquet")))
  expect_true(any(!is.na(div$paired_rate)))            # exactly 2 CDR3 chains -> computable
})

test_that("built-in VDJ panels render on an auto-detected AIRR project", {
  data <- scroll:::.scroll_load(vdj_airr_test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  pan <- function(id) Filter(function(p) identical(p$id, id),
                             scroll:::.scroll_assemble_panels(data$manifest))[[1]]

  shiny::testServer(pan("gene_usage")$server, args = list(data = data), {
    session$setInputs(segment = "TRBV", view = "Frequency", group = "")
    expect_s3_class(plot_r(), "ggplot")
  })
  shiny::testServer(pan("cdr3_length")$server, args = list(data = data), {
    session$setInputs(chain = "beta", style = "Density", colorby = "")
    expect_s3_class(plot_r(), "ggplot")
  })
  shiny::testServer(pan("diversity")$server, args = list(data = data), {
    session$setInputs(view = "Diversity metric", by = "group", metric = "shannon")
    expect_s3_class(plot_r(), "ggplot")
  })
})
