# VDJ (TCR/BCR repertoire) support: the baked store, the manifest block, panel
# gating via the `when` predicate, and each panel's plot rendering.

test_that("vdj_spec captures TCR defaults and requires group_col", {
  s <- vdj_spec("TCR", group_col = "celltype")
  expect_s3_class(s, "scroll_vdj_spec")
  expect_identical(s$chain_type, "TCR")
  expect_true("TRBV" %in% names(s$segments))
  expect_identical(unname(s$cdr3["beta"]), "cdr3_beta")
  expect_error(vdj_spec("TCR"), "group_col")
})

test_that("scroll_build with vdj = bakes the repertoire store + manifest block", {
  dir <- vdj_test_project()
  rep <- file.path(dir, "repertoire")
  expect_true(file.exists(file.path(rep, "rep_cells.parquet")))
  expect_true(file.exists(file.path(rep, "diversity.parquet")))

  rc <- as.data.frame(arrow::read_parquet(file.path(rep, "rep_cells.parquet")))
  expect_true(all(c("clone_id", "group", "clone_count", "TRBV", "TRAV") %in% names(rc)))
  expect_true("cdr3_beta_len" %in% names(rc))
  expect_true(all(rc$clone_count >= 1))

  man <- scroll_manifest(dir)
  expect_false(is.null(man$vdj))
  expect_identical(man$vdj$chain_type, "TCR")
  expect_true(man$vdj$n_clones >= 1)
  expect_true(man$vdj$has_antigen)
})

test_that("VDJ panels surface only for projects with a vdj manifest block", {
  vdj_ids <- c("clone_overview", "gene_usage", "cdr3_length", "diversity")

  vdj_man <- scroll_manifest(vdj_test_project())
  keep <- vapply(scroll:::.scroll_assemble_panels(vdj_man), `[[`, "", "id")
  expect_true(all(vdj_ids %in% keep))

  rna_man <- scroll_manifest(test_project())          # no vdj block
  drop <- vapply(scroll:::.scroll_assemble_panels(rna_man), `[[`, "", "id")
  expect_false(any(vdj_ids %in% drop))
})

test_that("all four VDJ panels assemble for a vdj project", {
  data <- scroll:::.scroll_load(vdj_test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  panels <- Filter(function(p) p$id %in%
                     c("clone_overview", "gene_usage", "cdr3_length", "diversity"),
                   scroll:::.scroll_assemble_panels(data$manifest))
  expect_length(panels, 4)
})

test_that("cdr3_length renders histogram and density styles", {
  data <- scroll:::.scroll_load(vdj_test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  p <- Filter(function(x) identical(x$id, "cdr3_length"),
              scroll:::.scroll_assemble_panels(data$manifest))[[1]]
  shiny::testServer(p$server, args = list(data = data), {
    session$setInputs(chain = "beta", style = "Histogram", colorby = "group")
    expect_error(force(output$plot), NA)
    session$setInputs(style = "Density", chain = "Combined")
    expect_error(force(output$plot), NA)
  })
})

test_that("clone_overview renders both views", {
  data <- scroll:::.scroll_load(vdj_test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  p <- Filter(function(x) identical(x$id, "clone_overview"),
              scroll:::.scroll_assemble_panels(data$manifest))[[1]]
  shiny::testServer(p$server, args = list(data = data), {
    session$setInputs(view = "Rank-abundance", group = "group")
    expect_error(force(output$plot), NA)
    session$setInputs(view = "Expansion composition")
    expect_error(force(output$plot), NA)
  })
})

test_that("diversity metric + gene-usage residual views render", {
  data <- scroll:::.scroll_load(vdj_test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  dv <- Filter(function(x) identical(x$id, "diversity"),
               scroll:::.scroll_assemble_panels(data$manifest))[[1]]
  shiny::testServer(dv$server, args = list(data = data), {
    session$setInputs(view = "Diversity metric", by = "group", metric = "shannon")
    expect_error(force(output$plot), NA)
  })
  gu <- Filter(function(x) identical(x$id, "gene_usage"),
               scroll:::.scroll_assemble_panels(data$manifest))[[1]]
  shiny::testServer(gu$server, args = list(data = data), {
    session$setInputs(segment = "TRBV", view = "Frequency")
    expect_error(force(output$plot), NA)
    session$setInputs(view = "Chi-square residuals")
    expect_error(force(output$plot), NA)
  })
})
