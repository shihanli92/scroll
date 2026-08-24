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

test_that("VDJ panels re-group by any baked categorical and expose CSV", {
  data <- scroll:::.scroll_load(vdj_test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  panel <- function(id) Filter(function(x) identical(x$id, id),
                               scroll:::.scroll_assemble_panels(data$manifest))[[1]]

  # clone_overview rank-abundance facets by the chosen group (here antigen, not group)
  shiny::testServer(panel("clone_overview")$server, args = list(data = data), {
    session$setInputs(view = "Rank-abundance", group = "antigen")
    tryCatch(force(output$plot), error = function(e) NULL)
    expect_false(is.null(output$csv))
  })
  # gene_usage frequency re-groups by antigen; the group control is Frequency-only
  gu <- panel("gene_usage")
  expect_true(any(vapply(environment(gu$ui)$controls, function(c)
    identical(c$id, "group") && !is.null(c$visible_when), logical(1))))
  shiny::testServer(gu$server, args = list(data = data), {
    session$setInputs(segment = "TRBV", view = "Frequency", group = "antigen")
    tryCatch(force(output$plot), error = function(e) NULL)
    expect_false(is.null(output$csv))
  })
})

test_that("clone-id vs group-by columns are split by cardinality (no high-card group-by)", {
  data <- scroll:::.scroll_load(vdj_test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  # clone ids = baked clone_id + the high-cardinality alternate; grouping = the rest
  cc <- scroll:::.scroll_vdj_clone_cols(data)
  sc <- scroll:::.scroll_vdj_split_cols(data)
  expect_true(all(c("clone_id", "clonotype_nt") %in% cc))
  expect_false(any(c("antigen", "tissue") %in% cc))
  expect_true(all(c("antigen", "tissue") %in% sc))
  expect_false("clonotype_nt" %in% sc)                 # high-card clone col NOT a group axis
  expect_true(length(intersect(cc, sc)) == 0)          # the two sets are disjoint
})

test_that("clone_overview / diversity / cdr3_length all take a Clone ID column", {
  data <- scroll:::.scroll_load(vdj_test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  panel <- function(id) Filter(function(x) identical(x$id, id),
                               scroll:::.scroll_assemble_panels(data$manifest))[[1]]

  shiny::testServer(panel("clone_overview")$server, args = list(data = data), {
    session$setInputs(view = "Rank-abundance", group = "group", clone_col = "clonotype_nt")
    expect_error(force(output$plot), NA)                 # finer nt definition renders
    expect_false(is.null(output$csv))
  })
  shiny::testServer(panel("diversity")$server, args = list(data = data), {
    session$setInputs(view = "Diversity metric", by = "group", metric = "shannon",
                      clone_col = "clonotype_nt")
    expect_error(force(output$plot), NA)                 # diversity recomputes under nt clones
  })
  shiny::testServer(panel("cdr3_length")$server, args = list(data = data), {
    session$setInputs(chain = "beta", style = "Density", colorby = "group",
                      clone_col = "clonotype_nt")
    expect_error(force(output$plot), NA)                 # dedup keys off the chosen clone id
  })
})

test_that("VDJ colour picker resolves palette vs manual, and group filter spans panels", {
  data <- scroll:::.scroll_load(vdj_test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  # .scroll_vdj_colours: a named palette by default; manual col_i inputs when "Manual"
  lv <- c("B8R", "gB")
  pal <- scroll:::.scroll_vdj_colours(list(palette = "Set2"), lv)
  expect_setequal(names(pal), lv)
  man <- scroll:::.scroll_vdj_colours(list(palette = "Manual", col_1 = "#111111",
                                           col_2 = "#222222"), lv)
  expect_equal(unname(man[lv]), c("#111111", "#222222"))     # col_i keyed to sorted levels

  panel <- function(id) Filter(function(x) identical(x$id, id),
                               scroll:::.scroll_assemble_panels(data$manifest))[[1]]
  # gene_usage: group filter + Manual palette renders the per-group pickers
  shiny::testServer(panel("gene_usage")$server, args = list(data = data), {
    session$setInputs(segment = "TRBV", view = "Frequency", group = "antigen",
                      group_levels = "gB", palette = "Manual")
    tryCatch(force(output$plot), error = function(e) NULL)
    if (requireNamespace("colourpicker", quietly = TRUE))
      expect_false(is.null(output$palette_manual))            # dynamic manual pickers
  })
  shiny::testServer(panel("cdr3_length")$server, args = list(data = data), {
    session$setInputs(chain = "beta", style = "Histogram", colorby = "antigen",
                      clone_col = "clone_id", group_levels = "gB", palette = "Set2")
    expect_error(force(output$plot), NA)                       # filter + palette
  })
  shiny::testServer(panel("diversity")$server, args = list(data = data), {
    session$setInputs(view = "Diversity metric", by = "antigen", metric = "shannon",
                      clone_col = "clone_id", group_levels = "gB", palette = "Okabe-Ito")
    expect_error(force(output$plot), NA)                       # filter + palette
  })
})

test_that("clone_overview group filter restricts to the selected levels", {
  data <- scroll:::.scroll_load(vdj_test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  lv <- scroll:::.scroll_vdj_col_levels(list(group = "antigen"), data, "group")
  expect_setequal(lv, c("gB", "B8R"))                    # levels of the selected group col

  p <- Filter(function(x) identical(x$id, "clone_overview"),
              scroll:::.scroll_assemble_panels(data$manifest))[[1]]
  shiny::testServer(p$server, args = list(data = data), {
    session$setInputs(view = "Expansion composition", clone_col = "clone_id",
                      group = "antigen", group_levels = character(0))
    all_src <- attr(plot_r(), "scroll_source")
    expect_setequal(unique(all_src$g), c("gB", "B8R"))   # both groups shown by default
    session$setInputs(group_levels = "gB")               # filter to one antigen
    one_src <- attr(plot_r(), "scroll_source")
    expect_setequal(unique(one_src$g), "gB")             # only the kept group survives
  })
})

test_that("diversity recomputes at runtime, matching the baked group-level table", {
  dir <- vdj_test_project()
  data <- scroll:::.scroll_load(dir)
  on.exit(scroll_disconnect(data$con), add = TRUE)
  rc <- as.data.frame(arrow::read_parquet(file.path(dir, "repertoire", "rep_cells.parquet")))
  baked <- as.data.frame(arrow::read_parquet(file.path(dir, "repertoire", "diversity.parquet")))
  baked_grp <- baked[baked$level_type == "group", , drop = FALSE]

  live <- scroll:::.scroll_vdj_diversity(rc, "clone_id", "group",
                                         c("cdr3_beta_len", "cdr3_alpha_len"))
  live <- live[order(live$level), ]; baked_grp <- baked_grp[order(baked_grp$level), ]
  expect_equal(live$shannon, baked_grp$shannon, tolerance = 1e-8)
  expect_equal(live$n_clones, baked_grp$n_clones)

  # runtime re-grouping by antigen (not a baked level_type) + paired_rate metric
  dv <- Filter(function(x) identical(x$id, "diversity"),
               scroll:::.scroll_assemble_panels(data$manifest))[[1]]
  shiny::testServer(dv$server, args = list(data = data), {
    session$setInputs(view = "Diversity metric", by = "antigen", metric = "paired_rate")
    tryCatch(force(output$plot), error = function(e) NULL)
    expect_false(is.null(output$csv))
  })
})
