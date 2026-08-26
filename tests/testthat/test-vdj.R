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
    session$setInputs(view = "Chi-square residuals", group = "group")
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
  # gene_usage re-groups by antigen; the group control is universal (applies to every
  # view now), while count/gene_order are gated to Frequency/Pairing.
  gu <- panel("gene_usage")
  ctl <- environment(gu$ui)$controls
  expect_true(any(vapply(ctl, function(c)
    identical(c$id, "group") && is.null(c$visible_when), logical(1))))   # NOT gated
  expect_true(any(vapply(ctl, function(c)
    identical(c$id, "gene_order") && !is.null(c$visible_when), logical(1))))
  shiny::testServer(gu$server, args = list(data = data), {
    session$setInputs(segment = "TRBV", view = "Frequency", group = "antigen")
    tryCatch(force(output$plot), error = function(e) NULL)
    expect_false(is.null(output$csv))
  })
})

test_that("download-scale slider value is read (and clamped) from the module input", {
  data <- scroll:::.scroll_load(vdj_test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  gu <- Filter(function(x) identical(x$id, "gene_usage"),
               scroll:::.scroll_assemble_panels(data$manifest))[[1]]
  shiny::testServer(gu$server, args = list(data = data), {
    expect_equal(scroll:::.scroll_dl_scale(session), 1)          # absent -> 1
    session$setInputs(dl_scale = 3); expect_equal(scroll:::.scroll_dl_scale(session), 3)
    session$setInputs(dl_scale = 9); expect_equal(scroll:::.scroll_dl_scale(session), 5)  # clamp
    session$setInputs(dl_scale = 0); expect_equal(scroll:::.scroll_dl_scale(session), 1)  # invalid
  })
})

test_that("export scale multiplies the exported figure dimensions (ggsave scale)", {
  png_w <- function(f) { r <- readBin(f, "raw", 24); sum(as.integer(r[17:20]) * 256^(3:0)) }
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(.data$x, .data$y)) +
    ggplot2::geom_point()
  f1 <- tempfile(fileext = ".png"); scroll:::.scroll_write_plot(f1, p, "png", 1)
  f2 <- tempfile(fileext = ".png"); scroll:::.scroll_write_plot(f2, p, "png", 2)
  expect_equal(png_w(f2), png_w(f1) * 2)                  # 2x scale -> 2x pixels wide
})

test_that("VDJ value plots have zero lower-end expansion on the continuous axis", {
  data <- scroll:::.scroll_load(vdj_test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  gu <- Filter(function(x) identical(x$id, "gene_usage"),
               scroll:::.scroll_assemble_panels(data$manifest))[[1]]
  shiny::testServer(gu$server, args = list(data = data), {
    session$setInputs(segment = "TRBV", view = "Frequency", group = "")
    ysc <- plot_r()$scales$get_scales("y")
    expect_false(is.null(ysc))
    expect_equal(ysc$expand[1], 0)                          # lower multiplicative expansion = 0
  })
})

test_that("ungrouped VDJ series colour follows the palette / Manual colour control", {
  data <- scroll:::.scroll_load(vdj_test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  # no group -> the level set collapses to one synthetic "all" level
  expect_equal(scroll:::.scroll_vdj_level_set(list(group = ""), data, "group"), "all")
  # palette mode: the first colour of the chosen palette
  expect_match(scroll:::.scroll_vdj_one_colour(list(group = "", palette = "Set2"), data, "group"),
               "^#")
  # Manual mode: the single picker (col_1) wins
  expect_equal(
    scroll:::.scroll_vdj_one_colour(list(group = "", palette = "Manual", col_1 = "#123456"),
                                    data, "group"),
    "#123456")
})

test_that("clone_overview / cdr3_length / diversity group-by is optional (None default)", {
  data <- scroll:::.scroll_load(vdj_test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  panel <- function(id) Filter(function(x) identical(x$id, id),
                               scroll:::.scroll_assemble_panels(data$manifest))[[1]]

  # clone_overview: no group -> a single un-faceted rank-abundance curve + one comp bar
  shiny::testServer(panel("clone_overview")$server, args = list(data = data), {
    session$setInputs(view = "Rank-abundance", group = "")
    expect_s3_class(plot_r(), "ggplot")
    session$setInputs(view = "Expansion composition", group = "")
    src <- attr(plot_r(), "scroll_source")
    expect_equal(length(unique(src$g)), 1L)                 # one pooled bar
  })
  # cdr3_length: no colour-by -> a single pooled distribution
  shiny::testServer(panel("cdr3_length")$server, args = list(data = data), {
    session$setInputs(chain = "beta", style = "Density", colorby = "")
    expect_s3_class(plot_r(), "ggplot")
    session$setInputs(style = "Histogram", colorby = "")
    expect_error(plot_r(), NA)
  })
  # diversity: no group -> one overall repertoire value
  shiny::testServer(panel("diversity")$server, args = list(data = data), {
    session$setInputs(view = "Diversity metric", by = "", metric = "shannon")
    expect_equal(nrow(attr(plot_r(), "scroll_source")), 1L)
  })
})

test_that("gene_usage group-by is optional: None pools into one ungrouped plot", {
  data <- scroll:::.scroll_load(vdj_test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  gu <- Filter(function(x) identical(x$id, "gene_usage"),
               scroll:::.scroll_assemble_panels(data$manifest))[[1]]
  shiny::testServer(gu$server, args = list(data = data), {
    # Frequency with no group -> a single overall usage barplot (geom_col)
    session$setInputs(segment = "TRBV", view = "Frequency", group = "")
    p <- plot_r()
    expect_s3_class(p, "ggplot")
    expect_true(any(vapply(p$layers, function(l) inherits(l$geom, "GeomCol"), logical(1))))
    # Pairing with no group -> a single un-faceted heatmap (one group value)
    session$setInputs(view = "Pairing", segment = "TRBV", segment_y = "TRBJ", group = "")
    expect_equal(length(unique(attr(plot_r(), "scroll_source")$grp)), 1L)
    # Chi-square still requires a group
    session$setInputs(view = "Chi-square residuals", group = "")
    expect_error(plot_r(), "Group by")
  })
})

test_that("gene_usage: pairing heatmap and runtime-recomputed chi-square group-by", {
  data <- scroll:::.scroll_load(vdj_test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  gu <- Filter(function(x) identical(x$id, "gene_usage"),
               scroll:::.scroll_assemble_panels(data$manifest))[[1]]

  # Pairing: x/y are two different segments; source has the tile grid + within-group freq
  shiny::testServer(gu$server, args = list(data = data), {
    session$setInputs(view = "Pairing", segment = "TRBV", segment_y = "TRBJ",
                      group = "group", count = "Cells", clone_col = "clone_id")
    src <- attr(plot_r(), "scroll_source")
    expect_true(all(c("gx", "gy", "grp", "freq") %in% names(src)))
    expect_true(all(src$freq >= 0 & src$freq <= 1))
    # absent pairings are completed to freq 0 (no NA gaps), and the grid is full:
    # one row per (gx, gy, group)
    expect_false(anyNA(src$freq))
    expect_true(any(src$freq == 0))
    expect_equal(nrow(src),
                 length(unique(src$gx)) * length(unique(src$gy)) * length(unique(src$grp)))
    # colour palette + quantile cut-offs render a fill scale
    session$setInputs(cpalette = "magma", cquant = c(0.05, 0.95))
    p <- plot_r()
    expect_s3_class(p, "ggplot")
    expect_true("fill" %in% p$scales$get_scales("fill")$aesthetics)
    # two different segments required
    session$setInputs(segment_y = "TRBV")
    expect_error(plot_r(), "two different segments")
  })

  # Chi-square recomputes at runtime, so a NON-baked group (antigen) re-groups the heatmap
  shiny::testServer(gu$server, args = list(data = data), {
    session$setInputs(view = "Chi-square residuals", segment = "TRBV", group = "antigen")
    src <- attr(plot_r(), "scroll_source")
    expect_true(all(c("group", "gene", "residual") %in% names(src)))
    expect_setequal(unique(src$group), c("gB", "B8R"))     # grouped by antigen, not celltype
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
    session$setInputs(by = "group", metric = "shannon", features = "clonotype_nt")
    expect_error(force(output$plot), NA)                 # diversity recomputes under nt clones
  })
  shiny::testServer(panel("cdr3_length")$server, args = list(data = data), {
    session$setInputs(chain = "beta", style = "Density", colorby = "group",
                      clone_col = "clonotype_nt")
    expect_error(force(output$plot), NA)                 # dedup keys off the chosen clone id
  })
})

test_that("diversity: multi-select of clone ids + genes facets one panel per feature", {
  data <- scroll:::.scroll_load(vdj_test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  dv <- Filter(function(x) identical(x$id, "diversity"),
               scroll:::.scroll_assemble_panels(data$manifest))[[1]]
  shiny::testServer(dv$server, args = list(data = data), {
    # diversity of a clone definition AND a V gene, grouped by donor
    session$setInputs(by = "group", metric = "shannon", features = c("clone_id", "TRBV"))
    src <- attr(plot_r(), "scroll_source")
    expect_true("feature" %in% names(src))
    expect_setequal(as.character(unique(src$feature)), c("clone_id", "TRBV"))
    # single feature -> no facet dimension (one feature)
    session$setInputs(features = "clone_id")
    expect_equal(length(unique(attr(plot_r(), "scroll_source")$feature)), 1L)
    # Combine -> the selected features become one composite feature (V-J pairing)
    session$setInputs(features = c("TRBV", "TRBJ"), combine = "Combine")
    src <- attr(plot_r(), "scroll_source")
    expect_equal(unique(as.character(src$feature)), "TRBV+TRBJ")
  })
})

test_that("gene-segment names sort in genomic (natural) order, not alphabetical", {
  x <- c("TRBV10-1", "TRBV2", "TRBV1", "TRBV11-2", "TRBV9", "TRBV11-1", "TRBV20-1")
  expect_equal(scroll:::.scroll_gene_natural_levels(x),
               c("TRBV1", "TRBV2", "TRBV9", "TRBV10-1", "TRBV11-1", "TRBV11-2", "TRBV20-1"))
  # dedups and is order-independent of input
  expect_equal(scroll:::.scroll_gene_natural_levels(c("TRAV13-1", "TRAV13-1", "TRAV3")),
               c("TRAV3", "TRAV13-1"))
})

test_that("gene_usage frequency can deduplicate expanded clones (count clones, not cells)", {
  data <- scroll:::.scroll_load(vdj_test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  panel <- function(id) Filter(function(x) identical(x$id, id),
                               scroll:::.scroll_assemble_panels(data$manifest))[[1]]
  gu <- panel("gene_usage")

  # the Count + Clone-ID controls exist and are Frequency-view only
  ctl <- environment(gu$ui)$controls
  expect_true(any(vapply(ctl, function(c)
    identical(c$id, "count") && !is.null(c$visible_when), logical(1))))

  cells_total <- NULL; clones_total <- NULL
  shiny::testServer(gu$server, args = list(data = data), {
    session$setInputs(segment = "TRBV", view = "Frequency", group = "group",
                      count = "Cells")
    cells_total <<- sum(attr(plot_r(), "scroll_source")$count)
    session$setInputs(count = "Clones (dedup)", clone_col = "clone_id")
    clones_total <<- sum(attr(plot_r(), "scroll_source")$count)
  })
  # cells >= clones, and with expansion present the dedup is a strict reduction
  expect_gt(cells_total, clones_total)
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

  # runtime re-grouping by antigen (not a baked level_type)
  dv <- Filter(function(x) identical(x$id, "diversity"),
               scroll:::.scroll_assemble_panels(data$manifest))[[1]]
  shiny::testServer(dv$server, args = list(data = data), {
    session$setInputs(by = "antigen", metric = "shannon")
    tryCatch(force(output$plot), error = function(e) NULL)
    expect_false(is.null(output$csv))
  })
})

test_that("VDJ panels honour the active cell filter (rep_cells cell key)", {
  data <- scroll:::.scroll_load(vdj_test_project())
  on.exit(scroll_disconnect(data$con))
  rc <- scroll:::.scroll_vdj_read(data, "rep_cells.parquet")
  expect_true("cell" %in% names(rc))            # baked cell key enables scoping
  # `cell` must not surface as a group-by or clone-id option (it's a barcode)
  expect_false("cell" %in% scroll:::.scroll_vdj_split_cols(data))
  expect_false("cell" %in% scroll:::.scroll_vdj_clone_cols(data))

  cp <- Filter(function(x) identical(x$id, "clone_overview"),
               scroll:::.scroll_assemble_panels(data$manifest))[[1]]
  n_clones <- function(sub) {
    out <- NULL
    shiny::testServer(cp$server, args = list(data = data, cells_r = shiny::reactive(sub)), {
      session$setInputs(view = "Rank-abundance", clone_col = "clone_id", group = "group",
                        group_levels = character(0), yscale = "Log",
                        palette = "Tableau 10", aspect = 1)
      out <<- nrow(attr(plot_r(), "scroll_source"))
    })
    out
  }
  half <- data$cells[seq_len(nrow(data$cells) %/% 2), , drop = FALSE]
  expect_gt(n_clones(data$cells), n_clones(half))    # fewer cells -> fewer clones
})
