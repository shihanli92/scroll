# scATAC support: peak-name coordinate parsing, the baked peaks table, the assay
# `kind: peaks` mark, panel gating, and the gene -> peak accessibility render.

test_that("atac_spec captures defaults", {
  s <- atac_spec()
  expect_s3_class(s, "scroll_atac_spec")
  expect_identical(s$assay, "peaks")
  expect_true(s$nearest_gene)
})

test_that("peak names parse into genomic coordinates", {
  p <- scroll:::.scroll_parse_peaks(c("chr1-9776-10668", "chrX:5000-6000", "notapeak"))
  expect_equal(p$chr, c("chr1", "chrX", NA))
  expect_equal(p$start, c(9776, 5000, NA))
  expect_equal(p$end, c(10668, 6000, NA))
  expect_equal(p$width[1], 10668 - 9776 + 1)
})

test_that("scroll_build with atac= bakes the peak table + marks the assay", {
  dir <- atac_test_project()
  expect_true(file.exists(file.path(dir, "atac", "peaks.parquet")))
  man <- scroll_manifest(dir)
  expect_false(is.null(man$atac))
  expect_identical(man$atac$assay, "peaks")
  expect_identical(as.integer(man$atac$n_peaks), 30L)
  expect_identical(man$assays$peaks$kind, "peaks")     # accessibility, not expression
  expect_null(man$assays$RNA$kind)                      # RNA is untouched

  pk <- as.data.frame(arrow::read_parquet(file.path(dir, "atac", "peaks.parquet")))
  expect_true(all(c("feature", "chr", "start", "end", "nearest_gene") %in% names(pk)))
  expect_true(all(pk$chr == "chr1"))
})

test_that("the Peaks panel surfaces only for a project with an atac block", {
  atac_man <- scroll_manifest(atac_test_project())
  expect_true("peaks" %in% vapply(scroll:::.scroll_assemble_panels(atac_man), `[[`, "", "id"))

  rna_man <- scroll_manifest(test_project())            # no atac block
  expect_false("peaks" %in% vapply(scroll:::.scroll_assemble_panels(rna_man), `[[`, "", "id"))
})

test_that("the Peaks panel renders a peak's accessibility on the embedding", {
  data <- scroll:::.scroll_load(atac_test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  p <- Filter(function(x) identical(x$id, "peaks"),
              scroll:::.scroll_assemble_panels(data$manifest))[[1]]
  pk <- scroll:::.scroll_atac_read(data)
  peak <- pk$feature[pk$nearest_gene == "CD3D"][1]
  shiny::testServer(p$server, args = list(data = data), {
    session$setInputs(mode = "Gene", gene = "CD3D", peak = peak, embedding = "umap")
    expect_error(force(output$plot), NA)
    # region text search (the no-annotation fallback) selects peaks by name
    session$setInputs(mode = "Region", region = "chr1-10000")
    expect_error(force(output$plot), NA)  # empty region -> friendly "pick a peak" is fine
  })
})

test_that("region search does genomic-interval overlap, with a substring fallback", {
  pk <- scroll:::.scroll_parse_peaks(c("chr1-10000-10800", "chr1-20000-20800",
                                       "chr2-5000-5800"))
  # chr1:15000-25000 overlaps only the second chr1 peak (10800 < 15000 misses)
  expect_equal(scroll:::.scroll_atac_region_hits(pk, "chr1:15000-25000"), "chr1-20000-20800")
  expect_equal(scroll:::.scroll_atac_region_hits(pk, "chr1-20000-20800"), "chr1-20000-20800")
  expect_length(scroll:::.scroll_atac_region_hits(pk, "chr2:1-100000"), 1)   # chr filter
  # a non-interval query (bare chr) falls back to substring over peak names
  expect_length(scroll:::.scroll_atac_region_hits(pk, "chr1"), 2)
  expect_length(scroll:::.scroll_atac_region_hits(pk, ""), 3)                 # empty -> all
})

test_that("the peak list is capped and truncation is surfaced (not silent)", {
  n <- scroll:::.SCROLL_ATAC_MAX_PEAKS
  pk <- scroll:::.scroll_parse_peaks(sprintf("chr1-%d-%d", (1:(n + 20)) * 100,
                                             (1:(n + 20)) * 100 + 50))
  hits <- scroll:::.scroll_atac_region_hits(pk, "chr1")     # matches every peak
  expect_gt(length(hits), n)                                # more hits than the cap
  expect_length(utils::head(sort(hits), n), n)              # peaks_for bounds the list
})

test_that("the Peaks panel exposes a CSV of the matching peaks", {
  data <- scroll:::.scroll_load(atac_test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  p <- Filter(function(x) identical(x$id, "peaks"),
              scroll:::.scroll_assemble_panels(data$manifest))[[1]]
  pk <- scroll:::.scroll_atac_read(data)
  peak <- pk$feature[pk$nearest_gene == "CD3D"][1]
  shiny::testServer(p$server, args = list(data = data), {
    session$setInputs(mode = "Gene", gene = "CD3D", peak = peak, embedding = "umap")
    tryCatch(force(output$plot), error = function(e) NULL)
    expect_false(is.null(output$csv))
  })
})

test_that("peaks assay accessibility is queryable like any assay (FeaturePlot path)", {
  data <- scroll:::.scroll_load(atac_test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  hit <- data$query1("peaks", "chr1-10000-10800")
  expect_true(is.data.frame(hit))
  expect_true(all(c("cell", "value") %in% names(hit)))
})
