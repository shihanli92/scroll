# The CDR3 logo panel: per-group amino-acid logos via ggseqlogo, gated on vdj + baked
# CDR3 sequences.

test_that("the bake stores CDR3 amino-acid sequences alongside lengths", {
  data <- scroll:::.scroll_load(vdj_test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  rc <- scroll:::.scroll_vdj_read(data, "rep_cells.parquet")
  chains <- unlist(data$manifest$vdj$cdr3_chains)
  expect_true(all(paste0("cdr3_", chains) %in% names(rc)))          # aa sequence columns
  expect_true(all(paste0("cdr3_", chains, "_len") %in% names(rc)))  # length columns still there
  # a sequence's length matches its baked *_len
  ok <- !is.na(rc[[paste0("cdr3_", chains[1])]])
  expect_equal(nchar(rc[[paste0("cdr3_", chains[1])]][ok]),
               rc[[paste0("cdr3_", chains[1], "_len")]][ok])
  # CDR3 sequence columns are not offered as grouping / clone-id axes
  expect_false(any(paste0("cdr3_", chains) %in% scroll:::.scroll_vdj_split_cols(data)))
  expect_false(any(paste0("cdr3_", chains) %in% scroll:::.scroll_vdj_clone_cols(data)))
})

test_that("cdr3_logo panel builds per-group logos (and a pooled one), honouring toggles", {
  skip_if_not_installed("ggseqlogo")
  data <- scroll:::.scroll_load(vdj_test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  ids <- vapply(scroll:::.scroll_assemble_panels(data$manifest), function(p) p$id, character(1))
  expect_equal(which(ids == "cdr3_logo"), which(ids == "clone_map") + 1L)   # after Clone map

  cl <- Filter(function(x) identical(x$id, "cdr3_logo"),
               scroll:::.scroll_assemble_panels(data$manifest))[[1]]
  shiny::testServer(cl$server, args = list(data = data), {
    session$setInputs(chain = "beta", group = "group",
                      count = "Clones (dedup)", units = "Bits")
    p <- plot_r()
    expect_s3_class(p, "ggplot")
    src <- attr(p, "scroll_source")
    expect_true(all(c("group", "cdr3") %in% names(src)))
    expect_true(all(nchar(src$cdr3) == nchar(src$cdr3[1])))   # a single length feeds the logo
    # units toggle + cell counting
    session$setInputs(units = "Probability", count = "Cells")
    expect_s3_class(plot_r(), "ggplot")
    # None -> a single pooled logo
    session$setInputs(group = "")
    expect_s3_class(plot_r(), "ggplot")
    # Manual colour scheme renders
    session$setInputs(group = "group", colscheme = "Manual",
                      chem_1 = "#ff0000", chem_2 = "#00ff00", chem_3 = "#000000",
                      chem_4 = "#0000ff", chem_5 = "#123456")
    expect_s3_class(plot_r(), "ggplot")
  })
})

test_that("cdr3_logo colour scheme: named passes through, Manual builds a col scheme", {
  skip_if_not_installed("ggseqlogo")
  expect_equal(scroll:::.scroll_logo_col_scheme(list(colscheme = "taylor")), "taylor")
  cs <- scroll:::.scroll_logo_col_scheme(list(colscheme = "Manual",
    chem_1 = "#ff0000", chem_2 = "#00ff00", chem_3 = "#000000",
    chem_4 = "#0000ff", chem_5 = "#123456"))
  d <- as.data.frame(cs)
  expect_equal(nrow(d), 20L)                    # all 20 amino acids coloured
  expect_equal(d$col[d$letter == "D"], "#ff0000")   # Acidic picker
  expect_equal(d$col[d$letter == "G"], "#123456")   # Polar picker
})
