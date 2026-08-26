# The Clone-map panel: interactive per-clone table over a greyed embedding, selected
# clones highlighted. Gated on manifest$vdj.

test_that("clone_map is gated on a repertoire and sits after the other VDJ panels", {
  data <- scroll:::.scroll_load(vdj_test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  ids <- vapply(scroll:::.scroll_assemble_panels(data$manifest), function(p) p$id, character(1))
  expect_true("clone_map" %in% ids)
  expect_equal(which(ids == "clone_map"), which(ids == "diversity") + 1L)   # after Diversity

  rna <- scroll:::.scroll_load(test_project())                # no vdj -> panel absent
  on.exit(scroll_disconnect(rna$con), add = TRUE)
  rna_ids <- vapply(scroll:::.scroll_assemble_panels(rna$manifest), function(p) p$id, character(1))
  expect_false("clone_map" %in% rna_ids)
})

test_that("clone table aggregates one row per clone, sorted by size", {
  data <- scroll:::.scroll_load(vdj_test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  rc <- scroll:::.scroll_vdj_read(data, "rep_cells.parquet")
  tab <- scroll:::.scroll_clone_table(rc)
  expect_equal(nrow(tab), scroll:::.scroll_ndistinct(rc$clone_id))   # one row per clone
  expect_true(all(c("clone_id", "size") %in% names(tab)))
  expect_true(all(diff(tab$size) <= 0))                              # size descending
  expect_equal(sum(tab$size), sum(!is.na(rc$clone_id)))             # sizes account for all cells
})

test_that("selecting clones highlights exactly their cells on the embedding", {
  data <- scroll:::.scroll_load(vdj_test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  cm <- Filter(function(x) identical(x$id, "clone_map"),
               scroll:::.scroll_assemble_panels(data$manifest))[[1]]
  shiny::testServer(cm$server, args = list(data = data), {
    session$setInputs(minsize = 1)
    tbl <- tbl_r()
    expect_true(nrow(tbl) >= 2)
    session$setInputs(table_rows_selected = c(1L, 2L))         # the two biggest clones
    sel <- selected_r()
    expect_equal(sel, as.character(tbl$clone_id[1:2]))
    p <- plot_r()
    expect_s3_class(p, "ggplot")
    rc <- scroll:::.scroll_vdj_read(data, "rep_cells.parquet")
    n_expected <- sum(!is.na(rc$clone_id) & rc$clone_id %in% sel)
    expect_equal(nrow(ggplot2::layer_data(p, 2L)), n_expected)  # highlight layer = selected cells
    expect_false(is.null(output$csv))
  })
})

test_that("clone-map highlight colours follow the palette / Manual pickers", {
  # palette mode: a named colour per selected clone
  c1 <- scroll:::.scroll_clone_colours(list(palette = "Set1"), c("A", "B"))
  expect_setequal(names(c1), c("A", "B"))
  expect_match(unname(c1[1]), "^#")
  # Manual mode: the per-clone pickers (col_i) win
  c2 <- scroll:::.scroll_clone_colours(
    list(palette = "Manual", col_1 = "#111111", col_2 = "#222222"), c("A", "B"))
  expect_equal(unname(c2[c("A", "B")]), c("#111111", "#222222"))
})

test_that("clone_map plot core greys the base and colours one hue per selected clone", {
  df <- data.frame(cell = paste0("c", 1:6), umap_1 = 1:6, umap_2 = 6:1,
                   clone = c("A", "A", "B", "B", "C", NA), stringsAsFactors = FALSE)
  p <- scroll:::.scroll_clone_map_plot(df, "umap", sel = c("A", "B"))
  expect_s3_class(p, "ggplot")
  hl <- ggplot2::layer_data(p, 2L)
  expect_equal(nrow(hl), 4L)                                   # the A + B cells (not C or NA)
  expect_equal(length(unique(hl$colour)), 2L)                  # one colour per selected clone

  # connect = TRUE adds a per-clone path layer under the points
  pc <- scroll:::.scroll_clone_map_plot(df, "umap", sel = c("A", "B"), connect = TRUE)
  expect_true(any(vapply(pc$layers, function(l)
    inherits(l$geom, "GeomPath") && !inherits(l$geom, "GeomLine"), logical(1))))
})
