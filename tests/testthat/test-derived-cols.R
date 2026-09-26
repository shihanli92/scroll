# Session columns: a Signature score added with "Add as column" becomes a column
# the other panels can use, for this session only (R/panel-helpers.R, R/app.R).

test_that("Add as column puts the score (and high/low groups) in every panel's reach", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  panels <- Filter(function(p) p$id %in% c("dimplot", "featureplot", "signature", "de"),
                   scroll:::.scroll_builtin_panels())
  n <- nrow(data$cells)
  shiny::testServer(function(input, output, session) {
    reg <- scroll:::.scroll_wire(input, output, session, data, panels)
  }, {
    session$setInputs(`signature-sig` = c("CD3D", "CD8A", "NKG7"), `signature-method` = "mean",
                      `signature-view` = "umap", `signature-compute` = 1)
    session$setInputs(`signature-col_name` = "my sig", `signature-col_groups` = TRUE,
                      `signature-col_cut` = 0.5, `signature-col_add` = 1)
    cols <- reg$derived()$cols
    expect_setequal(names(cols), c("my_sig", "my_sig_group"))   # name sanitized
    expect_length(cols$my_sig$values, n)
    expect_identical(cols$my_sig$meta$type, "numeric")
    expect_setequal(unique(cols$my_sig_group$values), c("high", "low"))
    expect_identical(cols$my_sig_group$values, ifelse(cols$my_sig$values >= 0.5, "high", "low"))

    # the active cells carry the columns; the plot cache key moves on
    expect_true(all(c("my_sig", "my_sig_group") %in% names(reg$cells())))
    k1 <- reg$cache_key()
    session$setInputs(`signature-col_groups` = FALSE, `signature-col_add` = 2)   # re-add
    expect_false(identical(reg$cache_key(), k1))

    # the column helpers list them (via this session's manifest)
    m <- data$manifest; attr(m, "scroll_derived") <- reg$derived
    expect_true("my_sig" %in% scroll:::.scroll_num_cols(m))
    expect_true("my_sig_group" %in% scroll:::.scroll_cat_cols(m))
    expect_true("my_sig" %in% unlist(scroll:::.scroll_colorby_choices(m)))
    d2 <- data; d2$manifest <- m
    expect_setequal(scroll:::.scroll_meta_levels(d2, "my_sig_group"), c("high", "low"))
    # ... while the shared manifest is untouched
    expect_false("my_sig" %in% scroll:::.scroll_num_cols(data$manifest))

    # a dataset column's name is refused
    session$setInputs(`signature-col_name` = "celltype", `signature-col_add` = 3)
    expect_false("celltype" %in% names(reg$derived()$cols))

    # Remove drops it again
    session$setInputs(`signature-col_drop` = c("my_sig", "my_sig_group"), `signature-col_remove` = 1)
    expect_length(reg$derived()$cols, 0)
    expect_false("my_sig" %in% names(reg$cells()))
  })
})

test_that("a filter-invariant score is added for every cell, even from a filtered view", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  panels <- Filter(function(p) p$id == "signature", scroll:::.scroll_builtin_panels())
  shiny::testServer(function(input, output, session) {
    reg <- scroll:::.scroll_wire(input, output, session, data, panels)
  }, {
    session$setInputs(scroll_subset_col = "celltype", scroll_subset_val = "T")
    shown <- nrow(reg$cells())
    expect_lt(shown, nrow(data$cells))
    session$setInputs(`signature-sig` = c("CD3D", "NKG7"), `signature-method` = "mean",
                      `signature-compute` = 1, `signature-col_name` = "s1", `signature-col_add` = 1)
    v <- reg$derived()$cols$s1$values
    expect_false(anyNA(v))                                   # scored over all cells
    full <- scroll:::.scroll_cell_mean(data$cells, c("CD3D", "NKG7"), data, "RNA")
    expect_equal(v, full)
    # Scaled is relative to the shown cells: others are NA
    session$setInputs(`signature-method` = "scaled", `signature-compute` = 2,
                      `signature-col_name` = "s2", `signature-col_add` = 2)
    expect_equal(sum(!is.na(reg$derived()$cols$s2$values)), shown)
  })
})
