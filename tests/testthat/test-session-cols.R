# New column from a level mapping (R/session-cols.R): a categorical column made in
# the app for this session, e.g. donor -> genotype.

test_that("the level-mapping parser reads groups, single levels and a default", {
  lv <- c("D1", "D2", "D3", "D4", "D5")
  p <- scroll:::.scroll_parse_level_map("[D1,D2]:2,[D3,D4]:1", lv)
  expect_equal(p$map, c(D1 = "2", D2 = "2", D3 = "1", D4 = "1"))
  expect_true(is.na(p$default))
  # spaces, new lines, quotes, a single level, a default and a stray bracket are fine
  p <- scroll:::.scroll_parse_level_map("[ D1, 'D2' ] : WT\nD3: KO\n*: other]", lv)
  expect_equal(p$map, c(D1 = "WT", D2 = "WT", D3 = "KO"))
  expect_equal(p$default, "other")
  # level names match case-insensitively when not exact
  expect_equal(scroll:::.scroll_parse_level_map("[d1]: a", lv)$map, c(D1 = "a"))
  # readable errors
  expect_error(scroll:::.scroll_parse_level_map("", lv), "Type a mapping")
  expect_error(scroll:::.scroll_parse_level_map("[D1, D9]: 2", lv), "Not a level.*D9")
  expect_error(scroll:::.scroll_parse_level_map("[D1]: a, [D1]: b", lv), "Mapped twice: D1")
  expect_error(scroll:::.scroll_parse_level_map("D1 D2", lv), "Couldn't read")
})

test_that("applying a mapping fills listed levels, the default, and keeps NA as NA", {
  p <- list(map = c(D1 = "2", D2 = "2"), default = NA_character_)
  expect_equal(scroll:::.scroll_apply_level_map(c("D1", "D3", NA, "D2"), p),
               c("2", NA, NA, "2"))
  p$default <- "other"
  expect_equal(scroll:::.scroll_apply_level_map(c("D1", "D3", NA), p), c("2", "other", NA))
})

test_that("the rail's New column adds a session column every panel can use", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  panels <- Filter(function(p) p$id %in% c("dimplot", "de"), scroll:::.scroll_builtin_panels())
  lv <- scroll:::.scroll_meta_levels(data, "celltype")
  expect_gte(length(lv), 2)
  shiny::testServer(function(input, output, session) {
    reg <- scroll:::.scroll_wire(input, output, session, data, panels)
  }, {
    map <- sprintf("[%s]: first, *: rest", lv[[1]])
    session$setInputs(scroll_newcol_from = "celltype", scroll_newcol_map = map,
                      scroll_newcol_name = "is first", scroll_newcol_add = 1)
    cols <- reg$derived()$cols
    expect_named(cols, "is_first")
    v <- cols$is_first$values
    expect_identical(v, ifelse(as.character(data$cells$celltype) == lv[[1]], "first", "rest"))
    expect_identical(cols$is_first$meta$type, "categorical")
    expect_true("is_first" %in% names(reg$cells()))
    m <- data$manifest; attr(m, "scroll_derived") <- reg$derived
    expect_true("is_first" %in% scroll:::.scroll_cat_cols(m))

    # a bad mapping or a dataset column's name adds nothing
    session$setInputs(scroll_newcol_name = "other", scroll_newcol_map = "[nope]: x",
                      scroll_newcol_add = 2)
    session$setInputs(scroll_newcol_name = "celltype", scroll_newcol_map = map,
                      scroll_newcol_add = 3)
    expect_named(reg$derived()$cols, "is_first")

    # remove
    session$setInputs(scroll_newcol_drop = "is_first", scroll_newcol_remove = 1)
    expect_length(reg$derived()$cols, 0)
  })
})

test_that("the New column section is in the rail unless config turns it off", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con), add = TRUE)
  html <- as.character(scroll:::.scroll_controls_ui(data))
  expect_match(html, "scroll_newcol_map", fixed = TRUE)
  data$config$session_columns <- FALSE
  html <- as.character(scroll:::.scroll_controls_ui(data))
  expect_false(grepl("scroll_newcol_map", html, fixed = TRUE))
})

test_that("New column can map a subset view's own column", {
  data <- scroll:::.scroll_load(subset_test_project())
  on.exit(scroll_disconnect(data$con))
  expect_true("tsub" %in% scroll:::.scroll_newcol_cols(data$manifest))
  panels <- Filter(function(p) p$id == "dimplot", scroll:::.scroll_builtin_panels())
  lv <- scroll:::.scroll_meta_levels(data, "tsub")
  shiny::testServer(function(input, output, session) {
    reg <- scroll:::.scroll_wire(input, output, session, data, panels, prewarm = FALSE)
  }, {
    session$setInputs(scroll_newcol_from = "tsub", scroll_newcol_name = "tlab",
                      scroll_newcol_map = sprintf("[%s]: a, *: b", lv[1]), scroll_newcol_add = 1)
    v <- reg$derived()$cols$tlab$values
    src <- as.character(data$cells$tsub)
    expect_true(all(is.na(v[is.na(src)])))              # outside the subset: no label
    expect_true(all(v[!is.na(src)] == ifelse(src[!is.na(src)] == lv[1], "a", "b")))
  })
})
