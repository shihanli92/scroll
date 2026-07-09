# The rework scaffolds a runnable Shiny app + an optional defaults config,
# not a closeread story.

test_that("scroll_scaffold_app writes app.R and a defaults config.yaml", {
  dir <- test_project()
  expect_true(file.exists(file.path(dir, "app.R")))
  expect_match(paste(readLines(file.path(dir, "app.R")), collapse = "\n"),
               "scroll_app", fixed = TRUE)

  cfg <- scroll_config(dir)
  expect_false(is.null(cfg))
  expect_equal(cfg$default_embedding, "umap")
  expect_equal(cfg$default_assay, "RNA")
  expect_gte(length(cfg$markers), 1)
})

test_that("scroll_config returns NULL when there is no config.yaml", {
  expect_null(scroll_config(tempfile()))
})

test_that(".scroll_pick_features matches markers case-insensitively (mouse)", {
  mouse <- c("Cd3d", "Cd8a", "Ms4a1", "Cd14", "Nkg7", "Gm1992", "Xkr4")
  picked <- scroll:::.scroll_pick_features(mouse)
  expect_true(all(c("Cd3d", "Cd8a", "Ms4a1") %in% picked))  # real markers, actual case
  expect_false("Xkr4" %in% picked)                          # not the junk head() fallback
})
