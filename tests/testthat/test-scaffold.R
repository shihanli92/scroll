# Guards for the scaffolded story: per-section stickies (no JS bridge) and the
# CSS-inlining workaround for the _extensions static-serving gap under
# server: shiny.

test_that("scaffolded story.qmd emits a persistent header + one sticky per section", {
  qmd <- readLines(file.path(test_project(), "story.qmd"))
  txt <- paste(qmd, collapse = "\n")
  expect_match(txt, "scroll_app_ui(dir = project_dir)", fixed = TRUE)
  expect_match(txt, "scroll_app_server(dir = project_dir)", fixed = TRUE)

  cfg <- scroll_config(test_project())
  # one cr-section per config section
  n_sections <- length(regmatches(txt, gregexpr("{.cr-section}", txt, fixed = TRUE))[[1]])
  expect_equal(n_sections, length(cfg$sections))
  for (s in cfg$sections) {
    expect_match(txt, sprintf("scroll_sticky_ui(\"%s\"", s$id), fixed = TRUE)
  }
  # the removed shared-sticky bridge must be gone
  expect_false(grepl("data-section", txt))
  expect_false(grepl("active_section", txt))
})

test_that(".scroll_closeread_css inlines the extension stylesheet when present", {
  dir <- test_project()
  ext <- file.path(dir, "_extensions", "qmd-lab", "closeread")
  dir.create(ext, recursive = TRUE, showWarnings = FALSE)
  writeLines(".cr-section { display: grid; }", file.path(ext, "closeread.css"))
  expect_match(scroll:::.scroll_closeread_css(dir), "display: grid", fixed = TRUE)
  expect_equal(scroll:::.scroll_closeread_css(tempfile()), "")
})
