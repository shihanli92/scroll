# Guards for the closeread integration details that live validation pinned down:
# the trigger syntax the bridge depends on, and the CSS-inlining workaround for
# the _extensions static-serving gap under server: shiny.

test_that("scaffolded story.qmd uses the closeread trigger syntax the bridge needs", {
  qmd <- readLines(file.path(test_project(), "story.qmd"))
  txt <- paste(qmd, collapse = "\n")
  # one persistent sticky holding the app UI
  expect_match(txt, "#cr-sticky .sticky", fixed = TRUE)
  expect_match(txt, "scroll_app_ui(\"main\"", fixed = TRUE)
  # every section id appears as a data-section span focusing the sticky
  cfg <- scroll_config(test_project())
  for (s in cfg$sections) {
    expect_match(txt, sprintf("data-section=\"%s\"", s$id), fixed = TRUE)
  }
  expect_match(txt, "@cr-sticky", fixed = TRUE)
})

test_that(".scroll_closeread_css returns content when the extension is present", {
  dir <- test_project()
  ext <- file.path(dir, "_extensions", "qmd-lab", "closeread")
  dir.create(ext, recursive = TRUE, showWarnings = FALSE)
  writeLines(".cr-section { display: grid; }", file.path(ext, "closeread.css"))
  css <- scroll:::.scroll_closeread_css(dir)
  expect_match(css, "display: grid", fixed = TRUE)
  # absent extension -> empty string, not an error
  expect_equal(scroll:::.scroll_closeread_css(tempfile()), "")
})
