# High-cardinality metadata: the manifest caches `levels` only up to `max_levels`
# and records `n_levels`; the runtime recomputes the level set on demand and gates
# contrasts off the count, so a level-uncached column stays fully usable.

test_that(".scroll_meta_levels returns cached levels, else recomputes from cells", {
  data <- scroll:::.scroll_load(test_project())
  expect_setequal(scroll:::.scroll_meta_levels(data, "celltype"), c("T", "B", "NK"))

  # Simulate an uncached (high-cardinality) column: drop the manifest cache.
  d2 <- data
  d2$manifest$meta$celltype$levels <- NULL
  expect_setequal(scroll:::.scroll_meta_levels(d2, "celltype"),
                  sort(unique(as.character(data$cells$celltype))))
})

test_that(".scroll_n_levels reads n_levels, falls back to cached levels length", {
  m <- list(meta = list(
    a = list(type = "categorical", n_levels = 500L),                 # uncached
    b = list(type = "categorical", levels = as.list(c("x", "y")))))  # old-style manifest
  expect_equal(scroll:::.scroll_n_levels(m, "a"), 500L)
  expect_equal(scroll:::.scroll_n_levels(m, "b"), 2L)
})

test_that(".scroll_has_contrast counts levels even when uncached (and back-compat)", {
  # uncached high-cardinality categorical still offers a contrast
  expect_true(scroll:::.scroll_has_contrast(
    list(meta = list(clone = list(type = "categorical", n_levels = 300L)))))
  # old-style manifest with only cached levels (no n_levels)
  expect_true(scroll:::.scroll_has_contrast(
    list(meta = list(grp = list(type = "categorical", levels = as.list(c("a", "b")))))))
  # a single-level column is no contrast
  expect_false(scroll:::.scroll_has_contrast(
    list(meta = list(one = list(type = "categorical", n_levels = 1L)))))
})

test_that(".scroll_manual_ui degrades to a note above the display cap", {
  ns <- function(x) x
  many <- paste0("lv", seq_len(scroll:::.SCROLL_MANUAL_CAP + 5L))
  ui <- scroll:::.scroll_manual_ui(ns, many)
  html <- as.character(ui)
  expect_match(html, "Manual colours are unavailable")
  expect_false(grepl("scroll-manual-grid", html))       # no per-level pickers rendered

  few <- c("a", "b", "c")
  expect_match(as.character(scroll:::.scroll_manual_ui(ns, few)), "scroll-manual-grid")
})
