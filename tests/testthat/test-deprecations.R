# Soft deprecation of exports being retired from the public API in 0.3.0: warn
# when a USER calls them, stay silent when the package calls them internally.

reset_seen <- function() {
  seen <- scroll:::.scroll_dep_seen
  rm(list = ls(seen, all.names = TRUE), envir = seen)
}

test_that(".scroll_soft_deprecate warns once for user code, never for package code", {
  reset_seen(); on.exit(reset_seen())
  user <- new.env(parent = globalenv())
  expect_warning(scroll:::.scroll_soft_deprecate("f", "Use g().", user),
                 "`f\\(\\)` is deprecated.*0\\.3\\.0.*Use g\\(\\)")
  expect_no_warning(scroll:::.scroll_soft_deprecate("f", "Use g().", user))  # once/session
  expect_no_warning(scroll:::.scroll_soft_deprecate("h", "x", asNamespace("scroll")))
})

test_that("a deprecated export warns when called from user code", {
  reset_seen(); on.exit(reset_seen())
  con <- scroll_connect(test_project()); on.exit(scroll_disconnect(con), add = TRUE)
  user <- new.env(parent = globalenv())
  assign("con", con, envir = user)
  expect_warning(evalq(scroll::scroll_query_cells(con, "RNA", integer()), user),
                 "scroll_query_cells.*scroll_de")
})

test_that("internal use of a deprecated export (scroll_de -> scroll_query_cells) is silent", {
  skip_if_not_installed("presto")
  reset_seen(); on.exit(reset_seen())
  data <- scroll:::.scroll_load(test_project()); on.exit(scroll_disconnect(data$con), add = TRUE)
  expect_no_warning(scroll_de(data, "RNA", "celltype", "T", "B"))
})
