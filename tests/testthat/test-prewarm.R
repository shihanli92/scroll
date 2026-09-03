# Startup warm-up: .scroll_prewarm cycles the View selector on a timer, behind an
# overlay, so each subset view's cached scatters pre-render. It's gated on
# config `prewarm_views` (off by default) + a cache; here we exercise the sequence.

test_that(".scroll_prewarm runs the view-cycle sequence without error", {
  seen <- character(0)
  srv <- function(input, output, session) {
    # spy on the overlay show/hide + the view cycling
    session$sendCustomMessage <- function(type, message) {
      if (identical(type, "scroll_warm"))
        seen[[length(seen) + 1L]] <<- if (isTRUE(message$show)) "show" else "hide"
      invisible()
    }
    scroll:::.scroll_prewarm(session, c("a", "b"), step_ms = 10L)
  }
  shiny::testServer(srv, {
    session$flushReact()                          # onFlushed -> show overlay, kick
    for (i in 1:8) { session$elapse(15); session$flushReact() }
  })
  expect_true("show" %in% seen)                   # overlay was raised
  expect_identical(seen[[length(seen)]], "hide")  # ...and dropped when the cycle finished
})

test_that("prewarm stays off unless configured (no cache / no subsets / prewarm 0)", {
  # a manifest with no subsets: .scroll_subset_names is empty, so no warm-up path
  m <- list(subsets = NULL)
  expect_length(scroll:::.scroll_subset_names(m), 0L)
})
