# Startup warm-up: .scroll_prewarm cycles the View selector on a timer, behind an
# overlay, so each subset view's cached scatters pre-render. It's gated on
# config `prewarm_views` (off by default) + a cache; here we exercise the sequence.

test_that(".scroll_prewarm advances on the view echo and hides only at the end", {
  seen <- character(0); texts <- character(0)
  srv <- function(input, output, session) {
    # spy on the overlay show/hide + the per-stage text
    session$sendCustomMessage <- function(type, message) {
      if (identical(type, "scroll_warm")) {
        seen[[length(seen) + 1L]] <<- if (isTRUE(message$show)) "show" else "hide"
        if (!is.null(message$text)) texts[[length(texts) + 1L]] <<- message$text
      }
      invisible()
    }
    scroll:::.scroll_prewarm(input, session, c("a", "b"), labels = c("Donor A", "Donor B"))
  }
  shiny::testServer(srv, {
    session$flushReact()                          # onFlushed -> step 1 (requests view "a")
    # simulate the client echoing each requested view; each echo advances one step
    for (v in c("a", "b", "")) { session$setInputs(scroll_view = v); session$flushReact() }
    session$flushReact()                          # let the finish() onFlushed fire
  })
  expect_true(any(grepl("Donor A", texts)))       # stage message uses the view labels
  expect_true("hide" %in% seen)
  expect_identical(seen[[length(seen)]], "hide")  # hide comes last (after the final view)
})

test_that(".scroll_prewarm watchdog hides if a step's echo never arrives", {
  seen <- character(0)
  srv <- function(input, output, session) {
    session$sendCustomMessage <- function(type, message) {
      if (identical(type, "scroll_warm"))
        seen[[length(seen) + 1L]] <<- if (isTRUE(message$show)) "show" else "hide"
      invisible()
    }
    scroll:::.scroll_prewarm(input, session, c("a", "b"), step_timeout = 0.2)
  }
  suppressWarnings(shiny::testServer(srv, {
    session$flushReact()                          # step 1 requests "a"; no echo follows
    Sys.sleep(0.35)                               # real time (watchdog uses Sys.time)
    session$elapse(1000); session$flushReact()    # watchdog observe fires, sees timeout
    session$flushReact()                          # finish()'s onFlushed sends hide
  }))
  expect_identical(seen[[length(seen)]], "hide")  # aborts + hides, never strands the overlay
})

test_that("prewarm stays off unless configured (no cache / no subsets / prewarm 0)", {
  # a manifest with no subsets: .scroll_subset_names is empty, so no warm-up path
  m <- list(subsets = NULL)
  expect_length(scroll:::.scroll_subset_names(m), 0L)
})

test_that(".scroll_prewarm_on gates on prewarm_views + subsets + cache_plots", {
  data <- scroll:::.scroll_load(subset_test_project())
  on.exit(scroll_disconnect(data$con))
  expect_gt(length(scroll:::.scroll_subset_names(data$manifest)), 0)   # fixture has a subset
  data$config$prewarm_views <- 0L
  expect_false(scroll:::.scroll_prewarm_on(data))     # off by default (0/absent)
  data$config$prewarm_views <- 3L
  expect_true(scroll:::.scroll_prewarm_on(data))      # subsets + prewarm -> on
  data$config$cache_plots <- FALSE
  expect_false(scroll:::.scroll_prewarm_on(data))     # cache disabled -> off
})

test_that(".scroll_page shows the warm-up overlay only when warming", {
  data <- scroll:::.scroll_load(subset_test_project())
  on.exit(scroll_disconnect(data$con))
  panels <- scroll:::.scroll_assemble_panels(data$manifest)
  html_on  <- as.character(scroll:::.scroll_page(data, "t", panels, warming = TRUE))
  html_off <- as.character(scroll:::.scroll_page(data, "t", panels, warming = FALSE))
  expect_true(grepl("scroll-warm-overlay", html_on))
  expect_false(grepl("scroll-warm-overlay", html_off))   # preview path: no stranded overlay
})
