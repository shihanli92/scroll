test_that(".scroll_check_zstd passes when arrow has the zstd codec", {
  skip_if_not(isTRUE(arrow::arrow_info()$capabilities[["zstd"]]),
              "this arrow build lacks zstd")
  expect_true(scroll:::.scroll_check_zstd())
})

test_that(".scroll_check_zstd stops with an actionable message when zstd is missing", {
  testthat::with_mocked_bindings(
    arrow_info = function() list(capabilities = c(zstd = FALSE)),
    {
      expect_error(scroll:::.scroll_check_zstd(), "zstd", ignore.case = TRUE)
      expect_error(scroll:::.scroll_check_zstd(), "install.packages")
    },
    .package = "arrow"
  )
})

test_that("scroll_connect surfaces the zstd guard before touching the store", {
  dir <- test_project()   # built with the real (zstd-capable) arrow
  testthat::with_mocked_bindings(
    arrow_info = function() list(capabilities = c(zstd = FALSE)),
    expect_error(scroll_connect(dir), "zstd", ignore.case = TRUE),
    .package = "arrow"
  )
})
