# A data/counts layer can come back as a DENSE base matrix (some v5 layers, small
# or ADT assays). The build must coerce it to sparse rather than reaching for @x.

test_that(".scroll_as_sparse coerces a dense base matrix to a dgCMatrix", {
  m <- matrix(c(0, 2, 0, 0, 3, 0), nrow = 2,
              dimnames = list(c("g1", "g2"), c("c1", "c2", "c3")))
  s <- scroll:::.scroll_as_sparse(m)
  expect_s4_class(s, "CsparseMatrix")
  expect_equal(length(s@x), 2L)                 # only the nonzeros
  expect_identical(scroll:::.scroll_as_sparse(s), s)   # already sparse -> passthrough
})

test_that("scroll_build handles an assay whose data layer is a dense matrix", {
  skip_if_not_installed("SeuratObject")
  obj <- make_test_object()
  dense <- as.matrix(SeuratObject::GetAssayData(obj, assay = "RNA", layer = "data"))
  expect_true(is.matrix(dense) && !isTRUE(methods::is(dense, "Matrix")))
  obj <- SeuratObject::SetAssayData(obj, assay = "RNA", layer = "data", new.data = dense)

  dir <- file.path(tempdir(), "scroll-dense-layer")
  expect_no_error(
    suppressMessages(scroll_build(obj, dir, quantize = FALSE, overwrite = TRUE)))

  data <- scroll:::.scroll_load(dir); on.exit(scroll_disconnect(data$con))
  truth <- as.numeric(dense["MS4A1", ])
  q <- data$query1("RNA", "MS4A1")
  got <- rep(0, ncol(dense)); got[q$cell] <- q$value
  expect_equal(got, truth, tolerance = 1e-5)
})
