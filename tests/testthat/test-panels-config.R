# The explicit panel override (config.yaml `panels:` / `exclude_panels:`): an
# allowlist that force-shows exactly those ids in order (bypassing the `when` gate),
# and a denylist applied afterwards.

test_that("include shows exactly the listed panels, in order, bypassing the gate", {
  scroll_reset_panels()
  # a single-region-style manifest: DE's contrast gate would normally drop it
  m <- list(scroll_version = "0.0",
            embeddings = list(umap = list(dims = 2)),
            assays = list(RNA = list(features = c("g1", "g2"))),
            meta = list(region = list(type = "categorical", levels = "only")),
            has_counts = FALSE)
  # gated: no de (single-level categorical)
  gated <- vapply(scroll:::.scroll_assemble_panels(m), `[[`, "", "id")
  expect_false("de" %in% gated)
  # include forces exactly these, in this order, gate bypassed
  ids <- vapply(scroll:::.scroll_assemble_panels(m, include = c("de", "dimplot")), `[[`, "", "id")
  expect_equal(ids, c("de", "dimplot"))
})

test_that("include renumbers contiguously by final position", {
  scroll_reset_panels()
  p <- scroll:::.scroll_assemble_panels(NULL, include = c("violin", "dimplot", "de"))
  expect_equal(vapply(p, `[[`, "", "id"), c("violin", "dimplot", "de"))
  expect_equal(vapply(p, `[[`, "", "num"), c("01", "02", "03"))
})

test_that("exclude drops the listed panels from the gated set", {
  scroll_reset_panels()
  full <- scroll_manifest(test_project())            # supports all 8 RNA panels
  base <- vapply(scroll:::.scroll_assemble_panels(full), `[[`, "", "id")
  expect_true(all(c("de", "pseudobulk") %in% base))
  kept <- vapply(scroll:::.scroll_assemble_panels(full, exclude = c("de", "pseudobulk")),
                 `[[`, "", "id")
  expect_false(any(c("de", "pseudobulk") %in% kept))
  expect_true("dimplot" %in% kept)
})

test_that("include and exclude compose (allowlist base, then remove)", {
  scroll_reset_panels()
  ids <- vapply(scroll:::.scroll_assemble_panels(
    NULL, include = c("dimplot", "featureplot", "de"), exclude = "de"), `[[`, "", "id")
  expect_equal(ids, c("dimplot", "featureplot"))
})

test_that("unknown panel ids warn but do not error", {
  scroll_reset_panels()
  expect_warning(scroll:::.scroll_assemble_panels(NULL, include = c("dimplot", "nope")),
                 "unknown panel id")
  expect_warning(scroll:::.scroll_assemble_panels(NULL, exclude = "notapanel"),
                 "unknown panel id")
  # the valid ids still resolve despite the warning
  ids <- suppressWarnings(
    vapply(scroll:::.scroll_assemble_panels(NULL, include = c("dimplot", "nope")), `[[`, "", "id"))
  expect_equal(ids, "dimplot")
})

test_that("scroll_build(panels=) round-trips through config.yaml", {
  obj <- make_test_object(60)
  dir <- file.path(tempdir(), "scroll-panels-cfg")
  suppressMessages(scroll_build(obj, dir, assays = "RNA",
                                panels = c("dimplot", "featureplot"),
                                exclude_panels = "violin", overwrite = TRUE))
  cfg <- scroll_config(dir)
  expect_equal(unlist(cfg$panels), c("dimplot", "featureplot"))
  expect_equal(unlist(cfg$exclude_panels), "violin")

  # the app assembles from those config keys
  data <- scroll:::.scroll_load(dir)
  on.exit(scroll_disconnect(data$con), add = TRUE)
  ids <- vapply(scroll:::.scroll_assemble_panels(
    data$manifest, include = data$config$panels, exclude = data$config$exclude_panels),
    `[[`, "", "id")
  expect_equal(ids, c("dimplot", "featureplot"))
})

test_that("scroll_build rejects a non-character panels argument", {
  obj <- make_test_object(40)
  expect_error(
    suppressMessages(scroll_build(obj, file.path(tempdir(), "bad-panels"),
                                  assays = "RNA", panels = 1:3, overwrite = TRUE)),
    "character vector")
})
