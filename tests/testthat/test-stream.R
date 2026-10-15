test_that("scroll_build_stream assembles two sources into one store", {
  dir <- file.path(tempdir(), "scroll-stream-ok")
  a <- make_test_object(n = 40, seed = 1)
  b <- make_test_object(n = 30, seed = 2)
  colnames(a) <- paste0("a_", colnames(a))   # disjoint barcodes across sources
  colnames(b) <- paste0("b_", colnames(b))
  reader <- function(s) if (s == "A") a else b

  suppressMessages(scroll_build_stream(
    dir, sources = c("A", "B"), reader = reader, assays = "RNA",
    embeddings = "umap", meta_cols = c("condition", "celltype"), overwrite = TRUE))

  man <- scroll_manifest(dir)
  expect_equal(man$n_cells, ncol(a) + ncol(b))

  con <- scroll_connect(dir); on.exit(scroll_disconnect(con))
  q <- scroll_query_feature(con, "RNA", "CD3D")     # spans both sources
  expect_s3_class(q, "data.frame")
  expect_gt(nrow(q), 0)
  # the two sources' parts are compacted into one feature-sorted file per assay
  expect_identical(list.files(file.path(dir, "expr", "RNA")), "part-0.parquet")
})

test_that("scroll_build_stream errors, naming the source, on a missing declared column", {
  dir <- file.path(tempdir(), "scroll-stream-bad")
  a <- make_test_object(n = 40, seed = 1)
  b <- make_test_object(n = 30, seed = 2)
  colnames(a) <- paste0("a_", colnames(a))
  colnames(b) <- paste0("b_", colnames(b))
  b$condition <- NULL                                # source B lacks a declared col
  reader <- function(s) if (s == "A") a else b

  expect_error(
    suppressMessages(scroll_build_stream(
      dir, sources = c("A", "B"), reader = reader, assays = "RNA",
      embeddings = "umap", meta_cols = c("condition", "celltype"), overwrite = TRUE)),
    "'B'.*condition")
})

test_that("scroll_build_stream(counts = TRUE) matches a single build and runs pseudobulk", {
  obj <- make_test_object()                          # 120 cells, RNA with a counts layer
  d1 <- file.path(tempdir(), "scroll-stream-counts-single")
  suppressMessages(scroll_build(obj, d1, assays = "RNA", embeddings = "umap",
    meta_cols = c("condition", "celltype"), quantize = FALSE, counts = TRUE,
    overwrite = TRUE))
  d2 <- file.path(tempdir(), "scroll-stream-counts")
  suppressMessages(scroll_build_stream(d2, list(1:60, 61:120), function(idx) obj[, idx],
    assays = "RNA", embeddings = "umap", meta_cols = c("condition", "celltype"),
    counts = TRUE, id_of = function(x) paste0("h", x[[1]]), overwrite = TRUE,
    verbose = FALSE))

  expect_true(isTRUE(scroll_manifest(d2)$has_counts))
  expect_length(list.files(file.path(d2, "counts", "RNA")), 2L)   # one part per source
  norm <- function(x) { x <- as.data.frame(x); x[order(x$feature, x$cell), ] |> `rownames<-`(NULL) }
  single <- norm(arrow::read_parquet(file.path(d1, "counts", "RNA.parquet")))
  stream <- norm(dplyr::collect(arrow::open_dataset(file.path(d2, "counts", "RNA"))))
  expect_equal(stream, single)                       # same global cell index, same counts

  skip_if_not_installed("edgeR"); skip_if_not_installed("limma")
  run <- function(d) {
    data <- scroll:::.scroll_load(d); on.exit(scroll_disconnect(data$con))
    scroll_pseudobulk_de(data, "RNA", aggregate_cols = "celltype", ident1 = "T",
                         ident2 = "B", replicate_col = "condition", min_cells = 5)
  }
  r1 <- run(d1); r2 <- run(d2)
  expect_equal(r2$gene, r1$gene)
  expect_equal(r2$logFC, r1$logFC, tolerance = 1e-8)
})

test_that("an append run must keep the project's counts setting", {
  obj <- make_test_object()
  dir <- file.path(tempdir(), "scroll-stream-counts-mismatch")
  reader <- function(idx) obj[, idx]
  args <- list(outdir = dir, reader = reader, assays = "RNA", embeddings = "umap",
               meta_cols = c("condition", "celltype"),
               id_of = function(x) paste0("h", x[[1]]), verbose = FALSE)
  suppressMessages(do.call(scroll_build_stream,
    c(args, list(sources = list(1:60), counts = TRUE, overwrite = TRUE))))
  expect_error(
    suppressMessages(do.call(scroll_build_stream,
      c(args, list(sources = list(1:60, 61:120), counts = FALSE)))),
    "built with counts = TRUE")
})

test_that("a streamed build's ranks and gene means match a single build's", {
  obj <- make_test_object()
  single <- tempfile("single"); streamed <- tempfile("streamed")
  suppressMessages(scroll_build(obj, single, assays = "RNA", embeddings = "umap",
                                meta_cols = "celltype", quantize = FALSE, ranks = TRUE,
                                verbose = FALSE))
  halves <- list(obj[, 1:60], obj[, 61:120])
  suppressMessages(scroll_build_stream(streamed, 1:2, reader = function(i) halves[[i]],
                                       assays = "RNA", embeddings = "umap",
                                       meta_cols = "celltype", ranks = TRUE, verbose = FALSE))
  rd <- function(d, f) as.data.frame(arrow::read_parquet(file.path(d, f)))
  a <- rd(single, "expr/RNA/part-0.parquet"); b <- rd(streamed, "expr/RNA/part-0.parquet")
  expect_identical(a$rank2[order(a$feature, a$cell)], b$rank2[order(b$feature, b$cell)])
  expect_equal(rd(streamed, "stats/RNA/genes.parquet")$mean, rd(single, "stats/RNA/genes.parquet")$mean)
  expect_identical(rd(streamed, "stats/RNA/cells.parquet"), rd(single, "stats/RNA/cells.parquet"))
  expect_true(isTRUE(scroll_manifest(streamed)$assays$RNA$ranks))
  # an append must keep the original ranks setting
  expect_error(suppressMessages(scroll_build_stream(streamed, 3, reader = function(i) halves[[1]],
                                                    assays = "RNA", meta_cols = "celltype",
                                                    ranks = FALSE, verbose = FALSE)),
               "ranks = TRUE")
})

stream_objects <- function() {
  objs <- list(a = make_test_object(30, 1), b = make_test_object(30, 2), c = make_test_object(30, 3))
  for (nm in names(objs)) objs[[nm]] <- SeuratObject::RenameCells(objs[[nm]], add.cell.id = nm)
  objs
}
# every stored value sits on its source cell's row of cells.parquet
expect_store_matches <- function(dir, objs) {
  cells <- as.data.frame(arrow::read_parquet(file.path(dir, "cells.parquet"), mmap = FALSE))
  src <- sub("_.*", "", cells$cell)
  con <- scroll_connect(dir); on.exit(scroll_disconnect(con))
  for (g in rownames(objs[[1]])) {
    q <- scroll_query_feature(con, "RNA", g)
    expect_false(anyDuplicated(q$cell) > 0)
    v <- numeric(nrow(cells)); v[q$cell] <- q$value
    truth <- unlist(lapply(unique(src), function(s) {
      o <- objs[[s]]
      if (g %in% rownames(o)) SeuratObject::GetAssayData(o, layer = "data")[g, ] else rep(0, ncol(o))
    }), use.names = FALSE)
    expect_equal(v, truth, tolerance = 1e-5)
  }
}

test_that("a streaming run that fails part-way keeps what it finished, and carries on", {
  dir <- tempfile("stream-fail")
  objs <- stream_objects()
  fail <- TRUE
  reader <- function(s) { if (s == "b" && fail) stop("disk hiccup"); objs[[s]] }
  run <- function() suppressMessages(scroll_build_stream(
    dir, c("a", "b", "c"), reader, meta_cols = "celltype", embeddings = "umap", verbose = FALSE))
  expect_error(run(), "disk hiccup")
  expect_true(file.exists(file.path(dir, ".stream", "a.rds")))        # a is committed
  fail <- FALSE
  run()
  cells <- arrow::read_parquet(file.path(dir, "cells.parquet"))
  expect_identical(as.vector(table(sub("_.*", "", cells$cell))), c(30L, 30L, 30L))
  expect_store_matches(dir, objs)
  expect_identical(readLines(file.path(dir, ".stream_sources.txt")), c("a", "b", "c"))
  expect_false(dir.exists(file.path(dir, ".stream")))
  expect_message(scroll_build_stream(dir, c("a", "b", "c"), reader, meta_cols = "celltype",
                                     embeddings = "umap", verbose = FALSE), "nothing new")
})

test_that("a source whose commit never landed is redone, not duplicated", {
  dir <- tempfile("stream-half")
  objs <- stream_objects()
  reader <- function(s) objs[[s]]
  suppressMessages(scroll_build_stream(dir, "a", reader, meta_cols = "celltype",
                                       embeddings = "umap", verbose = FALSE))
  # an append of b died after its parts were written but before its record was
  # renamed into place: leftovers in expr/ and .stream/, no b.rds
  file.copy(file.path(dir, "expr", "RNA", "part-0.parquet"),
            file.path(dir, "expr", "RNA", "b-part-30.parquet"))
  expect_false(dir.exists(file.path(dir, ".stream")))   # a finished run cleans up
  dir.create(file.path(dir, ".stream"))
  file.copy(file.path(dir, "cells.parquet"), file.path(dir, ".stream", "b.cells.parquet"))
  saveRDS(list(), file.path(dir, ".stream", "b.rds.tmp"))
  suppressMessages(scroll_build_stream(dir, c("a", "b", "c"), reader, meta_cols = "celltype",
                                       embeddings = "umap", verbose = FALSE))
  expect_identical(scroll_manifest(dir)$n_cells, 90L)
  expect_store_matches(dir, objs)
})

test_that("the manifest lists every source's features, not just the last one's", {
  dir <- tempfile("stream-union")
  objs <- stream_objects()
  drop <- rownames(objs$c)[1:3]
  objs$c <- subset(objs$c, features = c(setdiff(rownames(objs$c), drop), rownames(objs$c[["ADT"]])))
  suppressMessages(scroll_build_stream(dir, c("a", "b", "c"), function(s) objs[[s]],
                                       meta_cols = "celltype", embeddings = "umap", verbose = FALSE))
  man <- scroll_manifest(dir)
  expect_true(all(drop %in% unlist(man$assays$RNA$features)))
  expect_identical(man$assays$RNA$n_features, length(rownames(objs$a)))
})

test_that("a project streamed by an older scroll is adopted on append", {
  dir <- tempfile("stream-legacy")
  objs <- c(stream_objects(), list(d = SeuratObject::RenameCells(make_test_object(20, 4), add.cell.id = "d")))
  reader <- function(s) objs[[s]]
  suppressMessages(scroll_build_stream(dir, c("a", "b"), reader, meta_cols = "celltype",
                                       embeddings = "umap", verbose = FALSE))
  unlink(file.path(dir, ".stream"), recursive = TRUE)   # what scroll <= 0.3.15 left
  # ... including an old interrupted run: c's part moved in and listed, never finished
  file.copy(file.path(dir, "expr", "RNA", "part-0.parquet"),
            file.path(dir, "expr", "RNA", "c-part-60.parquet"))
  cat("c\n", file = file.path(dir, ".stream_sources.txt"), append = TRUE)
  suppressMessages(scroll_build_stream(dir, c("a", "b", "c", "d"), reader, meta_cols = "celltype",
                                       embeddings = "umap", verbose = FALSE))
  expect_identical(scroll_manifest(dir)$n_cells, 110L)
  expect_store_matches(dir, objs)
})

test_that("compaction is batched and an interrupted one is finished or rolled back", {
  o <- tempfile("compact"); a <- file.path(o, "expr", "RNA"); dir.create(a, recursive = TRUE)
  p1 <- data.frame(feature = c("B", "a", "b"), cell = c(3L, 1L, 2L), value = c(1, 2, 3))
  p2 <- data.frame(feature = c("a", "c"), cell = 4:5, value = c(4, 5))
  arrow::write_parquet(p1, file.path(a, "x-part-0.parquet"))
  arrow::write_parquet(p2, file.path(a, "y-part-3.parquet"))
  scroll:::.scroll_compact_store(o, quantize = FALSE, batch_rows = 1)   # one feature per batch
  expect_identical(list.files(a, all.files = TRUE, no.. = TRUE), "part-0.parquet")
  out <- as.data.frame(arrow::read_parquet(file.path(a, "part-0.parquet")))
  all <- rbind(p1, p2)
  expect_identical(nrow(out), nrow(all))
  expect_equal(out[order(out$feature, out$cell), c("feature", "cell", "value")],
               all[order(all$feature, all$cell), ], ignore_attr = TRUE)

  # a complete temp file with its marker is swapped in; one without is discarded
  arrow::write_parquet(p2, file.path(a, ".compact.parquet.tmp"))
  file.create(file.path(a, ".compact-done"))
  scroll:::.scroll_compact_recover(a)
  expect_identical(nrow(arrow::read_parquet(file.path(a, "part-0.parquet"))), 2L)
  arrow::write_parquet(p1, file.path(a, ".compact.parquet.tmp"))
  scroll:::.scroll_compact_recover(a)
  expect_identical(list.files(a, all.files = TRUE, no.. = TRUE), "part-0.parquet")
  expect_identical(nrow(arrow::read_parquet(file.path(a, "part-0.parquet"))), 2L)
})
