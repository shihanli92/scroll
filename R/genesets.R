# Gene sets and preranked GSEA.
#
# Gene sets come from two places. Sets saved into a project
# (genesets/sets.parquet: collection, set, gene) are written by
# scroll_add_genesets() -- from MSigDB via msigdbr, and/or from .gmt files -- and need
# nothing at runtime. When msigdbr is installed, the GSEA panel can also fetch the
# curated MSigDB collections live. scroll_gsea() runs fgsea on a DE / Pseudobulk
# result, ranking genes by an unrounded metric (the result's "ranking" attribute).

# ---- sharing DE results across panels ----------------------------------------

# Publish a DE-type result to this session's registry (data$results, from
# .scroll_wire) so the GSEA panel can use it. No-op outside an app session.
.scroll_publish_result <- function(data, key, df, label) {
  if (!is.function(data$results) || is.null(df)) return(invisible(NULL))
  cur <- shiny::isolate(data$results())
  cur[[key]] <- list(df = df, ranking = attr(df, "ranking"), label = label,
                     n_genes = nrow(df), time = Sys.time())
  data$results(cur)
  invisible(NULL)
}

# "celltype: T vs rest (<how>)" for a contrast.
.scroll_contrast_text <- function(group, ident1, ident2, how) {
  side <- function(x) if (length(x)) paste(x, collapse = ", ") else "rest"
  sprintf("%s: %s vs %s (%s)", paste(group, collapse = " | "), side(ident1), side(ident2), how)
}

# ---- MSigDB collections ----------------------------------------------------------

# The curated collections offered by name, mapped to msigdbr's collection /
# subcollection.
.SCROLL_MSIGDB <- list(
  H            = list(label = "Hallmark", collection = "H", sub = NULL),
  REACTOME     = list(label = "Reactome", collection = "C2", sub = "CP:REACTOME"),
  KEGG         = list(label = "KEGG (legacy)", collection = "C2", sub = "CP:KEGG_LEGACY"),
  WIKIPATHWAYS = list(label = "WikiPathways", collection = "C2", sub = "CP:WIKIPATHWAYS"),
  `GO:BP`      = list(label = "GO biological process", collection = "C5", sub = "GO:BP"),
  `GO:MF`      = list(label = "GO molecular function", collection = "C5", sub = "GO:MF"),
  `GO:CC`      = list(label = "GO cellular component", collection = "C5", sub = "GO:CC"),
  IMMUNESIGDB  = list(label = "ImmuneSigDB", collection = "C7", sub = "IMMUNESIGDB"),
  CELLTYPE     = list(label = "Cell type signatures", collection = "C8", sub = NULL))

# Guess the species from gene symbols: mouse symbols are Title-case (Cd8a),
# human ones upper-case (CD8A).
.scroll_guess_species <- function(features) {
  f <- features[grepl("^[A-Za-z][A-Za-z0-9]+$", features)]
  if (!length(f)) return("Homo sapiens")
  mouse <- mean(grepl("^[A-Z][a-z0-9]+$", f) & grepl("[a-z]", f))
  if (mouse > 0.5) "Mus musculus" else "Homo sapiens"
}

.scroll_project_species <- function(data) {
  data$config$species %||% data$manifest$genesets$species %||%
    .scroll_guess_species(.scroll_features_of(data$manifest, data$manifest$default_assay))
}

# One curated MSigDB collection as a named list (set -> genes), via msigdbr.
.scroll_msigdb_sets <- function(name, species) {
  spec <- .SCROLL_MSIGDB[[name]]
  if (is.null(spec)) stop("Unknown MSigDB collection '", name, "'. Choose from: ",
                          paste(names(.SCROLL_MSIGDB), collapse = ", "), call. = FALSE)
  if (!requireNamespace("msigdbr", quietly = TRUE))
    stop("MSigDB gene sets need the 'msigdbr' package: install.packages('msigdbr').",
         call. = FALSE)
  args <- list(species = species, collection = spec$collection)
  if (!is.null(spec$sub)) args$subcollection <- spec$sub
  if (!"collection" %in% names(formals(msigdbr::msigdbr)))     # msigdbr < 10
    names(args) <- sub("^collection$", "category", sub("^subcollection$", "subcategory", names(args)))
  tab <- do.call(msigdbr::msigdbr, args)
  split(as.character(tab$gene_symbol), as.character(tab$gs_name))
}

# Live MSigDB sets, memoized for the process (a msigdbr call takes a second or two).
.scroll_msigdb_cache <- new.env(parent = emptyenv())
.scroll_msigdb_live <- function(name, species) {
  key <- paste(species, name, sep = "|")
  if (is.null(.scroll_msigdb_cache[[key]]))
    .scroll_msigdb_cache[[key]] <- .scroll_msigdb_sets(name, species)
  .scroll_msigdb_cache[[key]]
}

# ---- .gmt files ------------------------------------------------------------------

# Read a .gmt file (name <tab> description <tab> genes...) into a named list.
.scroll_read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  parts <- strsplit(lines, "\t", fixed = TRUE)
  sets <- lapply(parts, function(p) unique(p[-(1:2)][nzchar(p[-(1:2)])]))
  stats::setNames(sets, vapply(parts, `[[`, "", 1L))
}

# ---- saving gene sets into a project -----------------------------------------------

#' Add gene sets to a built scroll project
#'
#' Saves gene-set collections into a project (`genesets/sets.parquet`) so the GSEA
#' panel can use them without any extra package at runtime. Needs no Seurat object,
#' so it also works on a deployed project. Collections already saved under the same
#' name are replaced; others are kept.
#'
#' @param dir A built scroll project directory.
#' @param msigdb Names of curated MSigDB collections to fetch with the msigdbr
#'   package: any of `"H"` (Hallmark), `"REACTOME"`, `"KEGG"`, `"WIKIPATHWAYS"`,
#'   `"GO:BP"`, `"GO:MF"`, `"GO:CC"`, `"IMMUNESIGDB"`, `"CELLTYPE"`.
#' @param gmt Paths to `.gmt` files. Each becomes a collection named after the file
#'   (or after the element name, when `gmt` is a named vector).
#' @param species Species for MSigDB (e.g. `"Homo sapiens"`, `"Mus musculus"`; mouse
#'   uses msigdbr's ortholog-mapped sets). `NULL` guesses from the gene symbols.
#' @param assay Assay whose genes the sets are trimmed to (default: the project's
#'   default assay). Genes absent from it are dropped, and so are sets left empty.
#' @return Invisibly, the names of the collections saved.
#' @seealso [scroll_gsea()]
#' @export
scroll_add_genesets <- function(dir, msigdb = NULL, gmt = NULL, species = NULL, assay = NULL) {
  man_path <- file.path(dir, "manifest.yaml")
  if (!file.exists(man_path)) stop("'", dir, "' is not a built scroll project.", call. = FALSE)
  m <- scroll_manifest(dir)
  assay <- assay %||% m$default_assay
  feats <- .scroll_features_of(m, assay)
  species <- species %||% m$genesets$species %||% .scroll_guess_species(feats)
  new <- list()
  for (nm in msigdb) new[[nm]] <- .scroll_msigdb_sets(nm, species)
  if (length(gmt)) {
    gnames <- if (!is.null(names(gmt)) && all(nzchar(names(gmt)))) names(gmt)
              else tools::file_path_sans_ext(basename(gmt))
    for (i in seq_along(gmt)) new[[gnames[[i]]]] <- .scroll_read_gmt(gmt[[i]])
  }
  if (!length(new)) stop("Name at least one MSigDB collection or .gmt file.", call. = FALSE)
  rows <- lapply(names(new), function(coll) {
    sets <- lapply(new[[coll]], function(g) intersect(g, feats))
    sets <- sets[lengths(sets) > 0]
    if (!length(sets)) return(NULL)
    data.frame(collection = coll, set = rep(names(sets), lengths(sets)),
               gene = unlist(sets, use.names = FALSE), stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows)
  if (is.null(tab)) stop("None of the gene sets share genes with assay '", assay, "'.", call. = FALSE)
  path <- file.path(dir, "genesets", "sets.parquet")
  if (file.exists(path)) {                                 # keep other saved collections
    old <- as.data.frame(arrow::read_parquet(path, mmap = FALSE))
    tab <- rbind(old[!old$collection %in% names(new), , drop = FALSE], tab)
  }
  dir.create(dirname(path), showWarnings = FALSE)
  arrow::write_parquet(tab, path, compression = "zstd")
  counts <- tapply(tab$set, tab$collection, function(s) length(unique(s)))
  m$genesets <- list(species = species, assay = assay,
                     collections = as.list(stats::setNames(as.integer(counts), names(counts))))
  yaml::write_yaml(m, man_path, unicode = TRUE)
  message(sprintf("scroll: saved %s (%s sets) to %s",
                  paste(names(new), collapse = ", "),
                  paste(as.integer(counts[names(new)]), collapse = ", "), path))
  invisible(names(new))
}

# ---- runtime access ------------------------------------------------------------------

# The collections the GSEA panel can offer: saved ones, then (with msigdbr) the
# curated MSigDB ones not saved. A named character vector: label -> id, where the id
# is "saved:<name>" or "msigdb:<name>".
.scroll_geneset_choices <- function(data) {
  saved <- names(data$manifest$genesets$collections)
  lab <- vapply(saved, function(n) .SCROLL_MSIGDB[[n]]$label %||% n, "", USE.NAMES = FALSE)
  out <- stats::setNames(paste0("saved:", saved), lab)
  if (requireNamespace("msigdbr", quietly = TRUE)) {
    live <- setdiff(names(.SCROLL_MSIGDB), saved)
    out <- c(out, stats::setNames(paste0("msigdb:", live),
                                  vapply(live, function(n) paste0(.SCROLL_MSIGDB[[n]]$label, " (MSigDB, live)"), "")))
  }
  out
}

# A collection's sets (named list) from its choice id.
.scroll_geneset_get <- function(data, id) {
  kind <- sub(":.*$", "", id); name <- sub("^[^:]+:", "", id)
  if (identical(kind, "saved")) {
    con <- data$con
    if (is.null(con$genesets)) {
      path <- file.path(data$dir, "genesets", "sets.parquet")
      con$genesets <- as.data.frame(arrow::read_parquet(path, mmap = FALSE))
    }
    t <- con$genesets[con$genesets$collection == name, , drop = FALSE]
    return(split(t$gene, t$set))
  }
  .scroll_msigdb_live(name, .scroll_project_species(data))
}

# ---- GSEA ----------------------------------------------------------------------------

# The rank metrics a result supports (label -> id), best first.
.scroll_rank_choices <- function(ranking) {
  n <- names(ranking)
  ch <- c("limma t" = "t", "avg_log2FC" = "avg_log2FC", "median logFC" = "median_logFC",
          "logFC" = "logFC", "AUC - 0.5" = "auc", "signed -log10 p" = "signed_logp")
  ok <- ch[ch %in% n | (ch == "signed_logp" & "p_val" %in% n)]
  if ("avg_log2FC" %in% n) ok <- ok[ok != "logFC"]      # presto's logFC is a mean-of-logs
  ok
}

# The ranked statistic for GSEA: the chosen metric per gene, with ties broken
# deterministically by the signed -log10 p (a nudge below half the smallest gap
# between distinct values, so the metric's order is never changed).
.scroll_rank_stats <- function(ranking, rank_by) {
  fc <- ranking$avg_log2FC %||% ranking$t %||% ranking$median_logFC %||% ranking$logFC
  slp <- if (!is.null(ranking$p_val) && !is.null(fc))
    sign(fc) * -log10(pmax(ranking$p_val, .Machine$double.xmin)) else NULL
  x <- switch(rank_by,
    auc = ranking$auc - 0.5,
    signed_logp = { if (is.null(slp)) stop("This result has no p-values to rank by.", call. = FALSE); slp },
    ranking[[rank_by]])
  if (is.null(x)) stop("This result has no '", rank_by, "' column to rank by.", call. = FALSE)
  keep <- !is.na(x) & !is.na(ranking$gene) & !duplicated(ranking$gene)
  x <- x[keep]; g <- ranking$gene[keep]
  second <- if (!is.null(slp)) slp[keep] else seq_along(x)
  gaps <- diff(sort(unique(x)))
  eps <- if (length(gaps)) min(gaps) / (2 * length(x) + 2) else 1e-12
  x <- x + eps * rank(second, ties.method = "first") / length(x)
  sort(stats::setNames(x, g), decreasing = TRUE)
}

#' Preranked GSEA on a differential-expression result
#'
#' Ranks genes from a [scroll_de()] or [scroll_pseudobulk_de()] result and runs
#' `fgsea::fgsea()` (multilevel) against a list of gene sets. Uses the result's
#' unrounded metrics when present, and breaks rank ties deterministically.
#'
#' @param de A result from [scroll_de()] or [scroll_pseudobulk_de()].
#' @param gene_sets A named list of gene vectors (e.g. from a `.gmt` file).
#' @param rank_by Ranking metric: `"t"` (Pseudobulk's limma t), `"avg_log2FC"`,
#'   `"median_logFC"`, `"logFC"`, `"auc"` (AUC - 0.5), or `"signed_logp"`
#'   (sign x -log10 p). `NULL` picks the first available of that order.
#' @param min_size,max_size Gene-set size limits, after keeping ranked genes only.
#' @param seed Seed for fgsea's sampling (the session's RNG is restored).
#' @return A data.frame (`pathway`, `size`, `ES`, `NES`, `pval`, `padj`,
#'   `leading_edge`) ordered by `padj`, with attributes `stats` (the ranked
#'   statistic) and `sets` (the gene sets tested).
#' @export
scroll_gsea <- function(de, gene_sets, rank_by = NULL, min_size = 15, max_size = 500, seed = 1) {
  if (!requireNamespace("fgsea", quietly = TRUE))
    stop("GSEA needs the 'fgsea' package (Bioconductor): BiocManager::install('fgsea').",
         call. = FALSE)
  ranking <- attr(de, "ranking") %||% de
  rank_by <- rank_by %||% unname(.scroll_rank_choices(ranking)[1])
  stats <- .scroll_rank_stats(ranking, rank_by)
  sets <- lapply(gene_sets, function(g) intersect(g, names(stats)))
  sets <- sets[lengths(sets) >= min_size & lengths(sets) <= max_size]
  if (!length(sets))
    stop(sprintf("No gene set has %d-%d of the ranked genes; widen the set-size range.",
                 min_size, max_size), call. = FALSE)
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    old <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
  } else on.exit(rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  set.seed(seed)
  res <- NULL                                  # fgsea draws a console progress bar: keep it quiet
  utils::capture.output(res <- fgsea::fgsea(pathways = sets, stats = stats, minSize = min_size,
                                            maxSize = max_size, eps = 0, nproc = 1))
  out <- data.frame(pathway = res$pathway, size = res$size, ES = res$ES, NES = res$NES,
                    pval = res$pval, padj = res$padj,
                    leading_edge = vapply(res$leadingEdge, paste, "", collapse = ","),
                    stringsAsFactors = FALSE)
  out <- out[order(out$padj, -abs(out$NES)), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "stats") <- stats
  attr(out, "sets") <- sets
  attr(out, "rank_by") <- rank_by
  out
}
