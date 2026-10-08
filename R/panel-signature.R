# --- Signature panel ----------------------------------------------------------

# Live per-cell gene-signature scoring from a typed/pasted gene list, shown on the
# embedding or as a violin by group. The scoring methods live in R/signature.R:
# Mean, Scaled (z over the shown cells), AddModuleScore (Seurat, exact on projects
# with build-time gene means), and the rank-based UCell / AUCell-style (projects
# built with within-cell ranks). "Add as column" hands the score -- optionally with
# high/low groups -- to every other panel for this session. An optional list of down
# genes is scored the same way and subtracted (.scroll_signed_score).

.SCROLL_SIG_METHODS <- c("Mean (log-norm)" = "mean", "Scaled (z-score, shown cells)" = "scaled",
                         "AddModuleScore (Seurat)" = "addmodulescore",
                         "UCell (rank-based)" = "ucell", "AUCell-style (rank-based)" = "aucell")

# The methods a project can run: the rank-based ones need within-cell ranks.
.scroll_sig_methods <- function(m) {
  if (any(vapply(m$assays, function(a) isTRUE(a$ranks), logical(1)))) .SCROLL_SIG_METHODS
  else .SCROLL_SIG_METHODS[!.SCROLL_SIG_METHODS %in% c("ucell", "aucell")]
}

signature_ui <- function(id, data) {
  ns <- NS(id)
  m <- data$manifest
  assays <- .scroll_assays_of(m); cats <- .scroll_cat_cols(m); reds <- .scroll_view_embeddings(m, NULL)
  n_feat <- m$assays[[m$default_assay]]$n_features %||% 20000
  when <- function(method) sprintf("input['%s'] == '%s'", ns("method"), method)
  bslib::layout_columns(
    col_widths = c(3, 9), class = "scroll-panel",
    div(
      class = "scroll-controls",
      .scroll_group("Signature",
        selectizeInput(ns("sig"), "Up genes", choices = NULL, multiple = TRUE,
                       options = list(placeholder = "Add genes, or paste a list...",
                                      maxOptions = 50, plugins = list("remove_button"))),
        .scroll_paste_handler(ns("sig"), ns("sig_paste")),
        # optional: scored with the same method and subtracted from the up score
        selectizeInput(ns("sig_down"), "Down genes (optional)", choices = NULL, multiple = TRUE,
                       options = list(placeholder = "Genes expected to be low",
                                      maxOptions = 50, plugins = list("remove_button"))),
        .scroll_paste_handler(ns("sig_down"), ns("sig_down_paste")),
        selectInput(ns("method"), "Score", .scroll_sig_methods(m)),
        # each method's own parameters, at the defaults of the tool it reproduces
        .scroll_details("Method settings", open = FALSE,
          numericInput(ns("w_neg"), "Down-gene weight", 1, min = 0, step = 0.1),
          .scroll_cond_panel(when("addmodulescore"),
            numericInput(ns("nbin"), "Expression bins", 24, min = 2, step = 1),
            numericInput(ns("ctrl"), "Control genes per gene", 100, min = 1, step = 1),
            numericInput(ns("seed"), "Seed", 1, step = 1)),
          .scroll_cond_panel(when("ucell"),
            numericInput(ns("max_rank"), "Max rank", 1500, min = 10, step = 50)),
          .scroll_cond_panel(when("aucell"),
            numericInput(ns("auc_max_rank"), "Max rank (top genes)", ceiling(0.05 * n_feat),
                         min = 5, step = 10)),
          .scroll_cond_panel(sprintf("['mean','scaled'].indexOf(input['%s']) >= 0", ns("method")),
            helpText("No other settings for this method."))),
        if (length(assays) > 1)
          selectInput(ns("assay"), "Assay", assays, selected = m$default_assay),
        actionButton(ns("compute"), "Calculate", class = "btn-primary", width = "100%"),
        uiOutput(ns("note"))),
      # hand the score to the other panels (this session only)
      uiOutput(ns("add_ui")),
      .scroll_group("View",
        selectInput(ns("view"), "Show as",
                    c("Feature UMAP" = "umap", "Violin by group" = "violin"))),
      .scroll_cond_panel(
        sprintf("input['%s'] == 'umap'", ns("view")),
        .scroll_group("Embedding",
          selectInput(ns("reduction"), "Reduction", reds,
                      selected = .scroll_default(data, "default_embedding", reds[[1]])))),
      .scroll_cond_panel(
        sprintf("input['%s'] == 'violin'", ns("view")),
        .scroll_group("Grouping",
          if (length(cats))
            selectInput(ns("group"), "Group by", stats::setNames(cats, cats), selected = cats[[1]])
          else helpText("No categorical column available for a violin.")))
    ),
    .scroll_plot_area(ns, csv = TRUE)
  )
}

# The look controls, shown in the panel's Style sheet (same input ids). The
# scatter options apply to the Feature UMAP view only.
signature_style_ui <- function(id, data) {
  ns <- NS(id)
  tagList(
    .scroll_cond_panel(
      sprintf("input['%s'] == 'umap'", ns("view")),
      selectInput(ns("palette"), "Palette", .scroll_continuous_palettes, selected = "grey-purple"),
      sliderInput(ns("clip"), "Color quantiles (%)", 0, 100, c(0, 100), 1),
      bslib::input_switch(ns("order"), "High-scoring cells on top", TRUE),
      bslib::input_switch(ns("legend"), "Legend", TRUE),
      bslib::input_switch(ns("raster"), "Rasterize (fast)", TRUE)),
    .scroll_aspect_input(ns))
}

signature_server <- function(id, data, cells_r = reactive(data$cells),
                             view_r = reactive(NULL), theme_r = reactive(NULL), style_r = reactive(NULL),
                             derived_rv = NULL) {
  moduleServer(id, function(input, output, session) {
    m <- data$manifest
    assay <- reactive(input$assay %||% m$default_assay)
    reds0 <- .scroll_view_embeddings(m, NULL)
    red_rv <- reactiveVal(.scroll_default(data, "default_embedding", reds0[[1]]))
    observeEvent(view_r(), {                             # reductions follow the active subset view
      reds <- .scroll_view_embeddings(m, view_r())
      sel <- if (is.null(view_r())) .scroll_default(data, "default_embedding", reds[[1]]) else reds[[1]]
      if (!sel %in% reds) sel <- reds[[1]]
      red_rv(sel); updateSelectInput(session, "reduction", choices = reds, selected = sel)
    }, ignoreNULL = FALSE, priority = 100)
    observeEvent(input$reduction, red_rv(input$reduction), ignoreInit = TRUE)
    .scroll_bind_view_cats(input, session, view_r, m, "group")
    # assay-aware gene list + pasted delimited gene list
    .scroll_bind_gene_box(input, session, m, assay, "sig", paste_id = "sig_paste")
    .scroll_bind_gene_box(input, session, m, assay, "sig_down", paste_id = "sig_down_paste")
    # Score only when Calculate is clicked (a score -- especially AddModuleScore's
    # control set -- shouldn't recompute on every keystroke). The view / embedding /
    # group / palette apply live from the cached score, so switching them never
    # re-queries.
    score_fn <- function(method, genes, a, prm, progress = NULL)
      function(cells) switch(method,
        mean   = .scroll_cell_mean(cells, genes, data, a),
        scaled = .scroll_signature_score(cells, genes, data$queryN(a, genes), "scaled"),
        addmodulescore = .scroll_module_score(cells, genes, data, a, nbin = prm$nbin,
                                              ctrl = prm$ctrl, seed = prm$seed,
                                              progress = progress),
        ucell  = .scroll_ucell_score(cells, genes, data, a, max_rank = prm$max_rank),
        aucell = .scroll_aucell_score(cells, genes, data, a, auc_max_rank = prm$auc_max_rank))
    num <- function(x, default, min = 1) {
      x <- suppressWarnings(as.numeric(x))
      if (length(x) != 1L || is.na(x) || x < min) default else x
    }
    score_r <- eventReactive(input$compute, {
      cells <- cells_r()
      feats <- .scroll_features_of(m, assay())
      genes <- intersect(input$sig, feats)
      down <- intersect(input$sig_down, feats)
      validate(need(length(genes) >= 1, "Add one or more up genes to build a signature."))
      both <- intersect(genes, down)
      validate(need(!length(both), sprintf("%s %s in both the up and the down list.",
        paste(both, collapse = ", "), if (length(both) == 1) "is" else "are")))
      method <- input$method %||% "mean"
      prm <- list(nbin = as.integer(num(input$nbin, 24, 2)), ctrl = as.integer(num(input$ctrl, 100)),
                  seed = as.integer(num(input$seed, 1, -Inf)), max_rank = num(input$max_rank, 1500, 10),
                  auc_max_rank = num(input$auc_max_rank, NULL, 5))
      w_neg <- num(input$w_neg, 1, 0)
      a <- assay()
      # up genes, then (if any) down genes with the same method, combined
      fn <- function(progress = NULL) function(cells) {
        half <- function(lo, hi) if (is.function(progress))
          function(f, d) progress(lo + (hi - lo) * f, d)
        if (!length(down)) return(score_fn(method, genes, a, prm, progress)(cells))
        up <- score_fn(method, genes, a, prm, half(0, 0.5))(cells)
        dn <- score_fn(method, down, a, prm, half(0.5, 1))(cells)
        out <- .scroll_signed_score(as.numeric(up), as.numeric(dn), method, w_neg)
        attr(out, "exact") <- attr(up, "exact")
        out
      }
      score <- tryCatch(
        withProgress(message = "Scoring signature", value = 0,
          fn(function(f, d) setProgress(value = f, detail = d))(cells)),
        error = function(e) e)
      validate(need(!inherits(score, "error"),
                    if (inherits(score, "error")) conditionMessage(score)))
      list(cells = cells, genes = genes, down = down, method = method, fn = fn(),
           exact = attr(score, "exact"),
           values = data.frame(cell = cells$cell, value = as.numeric(score),
                               stringsAsFactors = FALSE))
    })
    output$note <- renderUI({
      req(input$compute > 0)
      d <- tryCatch(score_r(), error = function(e) NULL)
      if (!is.null(d) && isFALSE(d$exact))
        helpText("Approximate: this project predates the build-time gene means, so the",
                 "control genes can differ slightly from Seurat's. Rebuild for an exact match.")
    })
    # ---- Add as column: the score (+ high/low groups) for the other panels ----
    # Session-only (the registry from .scroll_wire). A filter-invariant score is
    # recomputed over every cell so the column is complete; Scaled is relative to the
    # shown cells, so other cells get NA.
    if (is.function(derived_rv)) {
      output$add_ui <- renderUI({
        req(input$compute > 0)
        d <- tryCatch(score_r(), error = function(e) NULL)
        req(d)
        v <- d$values$value; rng <- range(v, na.rm = TRUE)
        step <- if (diff(rng) > 0) signif(diff(rng) / 100, 2) else 0.01
        .scroll_group("Use in other panels",
          textInput(session$ns("col_name"), "Column name", .scroll_sig_col_name(d)),
          checkboxInput(session$ns("col_groups"), "Also add high / low groups", FALSE),
          .scroll_cond_panel(sprintf("input['%s']", session$ns("col_groups")),
            sliderInput(session$ns("col_cut"), "High from (score)", rng[1], rng[2],
                        round(stats::median(v, na.rm = TRUE), 4), step = step)),
          if (identical(d$method, "scaled"))
            helpText("Scaled scores are relative to the shown cells; other cells get NA."),
          actionButton(session$ns("col_add"), "Add as column", class = "btn-sm btn-outline-primary",
                       width = "100%"),
          uiOutput(session$ns("col_list")))
      })
      output$col_list <- renderUI({
        added <- names(derived_rv()$cols)
        if (!length(added)) return(NULL)
        tagList(
          helpText(paste("Added this session:", paste(added, collapse = ", "))),
          selectizeInput(session$ns("col_drop"), NULL, added, multiple = TRUE,
                         options = list(placeholder = "Columns to remove...")),
          actionLink(session$ns("col_remove"), "Remove selected"))
      })
      observeEvent(input$col_add, {
        d <- score_r()
        nm <- .scroll_sig_clean_name(input$col_name)
        grp <- paste0(nm, "_group")
        taken <- intersect(c(nm, if (isTRUE(input$col_groups)) grp), names(data$cells))
        if (!nzchar(nm) || length(taken)) {
          showNotification(if (!nzchar(nm)) "Give the column a name."
                           else paste0("'", taken[[1]], "' is already a dataset column; pick another name."),
                           type = "error")
          return()
        }
        full <- if (identical(d$method, "scaled")) {
          v <- rep(NA_real_, nrow(data$cells))
          v[match(d$cells$.gidx, data$cells$.gidx)] <- d$values$value
          v
        } else if (nrow(d$cells) == nrow(data$cells)) d$values$value
        else withProgress(message = "Scoring every cell", value = 0.5,
                          as.numeric(d$fn(data$cells)))
        cur <- derived_rv(); cols <- cur$cols
        cols[[nm]] <- list(values = full, meta = .scroll_meta_entry(full))
        if (isTRUE(input$col_groups)) {
          cut <- input$col_cut %||% stats::median(full, na.rm = TRUE)
          g <- ifelse(full >= cut, "high", "low")
          cols[[grp]] <- list(values = g, meta = .scroll_meta_entry(g))
        }
        derived_rv(list(cols = cols, version = cur$version + 1L))
        showNotification(paste0("Added ", paste(c(nm, if (isTRUE(input$col_groups)) grp), collapse = " and "),
                                ": pick it in the other panels' column menus."), duration = 4)
      })
      observeEvent(input$col_remove, {
        drop <- input$col_drop
        req(length(drop))
        cur <- derived_rv()
        derived_rv(list(cols = cur$cols[setdiff(names(cur$cols), drop)], version = cur$version + 1L))
      })
    }

    cosmetic_r <- .scroll_cosmetic(reactive(
      list(theme = theme_r(), palette = input$palette, point_size = input$size,
           order = isTRUE(input$order), legend = isTRUE(input$legend),
           clip = input$clip / 100, aspect = input$aspect)))
    build <- function(raster) {
      validate(need(input$compute > 0, "Add a gene signature, then click Calculate."))
      d <- score_r(); st <- cosmetic_r(); st$raster <- raster
      lab <- .scroll_sig_label(d$genes, d$down)
      if (identical(input$view %||% "umap", "violin")) {
        grp <- input$group
        validate(need(!is.null(grp) && nzchar(grp), "Pick a categorical column to group the violin."))
        return(view_violin(d$cells, list(group_by = grp, feature = lab), d$values, st))
      }
      view_feature_plot(d$cells, list(embedding = red_rv(), feature = lab), d$values, st)
    }
    plot_r <- .scroll_lazy_plot(input, function() {
      n <- if (isTRUE(input$compute > 0)) nrow(score_r()$cells) else 0L
      build(.scroll_use_raster(input$raster, n))
    }, style_r)
    export_r <- reactive(.scroll_apply_style(build(FALSE), style_r()))
    output$plot <- renderPlot(plot_r())
    csv_r <- reactive({
      req(input$compute > 0)
      d <- score_r()
      data.frame(cell = d$cells$cell, signature_score = d$values$value, stringsAsFactors = FALSE)
    })
    .scroll_plot_downloads(output, export_r, id, csv_r = csv_r)
  })
}

# The plot label: "Signature (5 genes)", or "Signature (5 up, 3 down)".
.scroll_sig_label <- function(up, down = character(0)) {
  if (length(down)) sprintf("Signature (%d up, %d down)", length(up), length(down))
  else sprintf("Signature (%d gene%s)", length(up), if (length(up) == 1) "" else "s")
}

# A default column name for a score: sig_<first gene>[_<method>].
.scroll_sig_col_name <- function(d) {
  .scroll_sig_clean_name(paste0("sig_", d$genes[[1]], if (!identical(d$method, "mean")) paste0("_", d$method)))
}
.scroll_sig_clean_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_.]+", "_", trimws(as.character(x %||% "")))
  substr(x, 1L, 60L)
}
