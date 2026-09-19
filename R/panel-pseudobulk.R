# --- Pseudobulk DE panel ------------------------------------------------------

# Combined-interaction levels observed in `cells` for the given columns.
.scroll_combo_choices <- function(cells, cols) {
  if (!length(cols)) return(character(0))
  combo <- .scroll_combo_levels(cells, cols)
  sort(unique(combo[!is.na(combo)]))
}

# Tidy limma design column names for display: "(Intercept)" -> "int",
# "groupgroup1" -> "group1", "replicate<level>" -> "<level>".
.scroll_tidy_design_cols <- function(nm) {
  nm <- sub("^\\(Intercept\\)$", "int", nm)
  nm <- sub("^groupgroup", "group", nm)
  sub("^replicate", "", nm)
}

# Small HTML table for a (samples x terms) design matrix. `footer` (a numeric
# vector aligned to the columns) is appended as a highlighted "contrast" row;
# `counts` (per-row cell counts) is shown as a trailing "cells" column.
.scroll_matrix_html <- function(D, footer = NULL, counts = NULL) {
  cell_col <- !is.null(counts)
  ccell <- function(x) if (cell_col) tags$td(class = "scroll-designmat-n", x)
  foot <- if (!is.null(footer))
    tags$tr(class = "scroll-designmat-contrast",
            tags$th("contrast"),
            lapply(footer, function(v) tags$td(format(v))), ccell(""))
  tags$div(class = "scroll-designmat-wrap",
    tags$table(class = "scroll-designmat",
      tags$tr(tags$th(""), lapply(colnames(D), function(c) tags$th(c)),
              if (cell_col) tags$th(class = "scroll-designmat-n", "cells")),
      lapply(seq_len(nrow(D)), function(i)
        tags$tr(tags$th(rownames(D)[i]),
                lapply(D[i, ], function(v) tags$td(format(v))),
                ccell(format(counts[i], big.mark = ",")))),
      foot))
}

# Model formula + a collapsible, capped design matrix for the previewed samples.
# The matrix is the exact model scroll_pseudobulk_de() will fit; a high-cardinality
# replicate is capped to a one-line summary so it never bloats the sidebar.
.scroll_pseudobulk_design_ui <- function(sm) {
  n1 <- sum(sm$samp$group == "group1"); n2 <- sum(sm$samp$group == "group2")
  small <- function(txt) div(class = "scroll-desc", style = "margin:4px 0 0; font-size:11px", txt)
  if (n1 < 2 || n2 < 2)
    return(small(sprintf("Model: needs >= 2 samples/group (have %d vs %d).", n1, n2)))
  ds <- .scroll_pseudobulk_design(sm$samp, sm$regime, "auto")
  # cells per sample (= per group x replicate) and per-group totals
  ncell <- as.integer(table(sm$mapping$psample)[sm$samp$psample]); ncell[is.na(ncell)] <- 0L
  g1c <- sum(ncell[sm$samp$group == "group1"]); g2c <- sum(ncell[sm$samp$group == "group2"])
  head <- small(tagList(tags$b("Model: "), ds$formula,
                        sprintf("  (%d vs %d samples)", n1, n2)))
  cellsln <- small(tagList(tags$b("Cells: "),
                    sprintf("group1 = %s, group2 = %s",
                            format(g1c, big.mark = ","), format(g2c, big.mark = ","))))
  # the tested effect: group1 - group2 (positive logFC = up in ident1)
  pos <- .scroll_tidy_design_cols(names(ds$contrast)[ds$contrast > 0])
  neg <- .scroll_tidy_design_cols(names(ds$contrast)[ds$contrast < 0])
  contr <- small(tagList(tags$b("Contrast: "), paste(pos, "-", neg),
                         tags$span(style = "color:var(--sc-faint)", "  (up = ident1)")))
  D <- ds$design
  if (nrow(D) > 16 || ncol(D) > 12)                  # keep the sidebar sane
    return(tagList(head, cellsln, contr,
                   small(sprintf("design matrix %d x %d - too large to show.",
                                 nrow(D), ncol(D)))))
  rownames(D) <- sm$samp$psample
  colnames(D) <- .scroll_tidy_design_cols(colnames(D))
  tagList(head, cellsln, contr,
          .scroll_details("Design matrix", open = FALSE,
                          .scroll_matrix_html(D, footer = as.numeric(ds$contrast),
                                              counts = ncell)))
}

pseudobulk_de_ui <- function(id, data) {
  ns <- NS(id)
  m <- data$manifest
  cats <- .scroll_cat_cols(m); assays <- .scroll_assays_of(m)
  if (!isTRUE(m$has_counts))
    return(.scroll_empty_panel(paste(
      "Pseudobulk DE needs a counts store. Rebuild the project with",
      "scroll_build(..., counts = TRUE).")))
  if (!length(cats)) return(.scroll_empty_panel("Needs a categorical metadata column (none found)."))
  bslib::layout_columns(
    col_widths = c(3, 9), class = "scroll-panel",
    div(
      class = "scroll-controls",
      .scroll_group("Contrast",
        selectizeInput(ns("aggregate_by"), "Aggregate by", stats::setNames(cats, cats),
                       multiple = TRUE, selected = cats[[1]],
                       options = list(plugins = list("remove_button"))),
        selectizeInput(ns("ident1"), "Group 1", choices = NULL, multiple = TRUE,
                       options = list(plugins = list("remove_button"))),
        selectizeInput(ns("ident2"), "vs.", choices = NULL, multiple = TRUE,
                       options = list(plugins = list("remove_button"),
                                      placeholder = "rest (all other cells)")),
        selectInput(ns("replicate"), "Replicate",
                    c("No replicate (pseudo)" = "no_replicate", stats::setNames(cats, cats))),
        if (length(assays) > 1) selectInput(ns("assay"), "Assay", assays, selected = m$default_assay)),
      .scroll_group("Preview",                       # live: one coloured medoid per sample
        plotOutput(ns("preview"), height = "170px"),
        uiOutput(ns("design"))),                      # model formula + design matrix
      .scroll_group("Pseudobulk",
        sliderInput(ns("mincells"), "Min cells / sample", 3, 200, 10, 1),
        numericInput(ns("npseudo"), "Pseudo-reps (if no replicate)", 3, min = 2, max = 10),
        numericInput(ns("cellsper"), "Max cells / pseudo-rep", 50, min = 5, max = 2000),
        div(class = "scroll-desc", style = "margin:2px 0 0; font-size:11px",
            "Pseudo-reps are a disjoint partition (no shared cells); a group needs ",
            "at least Pseudo-reps x Min cells cells.")),
      .scroll_group("Stability",
        numericInput(ns("runs"), "Runs (re-sample; pseudo only)", 1, min = 1, max = 100),
        sliderInput(ns("stabcut"), "Consistency cutoff", 0, 1, 0.8, 0.05)),
      .scroll_group("Table",
        sliderInput(ns("topn"), "Show top", 10, 300, 50, 10),
        bslib::input_switch(ns("export_all"), "Export all genes (CSV)", FALSE)),
      .scroll_group("Volcano / stability",
        sliderInput(ns("lfc"), "logFC cutoff", 0, 3, 1, 0.1),
        numericInput(ns("padj"), "Adj. p cutoff", 0.05, min = 0, max = 1, step = 0.01),
        sliderInput(ns("labeln"), "Label top", 0, 40, 15, 1),
        .scroll_aspect_input(ns)),
      # input_task_button: the button shows a spinner + "Computing..." and disables
      # the instant it is clicked (immediate feedback during the blocking compute).
      bslib::input_task_button(ns("compute"), "Compute pseudobulk DE", type = "primary",
                               label_busy = "Computing...", class = "w-100")
    ),
    div(class = "scroll-plot",
        uiOutput(ns("note")),
        bslib::navset_tab(
          bslib::nav_panel("Table",
            div(class = "scroll-plot-bar", .scroll_dl_button(ns("csv"), "CSV")),
            div(class = "scroll-table",
                .scroll_spin(if (.scroll_has_dt()) DT::dataTableOutput(ns("table"))
                             else tableOutput(ns("table"))))),
          bslib::nav_panel("Plot",
            div(class = "scroll-plot-bar",
                .scroll_size_slider(ns),
                .scroll_dl_button(ns("png"), "PNG"), .scroll_dl_button(ns("pdf"), "PDF")),
            .scroll_spin(plotOutput(ns("plot"), height = .SCROLL_PLOT_H)))))
  )
}

pseudobulk_de_server <- function(id, data, cells_r = reactive(data$cells),
                                 view_r = reactive(NULL), theme_r = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    m <- data$manifest
    assay <- reactive(input$assay %||% m$default_assay)

    # Aggregate-by / Replicate follow the active view: subset-scoped categorical
    # columns (e.g. a subset's own clustering resolutions) surface only in their
    # view, and global columns stay available everywhere. The UI builds these from
    # the whole-dataset columns; refresh them when the view changes.
    observeEvent(view_r(), {
      cats_v <- .scroll_cat_cols(m, view_r())
      cur_agg <- intersect(isolate(input$aggregate_by), cats_v)
      if (!length(cur_agg) && length(cats_v)) cur_agg <- cats_v[[1]]
      updateSelectizeInput(session, "aggregate_by",
                           choices = stats::setNames(cats_v, cats_v), selected = cur_agg)
      cur_rep <- isolate(input$replicate)
      sel_rep <- if (!is.null(cur_rep) && (identical(cur_rep, "no_replicate") || cur_rep %in% cats_v))
                   cur_rep else "no_replicate"
      updateSelectInput(session, "replicate",
                        choices = c("No replicate (pseudo)" = "no_replicate",
                                    stats::setNames(cats_v, cats_v)),
                        selected = sel_rep)
    }, ignoreNULL = FALSE)

    # combined levels depend on the chosen aggregate-by columns + active subset
    observeEvent(list(input$aggregate_by, cells_r()), {
      combos <- .scroll_combo_choices(cells_r(), input$aggregate_by)
      updateSelectizeInput(session, "ident1", choices = combos,
                           selected = intersect(isolate(input$ident1), combos))
      updateSelectizeInput(session, "ident2", choices = combos,
                           selected = intersect(isolate(input$ident2), combos))
    })

    # Live sample preview (NOT gated by Compute): grey embedding outline + one
    # coloured point per pseudobulk sample at its medoid (red = group1, blue =
    # group2), so users see how cells compact into samples as they pick controls.
    preview_in <- .scroll_cosmetic(reactive(list(
      cells = cells_r(), agg = input$aggregate_by, ident1 = input$ident1,
      ident2 = input$ident2, rep = input$replicate, mincells = input$mincells,
      npseudo = input$npseudo, cellsper = input$cellsper,
      emb = .scroll_preview_embedding(data, view_r()))))   # default for whole dataset, else the view's

    # The exact samples the compute will build (metadata-only; no counts query), so
    # the preview medoids AND the design matrix match scroll_pseudobulk_de.
    pb_samples <- reactive({
      p <- preview_in()
      if (!length(p$agg) || !length(p$ident1)) return(NULL)
      rep_col <- if (!is.null(p$rep) && nzchar(p$rep)) p$rep else "no_replicate"
      tryCatch(
        .scroll_pseudobulk_samples(p$cells, seq_len(nrow(p$cells)), p$agg, p$ident1,
          if (length(p$ident2)) p$ident2 else NULL, rep_col,
          p$mincells %||% 10, p$npseudo %||% 3, p$cellsper %||% 50, seed = 1L),
        error = function(e) list(err = conditionMessage(e)))
    })

    # medoid per actual pseudobulk sample (one per group x replicate, or per
    # pseudo-partition bin) -- so the preview counts match what Compute uses.
    output$preview <- renderPlot({
      p <- preview_in(); req(length(p$agg) > 0, p$emb)
      sm <- pb_samples()
      validate(need(!is.null(sm), "Pick Group 1 to preview the samples."),
               need(is.null(sm$err), sm$err))
      samp_vec <- rep(NA_character_, nrow(p$cells)); grp_vec <- samp_vec
      samp_vec[sm$mapping$cell] <- sm$mapping$psample
      grp_vec[sm$mapping$cell]  <- sm$mapping$group
      view_contrast_medoids(p$cells, p$emb, samp_vec, grp_vec)
    })

    # the model that will be fit: formula + a collapsible, capped design matrix
    output$design <- renderUI({
      sm <- pb_samples(); if (is.null(sm) || !is.null(sm$err)) return(NULL)
      .scroll_pseudobulk_design_ui(sm)
    })

    # stability = re-run the pseudo-replication K times (pseudo mode only). The
    # K runs are computed once here; de_df() aggregates them cheaply when the
    # lfc/padj cutoffs change.
    result <- eventReactive(input$compute, {
      req(input$aggregate_by, input$ident1)
      args <- list(
        data, assay(), aggregate_cols = input$aggregate_by, ident1 = input$ident1,
        ident2 = if (length(input$ident2)) input$ident2 else NULL,
        replicate_col = input$replicate, min_cells = input$mincells,
        n_pseudo = input$npseudo, cells_per_pseudo = input$cellsper, cells = cells_r())
      stability <- isTRUE(input$runs > 1) && identical(input$replicate, "no_replicate")
      withProgress(message = "Pseudobulk DE", value = 0,
        tryCatch(
          if (stability) {
            setProgress(0.1, detail = "Running stability resamples...")
            list(runs = do.call(.scroll_pseudobulk_runs, c(args, list(runs = input$runs))))
          } else list(ok = do.call(scroll_pseudobulk_de,
                        c(args, list(progress = function(f, d) setProgress(value = f, detail = d))))),
          error = function(e) list(err = conditionMessage(e))))
    })

    de_df <- reactive({
      validate(need(input$compute > 0, "Pick a contrast and click Compute pseudobulk DE."))
      r <- result()
      validate(need(is.null(r$err), r$err))
      if (!is.null(r$runs)) .scroll_stability_aggregate(r$runs, input$lfc, input$padj)
      else r$ok
    })

    output$note <- renderUI({
      req(input$compute > 0); r <- result(); if (!is.null(r$err)) return(NULL)
      amber <- function(txt) div(class = "scroll-desc", style = "margin:0 0 8px; color:#B45309", txt)
      note  <- function(txt) div(class = "scroll-desc", style = "margin:0 0 8px", txt)
      if (!is.null(r$runs))
        return(amber(sprintf(paste("Stability across %d pseudo-replicate runs - sel_freq is a",
                         "robustness heuristic, not a p-value; prefer real replicates."),
                   length(r$runs))))
      res <- r$ok; n <- attr(res, "n_samples")
      na <- attr(res, "dropped_na") %||% 0L
      na_txt <- if (isTRUE(na > 0)) sprintf(" %d cells with a missing replicate excluded.", na) else ""
      switch(attr(res, "regime") %||% "pseudo",
        "real-paired" = note(sprintf(
          "Paired design (~ replicate + group): %d vs %d replicate samples.%s",
          n[["group1"]], n[["group2"]], na_txt)),
        "real-unpaired" = note(sprintf(
          "Real replicates (~ group): %d vs %d samples.%s",
          n[["group1"]], n[["group2"]], na_txt)),
        "mixed" = amber(sprintf(paste("Some groups lacked replicates; pseudo-replicates used",
          "there (%d vs %d samples) - treat with caution.%s"),
          n[["group1"]], n[["group2"]], na_txt)),
        amber(sprintf(paste("Pseudo-replicates used (no biological replicates): %d vs %d",
          "samples - treat p-values with caution."), n[["group1"]], n[["group2"]])))
    })

    # display columns differ between single-run and stability results
    table_rows <- function() {
      d <- de_df()
      if ("sel_freq" %in% names(d))
        d <- d[, c("gene", "sel_freq", "median_logFC", "sign_agree", "median_padj", "n_tested")]
      utils::head(d, input$topn)
    }
    if (.scroll_has_dt())
      output$table <- DT::renderDataTable({
        d <- table_rows()
        DT::formatSignif(
          DT::datatable(d, rownames = FALSE, options = list(pageLength = 15, dom = "tip")),
          columns = intersect(c("p_val", "p_val_adj", "median_padj"), names(d)), digits = 3)
      })
    else
      output$table <- renderTable(table_rows())

    plot_r <- reactive({
      d <- de_df()
      if ("sel_freq" %in% names(d))
        view_stability(d, params = list(lfc = input$lfc, cut = input$stabcut,
                                        label_n = input$labeln), state = list(aspect = input$aspect, theme = theme_r()))
      else
        view_volcano(d, params = list(lfc = input$lfc, padj = input$padj,
                                      label_n = input$labeln), state = list(aspect = input$aspect, theme = theme_r()))
    })
    output$plot <- renderPlot(plot_r())
    .scroll_plot_downloads(output, plot_r, id)
    output$csv <- .scroll_csv_handler(
      reactive(if (isTRUE(input$export_all)) de_df() else table_rows()),
      paste0("scroll_", id, ".csv"))
  })
}

