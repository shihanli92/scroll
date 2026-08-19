# --- Pseudobulk DE panel ------------------------------------------------------

# Combined-interaction levels observed in `cells` for the given columns.
.scroll_combo_choices <- function(cells, cols) {
  if (!length(cols)) return(character(0))
  combo <- .scroll_combo_levels(cells, cols)
  sort(unique(combo[!is.na(combo)]))
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
        plotOutput(ns("preview"), height = "170px")),
      .scroll_group("Pseudobulk",
        sliderInput(ns("mincells"), "Min cells / sample", 3, 200, 10, 1),
        numericInput(ns("npseudo"), "Pseudo-reps (if no replicate)", 3, min = 2, max = 10),
        numericInput(ns("cellsper"), "Cells / pseudo-rep", 50, min = 5, max = 2000)),
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
                .scroll_dl_button(ns("png"), "PNG"), .scroll_dl_button(ns("pdf"), "PDF")),
            .scroll_spin(plotOutput(ns("plot"), height = "520px")))))
  )
}

pseudobulk_de_server <- function(id, data, cells_r = reactive(data$cells),
                                 view_r = reactive(NULL)) {
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
      ident2 = input$ident2, rep = input$replicate,
      emb = .scroll_preview_embedding(data, view_r()))))   # default for whole dataset, else the view's
    output$preview <- renderPlot({
      p <- preview_in(); req(length(p$agg) > 0, p$emb)
      combo <- .scroll_combo_levels(p$cells, p$agg)
      grp <- .scroll_combo_group(combo, p$ident1, p$ident2)
      samp <- if (!is.null(p$rep) && nzchar(p$rep) && p$rep != "no_replicate" &&
                  p$rep %in% names(p$cells))
                paste(combo, as.character(p$cells[[p$rep]]), sep = " :: ") else combo
      view_contrast_medoids(p$cells, p$emb, samp, grp)
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
      tryCatch(
        if (stability)
          list(runs = do.call(.scroll_pseudobulk_runs, c(args, list(runs = input$runs))))
        else list(ok = do.call(scroll_pseudobulk_de, args)),
        error = function(e) list(err = conditionMessage(e)))
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
      mk <- function(txt) div(class = "scroll-desc", style = "margin:0 0 8px; color:#B45309", txt)
      if (!is.null(r$runs))
        mk(sprintf(paste("Stability across %d pseudo-replicate runs - sel_freq is a",
                         "robustness heuristic, not a p-value; prefer real replicates."),
                   length(r$runs)))
      else if (isTRUE(attr(r$ok, "pseudo"))) {
        n <- attr(r$ok, "n_samples")
        mk(sprintf(paste("Pseudo-replicates used (no biological replicates): %d vs %d",
                         "samples - treat p-values with caution."), n[["group1"]], n[["group2"]]))
      } else NULL
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
                                        label_n = input$labeln), state = list(aspect = input$aspect))
      else
        view_volcano(d, params = list(lfc = input$lfc, padj = input$padj,
                                      label_n = input$labeln), state = list(aspect = input$aspect))
    })
    output$plot <- renderPlot(plot_r())
    .scroll_plot_downloads(output, plot_r, id)
    output$csv <- .scroll_csv_handler(
      reactive(if (isTRUE(input$export_all)) de_df() else table_rows()),
      paste0("scroll_", id, ".csv"))
  })
}

