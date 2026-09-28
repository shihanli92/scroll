# --- DE panel -----------------------------------------------------------------

# One contrast, computed once via presto, shown two ways: a ranked marker
# Table and a Volcano. A single `result` reactive feeds both tabs; the table's
# min-%-expressing and top-N are display filters (so the volcano keeps every
# gene), and the volcano thresholds only restyle the plot.
de_ui <- function(id, data) {
  ns <- NS(id)
  m <- data$manifest
  cats <- .scroll_cat_cols(m); assays <- .scroll_assays_of(m)
  if (!length(cats)) return(.scroll_empty_panel("Needs a categorical metadata column (none found)."))
  bslib::layout_columns(
    col_widths = c(3, 9), class = "scroll-panel",
    div(
      class = "scroll-controls",
      .scroll_group("Contrast",
        selectizeInput(ns("group"), "Group by", stats::setNames(cats, cats), multiple = TRUE,
                       selected = cats[[1]], options = list(plugins = list("remove_button"))),
        selectizeInput(ns("ident1"), "Group 1", choices = NULL, multiple = TRUE,
                       options = list(plugins = list("remove_button"))),
        selectizeInput(ns("ident2"), "vs.", choices = NULL, multiple = TRUE,
                       options = list(plugins = list("remove_button"),
                                      placeholder = "rest (all other cells)")),
        if (length(assays) > 1) selectInput(ns("assay"), "Assay", assays, selected = m$default_assay),
        numericInput(ns("maxcells"), "Max cells / group (0 = all)", 0, min = 0, step = 500)),
      .scroll_group("Preview",                       # live: red = group1, blue = group2/rest
        plotOutput(ns("preview"), height = "170px")),
      .scroll_group("Table",
        sliderInput(ns("minpct"), "Min % expressing", 0, 50, 10, 1),
        sliderInput(ns("topn"), "Show top", 10, 300, 50, 10),
        .scroll_export_input(ns)),
      .scroll_group("Volcano",
        sliderInput(ns("lfc"), "avg_log2FC cutoff", 0, 3, 1, 0.1),
        numericInput(ns("padj"), "Adj. p cutoff", 0.05, min = 0, max = 1, step = 0.01)),
      # input_task_button (not actionButton) so the button itself shows a spinner +
      # "Computing..." and disables the instant it is clicked — immediate feedback
      # during the (blocking) DE compute, independent of the output-area spinner
      # (which needs the optional shinycssloaders and sits away from the click).
      bslib::input_task_button(ns("compute"), "Compute DE", type = "primary",
                               label_busy = "Computing DE...", class = "w-100")
    ),
    div(class = "scroll-plot",
        bslib::navset_tab(
          bslib::nav_panel("Table",
            div(class = "scroll-plot-bar", .scroll_dl_button(ns("csv"), "CSV")),
            div(class = "scroll-table",
                .scroll_spin(if (.scroll_has_dt()) DT::dataTableOutput(ns("table"))
                             else tableOutput(ns("table"))))),
          bslib::nav_panel("Volcano",
            div(class = "scroll-plot-bar",
                .scroll_style_button(ns), .scroll_size_slider(ns),
                .scroll_dl_button(ns("png"), "PNG"),
                .scroll_dl_button(ns("pdf"), "PDF")),
            .scroll_spin(plotOutput(ns("plot"), height = .SCROLL_PLOT_H)))))
  )
}

.scroll_has_dt <- function() requireNamespace("DT", quietly = TRUE)

# What the Table tab's CSV button writes: the rows on screen, only the genes that
# pass the panel's cutoffs, or every gene tested.
.scroll_export_input <- function(ns)
  selectInput(ns("export"), "CSV export",
              c("Shown rows" = "shown", "Genes passing cutoffs" = "sig", "All genes" = "all"))

# Genes passing a DE result's cutoffs, the same rule the volcano / stability plot
# colours by: adj. p < padj and |fold change| >= lfc (plus a min fraction
# expressing when the result has pct.1 / pct.2). A stability result (runs > 1)
# passes on selection frequency >= cut and |median logFC| >= lfc.
.scroll_de_passing <- function(df, lfc = 1, padj = 0.05, fc_col = "logFC", min_pct = 0, cut = 0.8) {
  if (is.null(df) || !nrow(df)) return(df)
  keep <- if ("sel_freq" %in% names(df)) {
    df$sel_freq >= cut & abs(df$median_logFC) >= lfc
  } else {
    if (!fc_col %in% names(df)) fc_col <- "logFC"
    ok <- df$p_val_adj < padj & abs(df[[fc_col]]) >= lfc
    if (all(c("pct.1", "pct.2") %in% names(df)))
      ok <- ok & pmax(df$pct.1, df$pct.2) >= min_pct / 100
    ok
  }
  df[!is.na(keep) & keep, , drop = FALSE]
}

# The volcano's look controls, shown in the panel's Style sheet (same input ids).
de_style_ui <- function(id, data) {
  ns <- NS(id)
  tagList(sliderInput(ns("labeln"), "Label top genes", 0, 40, 15, 1),
          .scroll_aspect_input(ns))
}

de_server <- function(id, data, cells_r = reactive(data$cells),
                      view_r = reactive(NULL), theme_r = reactive(NULL), style_r = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    m <- data$manifest
    assay <- reactive(input$assay %||% m$default_assay)

    # Group-by may combine several categorical columns into interaction levels
    # (e.g. "genotype | timepoint"); the idents then pick which combos form each
    # side. The column choices track the active subset view (scoped columns).
    observeEvent(.scroll_cols_event(view_r, m), {   # + session columns
      cats <- .scroll_cat_cols(m, view_r())
      sel <- intersect(input$group, cats); if (!length(sel)) sel <- cats[[1]]
      updateSelectizeInput(session, "group", choices = stats::setNames(cats, cats), selected = sel)
    }, ignoreNULL = FALSE)

    observeEvent(list(input$group, cells_r()), {
      req(length(input$group) > 0)
      combos <- .scroll_combo_choices(cells_r(), input$group)     # interaction levels
      cur1 <- intersect(isolate(input$ident1), combos)
      if (!length(cur1) && length(combos)) cur1 <- combos[[1]]    # default a ready contrast
      updateSelectizeInput(session, "ident1", choices = combos, selected = cur1)
      updateSelectizeInput(session, "ident2", choices = combos,
                           selected = intersect(isolate(input$ident2), combos))
    })

    # Live contrast preview (NOT gated by Compute): a mini-UMAP showing which cells
    # each side selects — red = group1, blue = group2/rest, grey = everything else.
    # Debounced (selectize multi-picks fire fast); embedding follows the active view.
    preview_in <- .scroll_cosmetic(reactive(list(
      cells = cells_r(), group = input$group, ident1 = input$ident1,
      ident2 = input$ident2, emb = .scroll_preview_embedding(data, view_r()))))
    output$preview <- renderPlot({
      p <- preview_in(); req(p$group, p$emb)
      view_contrast_preview(p$cells, p$emb,
                            .scroll_contrast_labels(p$cells, p$group, p$ident1, p$ident2))
    })

    # min_pct = 0 so the volcano keeps every gene; the table applies its own
    # %-expressing filter below. A failed presto run (e.g. a contrast with too
    # few cells) is caught and surfaced as a friendly inline message, not a raw
    # Shiny error.
    result <- eventReactive(input$compute, {
      req(input$ident1)
      withProgress(message = "Differential expression", value = 0,
        tryCatch(
          list(ok = scroll_de(data, assay(), input$group, input$ident1,
                              if (length(input$ident2)) input$ident2 else NULL,
                              min_pct = 0, cells = cells_r(),
                              max_cells = if (isTRUE(input$maxcells > 0)) input$maxcells else NULL,
                              progress = function(f, d) setProgress(value = f, detail = d))),
          error = function(e) list(err = conditionMessage(e))))
    })

    de_df <- reactive({
      validate(need(input$compute > 0, "Pick a contrast and click Compute DE."))
      r <- result()
      validate(need(is.null(r$err), r$err))
      r$ok
    })

    # share the latest result with the GSEA panel (this session)
    observeEvent(result(), {
      r <- result()
      if (!is.null(r$err)) return()
      .scroll_publish_result(data, "de", r$ok, .scroll_contrast_text(
        input$group, input$ident1, input$ident2, sprintf("Wilcoxon, %s", assay())))
    })

    table_rows <- function() {
      res <- de_df()
      res <- res[pmax(res$pct.1, res$pct.2) >= input$minpct / 100, , drop = FALSE]
      utils::head(res, input$topn)
    }
    if (.scroll_has_dt())
      output$table <- DT::renderDataTable(
        DT::formatSignif(
          DT::datatable(table_rows(), rownames = FALSE, options = list(pageLength = 15, dom = "tip")),
          columns = c("p_val", "p_val_adj"), digits = 3))
    else
      output$table <- renderTable(table_rows())

    volcano_r <- reactive({
      .scroll_style_plot(view_volcano(de_df(),
                   params = list(lfc = input$lfc, padj = input$padj, label_n = input$labeln,
                                 fc_col = "avg_log2FC", fc_label = "avg_log2FC"),
                   state = list(aspect = input$aspect, theme = theme_r())), style_r)
    })
    output$plot <- renderPlot(volcano_r())
    .scroll_plot_downloads(output, volcano_r, id)
    export_r <- reactive(switch(input$export %||% "shown",
      all = de_df(),
      sig = .scroll_de_passing(de_df(), lfc = input$lfc %||% 1, padj = input$padj %||% 0.05,
                               fc_col = "avg_log2FC", min_pct = input$minpct %||% 0),
      table_rows()))
    output$csv <- .scroll_csv_handler(export_r, paste0("scroll_", id, ".csv"))
  })
}

