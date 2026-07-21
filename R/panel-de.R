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
        selectInput(ns("group"), "Group by", stats::setNames(cats, cats), selected = cats[[1]]),
        selectizeInput(ns("ident1"), "Group 1", choices = NULL, multiple = TRUE,
                       options = list(plugins = list("remove_button"))),
        selectizeInput(ns("ident2"), "vs.", choices = NULL, multiple = TRUE,
                       options = list(plugins = list("remove_button"),
                                      placeholder = "rest (all other cells)")),
        if (length(assays) > 1) selectInput(ns("assay"), "Assay", assays, selected = m$default_assay)),
      .scroll_group("Table",
        sliderInput(ns("minpct"), "Min % expressing", 0, 50, 10, 1),
        sliderInput(ns("topn"), "Show top", 10, 300, 50, 10),
        bslib::input_switch(ns("export_all"), "Export all genes (CSV)", FALSE)),
      .scroll_group("Volcano",
        sliderInput(ns("lfc"), "logFC cutoff", 0, 3, 1, 0.1),
        numericInput(ns("padj"), "Adj. p cutoff", 0.05, min = 0, max = 1, step = 0.01),
        sliderInput(ns("labeln"), "Label top", 0, 40, 15, 1),
        .scroll_aspect_input(ns)),
      actionButton(ns("compute"), "Compute DE", class = "btn-primary", width = "100%")
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
                .scroll_dl_button(ns("png"), "PNG"),
                .scroll_dl_button(ns("pdf"), "PDF")),
            .scroll_spin(plotOutput(ns("plot"), height = "520px")))))
  )
}

.scroll_has_dt <- function() requireNamespace("DT", quietly = TRUE)

de_server <- function(id, data, cells_r = reactive(data$cells),
                      view_r = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    m <- data$manifest
    assay <- reactive(input$assay %||% m$default_assay)
    .scroll_bind_view_cats(input, session, view_r, m, "group")

    observeEvent(input$group, {
      lv <- .scroll_meta_levels(data, input$group)
      updateSelectizeInput(session, "ident1", choices = lv, selected = lv[[1]])
      updateSelectizeInput(session, "ident2", choices = lv, selected = character(0))
    })

    # min_pct = 0 so the volcano keeps every gene; the table applies its own
    # %-expressing filter below. A failed presto run (e.g. a contrast with too
    # few cells) is caught and surfaced as a friendly inline message, not a raw
    # Shiny error.
    result <- eventReactive(input$compute, {
      req(input$ident1)
      tryCatch(
        list(ok = scroll_de(data, assay(), input$group, input$ident1,
                            if (length(input$ident2)) input$ident2 else NULL,
                            min_pct = 0, cells = cells_r())),
        error = function(e) list(err = conditionMessage(e)))
    })

    de_df <- reactive({
      validate(need(input$compute > 0, "Pick a contrast and click Compute DE."))
      r <- result()
      validate(need(is.null(r$err), r$err))
      r$ok
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
      view_volcano(de_df(),
                   params = list(lfc = input$lfc, padj = input$padj, label_n = input$labeln),
                   state = list(aspect = input$aspect))
    })
    output$plot <- renderPlot(volcano_r())
    .scroll_plot_downloads(output, volcano_r, id)
    output$csv <- .scroll_csv_handler(
      reactive(if (isTRUE(input$export_all)) de_df() else table_rows()),
      paste0("scroll_", id, ".csv"))
  })
}

