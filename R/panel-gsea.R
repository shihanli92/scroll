# --- GSEA panel -----------------------------------------------------------------

# Preranked gene-set enrichment on the latest DE or Pseudobulk result of this
# session (published to data$results by those panels; see R/genesets.R). No contrast
# controls of its own: run a contrast there, then pick it here. Gene sets are the
# ones saved in the project (scroll_add_genesets) and, with msigdbr installed, the
# curated MSigDB collections fetched live.

gsea_ui <- function(id, data) {
  ns <- NS(id)
  m <- data$manifest
  sources <- c("DE (Wilcoxon)" = "de",
               if (isTRUE(m$has_counts)) c("Pseudobulk DE" = "pseudobulk"))
  colls <- .scroll_geneset_choices(data)
  bslib::layout_columns(
    col_widths = c(3, 9), class = "scroll-panel",
    div(
      class = "scroll-controls",
      .scroll_group("Ranked genes",
        radioButtons(ns("source"), "From", sources, inline = TRUE),
        uiOutput(ns("source_note")),
        selectInput(ns("rank_by"), "Rank by", character(0))),
      .scroll_group("Gene sets",
        if (length(colls)) selectInput(ns("collection"), "Collection", colls)
        else helpText("No gene sets: save some with scroll_add_genesets(), or install",
                      "the msigdbr package for MSigDB collections."),
        bslib::layout_columns(col_widths = c(6, 6), gap = "8px",
          numericInput(ns("min_size"), "Min size", 15, min = 1, step = 5),
          numericInput(ns("max_size"), "Max size", 500, min = 5, step = 50))),
      bslib::input_task_button(ns("run"), "Run GSEA", type = "primary",
                               label_busy = "Running GSEA...", class = "w-100"),
      uiOutput(ns("stale")),
      .scroll_group("Results",
        numericInput(ns("padj"), "Adj. p cutoff", 0.05, min = 0, max = 1, step = 0.01),
        radioButtons(ns("direction"), "Direction", c("Both" = "both", "Up" = "up", "Down" = "down"),
                     inline = TRUE),
        sliderInput(ns("topn"), "Show top", 5, 50, 15, 1)),
      .scroll_group("Plot",
        selectInput(ns("show"), "Show", c("Top pathways" = "top", "Enrichment plot" = "enrich")),
        .scroll_cond_panel(sprintf("input['%s'] == 'enrich'", ns("show")),
          selectizeInput(ns("pathway"), "Pathway", choices = NULL,
                         options = list(placeholder = "Pick a pathway (or click a table row)"))))
    ),
    div(class = "scroll-plot",
        bslib::navset_tab(
          bslib::nav_panel("Table",
            div(class = "scroll-plot-bar", .scroll_dl_button(ns("csv"), "CSV")),
            div(class = "scroll-table",
                .scroll_spin(if (.scroll_has_dt()) DT::dataTableOutput(ns("table"))
                             else tableOutput(ns("table"))))),
          bslib::nav_panel("Plot",
            div(class = "scroll-plot-bar",
                .scroll_style_button(ns), .scroll_size_slider(ns),
                .scroll_dl_button(ns("png"), "PNG"), .scroll_dl_button(ns("pdf"), "PDF")),
            .scroll_spin(plotOutput(ns("plot"), height = .SCROLL_PLOT_H)))))
  )
}

gsea_style_ui <- function(id, data) {
  ns <- NS(id)
  tagList(selectInput(ns("palette"), "Palette (top pathways)", .scroll_continuous_palettes,
                      selected = "viridis"))
}

# A readable pathway name: drop the collection prefix, underscores -> spaces.
.scroll_pathway_label <- function(x, width = 45) {
  x <- sub("^(HALLMARK|GOBP|GOMF|GOCC|REACTOME|KEGG|WP|BIOCARTA|PID)_", "", x)
  x <- gsub("_", " ", x)
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), "", USE.NAMES = FALSE)
}

# The top pathways up and down as NES bars, coloured by -log10 padj.
.scroll_gsea_top_plot <- function(res, n = 15, padj = 0.05, direction = "both",
                                  palette = "viridis", theme = NULL) {
  hit <- res[!is.na(res$padj) & res$padj <= padj, , drop = FALSE]
  if (!nrow(hit)) stop(sprintf("No pathway passes adj. p <= %g; raise the cutoff.", padj), call. = FALSE)
  up <- hit[hit$NES > 0, , drop = FALSE]; dn <- hit[hit$NES < 0, , drop = FALSE]
  up <- utils::head(up[order(-up$NES), , drop = FALSE], n)
  dn <- utils::head(dn[order(dn$NES), , drop = FALSE], n)
  df <- switch(direction, up = up, down = dn, rbind(up, dn))
  if (!nrow(df)) stop("No pathway in that direction passes the cutoff.", call. = FALSE)
  df$label <- factor(.scroll_pathway_label(df$pathway), levels = unique(.scroll_pathway_label(df$pathway[order(df$NES)])))
  df$score <- -log10(pmax(df$padj, .Machine$double.xmin))
  ggplot2::ggplot(df, ggplot2::aes(.data$NES, .data$label, fill = .data$score)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
    .scroll_continuous_scale(palette, name = "-log10 adj. p", aesthetic = "fill") +
    ggplot2::labs(x = "Normalized enrichment score (NES)", y = NULL,
                  title = "Top enriched gene sets") +
    .scroll_base_theme() + .scroll_ggtheme(theme)
}

# The running enrichment score for one pathway (fgsea's plotEnrichment).
.scroll_gsea_enrichment_plot <- function(res, pathway, theme = NULL) {
  sets <- attr(res, "sets"); stats <- attr(res, "stats")
  if (is.null(pathway) || !pathway %in% names(sets)) stop("Pick a pathway to plot.", call. = FALSE)
  row <- res[res$pathway == pathway, , drop = FALSE]
  fgsea::plotEnrichment(sets[[pathway]], stats) +
    ggplot2::labs(title = .scroll_pathway_label(pathway, 70),
                  subtitle = sprintf("NES %.2f - adj. p %.2g - %d genes", row$NES, row$padj, row$size),
                  x = sprintf("Rank (by %s)", attr(res, "rank_by")), y = "Enrichment score") +
    .scroll_base_theme() + .scroll_ggtheme(theme)
}

gsea_server <- function(id, data, cells_r = reactive(data$cells),
                        view_r = reactive(NULL), theme_r = reactive(NULL), style_r = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    results <- if (is.function(data$results)) data$results else reactiveVal(list())
    src <- reactive(results()[[input$source %||% "de"]])

    output$source_note <- renderUI({
      s <- src()
      if (is.null(s))
        return(helpText(if (identical(input$source, "pseudobulk"))
                          "Run a contrast in the Pseudobulk DE panel first."
                        else "Run a contrast in the DE panel first."))
      helpText(sprintf("%s (%s genes, %s)", s$label, format(s$n_genes, big.mark = ","),
                       format(s$time, "%H:%M")))
    })
    # rank metrics offered by the chosen result (limma t first for Pseudobulk)
    observeEvent(src(), {
      ch <- .scroll_rank_choices(src()$ranking %||% src()$df)
      cur <- isolate(input$rank_by)
      updateSelectInput(session, "rank_by", choices = ch,
                        selected = if (!is.null(cur) && cur %in% ch) cur else ch[1])
    }, ignoreNULL = FALSE)

    gsea_r <- eventReactive(input$run, {
      s <- src()
      validate(need(!is.null(s), "Run a contrast in the DE or Pseudobulk panel first, then Run GSEA."))
      validate(need(length(input$collection) && nzchar(input$collection), "Pick a gene-set collection."))
      withProgress(message = "GSEA", value = 0.2, {
        out <- tryCatch({
          setProgress(0.2, detail = "Loading gene sets...")
          sets <- .scroll_geneset_get(data, input$collection)
          setProgress(0.5, detail = sprintf("Testing %s sets...", format(length(sets), big.mark = ",")))
          de <- s$df; attr(de, "ranking") <- s$ranking
          scroll_gsea(de, sets, rank_by = input$rank_by,
                      min_size = max(1, input$min_size %||% 15), max_size = max(5, input$max_size %||% 500))
        }, error = function(e) e)
        validate(need(!inherits(out, "error"), if (inherits(out, "error")) conditionMessage(out)))
        list(res = out, time = s$time, source = input$source, label = s$label)
      })
    })
    output$stale <- renderUI({
      req(input$run > 0)
      g <- tryCatch(gsea_r(), error = function(e) NULL); s <- src()
      if (!is.null(g) && !is.null(s) && identical(g$source, input$source) && !identical(g$time, s$time))
        helpText("The ranked result changed since this run; Run GSEA again to update.")
    })

    res_r <- reactive({
      validate(need(input$run > 0, "Run a contrast in the DE / Pseudobulk panel, then click Run GSEA."))
      gsea_r()$res
    })
    table_r <- reactive({
      r <- res_r()
      r <- r[!is.na(r$padj) & r$padj <= (input$padj %||% 0.05), , drop = FALSE]
      if (identical(input$direction, "up")) r <- r[r$NES > 0, , drop = FALSE]
      if (identical(input$direction, "down")) r <- r[r$NES < 0, , drop = FALSE]
      r
    })
    observeEvent(res_r(), {                               # pathway picker: every tested set
      r <- res_r()
      updateSelectizeInput(session, "pathway", choices = r$pathway, selected = r$pathway[1],
                           server = TRUE)
    })
    shown <- function() {
      r <- table_r()
      le <- strsplit(r$leading_edge, ",", fixed = TRUE)
      data.frame(pathway = r$pathway, size = r$size, NES = round(r$NES, 2),
                 padj = signif(r$padj, 3),
                 leading_edge = vapply(le, function(g) paste0(paste(utils::head(g, 8), collapse = ", "),
                                                              if (length(g) > 8) sprintf(" (+%d)", length(g) - 8) else ""), ""),
                 stringsAsFactors = FALSE)
    }
    if (.scroll_has_dt()) {
      output$table <- DT::renderDataTable(
        DT::datatable(shown(), rownames = FALSE, selection = "single",
                      options = list(pageLength = 15, dom = "tip")))
      observeEvent(input$table_rows_selected, {           # a clicked row -> enrichment plot
        pw <- table_r()$pathway[input$table_rows_selected]
        updateSelectizeInput(session, "pathway", choices = res_r()$pathway, selected = pw, server = TRUE)
        updateSelectInput(session, "show", selected = "enrich")
      })
    } else output$table <- renderTable(shown())

    plot_r <- reactive({
      r <- res_r()
      p <- tryCatch(
        if (identical(input$show, "enrich")) .scroll_gsea_enrichment_plot(r, input$pathway, theme_r())
        else .scroll_gsea_top_plot(r, input$topn %||% 15, input$padj %||% 0.05,
                                   input$direction %||% "both", input$palette, theme_r()),
        error = function(e) e)
      validate(need(!inherits(p, "error"), if (inherits(p, "error")) conditionMessage(p)))
      .scroll_style_plot(p, style_r)
    })
    output$plot <- renderPlot(plot_r())
    .scroll_plot_downloads(output, plot_r, id)
    output$csv <- .scroll_csv_handler(reactive(res_r()), paste0("scroll_", id, ".csv"))
  })
}
