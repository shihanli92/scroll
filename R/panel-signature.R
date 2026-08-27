# --- Signature panel ----------------------------------------------------------

# Live per-cell gene-signature (module) scoring from a typed/pasted gene list, shown on
# the embedding or as a violin by group. Scoring methods (see R/signature.R): "mean"
# (mean of the genes' log-norm expression), "scaled" (per-gene z over the shown cells),
# and "addmodulescore" (Seurat AddModuleScore -- signature mean minus an expression-
# matched control mean, using a once-cached full-store scan of per-gene means).

signature_ui <- function(id, data) {
  ns <- NS(id)
  m <- data$manifest
  assays <- .scroll_assays_of(m); cats <- .scroll_cat_cols(m); reds <- .scroll_view_embeddings(m, NULL)
  bslib::layout_columns(
    col_widths = c(3, 9), class = "scroll-panel",
    div(
      class = "scroll-controls",
      .scroll_group("Signature",
        selectizeInput(ns("sig"), "Genes", choices = NULL, multiple = TRUE,
                       options = list(placeholder = "Add genes, or paste a list...",
                                      maxOptions = 50, plugins = list("remove_button"))),
        .scroll_paste_handler(ns("sig"), ns("sig_paste")),
        selectInput(ns("method"), "Score",
                    c("Mean (log-norm)" = "mean", "Scaled (z-score, shown cells)" = "scaled",
                      "AddModuleScore (Seurat)" = "addmodulescore")),
        if (length(assays) > 1)
          selectInput(ns("assay"), "Assay", assays, selected = m$default_assay)),
      .scroll_group("View",
        selectInput(ns("view"), "Show as",
                    c("Feature UMAP" = "umap", "Violin by group" = "violin"))),
      conditionalPanel(
        sprintf("input['%s'] == 'umap'", ns("view")),
        .scroll_group("Embedding",
          selectInput(ns("reduction"), "Reduction", reds,
                      selected = .scroll_default(data, "default_embedding", reds[[1]]))),
        .scroll_group("Appearance",
          selectInput(ns("palette"), "Palette", .scroll_continuous_palettes, selected = "grey-purple"),
          sliderInput(ns("size"), "Point size", 0.1, 5, 0.7, 0.1),
          sliderInput(ns("clip"), "Color quantiles (%)", 0, 100, c(0, 100), 1),
          bslib::input_switch(ns("order"), "High-scoring cells on top", TRUE),
          bslib::input_switch(ns("legend"), "Legend", TRUE),
          bslib::input_switch(ns("raster"), "Rasterize (fast)", TRUE))),
      conditionalPanel(
        sprintf("input['%s'] == 'violin'", ns("view")),
        .scroll_group("Grouping",
          if (length(cats))
            selectInput(ns("group"), "Group by", stats::setNames(cats, cats), selected = cats[[1]])
          else helpText("No categorical column available for a violin."))),
      .scroll_group("Layout", .scroll_aspect_input(ns))
    ),
    .scroll_plot_area(ns, "560px", csv = TRUE)
  )
}

signature_server <- function(id, data, cells_r = reactive(data$cells),
                             view_r = reactive(NULL), theme_r = reactive(NULL)) {
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
    # keep the gene list valid for the active assay
    observeEvent(assay(), {
      feats <- .scroll_features_of(m, assay())
      updateSelectizeInput(session, "sig", choices = feats, server = TRUE,
                           selected = intersect(isolate(input$sig), feats))
    })
    # a delimited list pasted into the gene box (see .scroll_paste_handler) is parsed,
    # case-insensitively matched, added to the selection, and unknown tokens reported.
    observeEvent(input$sig_paste, {
      feats <- .scroll_features_of(m, assay())
      parsed <- .scroll_parse_gene_list(input$sig_paste, feats)
      req(length(parsed$ok) > 0 || length(parsed$missing) > 0)
      sel <- unique(c(isolate(input$sig), parsed$ok))
      updateSelectizeInput(session, "sig", choices = feats, server = TRUE, selected = sel)
      if (length(parsed$missing))
        showNotification(paste("Not in this assay:", paste(parsed$missing, collapse = ", ")),
                         type = "warning", duration = 6)
    })
    # DATA reactive: query the signature genes and build the per-cell score once (both
    # views reuse it). Method is a DATA input; palette/size/etc. are cosmetic.
    data_r <- reactive({
      cells <- cells_r()
      genes <- intersect(input$sig, .scroll_features_of(m, assay()))
      validate(need(length(genes) >= 1, "Add one or more genes to build a signature."))
      method <- input$method %||% "mean"
      # AddModuleScore needs each gene's mean over all cells: computed once (a full-store
      # scan, cached on the handle), then reused. The other methods touch only the genes.
      score <- if (identical(method, "addmodulescore"))
                 .scroll_module_score(cells, genes, data, assay())
               else
                 .scroll_signature_score(cells, genes, data$queryN(assay(), genes), method)
      list(cells = cells, genes = genes, n = nrow(cells), view = input$view %||% "umap",
           embedding = red_rv(), group_by = input$group,
           values = data.frame(cell = cells$cell, value = score, stringsAsFactors = FALSE))
    })
    cosmetic_r <- .scroll_cosmetic(reactive(
      list(theme = theme_r(), palette = input$palette, point_size = input$size,
           order = isTRUE(input$order), legend = isTRUE(input$legend),
           clip = input$clip / 100, aspect = input$aspect)))
    build <- function(raster) {
      d <- data_r(); st <- cosmetic_r(); st$raster <- raster
      lab <- sprintf("Signature (%d gene%s)", length(d$genes), if (length(d$genes) == 1) "" else "s")
      if (identical(d$view, "violin")) {
        validate(need(!is.null(d$group_by) && nzchar(d$group_by),
                      "Pick a categorical column to group the violin."))
        return(view_violin(d$cells, list(group_by = d$group_by, feature = lab), d$values, st))
      }
      view_feature_plot(d$cells, list(embedding = d$embedding, feature = lab), d$values, st)
    }
    plot_r   <- .scroll_lazy_plot(input, function() build(.scroll_use_raster(input$raster, data_r()$n)))
    export_r <- reactive(build(FALSE))
    output$plot <- renderPlot(plot_r())
    csv_r <- reactive({
      d <- data_r()
      data.frame(cell = d$cells$cell, signature_score = d$values$value, stringsAsFactors = FALSE)
    })
    .scroll_plot_downloads(output, export_r, id, csv_r = csv_r)
  })
}
