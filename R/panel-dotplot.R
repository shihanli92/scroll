# --- DotPlot panel ------------------------------------------------------------

# A <script> that lets a multi-select selectize accept a pasted, delimited gene
# list in the SAME box: when text containing a separator is pasted into the
# control, split it and hand it to Shiny (input `paste_id`) instead of dropping it
# in as one token. Single tokens paste normally through selectize. Delegated on
# document (capture) so it survives selectize re-rendering.
.scroll_paste_handler <- function(select_id, paste_id) {
  shiny::tags$script(shiny::HTML(sprintf("
(function(){
  document.addEventListener('paste', function(e){
    var sel = document.getElementById('%s'); if(!sel) return;
    var ctrl = sel.parentNode.querySelector('.selectize-control');
    if(!ctrl || !ctrl.contains(e.target)) return;
    var txt = (e.clipboardData || window.clipboardData).getData('text');
    if(!txt || !/[\\s,;]/.test(txt)) return;   // single token: let selectize handle it
    e.preventDefault();
    if(window.Shiny) Shiny.setInputValue('%s', txt, {priority:'event'});
  }, true);
})();", select_id, paste_id)))
}

# Parse a pasted gene list (space/comma/tab/newline separated) and match it to a
# feature set: exact first, then case-insensitively. Returns the matched features
# in pasted order (`ok`) plus the unmatched tokens (`missing`).
.scroll_parse_gene_list <- function(text, feats) {
  toks <- strsplit(text %||% "", "[[:space:],;]+", perl = TRUE)[[1]]
  toks <- unique(toks[nzchar(toks)])
  if (!length(toks)) return(list(ok = character(0), missing = character(0)))
  lut <- stats::setNames(feats, toupper(feats))
  matched <- ifelse(toks %in% feats, toks, unname(lut[toupper(toks)]))
  list(ok = unique(matched[!is.na(matched)]), missing = toks[is.na(matched)])
}

# Soft guidance only (does not block): a dot plot with very many gene rows is slow to
# compute (one partition read per gene, plus O(n^2) row hclust) and hard to read.
.SCROLL_DOTPLOT_WARN <- 80L

dotplot_ui <- function(id, data) {
  ns <- NS(id)
  m <- data$manifest
  assays <- .scroll_assays_of(m); cats <- .scroll_cat_cols(m)
  if (!length(cats)) return(.scroll_empty_panel("Needs a categorical metadata column (none found)."))
  bslib::layout_columns(
    col_widths = c(3, 9), class = "scroll-panel",
    div(
      class = "scroll-controls",
      .scroll_group("Genes",
        selectizeInput(ns("markers"), "Marker genes", choices = NULL, multiple = TRUE,
                       options = list(placeholder = "Add genes, or paste a list...",
                                      maxOptions = 50, plugins = list("remove_button"))),
        .scroll_paste_handler(ns("markers"), ns("marker_paste"))),
      .scroll_group("Grouping",
        selectInput(ns("group"), "Group by", stats::setNames(cats, cats), selected = cats[[1]]),
        if (length(assays) > 1)
          selectInput(ns("assay"), "Assay", assays, selected = m$default_assay)),
      .scroll_group("Appearance",
        # same genes x groups data as dots or as a tile heatmap
        selectInput(ns("display"), "Display", c("Dots" = "dots", "Heatmap (tiles)" = "tiles")),
        bslib::input_switch(ns("scale"), "Scale expression (z-score)", TRUE),
        conditionalPanel("input['display'] == 'tiles'", ns = ns,       # label rows past ~60 genes
          selectizeInput(ns("mark"), "Label genes", choices = NULL, multiple = TRUE,
                         options = list(placeholder = "Pick genes to label with leader lines",
                                        plugins = list("remove_button"))))),
      .scroll_group("Layout",
        selectInput(ns("cluster"), "Cluster (hclust)",
                    c("Off" = "off", "Rows" = "rows", "Columns" = "columns", "Both" = "both")))
    ),
    .scroll_plot_area(ns, csv = TRUE)
  )
}

# The look controls, shown in the panel's Style sheet (same input ids).
dotplot_style_ui <- function(id, data) {
  ns <- NS(id)
  tagList(
    selectInput(ns("palette"), "Palette", .scroll_continuous_palettes, selected = "magma"),
    conditionalPanel("input['display'] == 'dots'", ns = ns,
      sliderInput(ns("dotrange"), "Dot size", 0, 10, c(1, 6), 0.5)),
    conditionalPanel("input['display'] == 'tiles' && input['scale']", ns = ns,
      sliderInput(ns("clip"), "Clip z at \u00b1", 0.5, 5, 2.5, 0.5)),
    .scroll_aspect_input(ns))
}

dotplot_server <- function(id, data, cells_r = reactive(data$cells),
                           view_r = reactive(NULL), theme_r = reactive(NULL), style_r = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    m <- data$manifest
    .scroll_bind_view_cats(input, session, view_r, m, "group")
    assay <- reactive(input$assay %||% m$default_assay)
    defaults <- intersect(unlist(data$config$markers), .scroll_features_of(m, m$default_assay))
    # assay-repopulate + pasted gene list + long-list warning, and mirror the pick
    # into the "Label genes" (heatmap-tiles marks) selectize.
    .scroll_bind_gene_box(input, session, m, assay, "markers", paste_id = "marker_paste",
                          defaults = defaults, warn_n = .SCROLL_DOTPLOT_WARN, mark_id = "mark",
                          warn_msg = "%d genes selected - the dot plot may be slow to compute and hard to read.")
    # DATA reactive: query + aggregation + hclust. `scale` and `cluster` change
    # the aggregation/clustering, so they are DATA inputs (not cosmetic); palette
    # and dot size are cosmetic. No rasterization (dots = features x groups).
    data_r <- reactive({
      req(input$group)
      feats <- input$markers
      validate(need(length(feats) > 0, "Add one or more marker genes to build the panel."))
      cells <- cells_r(); el <- data$queryN(assay(), feats)
      display <- input$display %||% "dots"
      assembly <- if (identical(display, "tiles"))
        .scroll_heatmap_assemble(cells, feats, input$group, el, stat = "mean",
          scale = if (isTRUE(input$scale)) "zscore" else "none", cluster = input$cluster %||% "off")
      else .scroll_dotplot_assemble(cells, feats, input$group, el,
          scale = isTRUE(input$scale), cluster = input$cluster)
      list(cells = cells, group_by = input$group, features = feats,
           display = display, expr_long = el, assembly = assembly)
    })
    cosmetic_r <- .scroll_cosmetic(reactive(
      list(theme = theme_r(), palette = input$palette, dot_size = input$dotrange,
           clip = input$clip %||% 2.5, mark_genes = input$mark, aspect = input$aspect)))
    plot_r <- .scroll_lazy_plot(input, function() {
      d <- data_r(); st <- c(cosmetic_r(), list(style = style_r()))
      if (identical(d$display, "tiles"))
        view_heatmap(d$cells, list(group_by = d$group_by, features = d$features),
                     NULL, st, assembly = d$assembly)
      else view_dotplot(d$cells, list(group_by = d$group_by, features = d$features),
                        NULL, st, assembly = d$assembly)
    }, style_r)
    output$plot <- renderPlot(plot_r())
    csv_r <- reactive({ d <- data_r()
      if (identical(d$display, "tiles"))
        .scroll_heatmap_source(d$cells, d$features, d$group_by, d$expr_long)
      else .scroll_dotplot_source(d$assembly) })
    .scroll_plot_downloads(output, plot_r, id, csv_r = csv_r)
  })
}

