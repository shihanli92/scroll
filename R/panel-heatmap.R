# --- Heatmap panel (single-cell) ----------------------------------------------
# Genes x individual CELLS. Cells are randomly subsampled to a safe cap
# (proportional per group, deterministic), rasterized, and faceted by group so the
# render stays bounded at any cell count. The aggregated genes x groups heatmap
# lives on the DotPlot panel (Display: Heatmap tiles).

heatmap_ui <- function(id, data) {
  ns <- NS(id)
  m <- data$manifest
  assays <- .scroll_assays_of(m); cats <- .scroll_cat_cols(m)
  if (!length(cats)) return(.scroll_empty_panel("Needs a categorical metadata column (none found)."))
  ordcols <- c(cats, .scroll_num_cols(m))                # columns cells can be ordered by
  bslib::layout_columns(
    col_widths = c(3, 9), class = "scroll-panel",
    div(
      class = "scroll-controls",
      # the assembly (query + subsample + PC1) is expensive, so it recomputes only
      # on Compute -- not on every gene/scale/order tweak. Cosmetics stay live.
      bslib::input_task_button(ns("compute"), "Compute heatmap", type = "primary"),
      .scroll_group("Genes",
        # maxOptions must exceed the default gene count, else a large server-side
        # pre-selection (e.g. 100 top-variable genes) fails to bind to the input.
        selectizeInput(ns("markers"), "Genes", choices = NULL, multiple = TRUE,
                       options = list(placeholder = "Add genes, or paste a list...",
                                      maxOptions = 500, plugins = list("remove_button"))),
        .scroll_paste_handler(ns("markers"), ns("marker_paste"))),
      .scroll_group("Grouping",
        # multi-select and OPTIONAL: group cells by the "a | b" interaction of
        # several columns, or leave empty for a single ungrouped block.
        selectizeInput(ns("group"), "Group by (optional)", stats::setNames(cats, cats),
                       selected = cats[[1]], multiple = TRUE,
                       options = list(placeholder = "Leave empty for one block",
                                      plugins = list("remove_button"))),
        if (length(assays) > 1)
          selectInput(ns("assay"), "Assay", assays, selected = m$default_assay)),
      .scroll_group("Cells",
        sliderInput(ns("cellcap"), "Max cells shown", 500, 20000, 5000, 500),
        # order cells within each block: grouped only, by PC1, or by a column
        selectInput(ns("cellorder"), "Cell order",
                    c("Grouped" = "group", "By PC1 (similar together)" = "pc1",
                      "By metadata column" = "column")),
        conditionalPanel("input.cellorder == 'column'", ns = ns,
          selectInput(ns("ordercol"), "Order by column",
                      stats::setNames(ordcols, ordcols), selected = ordcols[[1]]))),
      .scroll_group("Values",
        bslib::input_switch(ns("scale"), "Scale per gene (z-score)", TRUE)),
      .scroll_group("Labels",
        # gene rows are unlabelled past ~60 genes; mark a chosen subset with
        # side labels + leader lines (ComplexHeatmap anno_mark style)
        selectizeInput(ns("mark"), "Label genes", choices = NULL, multiple = TRUE,
                       options = list(placeholder = "Pick genes to label with leader lines",
                                      plugins = list("remove_button")))),
      .scroll_group("Layout",
        # cells are never clustered (O(n^2)); only the gene rows get a dendrogram
        selectInput(ns("cluster"), "Cluster genes", c("Off" = "off", "Rows" = "rows")))
    ),
    .scroll_plot_area(ns, csv = TRUE)
  )
}

# The look controls, shown in the panel's Style sheet (same input ids).
heatmap_style_ui <- function(id, data) {
  ns <- NS(id)
  tagList(
    selectInput(ns("palette"), "Palette", .scroll_continuous_palettes, selected = "RdBu"),
    conditionalPanel("input['scale']", ns = ns,
      sliderInput(ns("clip"), "Clip z at \u00b1", 0.5, 5, 2.5, 0.5)),
    bslib::input_switch(ns("legend"), "Legend", TRUE),
    .scroll_aspect_input(ns))
}

heatmap_server <- function(id, data, cells_r = reactive(data$cells),
                           view_r = reactive(NULL), theme_r = reactive(NULL), style_r = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    m <- data$manifest
    # multi-select, view-aware group-by; optional (may stay empty for one block)
    .scroll_bind_view_cats(input, session, view_r, m, "group",
                           multiple = TRUE, allow_empty = TRUE)
    assay <- reactive(input$assay %||% m$default_assay)
    # default gene set: a project's `heatmap_markers:` (e.g. top variable genes)
    # when set, else the curated `markers:` shared with DotPlot.
    defaults <- intersect(unlist(data$config$heatmap_markers %||% data$config$markers),
                          .scroll_features_of(m, m$default_assay))
    .scroll_bind_gene_box(input, session, m, assay, "markers", paste_id = "marker_paste",
                          defaults = defaults, warn_n = 150L, mark_id = "mark",
                          warn_msg = "%d genes selected - the heatmap may be slow to compute.")
    # DATA reactive: query + subsample + assemble. scale/cluster/cap/order change the
    # matrix or ordering, so they are DATA inputs; palette/clip/legend are cosmetic.
    # Recompute trigger: bumped by Compute, and once automatically after the gene
    # selectize first populates (so the panel self-renders on first view without a
    # click, but never on later gene/scale/order tweaks).
    kick <- reactiveVal(0L)
    observeEvent(input$compute, kick(isolate(kick()) + 1L), ignoreInit = TRUE)
    # auto-render once, the first time the gene selectize actually has genes (the
    # server-side selectize arrives empty then filled, so guard on kick == 0 rather
    # than once = TRUE, which an empty first value would consume). Later gene edits
    # do NOT recompute -- they wait for Compute.
    observeEvent(input$markers,
                 if (isolate(kick()) == 0L && length(input$markers)) kick(1L),
                 ignoreNULL = FALSE)
    # Only the EXPENSIVE assembly (query + subsample + PC1) is gated on Compute (and
    # on a fresh cells_r() -- an app-bar filter/view change is a legitimate
    # re-trigger). eventReactive isolates its body, so the data controls (gene set,
    # scale/cluster/cap/order) are snapshotted here and re-run only on Compute.
    data_r <- eventReactive(list(kick(), cells_r()), {
      feats <- input$markers                             # group-by is optional (one block if empty)
      validate(need(length(feats) > 0, "Add one or more genes, then click Compute heatmap."))
      cells <- cells_r(); el <- data$queryN(assay(), feats)
      assembly <- .scroll_heatmap_cells_assemble(cells, feats, input$group, el,
        scale = if (isTRUE(input$scale)) "zscore" else "none",
        cluster = input$cluster %||% "off", cap = input$cellcap %||% 5000,
        cell_order = input$cellorder %||% "group", order_col = input$ordercol)
      list(cells = cells, group_by = input$group, features = feats,
           expr_long = el, assembly = assembly)
    }, ignoreNULL = FALSE)
    # Cheap cosmetics stay live: they only re-skin the memoized assembly (no
    # re-query, no subsample/PC1), so palette/clip/legend/label/aspect redraw at once.
    cosmetic_r <- .scroll_cosmetic(reactive(
      list(theme = theme_r(), palette = input$palette, clip = input$clip %||% 2.5,
           mark_genes = input$mark, legend = isTRUE(input$legend), aspect = input$aspect)))
    plot_r <- .scroll_lazy_plot(input, function() {
      d <- data_r()
      view_heatmap(d$cells, list(group_by = d$group_by, features = d$features),
                   NULL, c(cosmetic_r(), list(style = style_r())), assembly = d$assembly)
    }, style_r)
    output$plot <- renderPlot(plot_r())
    csv_r <- reactive({ d <- data_r()
      .scroll_heatmap_source(d$cells, d$features, d$group_by, d$expr_long) })
    .scroll_plot_downloads(output, plot_r, id, csv_r = csv_r)
  })
}
