# The runtime: a polished bslib explorer with one scrolling section per analysis
# type, each with its own controls. The app never loads the Seurat object -- it
# reads the built artifacts (cells.parquet once globally; expression one feature
# at a time via duckdb), so runtime RAM stays flat.

# --- data handle --------------------------------------------------------------

# Load the artifacts once (shared across sessions of one app process).
.scroll_load <- function(dir) {
  con <- scroll_connect(dir)
  manifest <- scroll_manifest(dir)
  d <- list(
    dir = normalizePath(dir),
    cells = as.data.frame(arrow::read_parquet(file.path(dir, "cells.parquet"))),
    manifest = manifest,
    config = scroll_config(dir),
    con = con
  )
  # bound query helpers that dequantize to normalized units
  d$query1 <- function(assay, feature) {
    hit <- scroll_query_feature(con, assay, feature)
    if (nrow(hit)) hit$value <- scroll_dequantize(hit$value, manifest, assay)
    hit
  }
  d$queryN <- function(assay, features) {
    hit <- scroll_query_features(con, assay, features)
    if (nrow(hit)) hit$value <- scroll_dequantize(hit$value, manifest, assay)
    hit
  }
  d
}

.scroll_nz <- function(x) if (is.null(x) || !length(x) || !nzchar(x)) NULL else x

# --- manifest-derived control choices -----------------------------------------

.scroll_cat_cols <- function(m) names(Filter(function(x) identical(x$type, "categorical"), m$meta))
.scroll_num_cols <- function(m) names(Filter(function(x) identical(x$type, "numeric"), m$meta))
.scroll_reductions <- function(m) names(m$embeddings)
.scroll_assays_of <- function(m) names(m$assays)
.scroll_features_of <- function(m, assay) unlist(m$assays[[assay]]$features)

.scroll_colorby_choices <- function(m) {
  ch <- list()
  cats <- .scroll_cat_cols(m); nums <- .scroll_num_cols(m)
  if (length(cats)) ch[["Cell annotations"]] <- as.list(stats::setNames(cats, cats))
  if (length(nums)) ch[["Numeric / QC"]] <- as.list(stats::setNames(nums, nums))
  ch
}

.scroll_default <- function(data, key, fallback) {
  data$config[[key]] %||% data$manifest[[key]] %||% fallback
}

# --- DimPlot panel ------------------------------------------------------------

dimplot_ui <- function(id, data) {
  ns <- NS(id)
  m <- data$manifest
  reductions <- .scroll_reductions(m)
  cats <- .scroll_cat_cols(m)
  first_cat <- if (length(cats)) cats[[1]] else .scroll_num_cols(m)[[1]]
  bslib::layout_columns(
    col_widths = c(3, 9), class = "scroll-panel",
    div(
      class = "scroll-controls",
      .scroll_group("Embedding",
        selectInput(ns("reduction"), "Reduction", reductions,
                    selected = .scroll_default(data, "default_embedding", reductions[[1]])),
        selectInput(ns("colorby"), "Color by", .scroll_colorby_choices(m), selected = first_cat)),
      .scroll_group("Groups",
        selectizeInput(ns("highlight"), "Highlight", choices = NULL, multiple = TRUE,
                       options = list(placeholder = "All groups"))),
      .scroll_group("Appearance",
        selectInput(ns("palette"), "Palette", names(.scroll_discrete_palettes)),
        uiOutput(ns("manual")),
        sliderInput(ns("size"), "Point size", 0.1, 5, 0.6, 0.1),
        sliderInput(ns("alpha"), "Opacity", 0.1, 1, 0.85, 0.05),
        bslib::input_switch(ns("labels"), "Cluster labels", TRUE),
        bslib::input_switch(ns("legend"), "Legend", TRUE)),
      .scroll_group("Layout",
        selectInput(ns("split"), "Split by",
                    c("None" = "", stats::setNames(cats, cats))),
        .scroll_aspect_input(ns))
    ),
    div(class = "scroll-plot", plotOutput(ns("plot"), height = "460px"))
  )
}

dimplot_server <- function(id, data) {
  moduleServer(id, function(input, output, session) {
    m <- data$manifest
    is_cat <- reactive(identical(m$meta[[input$colorby]]$type, "categorical"))
    levels_of <- reactive(if (is_cat()) unlist(m$meta[[input$colorby]]$levels) else character(0))

    observeEvent(input$colorby, {
      if (is_cat()) {
        updateSelectInput(session, "palette",
                          choices = c(names(.scroll_discrete_palettes), "Manual"), selected = "Tableau 10")
        updateSelectizeInput(session, "highlight", choices = levels_of(), selected = character(0))
      } else {
        updateSelectInput(session, "palette",
                          choices = .scroll_continuous_palettes, selected = "viridis")
        updateSelectizeInput(session, "highlight", choices = character(0), selected = character(0))
      }
    })

    # per-group color pickers, shown only when the palette is "Manual"
    output$manual <- renderUI({
      if (!is_cat() || !identical(input$palette, "Manual")) return(NULL)
      lv <- levels_of()
      defaults <- .scroll_discrete_colors(lv, "Tableau 10")
      shiny::tagList(lapply(seq_along(lv), function(i)
        colourpicker::colourInput(session$ns(paste0("col_", i)), lv[i], value = defaults[[lv[i]]])))
    })
    manual_colors <- reactive({
      if (!is_cat() || !identical(input$palette, "Manual")) return(NULL)
      lv <- levels_of()
      vals <- lapply(seq_along(lv), function(i) input[[paste0("col_", i)]])
      names(vals) <- lv
      vals <- vals[!vapply(vals, is.null, logical(1))]
      if (length(vals)) unlist(vals) else NULL
    })

    output$plot <- renderPlot({
      req(input$reduction, input$colorby)
      params <- list(embedding = input$reduction, color_by = input$colorby)
      state <- list(palette = input$palette, point_size = input$size, alpha = input$alpha,
                    show_labels = isTRUE(input$labels), legend = isTRUE(input$legend),
                    split_by = .scroll_nz(input$split), aspect = input$aspect,
                    highlight = input$highlight, manual_colors = manual_colors())
      view_umap_colorby(data$cells, params, state)
    })
  })
}

# --- FeaturePlot panel --------------------------------------------------------

featureplot_ui <- function(id, data) {
  ns <- NS(id)
  m <- data$manifest
  assays <- .scroll_assays_of(m); cats <- .scroll_cat_cols(m)
  bslib::layout_columns(
    col_widths = c(3, 9), class = "scroll-panel",
    div(
      class = "scroll-controls",
      .scroll_group("Feature",
        selectizeInput(ns("feature"), "Gene", choices = NULL, multiple = FALSE,
                       options = list(placeholder = "Search a gene...", maxOptions = 50)),
        if (length(assays) > 1)
          selectInput(ns("assay"), "Assay", assays, selected = m$default_assay)),
      .scroll_group("Embedding",
        selectInput(ns("reduction"), "Reduction", .scroll_reductions(m),
                    selected = .scroll_default(data, "default_embedding", .scroll_reductions(m)[[1]]))),
      .scroll_group("Appearance",
        selectInput(ns("palette"), "Palette", .scroll_continuous_palettes, selected = "grey-purple"),
        sliderInput(ns("size"), "Point size", 0.1, 5, 0.7, 0.1),
        sliderInput(ns("clip"), "Color quantiles (%)", 0, 100, c(0, 100), 1),
        bslib::input_switch(ns("order"), "Expressing cells on top", TRUE),
        bslib::input_switch(ns("legend"), "Legend", TRUE)),
      .scroll_group("Layout",
        selectInput(ns("split"), "Split by", c("None" = "", stats::setNames(cats, cats))),
        .scroll_aspect_input(ns))
    ),
    div(class = "scroll-plot", plotOutput(ns("plot"), height = "460px"))
  )
}

featureplot_server <- function(id, data) {
  moduleServer(id, function(input, output, session) {
    m <- data$manifest
    assay <- reactive(input$assay %||% m$default_assay)
    observe(updateSelectizeInput(session, "feature", choices = .scroll_features_of(m, assay()),
                                 server = TRUE, selected = isolate(input$feature)))
    output$plot <- renderPlot({
      req(input$reduction)
      feat <- .scroll_nz(input$feature)
      validate(need(!is.null(feat), "Search for a gene to plot its expression."))
      params <- list(embedding = input$reduction, feature = feat)
      state <- list(palette = input$palette, point_size = input$size,
                    order = isTRUE(input$order), legend = isTRUE(input$legend),
                    clip = input$clip / 100, split_by = .scroll_nz(input$split),
                    aspect = input$aspect)
      view_feature_plot(data$cells, params, data$query1(assay(), feat), state)
    })
  })
}

# --- DotPlot panel ------------------------------------------------------------

dotplot_ui <- function(id, data) {
  ns <- NS(id)
  m <- data$manifest
  assays <- .scroll_assays_of(m); cats <- .scroll_cat_cols(m)
  bslib::layout_columns(
    col_widths = c(3, 9), class = "scroll-panel",
    div(
      class = "scroll-controls",
      .scroll_group("Genes",
        selectizeInput(ns("markers"), "Marker genes", choices = NULL, multiple = TRUE,
                       options = list(placeholder = "Add genes...", maxOptions = 50,
                                      plugins = list("remove_button")))),
      .scroll_group("Grouping",
        selectInput(ns("group"), "Group by", stats::setNames(cats, cats), selected = cats[[1]]),
        if (length(assays) > 1)
          selectInput(ns("assay"), "Assay", assays, selected = m$default_assay)),
      .scroll_group("Appearance",
        bslib::input_switch(ns("scale"), "Scale expression (z-score)", TRUE),
        selectInput(ns("palette"), "Palette", .scroll_continuous_palettes, selected = "magma"),
        sliderInput(ns("dotrange"), "Dot size", 0, 10, c(1, 6), 0.5)),
      .scroll_group("Layout",
        selectInput(ns("cluster"), "Cluster (hclust)",
                    c("Off" = "off", "Rows" = "rows", "Columns" = "columns", "Both" = "both")),
        .scroll_aspect_input(ns))
    ),
    div(class = "scroll-plot", plotOutput(ns("plot"), height = "520px"))
  )
}

dotplot_server <- function(id, data) {
  moduleServer(id, function(input, output, session) {
    m <- data$manifest
    assay <- reactive(input$assay %||% m$default_assay)
    defaults <- intersect(unlist(data$config$markers), .scroll_features_of(m, m$default_assay))
    observe(updateSelectizeInput(session, "markers", choices = .scroll_features_of(m, assay()),
                                 server = TRUE, selected = isolate(input$markers) %||% defaults))
    output$plot <- renderPlot({
      req(input$group)
      feats <- input$markers
      validate(need(length(feats) > 0, "Add one or more marker genes to build the panel."))
      params <- list(group_by = input$group, features = feats)
      state <- list(scale = isTRUE(input$scale), palette = input$palette,
                    dot_size = input$dotrange, cluster = input$cluster,
                    aspect = input$aspect)
      view_dotplot(data$cells, params, data$queryN(assay(), feats), state)
    })
  })
}

# --- Violin panel -------------------------------------------------------------

violin_ui <- function(id, data) {
  ns <- NS(id)
  m <- data$manifest
  assays <- .scroll_assays_of(m); cats <- .scroll_cat_cols(m)
  bslib::layout_columns(
    col_widths = c(3, 9), class = "scroll-panel",
    div(
      class = "scroll-controls",
      .scroll_group("Feature",
        selectizeInput(ns("feature"), "Gene", choices = NULL, multiple = FALSE,
                       options = list(placeholder = "Search a gene...", maxOptions = 50)),
        if (length(assays) > 1)
          selectInput(ns("assay"), "Assay", assays, selected = m$default_assay)),
      .scroll_group("Grouping",
        selectInput(ns("group"), "Group by", stats::setNames(cats, cats), selected = cats[[1]])),
      .scroll_group("Appearance",
        selectInput(ns("palette"), "Palette", names(.scroll_discrete_palettes)),
        bslib::input_switch(ns("jitter"), "Show points", FALSE),
        bslib::input_switch(ns("legend"), "Legend", FALSE)),
      .scroll_group("Layout", .scroll_aspect_input(ns))
    ),
    div(class = "scroll-plot", plotOutput(ns("plot"), height = "460px"))
  )
}

violin_server <- function(id, data) {
  moduleServer(id, function(input, output, session) {
    m <- data$manifest
    assay <- reactive(input$assay %||% m$default_assay)
    observe(updateSelectizeInput(session, "feature", choices = .scroll_features_of(m, assay()),
                                 server = TRUE, selected = isolate(input$feature)))
    output$plot <- renderPlot({
      req(input$group)
      feat <- .scroll_nz(input$feature)
      validate(need(!is.null(feat), "Search for a gene to plot its distribution."))
      params <- list(feature = feat, group_by = input$group)
      state <- list(palette = input$palette, jitter = isTRUE(input$jitter),
                    legend = isTRUE(input$legend), aspect = input$aspect)
      view_violin(data$cells, params, data$query1(assay(), feat), state)
    })
  })
}

# --- Proportions panel --------------------------------------------------------

proportions_ui <- function(id, data) {
  ns <- NS(id)
  cats <- .scroll_cat_cols(data$manifest)
  x_default <- if (length(cats) >= 2) cats[[2]] else cats[[1]]
  bslib::layout_columns(
    col_widths = c(3, 9), class = "scroll-panel",
    div(
      class = "scroll-controls",
      .scroll_group("Composition",
        selectInput(ns("group"), "Group by (x)", stats::setNames(cats, cats), selected = x_default),
        selectInput(ns("fill"), "Fill by", stats::setNames(cats, cats), selected = cats[[1]])),
      .scroll_group("Appearance",
        selectInput(ns("palette"), "Palette", names(.scroll_discrete_palettes)),
        bslib::input_switch(ns("normalize"), "Normalize to 100%", TRUE),
        bslib::input_switch(ns("legend"), "Legend", TRUE)),
      .scroll_group("Layout", .scroll_aspect_input(ns))
    ),
    div(class = "scroll-plot", plotOutput(ns("plot"), height = "460px"))
  )
}

proportions_server <- function(id, data) {
  moduleServer(id, function(input, output, session) {
    output$plot <- renderPlot({
      req(input$group, input$fill)
      params <- list(group_by = input$group, fill_by = input$fill)
      state <- list(palette = input$palette, normalize = isTRUE(input$normalize),
                    legend = isTRUE(input$legend), aspect = input$aspect)
      view_proportions(data$cells, params, state)
    })
  })
}

# --- section registry + shell -------------------------------------------------

# Panels available in v1, in scroll order. Each entry: display meta + its
# module's ui/server (the register_panel() extension point in embryo).
.scroll_panels <- function() list(
  list(id = "dimplot", num = "01", label = "DimPlot",
       title = "Cells by annotation",
       desc = "The embedding coloured by any cell metadata column.",
       ui = dimplot_ui, server = dimplot_server),
  list(id = "featureplot", num = "02", label = "FeaturePlot",
       title = "Gene expression",
       desc = "The embedding coloured by a gene's expression.",
       ui = featureplot_ui, server = featureplot_server),
  list(id = "dotplot", num = "03", label = "DotPlot",
       title = "Marker panel",
       desc = "Mean expression and fraction expressing across groups.",
       ui = dotplot_ui, server = dotplot_server),
  list(id = "violin", num = "04", label = "Violin",
       title = "Expression distribution",
       desc = "A gene's per-group expression distribution.",
       ui = violin_ui, server = violin_server),
  list(id = "proportions", num = "05", label = "Proportions",
       title = "Composition",
       desc = "Stacked composition of one annotation within another.",
       ui = proportions_ui, server = proportions_server)
)

.scroll_group <- function(title, ...) {
  div(class = "scroll-cgroup", div(class = "scroll-cgroup-h", title), ...)
}

# Shared aspect-ratio control, applied uniformly as a rendered-height multiplier
# (works for every panel including the aplot dendrogram composite, where
# theme(aspect.ratio) would detach the trees). Every panel adds the input and
# sets its renderPlot height via .scroll_plot_height(); future panels get it for
# free by doing the same.
.scroll_aspect_input <- function(ns) sliderInput(ns("aspect"), "Aspect ratio", 0.4, 3, 1, 0.1)

.scroll_stat <- function(value, label)
  div(class = "scroll-stat", span(class = "scroll-stat-v", value),
      span(class = "scroll-stat-l", label))

.scroll_appbar <- function(data, title) {
  m <- data$manifest
  assay <- m$default_assay
  brand <- list(span(class = "scroll-logo", "scroll"))
  if (!is.null(title))
    brand <- c(brand, list(span(class = "scroll-slash", "/"),
                           span(class = "scroll-dataset", title)))
  div(
    class = "scroll-appbar",
    div(class = "scroll-brand", brand),
    div(class = "scroll-stats",
        .scroll_stat(format(m$n_cells, big.mark = ","), "cells"),
        .scroll_stat(format(m$assays[[assay]]$n_features, big.mark = ","), "genes"),
        .scroll_stat(paste(.scroll_assays_of(m), collapse = ", "), "assays"),
        .scroll_stat(paste(.scroll_reductions(m), collapse = ", "), "reductions"))
  )
}

.scroll_rail <- function(panels) {
  tags$nav(
    class = "scroll-rail",
    lapply(panels, function(s) tags$a(
      class = "scroll-rail-item", href = paste0("#", s$id),
      span(class = "scroll-rail-num", s$num), span(s$label))),
    div(class = "scroll-rail-foot", "auto-generated from manifest.yaml")
  )
}

.scroll_section_card <- function(sec, data) {
  bslib::card(
    id = sec$id, class = "scroll-section", full_screen = FALSE,
    bslib::card_header(
      div(class = "scroll-eyebrow",
          span(class = "scroll-num", sec$num), span(class = "scroll-kicker", sec$label)),
      tags$h2(class = "scroll-title", sec$title),
      tags$p(class = "scroll-desc", sec$desc)),
    sec$ui(sec$id, data)
  )
}

.scroll_page <- function(data, title, panels) {
  bslib::page_fluid(
    theme = .scroll_theme(),
    tags$head(tags$style(HTML(.scroll_css())), tags$script(HTML(.scroll_spy_js()))),
    .scroll_appbar(data, title),
    div(
      class = "scroll-layout",
      .scroll_rail(panels),
      div(class = "scroll-content",
          lapply(panels, function(s) .scroll_section_card(s, data)))
    )
  )
}

.scroll_theme <- function() {
  bslib::bs_theme(
    version = 5,
    bg = "#FBFCFD", fg = "#111826", primary = "#2563A8",
    "border-color" = "#E4E8EE"
  )
}

#' Build the scroll explorer app
#'
#' Returns a `shiny::shinyApp` that reads a built project directory. Deploy the
#' scaffolded `app.R` (which calls this) to a Shiny Server, or run it locally
#' with [scroll_serve()].
#'
#' @param dir A built scroll project directory.
#' @return A `shiny.appobj`.
#' @export
scroll_app <- function(dir = ".") {
  data <- .scroll_load(dir)
  panels <- .scroll_panels()
  title <- .scroll_nz(data$config$title)

  ui <- .scroll_page(data, title, panels)
  server <- function(input, output, session) {
    for (sec in panels) sec$server(sec$id, data)
  }
  shiny::shinyApp(ui, server)
}

#' Run the scroll explorer locally
#'
#' @param dir A built scroll project directory.
#' @param ... Passed to [shiny::runApp()] (e.g. `port`, `host`, `launch.browser`).
#' @export
scroll_serve <- function(dir = ".", ...) {
  shiny::runApp(scroll_app(dir), ...)
}
