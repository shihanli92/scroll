# --- Biaxial panel ------------------------------------------------------------

# Prefer hashtag / antibody CLR columns for the default axes, else the first few
# numerics — so an HTO or CITE-seq dataset opens on its biaxial signal plots.
.scroll_default_biaxial <- function(nums) {
  hit <- grep("^hto_|hashtag|adt_|_adt$|^ab_", nums, value = TRUE, ignore.case = TRUE)
  if (length(hit) >= 2) hit else utils::head(nums, 3)
}

biaxial_ui <- function(id, data) {
  ns <- NS(id)
  m <- data$manifest
  nums <- .scroll_num_cols(m); cats <- .scroll_cat_cols(m)
  if (length(nums) < 2)
    return(.scroll_empty_panel("Needs at least two numeric metadata columns (none found)."))
  if (!length(cats))
    return(.scroll_empty_panel("Needs a categorical column to colour by (none found)."))
  color_default <- if ("hto" %in% cats) "hto" else cats[[1]]
  bslib::layout_columns(
    col_widths = c(3, 9), class = "scroll-panel",
    div(
      class = "scroll-controls",
      .scroll_group("Axes",
        selectizeInput(ns("features"), "Numeric columns", choices = nums,
                       selected = .scroll_default_biaxial(nums), multiple = TRUE,
                       options = list(placeholder = "Pick 2+ numeric columns"))),
      .scroll_group("Colour",
        selectInput(ns("colorby"), "Colour by", stats::setNames(cats, cats),
                    selected = color_default)),
      .scroll_group("Appearance",
        selectInput(ns("palette"), "Palette", .scroll_cat_palettes()),
        uiOutput(ns("manual")),
        sliderInput(ns("size"), "Point size", 0.1, 3, 0.5, 0.1),
        sliderInput(ns("alpha"), "Opacity", 0.1, 1, 0.6, 0.05),
        bslib::input_switch(ns("legend"), "Legend", TRUE),
        bslib::input_switch(ns("raster"), "Rasterize (fast)", TRUE)),
      .scroll_group("Layout", .scroll_aspect_input(ns))
    ),
    .scroll_plot_area(ns, "460px")
  )
}

biaxial_server <- function(id, data, cells_r = reactive(data$cells),
                           view_r = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    m <- data$manifest
    .scroll_bind_view_cats(input, session, view_r, m, "colorby")
    # per-level color pickers when the palette is "Manual" (levels of Colour by)
    lvl_r <- reactive({ req(input$colorby); unlist(m$meta[[input$colorby]]$levels) })
    output$manual <- renderUI(
      if (identical(input$palette, "Manual")) .scroll_manual_ui(session$ns, lvl_r()))
    manual_colors <- reactive(
      if (identical(input$palette, "Manual")) .scroll_manual_colors(input, lvl_r()))
    # DATA reactive: build the (potentially large) pairwise long df once; cosmetic
    # drags no longer re-expand it. `params` is carried for the colour label.
    data_r <- reactive({
      req(input$colorby)
      feats <- input$features
      validate(need(length(feats) >= 2, "Pick at least two numeric columns."))
      params <- list(features = feats, color_by = input$colorby)
      list(params = params, df = .scroll_biaxial_df(cells_r(), params))
    })
    cosmetic_r <- .scroll_cosmetic(reactive(
      list(palette = input$palette, point_size = input$size,
           alpha = input$alpha, legend = isTRUE(input$legend), aspect = input$aspect,
           manual_colors = manual_colors())))
    build <- function(raster) {
      d <- data_r(); st <- cosmetic_r(); st$raster <- raster
      view_biaxial(NULL, d$params, st, df = d$df)
    }
    plot_r   <- reactive(build(.scroll_use_raster(input$raster, nrow(data_r()$df))))
    export_r <- reactive(build(FALSE))
    output$plot <- renderPlot(plot_r())
    .scroll_plot_downloads(output, export_r, id)
  })
}

