# The Clone-map panel: an interactive per-clone table over a greyed UMAP, where the
# clones selected in the table are highlighted (one colour per clone). A bespoke
# ui/server pair (not the declarative builder) because it needs two coordinated
# outputs — the table selection drives the plot. Gated on a repertoire (manifest$vdj).

# ---- data -------------------------------------------------------------------

# One row per clone from the per-cell repertoire table: size (cells) + a representative
# row for the descriptive columns (V/J genes, CDR3 lengths, group are constant within a
# clonotype by definition, so the first row is a faithful representative). Sorted by size.
.scroll_clone_table <- function(rc) {
  if (is.null(rc) || !"clone_id" %in% names(rc)) return(NULL)
  rc <- rc[!is.na(rc$clone_id), , drop = FALSE]
  if (!nrow(rc)) return(NULL)
  size <- as.data.frame(dplyr::count(rc, .data$clone_id, name = "size"))
  rep1 <- rc[!duplicated(rc$clone_id), , drop = FALSE]           # representative per clone
  tab <- dplyr::left_join(size, rep1, by = "clone_id")
  tab[order(-tab$size, tab$clone_id), , drop = FALSE]
}

# The columns to show in the table (clone id + size + whichever descriptors are baked).
.scroll_clone_disp_cols <- function(tab, segments) {
  keep <- c("clone_id", "size", "group", "antigen", unlist(segments),
            grep("_len$", names(tab), value = TRUE))
  intersect(keep, names(tab))
}

# ---- plot core --------------------------------------------------------------

# Greyed embedding with the selected clones' cells drawn on top, one colour per clone.
# Pure, so it renders identically on screen and on export. `sel` = selected clone ids.
.scroll_clone_map_plot <- function(df, emb, sel, palette = "Tableau 10", size = 0.8, aspect = 1) {
  df$.x <- df[[paste0(emb, "_1")]]; df$.y <- df[[paste0(emb, "_2")]]
  df$.hl <- ifelse(!is.na(df$clone) & df$clone %in% sel, as.character(df$clone), NA_character_)
  base <- df[is.na(df$.hl), , drop = FALSE]
  hi   <- df[!is.na(df$.hl), , drop = FALSE]
  p <- ggplot2::ggplot(df, ggplot2::aes(.data$.x, .data$.y)) +
    .scroll_point_layer(data = base, size = size, colour = "grey85",
                        raster = nrow(base) > .scroll_raster_threshold)
  if (nrow(hi)) {
    cols <- .scroll_discrete_colors(sel, palette)                # one colour per clone
    hi$.hl <- factor(hi$.hl, levels = sel)
    p <- p + ggplot2::geom_point(data = hi, ggplot2::aes(colour = .data$.hl),
                                 size = size + 0.7) +
      ggplot2::scale_colour_manual(values = cols, name = "Clone",
                                   guide = if (length(sel) > 24) "none" else "legend")
  }
  ttl <- if (length(sel)) paste(length(sel), "clone(s) highlighted",
                                sprintf("(%d cells)", nrow(hi)))
         else "Select clones in the table to highlight them"
  p <- p + ggplot2::labs(x = paste0(emb, "_1"), y = paste0(emb, "_2"), title = ttl) +
    .scroll_base_theme(legend = TRUE, axis_text = FALSE)
  # aspect 1 (default) keeps the embedding undistorted (coord_equal); otherwise honour
  # the requested panel aspect ratio.
  if (is.null(aspect) || abs(aspect - 1) < 1e-6) p + ggplot2::coord_equal()
  else p + ggplot2::theme(aspect.ratio = aspect)
}

# ---- the panel --------------------------------------------------------------

clone_map_ui <- function(id, data) {
  ns <- shiny::NS(id)
  m <- data$manifest
  if (is.null(m$vdj)) return(.scroll_empty_panel("No repertoire in this project."))
  emb <- .scroll_global_embeddings(m)
  if (!length(emb)) return(.scroll_empty_panel("No embedding to draw clones on."))
  ctl <- list(
    if (length(emb) > 1) selectInput(ns("embedding"), "Embedding", stats::setNames(emb, emb),
                                     selected = .scroll_default(data, "default_embedding", emb[[1]])),
    sliderInput(ns("minsize"), "Min clone size", 1, 50, 1, 1),
    selectInput(ns("palette"), "Highlight palette", names(.scroll_discrete_palettes)),
    sliderInput(ns("size"), "Point size", 0.2, 3, 0.8, 0.1),
    sliderInput(ns("aspect"), "Aspect ratio", 0.4, 3, 1, 0.1),
    actionButton(ns("clear"), "Clear selection", class = "btn-sm btn-outline-secondary"),
    tags$p(class = "scroll-desc", "Select clones in the table to colour their cells."))
  # table (top) + greyed UMAP (below); selection in the table drives the highlight
  table_ui <- if (.scroll_has_dt())
      .scroll_spin(DT::dataTableOutput(ns("table")))
    else
      selectizeInput(ns("clones"), "Clones", choices = NULL, multiple = TRUE,
                     width = "100%", options = list(placeholder = "pick clone ids",
                                                    plugins = list("remove_button")))
  bslib::layout_columns(
    col_widths = c(3, 9), class = "scroll-panel",
    div(class = "scroll-controls", do.call(.scroll_group, c(list("Controls"), ctl))),
    div(class = "scroll-plot",
        div(class = "scroll-table", table_ui),
        div(class = "scroll-plot-bar", .scroll_size_slider(ns),
            .scroll_dl_button(ns("png"), "PNG"), .scroll_dl_button(ns("pdf"), "PDF"),
            .scroll_dl_button(ns("csv"), "CSV")),
        .scroll_spin(plotOutput(ns("plot"), height = "480px"))))
}

clone_map_server <- function(id, data, cells_r = shiny::reactive(data$cells),
                             view_r = shiny::reactive(NULL), theme_r = shiny::reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    m <- data$manifest
    segs <- m$vdj$segments
    # scoped per-cell repertoire (honours the global filter / subset view), read once
    rc_r <- reactive(.scroll_vdj_scope(.scroll_vdj_read(data, "rep_cells.parquet"), cells_r()))
    # per-clone table, size-filtered
    tbl_r <- reactive({
      t <- .scroll_clone_table(rc_r())
      if (is.null(t)) return(NULL)
      t[t$size >= (input$minsize %||% 1), , drop = FALSE]
    })

    if (.scroll_has_dt()) {
      output$table <- DT::renderDataTable({
        t <- tbl_r(); validate(need(!is.null(t) && nrow(t), "No clones for this selection."))
        disp <- .scroll_clone_disp_cols(t, segs)
        dd <- t[, disp, drop = FALSE]
        # Every column gets a plain TEXT filter: numeric columns are stringified (else DT
        # renders a noUiSlider range filter that throws "$x.noUiSlider is not a function"
        # and flakes the session), with type="num" kept so they still sort numerically.
        # Factor columns are avoided too (they break DT's server-side filtering).
        num_idx <- unname(which(vapply(dd, is.numeric, logical(1)))) - 1L  # 0-based
        for (c in names(dd)[num_idx + 1L]) dd[[c]] <- as.character(dd[[c]])
        DT::datatable(dd, rownames = FALSE, selection = "multiple",
                      filter = "top",                       # a search input under each column
                      class = "compact stripe hover nowrap scroll-clone-dt",  # condensed rows
                      options = list(pageLength = 10, dom = "ftip", scrollX = TRUE,
                        order = list(list(1, "desc")),                   # size desc
                        columnDefs = list(
                          list(targets = num_idx, type = "num"),         # numeric sort on strings
                          # truncate the long composite clone id; full id on hover
                          list(targets = 0, render = DT::JS(
                            "function(d,t){return t==='display'&&d&&d.length>26 ?",
                            "'<span title=\"'+d+'\">'+d.substr(0,24)+'\\u2026'+'</span>' : d;}")))))
      }, server = TRUE)
      selected_r <- reactive({
        idx <- input$table_rows_selected; t <- tbl_r()
        if (!length(idx) || is.null(t)) character(0) else as.character(t$clone_id[idx])
      })
      observeEvent(input$clear, DT::selectRows(DT::dataTableProxy("table"), NULL))
    } else {
      observe({
        t <- tbl_r()
        updateSelectizeInput(session, "clones", server = TRUE,
                             choices = if (is.null(t)) character(0) else as.character(t$clone_id),
                             selected = isolate(input$clones))
      })
      selected_r <- reactive(as.character(input$clones %||% character(0)))
      observeEvent(input$clear, updateSelectizeInput(session, "clones", selected = character(0)))
    }

    plot_r <- reactive({
      cells <- cells_r()
      embs <- .scroll_global_embeddings(m)
      emb <- if (!is.null(input$embedding) && input$embedding %in% embs) input$embedding
             else .scroll_default(data, "default_embedding", embs[[1]])
      if (!emb %in% embs) emb <- embs[[1]]
      df <- .scroll_embedding_xy(cells, emb)
      validate(need(nrow(df), "No cells with coordinates."))
      rc <- rc_r()
      df$clone <- if (!is.null(rc) && "cell" %in% names(rc))
                    stats::setNames(rc$clone_id, rc$cell)[df$cell] else NA_character_
      .scroll_clone_map_plot(df, emb, selected_r(),
                             input$palette %||% "Tableau 10", input$size %||% 0.8,
                             input$aspect %||% 1) +
        .scroll_ggtheme(theme_r())
    })
    output$plot <- renderPlot(plot_r())
    .scroll_plot_downloads(output, plot_r, id)
    output$csv <- .scroll_csv_handler(reactive({
      t <- tbl_r(); if (is.null(t)) t else t[, .scroll_clone_disp_cols(t, segs), drop = FALSE]
    }), paste0("scroll_", id, ".csv"))
  })
}

# Assembled as a built-in gated on manifest$vdj (placed after the other VDJ panels).
.scroll_clone_map_panels <- function()
  list(list(id = "clone_map", label = "Clone map", title = "Clones on the embedding",
            desc = "Pick clones in the table to highlight their cells on the greyed embedding.",
            when = function(m) !is.null(m$vdj),
            ui = clone_map_ui, server = clone_map_server))
