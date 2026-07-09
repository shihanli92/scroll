# Two custom scroll panels for the mouse RNA + hashtag dataset:
#
#   qc_filter_*  — RNA QC sliders + histograms; publishes the set of QC-passing
#                  barcodes for the gating panel to consume.
#   hto_gate_*   — interactive manual hashtag gating: draw a lasso around a
#                  population on a biaxial hashtag-CLR plot, label it, and every
#                  cell inside becomes that population. Exports the barcodes to
#                  keep (in a kept gate AND passing QC) as CSV.
#
# Both follow the register_panel() contract (ui(id, data) / server(id, data,
# cells_r)). They are separate sections but cooperate: the export must reflect
# both steps, so they share one reactive through session$userData (which Shiny
# shares across every module in a session) — no changes to scroll's core.
#
# Hand-drawn gating is deliberate: with only 3 hashtags and a messy background,
# neither automated demux nor per-axis cutoffs separate the populations cleanly.

suppressMessages({
  library(shiny); library(bslib); library(plotly); library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

# The hashtag metadata columns the build wrote (hto_Hashtag-1, ...).
hto_cols_of <- function(manifest) grep("^hto_", names(manifest$meta), value = TRUE)

# Numeric [min, max] for a metadata column, from the manifest.
num_range <- function(manifest, col) {
  r <- manifest$meta[[col]]$range
  c(min = as.numeric(r$min), max = as.numeric(r$max))
}

# Seed the per-hashtag singlet gates from a saved keep-barcodes CSV (the app's own
# export: cell, sample, assignment). Returns a named list hashtag-label -> ids for
# the singlet labels; missing/invalid file -> all-empty. Doublet/Negative are not
# in the export, so they are re-derived by the prefilter.
load_saved_singlets <- function(path, labs) {
  out <- stats::setNames(vector("list", length(labs)), labs)
  if (is.null(path) || !file.exists(path)) return(out)
  df <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (!all(c("cell", "assignment") %in% names(df))) return(out)
  for (l in labs) out[[l]] <- as.character(df$cell[df$assignment == l])
  out
}

# Load a saved keep-barcodes CSV as a reference assignment (named vector
# cell -> assignment). Used to colour the gating plots by a previous demux.
load_ref_assignment <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  df <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (!all(c("cell", "assignment") %in% names(df))) return(NULL)
  stats::setNames(as.character(df$assignment), df$cell)
}

# Rough prefilter classification from a single CLR cutoff, so the biaxial plots
# are pre-coloured and the singlet/doublet structure is visible before gating. A
# cell positive (CLR >= thr) for one hashtag is that Singlet, for >=2 a Doublet,
# for none Negative. Pure -> unit-testable.
prefilter_classify <- function(cells, hto_cols, thr) {
  labs <- sub("^hto_", "", hto_cols)
  pos <- vapply(hto_cols, function(h) cells[[h]] >= thr, logical(nrow(cells)))
  if (is.null(dim(pos))) pos <- matrix(pos, ncol = length(hto_cols))
  n <- rowSums(pos)
  ifelse(n == 0, "Negative",
         ifelse(n >= 2, "Doublet", labs[max.col(pos, ties.method = "first")]))
}

# Base ("provisional") class for each cell: the saved CSV assignment where the
# cell is in it, else the prefilter class. This is what colours the biaxial plots
# and is the fallback the manual gates override.
provisional_class <- function(cells, hto_cols, thr, ref = NULL) {
  cls <- prefilter_classify(cells, hto_cols, thr)
  if (!is.null(ref)) {
    idx <- match(cells$cell, names(ref)); has <- !is.na(idx)
    cls[has] <- ref[idx[has]]
  }
  cls
}

# Classify cells from the per-hashtag single-positive gates plus any manually
# gated doublets. `singlets` is a named list hashtag-label -> cell ids gated
# single-positive for that hashtag; `doublet_ids` are cells hand-gated as doublets.
# A cell is a Doublet if manually gated doublet OR in >=2 singlet gates; a Singlet
# (its hashtag) if in exactly one gate; otherwise it falls back to `provisional`
# (the prefilter class) if given, else Negative. So manual gates override the
# prefilter. Pure -> unit-testable.
assign_hto <- function(cell_ids, singlets, doublet_ids = character(0), provisional = NULL) {
  labs <- names(singlets)
  memb <- vapply(labs, function(l) cell_ids %in% (singlets[[l]] %||% character(0)),
                 logical(length(cell_ids)))
  if (is.null(dim(memb))) memb <- matrix(memb, ncol = length(labs))
  n <- rowSums(memb)
  out <- ifelse(cell_ids %in% doublet_ids | n >= 2, "Doublet",
                ifelse(n == 1, labs[max.col(memb, ties.method = "first")], NA_character_))
  ung <- is.na(out)                                # ungated -> prefilter, else Negative
  out[ung] <- if (is.null(provisional)) "Negative" else provisional[ung]
  out
}

# Session-shared reactives, held in session$userData so the panels (separate
# modules) can cooperate: the QC panel publishes its passing barcodes, and the
# gating panel publishes its gate state, for the results panel to consume.
shared_qc_pass <- function(session) {
  if (is.null(session$userData$qc_pass))
    session$userData$qc_pass <- shiny::reactiveVal(NULL)
  session$userData$qc_pass
}
shared_hto_gates <- function(session) {
  if (is.null(session$userData$hto_gates))
    session$userData$hto_gates <-
      shiny::reactiveVal(list(singlets = list(), doublet = character(0), thr = 1))
  session$userData$hto_gates
}

# ============================ QC filter section ==============================

qc_filter_ui <- function(id, data) {
  ns <- NS(id)
  m  <- data$manifest
  nc <- num_range(m, "nCount_RNA"); nf <- num_range(m, "nFeature_RNA")
  mt <- num_range(m, "percent.mt")
  bslib::layout_columns(
    col_widths = c(3, 9), class = "scroll-panel",
    div(
      class = "scroll-controls",
      div(class = "scroll-cgroup", div(class = "scroll-cgroup-h", "RNA QC thresholds"),
          sliderInput(ns("ncount"), "nCount_RNA", min = floor(nc["min"]),
                      max = ceiling(nc["max"]), value = c(floor(nc["min"]), ceiling(nc["max"]))),
          sliderInput(ns("nfeature"), "nFeature_RNA", min = floor(nf["min"]),
                      max = ceiling(nf["max"]),
                      value = c(floor(nf["min"]), min(6000, ceiling(nf["max"])))),
          sliderInput(ns("pmt"), "percent.mt (max)", min = 0,
                      max = ceiling(mt["max"]), value = min(10, ceiling(mt["max"])))),
      div(class = "scroll-cgroup", div(class = "scroll-cgroup-h", "Pass"),
          div(class = "scroll-stat-v", style = "font-size:20px", textOutput(ns("passcount"))))
    ),
    div(class = "scroll-plot", plotOutput(ns("hist"), height = "320px"))
  )
}

qc_filter_server <- function(id, data, cells_r = reactive(data$cells)) {
  moduleServer(id, function(input, output, session) {
    qc_pass <- shared_qc_pass(session)

    passing <- reactive({
      cells <- cells_r()
      req(nrow(cells) > 0, input$ncount, input$pmt)
      cells$nCount_RNA   >= input$ncount[1]   & cells$nCount_RNA   <= input$ncount[2] &
      cells$nFeature_RNA >= input$nfeature[1] & cells$nFeature_RNA <= input$nfeature[2] &
      cells$percent.mt   <= input$pmt
    })

    # publish the QC-passing barcodes for the gating panel's export
    observe(qc_pass(cells_r()$cell[passing()]))

    output$passcount <- renderText({
      cells <- cells_r()
      sprintf("%s of %s", format(sum(passing()), big.mark = ","),
              format(nrow(cells), big.mark = ","))
    })

    output$hist <- renderPlot({
      cells <- cells_r()
      op <- graphics::par(mfrow = c(1, 3), mar = c(4, 4, 2, 1)); on.exit(graphics::par(op))
      graphics::hist(cells$nCount_RNA, breaks = 60, main = "nCount_RNA", xlab = "", col = "grey85", border = NA)
      graphics::abline(v = input$ncount, col = "#D64550", lwd = 2)
      graphics::hist(cells$nFeature_RNA, breaks = 60, main = "nFeature_RNA", xlab = "", col = "grey85", border = NA)
      graphics::abline(v = input$nfeature, col = "#D64550", lwd = 2)
      graphics::hist(cells$percent.mt, breaks = 60, main = "percent.mt", xlab = "", col = "grey85", border = NA)
      graphics::abline(v = input$pmt, col = "#D64550", lwd = 2)
    })
  })
}

# ============================ HTO gating section =============================

hto_gate_ui <- function(id, data) {
  ns <- NS(id)
  hto  <- hto_cols_of(data$manifest)
  labs <- sub("^hto_", "", hto)
  hto_max <- max(vapply(hto, function(h) num_range(data$manifest, h)["max"], numeric(1)))
  pairs <- utils::combn(hto, 2, simplify = FALSE)          # every hashtag pair
  plots <- lapply(seq_along(pairs), function(i) {
    pr <- sub("^hto_", "", pairs[[i]])
    div(class = "scroll-biax",
        tags$div(class = "scroll-kicker", paste(pr[1], "×", pr[2])),
        plotlyOutput(ns(paste0("gate", i)), height = "300px"))
  })
  keep_choices <- c(labs, "Doublet", "Negative")
  bslib::layout_columns(
    col_widths = c(3, 9), class = "scroll-panel",
    div(
      class = "scroll-controls",
      div(class = "scroll-cgroup", div(class = "scroll-cgroup-h", "Prefilter"),
          tags$p(class = "scroll-desc", style = "margin:0 0 8px",
                 paste("Plots are coloured by the saved demux (data/hto_keep_barcodes.csv)",
                       "where available, else by this rough CLR cutoff; your gates override both.")),
          sliderInput(ns("prethr"), "Positive if CLR ≥",
                      min = 0, max = ceiling(hto_max), value = 1, step = 0.05)),
      div(class = "scroll-cgroup", div(class = "scroll-cgroup-h", "Gates"),
          tags$p(class = "scroll-desc", style = "margin:0 0 8px",
                 paste("Lasso a single-positive cloud and assign it to its hashtag;",
                       "negatives are inferred. Doublets are inferred from double gating,",
                       "or lasso a doublet cloud and assign it to Doublet.")),
          selectInput(ns("target"), "Assign lasso to",
                      c(stats::setNames(hto, labs), Doublet = "Doublet"), width = "100%"),
          actionButton(ns("set"), "Set gate from lasso", class = "btn-primary", width = "100%"),
          actionButton(ns("clear"), "Clear all gates", width = "100%", style = "margin-top:6px"),
          div(style = "margin-top:10px", tableOutput(ns("gatestatus")))),
      div(class = "scroll-cgroup", div(class = "scroll-cgroup-h", "Export"),
          checkboxGroupInput(ns("keep"), "Keep classes", choices = keep_choices, selected = labs),
          downloadButton(ns("csv"), "Download kept barcodes (CSV)",
                         class = "btn-primary", style = "width:100%"))
    ),
    div(class = "scroll-plot",
        bslib::navset_tab(
          bslib::nav_panel(
            "Gates (lasso)",
            tags$p(class = "scroll-desc", style = "margin:6px 2px",
                   paste("Pick the lasso tool on any plot, draw around one hashtag's",
                         "single-positive cloud, choose that hashtag on the left, and",
                         "click “Set single-positive gate”. Repeat for each hashtag.")),
            plots),
          bslib::nav_panel("Summary", tableOutput(ns("summary")))))
  )
}

hto_gate_server <- function(id, data, cells_r = reactive(data$cells),
                            saved_csv = NULL,   # pass a path to seed gates from a saved export
                            ref_csv = "../../data/hto_keep_barcodes.csv") {  # colour plots by this demux
  moduleServer(id, function(input, output, session) {
    qc_pass <- shared_qc_pass(session)
    hto   <- hto_cols_of(data$manifest)
    labs  <- sub("^hto_", "", hto)
    pairs <- utils::combn(hto, 2, simplify = FALSE)
    ref    <- load_ref_assignment(ref_csv)                 # saved demux for colouring
    seeded <- load_saved_singlets(saved_csv, labs)         # saved selection, if any
    singlets <- reactiveVal(seeded)        # label -> ids (seeded from CSV on startup)
    doublet  <- reactiveVal(character(0))  # cells hand-gated as doublets
    last_sel <- reactiveVal(NULL)          # most recent lasso selection, any plot

    if (sum(lengths(seeded)) > 0)
      showNotification(sprintf("Loaded %s saved singlet assignments from %s.",
                               format(sum(lengths(seeded)), big.mark = ","),
                               basename(saved_csv)), type = "message", duration = 6)

    # publish the gate state + prefilter cutoff for the results panel
    hto_gates <- shared_hto_gates(session)
    observe(hto_gates(list(singlets = singlets(), doublet = doublet(), thr = input$prethr)))

    # cells (active subset) classified: manual gates override the base class
    # (saved-CSV demux where present, else the prefilter)
    gated <- reactive({
      cells <- cells_r(); req(nrow(cells) > 0, input$prethr)
      prov <- provisional_class(cells, hto, input$prethr, ref)
      cells$assignment <- assign_hto(cells$cell, singlets(), doublet(), prov)
      cells
    })

    kept <- reactive({
      a <- gated(); pass_ids <- qc_pass()
      qc_ok <- if (is.null(pass_ids)) TRUE else a$cell %in% pass_ids
      a[qc_ok & a$assignment %in% input$keep, , drop = FALSE]
    })

    # one biaxial per hashtag pair; each feeds the shared "last selection"
    lvl <- c(stats::setNames(scales::hue_pal()(length(labs)), labs),
             Doublet = "#D64550", Negative = "grey75")
    for (i in seq_along(pairs)) local({
      ii  <- i
      pr  <- pairs[[ii]]
      sid <- session$ns(paste0("gate", ii))
      output[[paste0("gate", ii)]] <- renderPlotly({
        a <- gated()
        plot_ly(source = sid) |>
          add_markers(x = a[[pr[1]]], y = a[[pr[2]]], customdata = a$cell,
                      color = factor(a$assignment, levels = names(lvl)),
                      colors = lvl, type = "scattergl",
                      marker = list(size = 3, opacity = 0.5), hoverinfo = "none") |>
          layout(dragmode = "lasso", showlegend = ii == 1,
                 xaxis = list(title = sub("^hto_", "", pr[1])),
                 yaxis = list(title = sub("^hto_", "", pr[2]))) |>
          event_register("plotly_selected")
      })
      observeEvent(event_data("plotly_selected", source = sid), {
        ev <- event_data("plotly_selected", source = sid)
        if (!is.null(ev) && nrow(ev)) last_sel(ev$customdata)
      })
    })

    observeEvent(input$set, {
      sel <- last_sel()
      validate(need(length(sel) > 0, "Draw a lasso selection first."))
      if (identical(input$target, "Doublet")) {
        doublet(union(doublet(), as.character(sel)))     # accumulate doublet regions
      } else {
        lbl <- sub("^hto_", "", input$target)
        sg <- singlets(); sg[[lbl]] <- as.character(sel); singlets(sg)
      }
    })
    observeEvent(input$clear, {
      singlets(stats::setNames(vector("list", length(labs)), labs))
      doublet(character(0))
    })

    output$gatestatus <- renderTable({
      sg <- singlets()
      data.frame(Gate = c(labs, "Doublet"),
                 Gated = c(vapply(labs, function(l) length(sg[[l]] %||% character(0)), integer(1)),
                           length(doublet())),
                 check.names = FALSE)
    }, digits = 0, align = "lr")

    output$summary <- renderTable({
      a <- gated(); pass_ids <- qc_pass()
      qc_ok <- if (is.null(pass_ids)) rep(TRUE, nrow(a)) else a$cell %in% pass_ids
      cls <- c(labs, "Doublet", "Negative")
      data.frame(
        Class = c(cls, "KEPT (export)"),
        Cells = c(vapply(cls, function(l) sum(a$assignment == l), integer(1)), nrow(kept())),
        `Pass QC` = c(vapply(cls, function(l) sum(qc_ok & a$assignment == l), integer(1)), NA_integer_),
        check.names = FALSE)
    }, digits = 0, align = "lrr", na = "")

    output$csv <- downloadHandler(
      filename = function() "hto_keep_barcodes.csv",
      content  = function(file)
        utils::write.csv(kept()[, c("cell", "sample", "assignment")], file, row.names = FALSE))
  })
}

# ====================== Final selection / ridgeplots =========================

# Bottom section: ridgeplots of each hashtag's CLR by the final class, and the
# filtered selection (singlets passing QC). Reads the gate state + QC-pass set the
# other two sections publish via session$userData, so it always reflects them.

result_ui <- function(id, data) {
  ns <- NS(id)
  bslib::layout_columns(
    col_widths = c(3, 9), class = "scroll-panel",
    div(
      class = "scroll-controls",
      div(class = "scroll-cgroup", div(class = "scroll-cgroup-h", "Show"),
          radioButtons(ns("show"), NULL,
                       c("Final kept (singlets, QC-pass)" = "kept",
                         "All cells (demux QC)" = "all"), selected = "kept"),
          tags$p(class = "scroll-desc", style = "margin-top:8px",
                 paste("Each ridge is a hashtag's CLR distribution per class.",
                       "Clean singlets sit high on their own hashtag, low on the others."))),
      div(class = "scroll-cgroup", div(class = "scroll-cgroup-h", "Final selection"),
          tableOutput(ns("summary")),
          downloadButton(ns("csv"), "Download filtered barcodes (CSV)",
                         class = "btn-primary", style = "width:100%"))
    ),
    div(class = "scroll-plot", plotOutput(ns("ridges"), height = "460px"))
  )
}

result_server <- function(id, data, cells_r = reactive(data$cells),
                          ref_csv = "../../data/hto_keep_barcodes.csv") {
  moduleServer(id, function(input, output, session) {
    qc_pass   <- shared_qc_pass(session)
    hto_gates <- shared_hto_gates(session)
    hto  <- hto_cols_of(data$manifest)
    labs <- sub("^hto_", "", hto)
    ref  <- load_ref_assignment(ref_csv)

    classified <- reactive({
      cells <- cells_r(); req(nrow(cells) > 0)
      g <- hto_gates()
      prov <- provisional_class(cells, hto, g$thr %||% 1, ref)
      cells$assignment <- factor(assign_hto(cells$cell, g$singlets, g$doublet, prov),
                                 levels = c(labs, "Doublet", "Negative"))
      pass <- qc_pass()
      cells$qc_pass <- if (is.null(pass)) TRUE else cells$cell %in% pass
      cells
    })
    kept <- reactive({
      a <- classified(); a[a$qc_pass & a$assignment %in% labs, , drop = FALSE]
    })

    output$ridges <- renderPlot({
      df <- if (identical(input$show, "kept")) kept() else classified()
      validate(need(nrow(df) > 0, "Set single-positive gates in the Hashtag gating section first."))
      long <- do.call(rbind, lapply(hto, function(h)
        data.frame(feature = sub("^hto_", "", h), clr = df[[h]],
                   assignment = df$assignment)))
      long <- long[!is.na(long$assignment), , drop = FALSE]
      ggplot(long, aes(x = clr, y = assignment, fill = assignment)) +
        ggridges::geom_density_ridges(scale = 1.8, alpha = 0.85,
                                      colour = "white", linewidth = 0.2) +
        facet_wrap(~ feature) +
        labs(x = "hashtag CLR", y = NULL) +
        theme_minimal(base_size = 13) +
        theme(legend.position = "none", panel.grid.minor = element_blank(),
              strip.text = element_text(face = "bold"))
    })

    output$summary <- renderTable({
      k <- kept(); asn <- as.character(k$assignment)
      data.frame(Class = c(labs, "Total kept"),
                 Cells = c(vapply(labs, function(l) sum(asn == l), integer(1)), nrow(k)))
    }, digits = 0, align = "lr")

    output$csv <- downloadHandler(
      filename = function() "hto_filtered_barcodes.csv",
      content  = function(file)
        utils::write.csv(kept()[, c("cell", "sample", "assignment")], file, row.names = FALSE))
  })
}
