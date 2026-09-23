# App-layer logic: the data handle and the DimPlot module server (no browser).

test_that(".scroll_load exposes cells, manifest, config and query helpers", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  expect_true(all(c("cells", "manifest", "config", "con", "query1", "queryN") %in% names(data)))
  hit <- data$query1("RNA", "MS4A1")
  expect_true(all(c("cell", "value") %in% names(hit)))
})

test_that("scroll_multi_app mounts several projects into one namespaced app", {
  p <- test_project()
  app <- scroll_multi_app(c("First" = p, "Second" = p))
  expect_s3_class(app, "shiny.appobj")     # both projects mount, namespaced (ds1-*/ds2-*)
  expect_error(scroll_multi_app(list()), "at least one project")
})

test_that("scroll_app returns a shiny app and handle cleanup is a safe no-op", {
  expect_s3_class(scroll_app(test_project()), "shiny.appobj")
  # scroll_disconnect just drops the handle's cached datasets (no live db
  # connection), so it is idempotent and silent — onStop can't crash the session.
  data <- scroll:::.scroll_load(test_project())
  expect_silent(scroll_disconnect(data$con))
  expect_silent(scroll_disconnect(data$con))     # second call is a no-op
})

test_that("dimplot_server renders a plot from its controls", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  shiny::testServer(scroll:::dimplot_server, args = list(data = data), {
    session$setInputs(reduction = "umap", colorby = "celltype", palette = "Tableau 10",
                      size = 0.6, alpha = 0.85, labels = TRUE, split = "")
    expect_false(is.null(output$plot))
    # switching to a numeric color-by swaps the palette family (no error)
    session$setInputs(colorby = "nCount_RNA")
    expect_false(is.null(output$plot))
  })
})

test_that("grouping panels degrade gracefully with no categorical metadata", {
  fake <- list(manifest = list(
    meta = list(nCount = list(type = "numeric", range = list(min = 0, max = 1))),
    assays = list(RNA = list(features = list("A", "B"), max = 1, n_features = 2)),
    embeddings = list(umap = list(dims = 2)),
    default_assay = "RNA", default_embedding = "umap"))
  for (ui in list(scroll:::dotplot_ui, scroll:::heatmap_ui, scroll:::violin_ui,
                  scroll:::proportions_ui, scroll:::de_ui)) {
    tag <- ui("p", fake)
    expect_s3_class(tag, "shiny.tag")
    expect_true(grepl("categorical metadata", as.character(tag)))
  }
})

test_that(".scroll_subset_cells filters to selected levels (and is a no-op otherwise)", {
  cells <- data.frame(cell = paste0("c", 1:6), grp = rep(c("A", "B", "C"), 2),
                      stringsAsFactors = FALSE)
  expect_equal(nrow(scroll:::.scroll_subset_cells(cells, "grp", c("A", "B"))), 4)
  expect_equal(unique(scroll:::.scroll_subset_cells(cells, "grp", "C")$grp), "C")
  # no column / no values / unknown column -> unchanged
  expect_identical(scroll:::.scroll_subset_cells(cells, NULL, NULL), cells)
  expect_identical(scroll:::.scroll_subset_cells(cells, "grp", character(0)), cells)
  expect_identical(scroll:::.scroll_subset_cells(cells, "missing", "A"), cells)
})

test_that("a panel server honors a subsetted cells_r (global filter)", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  ct <- sort(unique(as.character(data$cells$celltype)))[[1]]
  sub <- reactive(scroll:::.scroll_subset_cells(data$cells, "celltype", ct))
  shiny::testServer(scroll:::dimplot_server, args = list(data = data, cells_r = sub), {
    session$setInputs(reduction = "umap", colorby = "celltype", palette = "Tableau 10",
                      size = 0.6, alpha = 0.85, labels = TRUE, legend = TRUE, split = "",
                      aspect = 1, highlight = character(0))
    # the plot reactive draws only the subset's cells
    expect_equal(nrow(plot_r()$data), sum(as.character(data$cells$celltype) == ct))
  })
})

test_that("violin_server and proportions_server render from controls", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  shiny::testServer(scroll:::violin_server, args = list(data = data), {
    session$setInputs(feature = "CD3D", group = "celltype", palette = "Tableau 10",
                      jitter = FALSE, legend = FALSE)
    expect_false(is.null(output$plot))
  })
  shiny::testServer(scroll:::proportions_server, args = list(data = data), {
    session$setInputs(group = "condition", fill = "celltype", palette = "Tableau 10",
                      normalize = TRUE, legend = TRUE)
    expect_false(is.null(output$plot))
  })
  # Heatmap panel is single-cell (genes x subsampled cells) only; assembly is
  # gated on Compute (auto-fires once), so data-input changes need a Compute click.
  shiny::testServer(scroll:::heatmap_server, args = list(data = data), {
    session$setInputs(markers = c("CD3D", "CD8A"), group = "celltype", scale = TRUE,
                      clip = 2.5, cluster = "off", cellcap = 40, cellorder = "group",
                      palette = "RdBu", legend = TRUE, aspect = 1, compute = 1)
    expect_true(inherits(plot_r(), c("ggplot", "aplot")))
    expect_s3_class(csv_r(), "data.frame")
    session$setInputs(cellorder = "pc1", compute = 2)    # recompute with PC1 ordering
    expect_true(inherits(plot_r(), c("ggplot", "aplot")))
    session$setInputs(group = character(0), compute = 3) # ungrouped (one block)
    expect_true(inherits(plot_r(), c("ggplot", "aplot")))
  })
  # DotPlot's Heatmap-tiles display renders the aggregated genes x groups heatmap
  shiny::testServer(scroll:::dotplot_server, args = list(data = data), {
    session$setInputs(markers = c("CD3D", "CD8A"), group = "celltype", display = "tiles",
                      scale = TRUE, clip = 2.5, cluster = "off", palette = "RdBu", aspect = 1)
    expect_true(inherits(plot_r(), c("ggplot", "aplot")))
    session$setInputs(display = "dots", dotrange = c(1, 6))
    expect_true(inherits(plot_r(), c("ggplot", "aplot")))
  })
})

test_that("categorical panels expose Manual per-level colour pickers", {
  data <- scroll:::.scroll_load(test_project())
  on.exit(scroll_disconnect(data$con))
  shiny::testServer(scroll:::violin_server, args = list(data = data), {
    session$setInputs(group = "celltype", feature = "CD3D", palette = "Manual",
                      jitter = FALSE, legend = FALSE, aspect = 1)
    session$flushReact()
    lv <- lvl_r()
    expect_gt(length(lv), 0)
    expect_false(is.null(output$manual))                 # colour pickers rendered
    session$setInputs(col_1 = "#010203")                 # set first level's colour
    session$flushReact()
    expect_equal(unname(manual_colors()[[lv[[1]]]]), "#010203")
  })
})

test_that("base CSS suppresses scrollbar-toggle resize loops (both axes)", {
  css <- scroll:::.scroll_css()
  # Reserve the vertical gutter so a toggling vertical bar can't change width,
  # and hide page-level horizontal overflow so a toggling horizontal bar can't
  # change height -- either would fire window 'resize' and re-render every plot.
  expect_match(css, "overflow-y:scroll")
  expect_match(css, "scrollbar-gutter:stable")
  expect_match(css, "overflow-x:hidden")
  # .scroll-plot pins min-width:0, and the plot output snaps to a whole pixel
  # (round(down, 100%, 1px)) so a fractional grid width can't drive a HiDPI
  # (devicePixelRatio 2) sub-pixel ResizeObserver re-render loop.
  expect_match(css, "\\.scroll-plot\\{[^}]*min-width:0")
  expect_match(css, "round\\(down, 100%, 1px\\)")
  # The horizontal clip must live on <html> only. overflow-x:hidden on <body>
  # forces its computed overflow-y to `auto`, making body a scroll container that
  # breaks position:sticky on the app bar + section rail (they scroll away).
  expect_no_match(css, "body\\{[^}]*overflow-x")
  expect_match(css, "\\.scroll-rail\\{position:sticky")
})

test_that("the panel layout is a container query: clamp side-by-side, stack when narrow", {
  css <- scroll:::.scroll_css()
  # a clamped width knob (usable min on laptops, capped on big monitors) ...
  expect_match(css, "--sc-ctl-w:clamp\\(")
  # ... driving a container query on each CARD (so a half-width two-up card adapts)
  expect_match(css, "container-type:inline-size")
  expect_match(css, "@container sc-card \\(min-width:720px\\)")
  expect_match(css, "grid-template-columns:var\\(--sc-ctl-w\\) minmax\\(0,1fr\\)")
  # narrow cards stack controls above a full-width plot
  expect_match(css, "@container sc-card \\(max-width:719\\.98px\\)")
})

test_that("two-up mode packs two panel cards per row (auto-fit, container-scoped)", {
  css <- scroll:::.scroll_css()
  js  <- scroll:::.scroll_spy_js()
  expect_match(css, "\\.scroll-content\\.two-up\\{grid-template-columns:repeat\\(auto-fit")
  expect_match(css, "\\.scroll-layout:has\\(\\.two-up\\)\\{--sc-content-max:2200px")
  expect_match(css, "container-name:sc-card")
  expect_match(js, "scrollToggleTwoUp")
  expect_match(js, "scroll:two-up")            # persisted in localStorage
  # the toggle button is emitted into the app bar
  data <- scroll:::.scroll_load(test_project()); on.exit(scroll_disconnect(data$con))
  bar <- as.character(scroll:::.scroll_appbar(data, "t"))
  expect_match(bar, "scroll-two-toggle")
  expect_match(bar, "scrollToggleTwoUp")
})

test_that("the app bar is title-led: cells readout + dataset popover, no stats strip", {
  data <- scroll:::.scroll_load(test_project()); on.exit(scroll_disconnect(data$con))
  bar <- as.character(scroll:::.scroll_appbar(data, "My dataset"))
  expect_match(bar, "scroll-dataset")                    # dataset title present
  expect_match(bar, "scroll-cells")                      # live cell readout beside the pills
  expect_match(bar, "scroll-appbar-actions")             # the right-hugging action cluster
  expect_match(bar, "scroll-info-dl")                    # dataset facts moved to a popover
  expect_match(bar, "fa-sliders")                        # real icons, not glyphs
  expect_no_match(bar, "scroll-stats")                   # the old stats strip is gone
  expect_no_match(bar, "scroll-version")                 # version badge left the bar
})

test_that("app-bar buttons carry a visible hover label saying what they do", {
  data <- scroll:::.scroll_load(test_project()); on.exit(scroll_disconnect(data$con))
  bar <- as.character(scroll:::.scroll_appbar(data, "t"))
  tips <- scroll:::.SCROLL_TIPS
  # each button starts with its label (entity-escaped '&' in the markup) ...
  for (k in c("info", "two_on", "ctl_hide"))
    expect_match(bar, gsub("&", "&amp;", tips[[k]], fixed = TRUE), fixed = TRUE)
  expect_equal(lengths(regmatches(bar, gregexpr("data-tip=", bar))), 3L)
  # ... and none keeps a native title (it would double up with the styled label)
  expect_no_match(bar, "<button[^>]*scroll-(info|two|ctl)-toggle[^>]*title=")
  # the JS knows every state's label and recomputes it from live state on hover/focus
  js <- scroll:::.scroll_spy_js()
  expect_match(js, "window.SCROLL_TIPS=", fixed = TRUE)
  for (k in names(tips)) expect_match(js, tips[[k]], fixed = TRUE)
  expect_match(js, "function scrollTip(btn)", fixed = TRUE)
  expect_match(js, "'pointerover','focusin'", fixed = TRUE)
  css <- scroll:::.scroll_css()
  expect_match(css, "[data-tip]::after{content:attr(data-tip)", fixed = TRUE)
  expect_match(css, ":focus-visible::after", fixed = TRUE)   # keyboard users get it too
})

test_that("rail items are draggable to reorder panels, persisted and resettable", {
  js  <- scroll:::.scroll_spy_js()
  expect_match(js, "dragstart"); expect_match(js, "scrollResetOrder")
  expect_match(js, "scroll-order:")            # localStorage key by panel-id set
  rail <- as.character(scroll:::.scroll_rail(
    list(list(id = "dimplot", num = "01", label = "DimPlot"),
         list(id = "featureplot", num = "02", label = "FeaturePlot"))))
  expect_match(rail, "draggable=\"true\"")
  expect_match(rail, "scroll-rail-reset")
  expect_match(rail, "scrollResetOrder")
})

test_that("config.yaml layout: overrides the layout CSS vars and the two-up default", {
  # numbers -> px, strings verbatim; nothing set -> empty (appended to base CSS)
  expect_identical(scroll:::.scroll_layout_css(list()), "")
  st <- scroll:::.scroll_layout_css(
    list(layout = list(content_max = 2000, rail_width = 180,
                       control_width = "clamp(220px, 22%, 300px)")))
  expect_match(st, "--sc-content-max:2000px")
  expect_match(st, "--sc-rail-w:180px")
  expect_match(st, "--sc-ctl-w:clamp\\(220px, 22%, 300px\\)")
  # two_up default reaches the content div + the app-bar toggle
  data <- scroll:::.scroll_load(test_project()); on.exit(scroll_disconnect(data$con))
  data$config$layout <- list(two_up = TRUE)
  body <- as.character(scroll:::.scroll_body(data, "t", scroll:::.scroll_builtin_panels()[1]))
  expect_match(body, "scroll-content two-up")
  bar <- as.character(scroll:::.scroll_appbar(data, "t"))
  expect_match(bar, "scroll-two-toggle is-on")
  expect_match(bar, "aria-pressed=\"true\"")
})

test_that("wide screens go full-bleed with a capped content width (no fixed max-width)", {
  css <- scroll:::.scroll_css()
  expect_match(css, "--sc-content-max:1600px")
  expect_match(css, "--sc-rail-w:160px")
  # content track is capped (not a fixed centered max-width block) and edges are hugged
  expect_match(css, "grid-template-columns:var\\(--sc-rail-w\\) minmax\\(0,var\\(--sc-content-max\\)\\)")
  expect_match(css, "justify-content:space-between")
  # the old fixed centered cap is gone
  expect_no_match(css, "max-width:1320px")
  expect_no_match(css, "max-width:1620px")
})

test_that("below 1400 the filters column is an off-canvas drawer with a live app-bar offset", {
  css <- scroll:::.scroll_css()
  js  <- scroll:::.scroll_spy_js()
  # the right rail becomes a fixed, transform-hidden drawer, revealed by .filters-open
  expect_match(css, "@media \\(max-width:1399\\.98px\\)")
  expect_match(css, "\\.scroll-layout\\.filters-open>\\.scroll-filters\\{transform:none")
  # sticky offsets are driven by the measured app-bar bottom (survives wrapping)
  expect_match(css, "var\\(--sc-appbar-h,70px\\)")
  # the toggle opens the drawer on narrow viewports (matchMedia branch) and the JS
  # keeps --sc-appbar-h in sync + closes the drawer on Escape / outside click
  expect_match(js, "matchMedia\\('\\(max-width:1399\\.98px\\)'\\)")
  expect_match(js, "filters-open")
  expect_match(js, "--sc-appbar-h")
  expect_match(js, "scrollCloseDrawers")
})

test_that("a per-card Controls toggle collapses the controls when a card is stacked", {
  css <- scroll:::.scroll_css()
  js  <- scroll:::.scroll_spy_js()
  # button is display:none until the card is narrow enough to stack
  expect_match(css, "\\.scroll-ctl-btn\\{display:none")
  expect_match(css, "\\.scroll-panel-card:has\\(\\.scroll-panel\\.bslib-grid\\) \\.scroll-ctl-btn\\{display:inline-flex")
  expect_match(css, "\\.scroll-ctl-hidden \\.scroll-panel\\.bslib-grid>\\.bslib-grid-item:first-child\\{display:none")
  expect_match(js, "scrollTogglePanelControls")
  # the button is emitted into every panel card header
  data <- scroll:::.scroll_load(test_project()); on.exit(scroll_disconnect(data$con))
  card <- as.character(scroll:::.scroll_panel_card(
    list(id = "dimplot", num = "01", label = "DimPlot", title = "t", desc = "d",
         ui = scroll:::dimplot_ui), data))
  expect_match(card, "scrollTogglePanelControls")
  expect_match(card, "scroll-ctl-btn")
})

test_that("plot toolbars carry a download-scale slider wired to the dl_scale input", {
  html <- as.character(scroll:::.scroll_plot_area(identity))
  expect_match(html, "scroll-size")                     # the slider is in the toolbar
  expect_match(html, "dl_scale")                        # targets the module's dl_scale input
  expect_match(html, "Shiny.setInputValue")             # writes the value to Shiny
  expect_match(html, "scroll-size-val")                 # live value readout
  expect_match(scroll:::.scroll_css(), "\\.scroll-size\\{")    # and styled
})

test_that("on-screen gate is a Shiny input binding (value ready before first flush)", {
  js <- scroll:::.scroll_lazy_js()
  expect_match(js, "Shiny\\.InputBinding")
  expect_match(js, "inputBindings\\.register")
  expect_match(js, "getValue")                          # supplies the initial on-screen value
  expect_match(js, "-onscreen")                         # id maps to the module's input$onscreen
  expect_match(js, "DOMContentLoaded")                  # registered before Shiny.initialize/bindAll
  expect_no_match(js, "shiny:connected")                # no longer deferred to post-connect
})
