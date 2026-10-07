# Regenerates the screenshots in the "A tour of the app" article
# (vignettes/articles/app-tour.Rmd) by driving a running scroll app in headless
# Chrome (chromote). Built against the public 10x "CD8+ T cells, 4 donors
# (GEX + TCR + antigen)" project, with gene sets added by scroll_add_genesets().
#
#   1. serve the project:  scroll::scroll_serve("cd8_merged", port = 7797)
#   2. from the repo root: Rscript vignettes/articles/tour-screenshots.R
#
# Set SCROLL_TOUR_URL to point at another address, and SCROLL_TOUR_ONLY to a
# comma-separated list of names (e.g. "tour-dotplot,tour-gsea") to redo only those.

url <- Sys.getenv("SCROLL_TOUR_URL", "http://127.0.0.1:7797")
out <- "vignettes/articles/figures"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

b <- chromote::ChromoteSession$new(width = 1440, height = 1000)
on.exit(b$close(), add = TRUE)

js <- function(code) {
  r <- b$Runtime$evaluate(sprintf("(async () => { %s })()", code),
                          awaitPromise = TRUE, returnByValue = TRUE, timeout_ = 180)
  if (!is.null(r$exceptionDetails)) stop(r$exceptionDetails$exception$description)
  r$result$value
}
pause <- function(s) Sys.sleep(s)
# wait until Shiny has nothing in flight (and the startup warm-up is over); with a
# selector, also until a plot has been drawn inside it and nothing there is still
# recalculating
idle <- function(selector = NULL, timeout = 120, settle = 1.5) {
  sel <- if (is.null(selector)) "false" else sprintf(
    "!!document.querySelector('%1$s .recalculating') ||
     !document.querySelector('%1$s .shiny-plot-output img, %1$s .scroll-table table')", selector)
  t0 <- Sys.time()
  quiet <- 0
  repeat {
    pause(0.5)
    busy <- js(sprintf("const h = document.documentElement.classList;
      return h.contains('shiny-busy') || h.contains('scroll-warming') || %s;", sel))
    quiet <- if (isTRUE(busy)) 0 else quiet + 0.5
    if (quiet >= settle) break
    if (difftime(Sys.time(), t0, units = "secs") > timeout) {
      message("idle timeout: ", selector, " ", js(sprintf("const o = document.querySelector('%s .shiny-plot-output');
        return o ? o.className + ' | ' + o.innerText.slice(0, 80) : 'no plot output';", selector)))
      break
    }
  }
}
# set a Shiny input through its widget so the page and the server agree (selectize
# options loaded from the server may not exist yet, so add them first)
set <- function(id, value) {
  v <- jsonlite::toJSON(as.list(value), auto_unbox = TRUE)
  js(sprintf("const el = document.getElementById('%1$s'), v = %2$s;
    if (el && el.selectize) {
      const sz = el.selectize;
      v.forEach(x => { if (!sz.options[x]) sz.addOption({value: x, label: x}); });
      sz.setValue(v.length === 1 && !el.multiple ? v[0] : v);
    } else if (el) { $(el).val(v.length === 1 ? v[0] : v).trigger('change'); }
    else { Shiny.setInputValue('%1$s', v.length === 1 ? v[0] : v); }
    return true;", id, v))
  pause(0.5)
}
click <- function(id) js(sprintf("document.getElementById('%s').click(); return true;", id))
# scroll an element on screen; the small nudge makes sure the panel's on-screen
# observer fires, since panels draw only while visible
show <- function(selector, block = "center") {
  js(sprintf("const el = document.querySelector('%s');
    el.scrollIntoView({block: '%s', behavior: 'instant'});
    window.scrollBy(0, 40); await new Promise(r => setTimeout(r, 300));
    window.scrollBy(0, -40); return true;", selector, block))
  pause(0.8)
}
only <- strsplit(Sys.getenv("SCROLL_TOUR_ONLY"), ",")[[1]]
step <- function(name, expr) if (!length(only) || name %in% only) expr
# capture an element: scroll it on screen (panels render only while visible), then
# clip the page at its document position
shot <- function(name, selector) {
  rect <- js(sprintf("const r = document.querySelector('%s').getBoundingClientRect();
    return [r.left, r.top + window.scrollY, r.width, r.height];", selector))
  f <- file.path(out, paste0(name, ".png"))
  img <- b$Page$captureScreenshot(format = "png", captureBeyondViewport = TRUE,
    clip = list(x = rect[[1]], y = rect[[2]], width = rect[[3]], height = rect[[4]], scale = 1))
  writeBin(jsonlite::base64_dec(img$data), f)
  message("wrote ", f)
}
# the whole visible window
screen <- function(name) {
  f <- file.path(out, paste0(name, ".png"))
  img <- b$Page$captureScreenshot(format = "png")
  writeBin(jsonlite::base64_dec(img$data), f)
  message("wrote ", f)
}
tab <- function(card, label)            # click a bslib nav tab inside a card
  js(sprintf("const a = [...document.querySelectorAll('#%s .nav-link')]
      .find(x => x.textContent.trim() === '%s'); a.click(); return true;", card, label))

b$Page$navigate(url)
pause(5); idle("#dimplot", timeout = 240)

# 1. the whole page: app bar, panel rail, first panel, filter rail
step("tour-overview", { js("window.scrollTo(0, 0); return true;"); pause(1); screen("tour-overview") })
# from here on, keep the sticky app bar from covering the panel being captured
js("const s = document.createElement('style'); s.id = 'tour-css';
    s.textContent = '.scroll-appbar{position:relative !important}'; document.head.appendChild(s); return true;")

# 2. DimPlot coloured by donor, with its Style sheet open
step("tour-dimplot", {
  show("#dimplot")
  set("dimplot-colorby", "donor"); idle("#dimplot")
  shot("tour-dimplot", "#dimplot")
})

# 3. FeaturePlot: several genes as a grid
step("tour-featureplot", {
  show("#featureplot")
  set("featureplot-feature", c("CD8A", "GZMB", "CCR7", "GNLY")); idle("#featureplot")
  shot("tour-featureplot", "#featureplot")
})

# 4. Signature: UCell on a cytotoxicity set
step("tour-signature", {
  show("#signature")
  set("signature-sig", c("GZMB", "PRF1", "NKG7", "GNLY", "GZMH"))
  set("signature-method", "ucell"); click("signature-compute"); idle("#signature")
  shot("tour-signature", "#signature")
})

# 5. DotPlot of the configured markers
step("tour-dotplot", {
  show("#dotplot"); idle("#dotplot")
  shot("tour-dotplot", "#dotplot")
})

# 6. Proportions: donor composition of clusters
step("tour-proportions", {
  show("#proportions")
  set("proportions-group", "seurat_clusters"); set("proportions-fill", "donor"); idle("#proportions")
  shot("tour-proportions", "#proportions")
})

# 7. DE: a cluster against the rest, volcano tab
step(if ("tour-gsea" %in% only) "tour-gsea" else "tour-de", {
  show("#de")
  set("de-group", "seurat_clusters"); pause(1); set("de-ident1", "3")
  set("de-maxcells", 3000)                   # subsample: keeps the app's memory modest
  click("de-compute"); idle("#de", timeout = 300)
  tab("de", "Volcano"); idle("#de")
  shot("tour-de", "#de")
})

# 8. GSEA on that DE result (Hallmark)
step("tour-gsea", {
  show("#gsea")
  set("gsea-source", "de"); pause(1)
  set("gsea-collection", "saved:H"); click("gsea-run"); idle("#gsea", timeout = 300)
  tab("gsea", "Plot"); idle("#gsea")
  shot("tour-gsea", "#gsea")
})

# 9. a repertoire panel
step("tour-repertoire", {
  show("#clone_overview"); idle("#clone_overview")
  shot("tour-repertoire", "#clone_overview")
})

# 10. the Style sheet on the DimPlot
step("tour-style", {
  js("document.getElementById('tour-css').remove(); return true;")
  show("#dimplot")
  js("document.querySelector('#dimplot .scroll-style-btn').click(); return true;"); idle("#dimplot")
  screen("tour-style")
  js("document.querySelector('.scroll-sheet.is-open .scroll-sheet-close, .scroll-sheet [aria-label=Close]')?.click(); return true;")
  js("document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape'})); return true;"); pause(1)
})

# 11. the View selector: switch to one donor's own UMAP
step("tour-view", {
  js("window.scrollTo(0, 0); return true;")
  set("scroll_view", "donor1"); idle("#dimplot", timeout = 180)
  screen("tour-view")
  set("scroll_view", ""); idle(timeout = 180)
})

# 12. the rail: Filters and New column
step("tour-rail", {
  js("[...document.querySelectorAll('.scroll-filters summary')]
      .find(s => s.textContent === 'New column').click(); return true;")
  set("scroll_newcol_map", "[Donor 1, Donor 2]: A, [Donor 3, Donor 4]: B")
  set("scroll_newcol_name", "batch")
  js("const r = document.querySelector('.scroll-filters'); r.scrollTop = 1e5; return true;"); pause(0.8)
  shot("tour-rail", ".scroll-filters")
})
