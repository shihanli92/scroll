# App chrome: the design system as inline CSS (self-contained -- no external
# fonts or assets, so it works offline on a Shiny Server) plus a small
# scroll-spy that highlights the active panel in the rail.

.scroll_css <- function() "
:root{
  --sc-ground:#FBFCFD; --sc-card:#FFFFFF; --sc-ink:#111826; --sc-muted:#5A6472;
  --sc-faint:#8A93A2; --sc-line:#E4E8EE; --sc-line-2:#EEF1F5;
  --sc-accent:#2563A8; --sc-wash:#EAF1F8; --sc-accent-deep:#1B4C86;
}
/* Kill scrollbar-driven resize loops. On a tall multi-panel page a scrollbar can
   toggle on/off as plots render: the VERTICAL bar changes content width and the
   HORIZONTAL bar changes content height, either of which fires window 'resize',
   which re-renders every fluid-width plotOutput, which nudges the size back --
   an endless self-triggering loop. Reserve the vertical gutter always, and never
   show a page-level horizontal bar (wide plots scroll inside .scroll-plot).
   NOTE: the horizontal clip lives on <html> only. Putting overflow-x:hidden on
   <body> too forces body's computed overflow-y to `auto`, turning body into its
   own scroll container -- which breaks position:sticky on the app bar and panel
   rail (they anchor to the scrolling body box instead of the viewport, so they
   scroll away). <html> already clips horizontal overflow for the whole page. */
html{overflow-x:hidden; overflow-y:scroll; scrollbar-gutter:stable;}
body{background:var(--sc-ground); color:var(--sc-ink);
  font-family:Inter,-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;
  -webkit-font-smoothing:antialiased;}
.container-fluid{padding:0;}

/* app bar */
.scroll-appbar{position:sticky; top:0; z-index:1000; display:flex; align-items:center;
  justify-content:space-between; gap:24px; padding:14px 28px;
  background:rgba(251,252,253,.85); backdrop-filter:saturate(1.4) blur(8px);
  border-bottom:1px solid var(--sc-line);}
.scroll-brand{display:flex; align-items:baseline; gap:8px; font-size:19px;}
.scroll-logo{font-weight:800; letter-spacing:-.02em;}
.scroll-slash{color:var(--sc-faint);}
.scroll-dataset{font-weight:600; color:var(--sc-muted);}
.scroll-stats{display:flex; gap:0;}
.scroll-stat{display:flex; flex-direction:column; padding:0 18px; border-left:1px solid var(--sc-line-2);}
.scroll-stat:first-child{border-left:0;}
.scroll-stat-v{font-weight:700; font-variant-numeric:tabular-nums; line-height:1.1;}
.scroll-stat-l{font-size:11px; text-transform:uppercase; letter-spacing:.08em; color:var(--sc-faint);}

/* global cell-subset filter (app bar pill) */
.scroll-subset{display:flex; align-items:center; gap:8px; padding:5px 6px 5px 12px;
  background:var(--sc-card); border:1px solid var(--sc-line); border-radius:999px;
  box-shadow:0 1px 2px rgba(17,24,38,.04);}
.scroll-subset-label{text-transform:uppercase; letter-spacing:.08em; font-size:11px;
  font-weight:700; color:var(--sc-faint);}
.scroll-subset .shiny-input-container{margin-bottom:0 !important; width:auto !important;}
.scroll-subset .form-select,.scroll-subset .selectize-input{min-height:32px; padding-top:3px;
  padding-bottom:3px; border-radius:8px; font-size:13px;}
.scroll-subset .selectize-control{margin-bottom:0;}
/* keep the selected label clear of the dropdown caret (which flips on open) */
.scroll-subset .selectize-input{padding-right:28px !important;}
.scroll-subset .selectize-control.single .selectize-input:after{right:12px;}
.scroll-subset .selectize-input>.item{overflow:hidden; text-overflow:ellipsis; white-space:nowrap;}

/* Manual per-level colour pickers: compact hex pills that wrap into rows */
.scroll-manual-grid{display:flex; flex-wrap:wrap; gap:6px 8px; margin:4px 0 2px;}
.scroll-swatch{display:flex; flex-direction:column; align-items:center; width:74px;}
.scroll-swatch .shiny-input-container,.scroll-swatch .form-group{width:74px !important;
  margin:0 !important; padding:0 !important;}
.scroll-swatch .shiny-colour-input{width:74px !important; height:24px !important;
  min-height:24px !important; padding:0 4px !important; border:1px solid var(--sc-line) !important;
  border-radius:6px !important; box-shadow:none !important; cursor:pointer;
  font-size:11px; font-family:ui-monospace,SFMono-Regular,Menlo,monospace; letter-spacing:-.02em;}
.scroll-swatch-label{font-size:9px; line-height:1.1; color:var(--sc-muted); margin-top:3px;
  max-width:74px; overflow:hidden; text-overflow:ellipsis; white-space:nowrap; text-align:center;}

/* layout */
.scroll-layout{display:grid; grid-template-columns:200px minmax(0,1fr); gap:32px;
  max-width:1320px; margin:0 auto; padding:28px;}
.scroll-rail{position:sticky; top:80px; align-self:start; display:flex; flex-direction:column; gap:4px;}
.scroll-rail-item{display:flex; align-items:center; gap:10px; padding:9px 12px; border-radius:9px;
  color:var(--sc-muted); text-decoration:none; font-weight:600; font-size:14px;
  border-left:2px solid transparent;}
.scroll-rail-item:hover{background:var(--sc-line-2); color:var(--sc-ink);}
.scroll-rail-item.active{background:var(--sc-wash); color:var(--sc-accent-deep); border-left-color:var(--sc-accent);}
.scroll-rail-num{font-family:ui-monospace,SFMono-Regular,Menlo,monospace; font-size:12px; color:var(--sc-faint);}
.scroll-rail-item.active .scroll-rail-num{color:var(--sc-accent);}
.scroll-rail-foot{margin-top:14px; padding:0 12px; font-size:11px; color:var(--sc-faint);
  font-family:ui-monospace,monospace;}

/* Multi-dataset (scroll_multi_app): pin the dataset tab strip at the very top and
   drop each dataset's sticky app bar + panel rail below it, so the dataset
   selector is always reachable (not just at the top of the scroll). */
.scroll-multi .nav-tabs{position:sticky; top:0; z-index:1001; margin:0;
  padding:6px 20px 0; background:rgba(251,252,253,.92);
  backdrop-filter:saturate(1.4) blur(8px); border-bottom:1px solid var(--sc-line);}
.scroll-multi .scroll-appbar{top:43px;}
.scroll-multi .scroll-rail{top:123px;}

/* panels */
.scroll-content{display:flex; flex-direction:column; gap:28px;}
.scroll-panel-card{scroll-margin-top:84px; border:1px solid var(--sc-line); border-radius:12px;
  background:var(--sc-card); box-shadow:0 6px 20px -12px rgba(17,24,38,.18); overflow:hidden;}
.scroll-panel-card>.card-header{background:var(--sc-card); border-bottom:1px solid var(--sc-line-2); padding:18px 22px;}
.scroll-eyebrow{display:flex; align-items:center; gap:10px;}
.scroll-num{font-family:ui-monospace,monospace; font-size:12px; font-weight:600; color:var(--sc-accent);
  background:var(--sc-wash); padding:2px 7px; border-radius:6px;}
.scroll-kicker{text-transform:uppercase; letter-spacing:.1em; font-size:12px; font-weight:700; color:var(--sc-muted);}
.scroll-title{margin:8px 0 2px; font-size:22px; font-weight:800; letter-spacing:-.02em;}
.scroll-desc{margin:0; color:var(--sc-muted); font-size:14px;}

/* controls */
.scroll-panel{padding:20px 22px;}
.scroll-controls{display:flex; flex-direction:column; gap:18px;}
.scroll-cgroup-h{text-transform:uppercase; letter-spacing:.08em; font-size:11px; font-weight:700;
  color:var(--sc-faint); padding-bottom:8px; margin-bottom:10px; border-bottom:1px solid var(--sc-line-2);}
.scroll-controls .form-label{font-size:13px; font-weight:600; margin-bottom:3px;}
.scroll-controls .form-group,.scroll-controls .shiny-input-container{margin-bottom:12px;}
/* min-width:0 is essential: .scroll-plot is a flex item, so its default
   min-width:auto resolves to min-content = the rendered plot image's width. That
   makes the container size to the image, the image size to the container, and on a
   fractional grid column at HiDPI (devicePixelRatio 2) the sub-pixel rounding
   oscillates -> Shiny's ResizeObserver re-renders endlessly at narrow widths.
   Pinning min-width:0 lets the wrapper track its grid cell; a wide plot scrolls. */
/* overflow-y MUST be hidden, not the default. Setting overflow-x:auto alone forces
   the computed overflow-y from visible to auto, so at HiDPI a plot rendered a
   sub-pixel too tall pops a VERTICAL scrollbar that steals ~15px of width -> the
   plot re-renders narrower -> it fits -> scrollbar vanishes -> loop (the 15px
   width oscillation). Keep horizontal scroll for wide plots, never vertical. */
.scroll-plot{border-radius:10px; background:var(--sc-card); overflow-x:auto; overflow-y:hidden; min-width:0;}
/* Snap the plot output to a whole CSS pixel. bslib's grid columns are fractional
   (e.g. 518.25px); at devicePixelRatio 2 that .25px is 0.5 device px, so Shiny
   renders the image at a rounded device size that displays back a hair different,
   ResizeObserver fires, and every plot re-renders forever (only at HiDPI + narrow
   widths). round(down, 100%, 1px) makes the measured width an integer so device
   pixels land exactly; wide plots still scroll via the overflow above. */
.scroll-plot .shiny-plot-output,
.scroll-plot .shiny-spinner-output-container{width:round(down, 100%, 1px); min-width:0;}
.scroll-plot-bar{display:flex; justify-content:flex-end; padding:0 2px 8px;}
/* .scroll-plot-hold reserves the plot height so panels keep a stable size; the
   lazy-rendering gate (.scroll_lazy_js + .scroll_lazy_plot) recomputes a panel
   only while it is on screen. */
.scroll-dl.btn{padding:3px 10px; font-size:12px; font-weight:600; color:var(--sc-muted);
  background:var(--sc-card); border:1px solid var(--sc-line); border-radius:8px;}
.scroll-dl.btn:hover{color:var(--sc-accent-deep); border-color:var(--sc-accent);
  background:var(--sc-wash);}
.scroll-dl.btn .fa,.scroll-dl.btn svg{margin-right:5px; opacity:.7;}

@media (max-width:900px){
  .scroll-layout{grid-template-columns:1fr; gap:16px;}
  .scroll-rail{position:static; flex-direction:row; overflow-x:auto; top:auto;}
  .scroll-rail-foot{display:none;}
  .scroll-stats{display:none;}
}
"

.scroll_spy_js <- function() "
(function(){
  function spy(){
    var items=document.querySelectorAll('.scroll-rail-item');
    var secs=document.querySelectorAll('.scroll-panel-card');
    if(!secs.length) return;
    var io=new IntersectionObserver(function(es){
      es.forEach(function(e){ if(e.isIntersecting){
        items.forEach(function(it){ it.classList.toggle('active', it.getAttribute('href')==='#'+e.target.id); });
      }});
    },{rootMargin:'-45% 0px -45% 0px'});
    secs.forEach(function(s){io.observe(s);});
  }
  if(document.readyState!=='loading') spy(); else document.addEventListener('DOMContentLoaded',spy);
})();
"

# Lazy panel rendering: report each panel card's on-screen state (viewport + a
# margin) to its module as input$onscreen. .scroll_lazy_plot uses that to recompute
# a panel only while it is on screen, so a View/filter change redraws just the
# panels in view; an off-screen one refreshes when scrolled to. The card id is the
# module id, so `<id>-onscreen` lands on the module's input$onscreen.
.scroll_lazy_js <- function() "
(function(){
  function lazy(){
    var cards=document.querySelectorAll('.scroll-panel-card');
    if(!cards.length || !window.Shiny) return;
    var io=new IntersectionObserver(function(es){
      es.forEach(function(e){ Shiny.setInputValue(e.target.id+'-onscreen', e.isIntersecting); });
    },{rootMargin:'300px 0px 300px 0px'});
    cards.forEach(function(c){io.observe(c);});
  }
  if(document.readyState!=='loading') lazy(); else document.addEventListener('DOMContentLoaded',lazy);
})();
"
