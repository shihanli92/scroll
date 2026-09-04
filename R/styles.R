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

/* Startup warm-up overlay (shown while the subset views pre-render) */
.scroll-warm-overlay{position:fixed; inset:0; z-index:9999; display:flex;
  align-items:center; justify-content:center; background:rgba(251,252,253,.86);
  backdrop-filter:blur(1px); opacity:1; transition:opacity .3s ease;}
.scroll-warm-overlay.scroll-warm-hiding{opacity:0;}   /* fade out before display:none */
.scroll-warm-box{display:flex; flex-direction:column; align-items:center; gap:12px;
  font-size:14px; color:#4b5563; font-weight:500; min-width:240px;}
.scroll-warm-row{display:flex; align-items:center; gap:10px;}
.scroll-warm-spin{width:16px; height:16px; border:2px solid #d7dce3;
  border-top-color:#2563A8; border-radius:50%; animation:scrollwarmspin .8s linear infinite;}
@keyframes scrollwarmspin{to{transform:rotate(360deg);}}
.scroll-warm-bar{width:100%; height:4px; border-radius:2px; background:#e4e8ee; overflow:hidden;}
.scroll-warm-fill{height:100%; width:0; background:#2563A8; border-radius:2px;
  transition:width .35s ease;}
/* lock page scroll while warming so scrolling can't flip panels on-screen and pollute
   the plot cache mid-cycle (the overlay covers the page; this stops wheel/keyboard scroll) */
body.scroll-warming{overflow:hidden;}

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
.scroll-layout.has-filters{grid-template-columns:200px minmax(0,1fr) 300px; max-width:1620px;}
/* right-hand global control rail: filters that narrow every panel */
.scroll-filters{position:sticky; top:80px; align-self:start; display:flex; flex-direction:column;
  gap:12px; max-height:calc(100vh - 100px); overflow-y:auto; padding-right:2px;}
.scroll-multi .scroll-filters{top:123px;}
/* collapse toggle: hide the whole control rail to give the plots full width */
.scroll-ctl-toggle{flex:none; width:32px; height:32px; padding:0; cursor:pointer; line-height:1;
  border:1px solid var(--sc-line); background:#fff; color:var(--sc-muted); border-radius:8px;
  font-size:15px;}
.scroll-ctl-toggle:hover{color:var(--sc-ink); border-color:var(--sc-accent);}
.scroll-ctl-toggle::before{content:'\\00BB';}                 /* > : hide the rail */
.scroll-ctl-toggle.is-collapsed::before{content:'\\00AB';}    /* < : show the rail */
.scroll-layout.controls-collapsed{grid-template-columns:200px minmax(0,1fr);}
.scroll-layout.controls-collapsed>.scroll-filters{display:none;}
/* collapsible <details> sections (Theme / Filters, and Theme's sub-sections) */
.scroll-ctl-details{border-bottom:1px solid var(--sc-line);}
.scroll-ctl-details:last-child{border-bottom:none;}
.scroll-ctl-summary{cursor:pointer; list-style:none; display:flex; align-items:center; gap:6px;
  font-weight:700; font-size:12px; letter-spacing:.05em; text-transform:uppercase;
  color:var(--sc-muted); padding:7px 0;}
.scroll-ctl-summary::-webkit-details-marker{display:none;}
.scroll-ctl-summary::before{content:'\\25B8'; font-size:9px; color:var(--sc-faint);
  transition:transform .12s;}
.scroll-ctl-details[open]>.scroll-ctl-summary::before{transform:rotate(90deg);}
.scroll-ctl-body{display:flex; flex-direction:column; gap:6px; padding:1px 0 8px;}
/* nested sub-sections (Text & fonts / Axes / Panel inside Theme) read lighter */
.scroll-ctl-body .scroll-ctl-details{border-bottom:none;}
.scroll-ctl-body .scroll-ctl-summary{text-transform:none; letter-spacing:0; font-size:12px;
  font-weight:600; color:var(--sc-ink); padding:5px 0;}
/* theme sub-section controls lay out two per row (the rail is wide enough) */
.scroll-ctl-body .scroll-ctl-body{display:grid; grid-template-columns:1fr 1fr;
  gap:6px 10px; align-items:end; padding:0 0 6px 4px;}
.scroll-reset-row{display:flex; justify-content:flex-end; margin-bottom:2px;}
.scroll-reset-row a{font-size:11px; color:var(--sc-accent); text-decoration:none;}
/* Apply / Reset row at the top of the Theme and Filters sections */
.scroll-apply-row{display:flex; align-items:center; gap:10px; margin:2px 0 8px;}
.scroll-apply-row .btn{padding:3px 16px; font-size:12px; font-weight:600;}
.scroll-apply-row a{font-size:11px; color:var(--sc-accent); text-decoration:none; margin-left:auto;}
/* compact controls: kill the default form-group margin + shrink the selectize/slider */
.scroll-filter{display:flex; flex-direction:column; min-width:0;}
.scroll-filter .form-group,.scroll-filter .shiny-input-container{margin-bottom:0;}
.scroll-filter .control-label,.scroll-filter .form-label,.scroll-filter label{
  font-size:11.5px; font-weight:600; margin-bottom:1px; color:var(--sc-ink);}
/* colour options as a small circular swatch with the label beside it */
.scroll-colour{flex-direction:row; align-items:center; gap:8px;}
.scroll-colour .control-label,.scroll-colour label{flex:1; margin:0;}
.scroll-colour .shiny-input-container{display:flex; align-items:center; gap:8px; width:100%;}
.scroll-colour .input-group{width:auto;}
.scroll-colour input.shiny-colour-input{width:22px; height:22px; min-height:22px; padding:0;
  border-radius:50%; border:1px solid var(--sc-line); font-size:0; box-shadow:none; cursor:pointer;}
.scroll-filters .selectize-input{min-height:0; padding:4px 8px; line-height:16px;}
.scroll-filters .selectize-input.items{padding-top:4px; padding-bottom:4px;}
.scroll-filters .selectize-dropdown{font-size:12px;}
.scroll-filter .irs{font-size:10px;}
.scroll-filter .irs--shiny{top:0; height:26px;}
.scroll-filter .irs-with-grid{height:34px;}
.scroll-filters-empty{font-size:11px; color:var(--sc-faint);}
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
.scroll-plot-bar{display:flex; justify-content:flex-end; align-items:center; gap:8px; padding:0 2px 8px;}
/* condensed clone-map table: smaller text + tighter cells, single-line rows */
table.scroll-clone-dt{font-size:12px;}
table.scroll-clone-dt td,table.scroll-clone-dt th{padding:3px 8px;}
/* smaller per-column + global search inputs */
#clone_map .dataTables_filter input,
table.scroll-clone-dt thead input{
  font-size:11px; padding:1px 6px; height:24px; min-height:0; line-height:1.3;}
table.scroll-clone-dt thead th{padding-top:2px; padding-bottom:2px;}
/* compact download-scale slider + live value readout, sitting with the downloads */
.scroll-size-wrap{display:inline-flex; align-items:center; gap:6px;}
.scroll-size{width:80px; height:16px; accent-color:var(--sc-accent); cursor:pointer;}
.scroll-size::-webkit-slider-thumb{cursor:pointer;}
.scroll-size-val{font-size:12px; font-weight:600; color:var(--sc-muted);
  font-variant-numeric:tabular-nums; min-width:26px; text-align:right;}
/* .scroll-plot-hold reserves the plot height so panels keep a stable size; the
   lazy-rendering gate (.scroll_lazy_js + .scroll_lazy_plot) recomputes a panel
   only while it is on screen. */
.scroll-dl.btn{padding:3px 10px; font-size:12px; font-weight:600; color:var(--sc-muted);
  background:var(--sc-card); border:1px solid var(--sc-line); border-radius:8px;}
.scroll-dl.btn:hover{color:var(--sc-accent-deep); border-color:var(--sc-accent);
  background:var(--sc-wash);}
.scroll-dl.btn .fa,.scroll-dl.btn svg{margin-right:5px; opacity:.7;}

@media (max-width:900px){
  .scroll-layout,.scroll-layout.has-filters,.scroll-layout.controls-collapsed{
    grid-template-columns:1fr; gap:16px;}
  .scroll-rail{position:static; flex-direction:row; overflow-x:auto; top:auto;}
  .scroll-rail-foot{display:none;}
  .scroll-stats{display:none;}
  .scroll-filters{position:static; top:auto; max-height:none;}
}
"

.scroll_spy_js <- function() "
window.scrollToggleControls=function(btn){
  var body=btn.closest('.scroll-appbar').parentElement;
  var lay=body?body.querySelector('.scroll-layout'):document.querySelector('.scroll-layout');
  if(!lay) return;
  var collapsed=lay.classList.toggle('controls-collapsed');
  btn.classList.toggle('is-collapsed', collapsed);
  btn.setAttribute('aria-expanded', String(!collapsed));
};
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
# a panel only while it is on screen, so a View/filter change redraws just the panels
# in view; an off-screen one refreshes when scrolled to. The card id is the module id,
# so `<id>-onscreen` lands on the module's input$onscreen.
#
# This is a custom Shiny INPUT BINDING (not setInputValue on shiny:connected) so that
# getValue() runs during Shiny's initial bindAll -- BEFORE the first output flush --
# giving every panel a correct on-screen value up front. The old approach left
# input$onscreen unset until the observer reported post-connect, so the first flush saw
# NULL (-> the `%||% TRUE` default) and every live panel rendered at once. Registering
# on DOMContentLoaded is in time: native DOMContentLoaded listeners run before jQuery's
# ready callbacks (where Shiny.initialize/bindAll live), and window.Shiny is defined by
# then (its <script> executes during parse).
.scroll_lazy_js <- function() "
(function(){
  var M=80;                                            // viewport margin (matches rootMargin)
  function onscreen(el){
    var r=el.getBoundingClientRect();
    var vh=window.innerHeight||document.documentElement.clientHeight;
    return r.bottom > -M && r.top < vh + M;            // vertical intersection (+ margin)
  }
  function reg(){
    if(window.__scrollOnReg) return;
    if(!window.Shiny || !Shiny.InputBinding || !Shiny.inputBindings) return;
    window.__scrollOnReg=true;
    var b=new Shiny.InputBinding();
    $.extend(b,{
      find: function(scope){ return $(scope).find('.scroll-panel-card'); },
      getId: function(el){ return el.id+'-onscreen'; },
      // read at bind time (pre-first-flush), so off-screen panels report false up front
      getValue: function(el){
        return (typeof el.__scrollOn==='boolean') ? el.__scrollOn : onscreen(el);
      },
      subscribe: function(el, cb){
        var io=new IntersectionObserver(function(es){
          es.forEach(function(e){ el.__scrollOn=e.isIntersecting; });
          cb();
        },{rootMargin:M+'px 0px '+M+'px 0px'});
        io.observe(el); el.__scrollIO=io;
      },
      unsubscribe: function(el){ if(el.__scrollIO){ el.__scrollIO.disconnect(); delete el.__scrollIO; } }
    });
    Shiny.inputBindings.register(b, 'scroll.onscreen');
  }
  if(document.readyState!=='loading') reg();           // injected after load: Shiny present
  document.addEventListener('DOMContentLoaded', reg);  // fires before Shiny.initialize/bindAll
})();
"

# Startup warm-up overlay: the server sends 'scroll_warm' {show:true/false} while it
# cycles the subset views to pre-render their cached scatters; the overlay covers the
# viewport so the brief view-cycling isn't visible and can't be interacted with.
.scroll_warm_js <- function() "
(function(){
  // Client failsafe: hide the overlay only after the server has been SILENT this long.
  // The server drives a step at least every step_timeout (20s) -- or aborts and sends
  // {show:false} -- so a healthy (even long) warm-up keeps resetting this and never
  // trips it; it fires only when the server never responds (dead/dropped session,
  // where the overlay -- painted in the initial HTML -- would otherwise sit forever).
  var SILENCE_MS = 30000;
  var fuse = null;
  function hide(ov){
    if(fuse){ clearTimeout(fuse); fuse=null; }
    document.body.classList.remove('scroll-warming');
    ov.classList.add('scroll-warm-hiding');           // fade, then remove from flow
    setTimeout(function(){ ov.style.display='none'; }, 340);
  }
  // (re)arm the silence fuse; each progress message resets the clock. No-op when there
  // is no overlay (no warm-up), so it never fires on a plain app.
  function arm(){
    if(!document.getElementById('scroll-warm-overlay')) return;
    if(fuse) clearTimeout(fuse);
    fuse = setTimeout(function(){
      var ov=document.getElementById('scroll-warm-overlay');
      if(ov && ov.style.display!=='none') hide(ov);
    }, SILENCE_MS);
  }
  function reg(){
    if(!window.Shiny || !Shiny.addCustomMessageHandler) return;
    Shiny.addCustomMessageHandler('scroll_warm', function(m){
      var ov=document.getElementById('scroll-warm-overlay'); if(!ov) return;
      arm();                                            // progress -> reset the silence fuse
      if(m && m.text){ var el=ov.querySelector('.scroll-warm-msg'); if(el) el.textContent=m.text; }
      if(m && m.frac!=null){ var f=ov.querySelector('.scroll-warm-fill');
        if(f) f.style.width=(Math.max(0,Math.min(1,m.frac))*100)+'%'; }
      if(m && m.show){ ov.classList.remove('scroll-warm-hiding'); ov.style.display='flex';
        document.body.classList.add('scroll-warming'); }
      else { hide(ov); }
    });
  }
  function onConnect(){ reg(); arm(); }                 // fresh window for the first message
  if(window.Shiny && Shiny.addCustomMessageHandler) onConnect();
  else if(window.jQuery) jQuery(document).on('shiny:connected', onConnect);
  else document.addEventListener('shiny:connected', onConnect);
  // lock scroll + arm the fuse from the initial paint (body exists only after parse)
  function initLock(){ if(document.getElementById('scroll-warm-overlay')){
    document.body.classList.add('scroll-warming'); arm(); } }
  if(document.readyState!=='loading') initLock();
  else document.addEventListener('DOMContentLoaded', initLock);
})();
"
