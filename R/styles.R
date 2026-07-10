# App chrome: the design system as inline CSS (self-contained -- no external
# fonts or assets, so it works offline on a Shiny Server) plus a small
# scroll-spy that highlights the active section in the rail.

.scroll_css <- function() "
:root{
  --sc-ground:#FBFCFD; --sc-card:#FFFFFF; --sc-ink:#111826; --sc-muted:#5A6472;
  --sc-faint:#8A93A2; --sc-line:#E4E8EE; --sc-line-2:#EEF1F5;
  --sc-accent:#2563A8; --sc-wash:#EAF1F8; --sc-accent-deep:#1B4C86;
}
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

/* sections */
.scroll-content{display:flex; flex-direction:column; gap:28px;}
.scroll-section{scroll-margin-top:84px; border:1px solid var(--sc-line); border-radius:12px;
  background:var(--sc-card); box-shadow:0 6px 20px -12px rgba(17,24,38,.18); overflow:hidden;}
.scroll-section>.card-header{background:var(--sc-card); border-bottom:1px solid var(--sc-line-2); padding:18px 22px;}
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
.scroll-plot{border-radius:10px; background:var(--sc-card); overflow-x:auto;}
.scroll-plot-bar{display:flex; justify-content:flex-end; padding:0 2px 8px;}
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
    var secs=document.querySelectorAll('.scroll-section');
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
