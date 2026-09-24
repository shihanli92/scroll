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
  justify-content:flex-start; gap:16px; padding:12px 28px; flex-wrap:wrap; row-gap:8px;
  background:color-mix(in srgb, var(--sc-ground) 85%, transparent);
  backdrop-filter:saturate(1.4) blur(8px);
  border-bottom:1px solid var(--sc-line); box-shadow:0 6px 16px -14px rgba(17,24,38,.25);}
/* the dataset title is the dominant element; the wordmark is a small muted prefix */
.scroll-brand{display:flex; align-items:baseline; gap:7px; min-width:0; flex:0 1 auto;}
.scroll-logo{font-size:14px; font-weight:700; color:var(--sc-muted); letter-spacing:0;}
.scroll-slash{color:var(--sc-faint); font-size:14px;}
.scroll-dataset{font-size:18px; font-weight:700; color:var(--sc-ink); letter-spacing:-.01em;
  min-width:0; overflow:hidden; text-overflow:ellipsis; white-space:nowrap;}   /* ellipsize a long title */
/* live cell readout -- the one metric that changes with the View/Subset pills */
.scroll-cells{font-size:13px; font-weight:600; color:var(--sc-muted);
  font-variant-numeric:tabular-nums; white-space:nowrap;}
/* the layout actions travel together and hug the right edge of whichever row they land on */
.scroll-appbar-actions{display:flex; align-items:center; gap:6px; margin-left:auto; flex:none;}
/* dataset-details popover body */
.scroll-info-dl{display:grid; grid-template-columns:auto 1fr; gap:2px 14px; margin:0; font-size:12px; max-width:320px;}
.scroll-info-dl dt{color:var(--sc-faint); font-weight:600;}
.scroll-info-dl dd{margin:0; color:var(--sc-ink); font-variant-numeric:tabular-nums;
  overflow-wrap:anywhere;}

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
/* Full-bleed shell: the rail (and, when present, the filters column) hug the screen
   edges while the content fills the middle up to --sc-content-max, so a wide monitor
   gives the plot more room instead of dead margins -- but the content stays capped so
   a scatter/UMAP never grows to an unusable width. The grid grows the capped content
   track into free space up to its limit, then justify-content spreads the leftover to
   the edges. Below the cap it behaves exactly like the old 1fr. */
:root{--sc-rail-w:160px; --sc-content-max:1600px;}
.scroll-layout{display:grid; gap:24px; padding:28px; justify-content:start;
  grid-template-columns:var(--sc-rail-w) minmax(0,var(--sc-content-max));}
.scroll-layout.has-filters{justify-content:space-between;
  grid-template-columns:var(--sc-rail-w) minmax(0,var(--sc-content-max)) 300px;}
/* right-hand global control rail: filters that narrow every panel. Sticky offset is
   driven by --sc-appbar-h (the app bar's measured bottom, set in JS), so it stays
   correct when the app bar wraps to two rows on a narrow screen; 70 is the fallback
   before the script runs, and it resolves to the old 123px in a multi-app tab strip. */
.scroll-filters{position:sticky; top:calc(var(--sc-appbar-h,70px) + 10px); align-self:start;
  display:flex; flex-direction:column;
  gap:12px; max-height:calc(100vh - var(--sc-appbar-h,70px) - 30px); overflow-y:auto; padding-right:2px;}
/* app-bar icon buttons (dataset info, two-up, controls) -- one shared pill, real
   Font Awesome icons instead of font glyphs so they match across platforms. State
   is shown by a pressed style (+ the title tooltip updates in JS), not a glyph swap. */
.scroll-ctl-toggle,.scroll-two-toggle,.scroll-info-toggle{
  flex:none; display:inline-flex; align-items:center; justify-content:center;
  width:32px; height:32px; padding:0; cursor:pointer; line-height:1; font-size:13px;
  border:1px solid var(--sc-line); background:var(--sc-card); color:var(--sc-muted); border-radius:8px;}
.scroll-ctl-toggle:hover,.scroll-two-toggle:hover,.scroll-info-toggle:hover{
  color:var(--sc-ink); border-color:var(--sc-accent);}
/* Hover/focus label under each app-bar button, saying what it does in its current
   state (text is set in JS on hover/focus, see scrollTip). Replaces the native title
   tooltip, which appears only after ~1s and is easy to miss. Anchored to the right
   edge because the cluster hugs the right side of the bar. */
.scroll-appbar-actions [data-tip],.scroll-plot-bar [data-tip]{position:relative;}
.scroll-appbar-actions [data-tip]::after,.scroll-plot-bar [data-tip]::after{content:attr(data-tip); position:absolute;
  top:calc(100% + 8px); right:0; z-index:1001; width:max-content; max-width:240px;
  padding:6px 9px; border-radius:6px; background:var(--sc-ink); color:var(--sc-card);
  font-size:12px; font-weight:500; line-height:1.35; letter-spacing:0; white-space:normal;
  text-align:left; box-shadow:0 4px 14px rgba(17,24,38,.18); pointer-events:none;
  opacity:0; transform:translateY(-3px); transition:opacity .12s ease, transform .12s ease;}
.scroll-appbar-actions [data-tip]:hover::after,
.scroll-appbar-actions [data-tip]:focus-visible::after,
.scroll-plot-bar [data-tip]:hover::after,
.scroll-plot-bar [data-tip]:focus-visible::after{opacity:1; transform:none;
  transition-delay:.25s;}
/* no hover label while the button's sheet is open (the sheet names itself) */
.scroll-style-btn.is-open::after{display:none;}
/* the info button's popover is open -> its label would sit on top of it */
.scroll-info-toggle[aria-describedby]::after{display:none;}
/* pressed / active: two-up on, rail collapsed, or the filters drawer open */
.scroll-two-toggle.is-on,.scroll-ctl-toggle.is-collapsed,.scroll-ctl-toggle.is-open{
  color:var(--sc-accent-deep); border-color:var(--sc-accent); background:var(--sc-wash);}
/* two-up only shown where two columns can actually fit */
@media (max-width:1549.98px){.scroll-two-toggle{display:none;}}   /* only where two columns can actually fit (with filters drawered) */
.scroll-layout.controls-collapsed{justify-content:start;
  grid-template-columns:var(--sc-rail-w) minmax(0,var(--sc-content-max));}
.scroll-layout.controls-collapsed>.scroll-filters{display:none;}
/* a multi-select with many chips (e.g. 100 default genes) scrolls instead of
   growing the control column to thousands of px (which left the plot stranded). */
.scroll-controls .selectize-input{max-height:170px; overflow-y:auto;}
/* a tall control stack scrolls within the viewport instead of forcing the card
   (and its plot column) taller than the plot -- keeps the plot filling the space. */
.scroll-controls{max-height:calc(100vh - 150px); overflow-y:auto; padding-right:6px;}
/* pseudobulk design matrix (compact, monospace; scrolls if wide) */
.scroll-designmat-wrap{overflow-x:auto; max-width:100%;}
.scroll-designmat{border-collapse:collapse; font-family:ui-monospace,SFMono-Regular,Menlo,monospace;
  font-size:10px; font-variant-numeric:tabular-nums;}
.scroll-designmat th,.scroll-designmat td{padding:1px 5px; text-align:center; white-space:nowrap;}
.scroll-designmat th{color:var(--sc-muted); font-weight:600;}
.scroll-designmat th:first-child{text-align:left; color:var(--sc-faint); font-weight:500;}
.scroll-designmat td{color:var(--sc-ink);}
.scroll-designmat-contrast th,.scroll-designmat-contrast td{border-top:1px solid var(--sc-line-2);
  color:var(--sc-accent); font-weight:700; padding-top:2px;}
.scroll-designmat-n{border-left:1px solid var(--sc-line-2); color:var(--sc-muted); text-align:right;}
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
.scroll-rail{position:sticky; top:calc(var(--sc-appbar-h,70px) + 10px); align-self:start; display:flex; flex-direction:column; gap:2px;}
.scroll-rail-item{display:flex; align-items:center; gap:8px; padding:7px 10px; border-radius:8px;
  color:var(--sc-muted); text-decoration:none; font-weight:600; font-size:13px;
  border-left:2px solid transparent;}
/* long/custom labels ellipsize instead of forcing the (now 160px) rail wider; the
   title attr keeps the full label on hover */
.scroll-rail-item>span:last-child{min-width:0; white-space:nowrap; overflow:hidden; text-overflow:ellipsis;}
/* drag-to-reorder affordance (see .scroll_spy_js): a grip of raised dots on the left
   signals the item is movable; the item dims while dragging and siblings slide (FLIP)
   as it is pulled out and re-inserted. */
.scroll-rail-item{cursor:grab;}
.scroll-rail-item::before{content:''; flex:none; align-self:center; width:8px; height:14px;
  opacity:.28; color:var(--sc-faint); background-repeat:repeat;
  background-image:radial-gradient(currentColor 42%, transparent 46%); background-size:4px 4px;}
.scroll-rail-item:hover::before{opacity:.55;}
.scroll-rail-item.is-dragging{opacity:.4; cursor:grabbing;}
.scroll-rail-reset{display:block; margin-top:6px; background:none; border:none; padding:0;
  color:var(--sc-accent); font:inherit; font-size:11px; cursor:pointer; text-decoration:underline;}
.scroll-rail-item:hover{background:var(--sc-line-2); color:var(--sc-ink);}
.scroll-rail-item.active{background:var(--sc-wash); color:var(--sc-accent-deep); border-left-color:var(--sc-accent);}
.scroll-rail-num{font-family:ui-monospace,SFMono-Regular,Menlo,monospace; font-size:12px; color:var(--sc-faint);}
.scroll-rail-item.active .scroll-rail-num{color:var(--sc-accent);}
.scroll-rail-foot{margin-top:14px; padding:0 12px; font-size:11px; color:var(--sc-faint);
  font-family:ui-monospace,monospace;}
.scroll-rail-ver{display:block; font-size:10px; color:var(--sc-faint);}   /* scroll version, moved out of the app bar */

/* Multi-dataset (scroll_multi_app): pin the dataset tab strip at the very top and
   drop each dataset's sticky app bar + panel rail below it, so the dataset
   selector is always reachable (not just at the top of the scroll). */
.scroll-multi .nav-tabs{position:sticky; top:0; z-index:1001; margin:0;
  padding:6px 20px 0; background:rgba(251,252,253,.92);
  backdrop-filter:saturate(1.4) blur(8px); border-bottom:1px solid var(--sc-line);}
.scroll-multi .scroll-appbar{top:43px;}
/* rail/filters/scroll-margin offsets now come from --sc-appbar-h (the JS measures
   the app bar's bottom, which in multi-app already includes the 43px tab strip). */

/* panels */
/* container-type lets each panel size its control/plot split from the width the
   card actually has (see the @container rules below). inline-size also means panel
   content can never widen the column, strengthening the HiDPI no-reflow invariant. */
/* one card per row by default; .two-up packs two per row on a wide screen (auto-fit
   self-degrades to one column when there isn't room). align-items:start keeps
   ragged-height row-mates from stretching. */
.scroll-content{display:grid; grid-template-columns:minmax(0,1fr); gap:28px; align-items:start;
  container-type:inline-size; container-name:sc-content;}
.scroll-content.two-up{grid-template-columns:repeat(auto-fit, minmax(min(640px,100%),1fr));}
/* Each CARD is the query container for its control/plot split, so a half-width card
   in two-up mode stacks its own controls independently (see the @container sc-card
   rules below). --sc-content-max is bumped when two-up is on so two cards fit. */
.scroll-layout:has(.two-up){--sc-content-max:2200px;}
/* Two-up needs the whole width: when it's on, the filters column becomes an off-canvas
   drawer even on a wide screen (otherwise its 300px track leaves no room for a second
   card, so two-up looked like it did nothing). The app-bar controls button opens it. */
.scroll-layout.has-filters:has(.two-up){justify-content:start;
  grid-template-columns:var(--sc-rail-w) minmax(0,var(--sc-content-max));}
.scroll-layout.has-filters:has(.two-up)>.scroll-filters,
.scroll-layout.controls-collapsed:has(.two-up)>.scroll-filters{
  position:fixed; top:var(--sc-appbar-h,70px); right:0; z-index:1050; display:flex;
  width:min(300px,90vw); height:calc(100dvh - var(--sc-appbar-h,70px)); max-height:none;
  padding:16px 16px 24px; background:var(--sc-card); border-left:1px solid var(--sc-line);
  box-shadow:-14px 0 34px -18px rgba(17,24,38,.4); overflow-y:auto; overscroll-behavior:contain;
  transform:translateX(105%); transition:transform .18s ease;}
.scroll-layout.filters-open:has(.two-up)>.scroll-filters{transform:none;}
.scroll-layout:has(.two-up) .scroll-filters .scroll-ctl-body .scroll-ctl-body{grid-template-columns:1fr;}
.scroll-panel-card{scroll-margin-top:calc(var(--sc-appbar-h,70px) + 14px); border:1px solid var(--sc-line); border-radius:12px;
  container-type:inline-size; container-name:sc-card; min-width:0;
  background:var(--sc-card); box-shadow:0 6px 20px -12px rgba(17,24,38,.18); overflow:hidden;}
.scroll-panel-card>.card-header{background:var(--sc-card); border-bottom:1px solid var(--sc-line-2); padding:18px 22px;}
.scroll-eyebrow{display:flex; align-items:center; gap:10px;}
.scroll-num{font-family:ui-monospace,monospace; font-size:12px; font-weight:600; color:var(--sc-accent);
  background:var(--sc-wash); padding:2px 7px; border-radius:6px;}
.scroll-kicker{text-transform:uppercase; letter-spacing:.1em; font-size:12px; font-weight:700; color:var(--sc-muted);}
/* cap the text measure so title/description don't stretch across a full-bleed card */
.scroll-title{margin:8px 0 2px; font-size:22px; font-weight:800; letter-spacing:-.02em; max-width:48ch;}
.scroll-desc{margin:0; color:var(--sc-muted); font-size:14px; max-width:74ch;}

/* controls */
.scroll-panel{padding:20px 22px;}
:root{--sc-ctl-w:clamp(232px, 22%, 320px);}   /* control-column width knob (narrower = tighter inputs + more plot) */
.scroll-controls{display:flex; flex-direction:column; gap:12px;}
.scroll-cgroup-h{text-transform:uppercase; letter-spacing:.08em; font-size:11px; font-weight:700;
  color:var(--sc-faint); padding-bottom:6px; margin-bottom:8px; border-bottom:1px solid var(--sc-line-2);}
.scroll-controls .form-label{font-size:13px; font-weight:600; margin-bottom:2px;}
.scroll-controls .form-group,.scroll-controls .shiny-input-container{margin-bottom:8px;}
.scroll-controls .bslib-input-switch{margin-bottom:6px;}       /* switches are short -> less gap */
/* Responsive panel split (container query on .scroll-content). Side-by-side only
   when the card is wide enough to keep the plot usable (>=720px = ~78 chrome + 240
   controls + 24 gap + 378 plot); narrower cards -- a minimized window, a portrait
   monitor -- STACK: controls become a capped multi-column box ABOVE a full-width
   plot, so the plot is never squeezed to an unusable sliver. One rule adapts to
   filters on/off, collapsed rail, or multi-app tabs. bslib lays c(3,9) out as a
   12-track grid with g-col-sm-* spans; both branches override that. */
@container sc-card (min-width:720px){
  .scroll-panel.bslib-grid{grid-template-columns:var(--sc-ctl-w) minmax(0,1fr) !important;}
  .scroll-panel.bslib-grid>.bslib-grid-item{grid-column:auto !important;}
}
@container sc-card (max-width:719.98px){
  .scroll-panel.bslib-grid{grid-template-columns:minmax(0,1fr) !important; gap:14px;}
  .scroll-panel.bslib-grid>.bslib-grid-item{grid-column:auto !important;}
  .scroll-controls{max-height:min(42vh,420px);
    display:grid; grid-template-columns:repeat(auto-fill,minmax(210px,1fr)); gap:0 18px; align-items:start;}
  .scroll-controls>.btn,.scroll-controls>button,.scroll-controls>.bslib-input-task-button{grid-column:1/-1;}
  .scroll-controls .selectize-input{max-height:120px;}
}
@container sc-card (max-width:479.98px){   /* phones: reclaim card chrome */
  .scroll-panel{padding:12px 10px;}
  .scroll-panel-card>.card-body{padding:8px;}
  .scroll-panel-card>.card-header{padding:14px 16px;}
}
/* per-card Controls toggle: hidden by default, revealed only when the card is
   stacked (narrow) so a reader can collapse the controls to just the plot. */
.scroll-ctl-btn{display:none; margin-left:auto; cursor:pointer; font:600 11px/1 inherit;
  padding:5px 10px; border:1px solid var(--sc-line); border-radius:7px; background:#fff;
  color:var(--sc-muted);}
.scroll-ctl-btn:hover{color:var(--sc-accent-deep); border-color:var(--sc-accent);}
.scroll-ctl-btn::after{content:' \\25BE'; opacity:.7;}                 /* down-triangle: open */
.scroll-ctl-btn[aria-expanded='false']::after{content:' \\25B8';}     /* right-triangle: collapsed */
@container sc-card (max-width:719.98px){
  .scroll-panel-card:has(.scroll-panel.bslib-grid) .scroll-ctl-btn{display:inline-flex;}
  .scroll-panel-card.scroll-ctl-hidden .scroll-panel.bslib-grid>.bslib-grid-item:first-child{display:none;}
}
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
/* Style button: opens this plot's side sheet (R/style.R) */
.scroll-style-btn{display:inline-flex; align-items:center; justify-content:center; width:28px;
  height:26px; padding:0; font-size:12px; cursor:pointer; color:var(--sc-muted);
  background:var(--sc-card); border:1px solid var(--sc-line); border-radius:8px;}
.scroll-style-btn:hover,.scroll-style-btn.is-open{color:var(--sc-accent-deep);
  border-color:var(--sc-accent); background:var(--sc-wash);}
/* Per-plot Style sheet. Fixed to the right edge under the app bar, overlaying the
   page (no reflow, so opening it never resizes -- and re-renders -- the plots).
   Closed = translated off-screen + visibility:hidden (+ inert), NOT display:none, so
   controls inside keep a real size (sliders measure their width at bind time). */
.scroll-sheet{position:fixed; top:var(--sc-appbar-h,70px); right:0; z-index:1060;
  width:min(360px,92vw); height:calc(100dvh - var(--sc-appbar-h,70px));
  display:flex; flex-direction:column; background:var(--sc-card);
  border-left:1px solid var(--sc-line); box-shadow:-14px 0 34px -18px rgba(17,24,38,.4);
  transform:translateX(105%); visibility:hidden;
  transition:transform .18s ease, visibility 0s linear .18s;}
.scroll-sheet.is-open{transform:none; visibility:visible; transition:transform .18s ease;}
.scroll-sheet-head{display:flex; align-items:center; gap:10px; padding:12px 16px;
  border-bottom:1px solid var(--sc-line);}
.scroll-sheet-titles{display:flex; flex-direction:column; min-width:0; flex:1;}
.scroll-sheet-kicker{text-transform:uppercase; letter-spacing:.1em; font-size:11px;
  font-weight:700; color:var(--sc-muted);}
.scroll-sheet-title{font-size:16px; font-weight:800; letter-spacing:-.01em;
  white-space:nowrap; overflow:hidden; text-overflow:ellipsis;}
.scroll-sheet-close{flex:none; width:30px; height:30px; border:1px solid var(--sc-line);
  border-radius:8px; background:var(--sc-card); color:var(--sc-muted); cursor:pointer;}
.scroll-sheet-close:hover{color:var(--sc-ink); border-color:var(--sc-accent);}
.scroll-sheet-body{flex:1; overflow-y:auto; overscroll-behavior:contain; padding:6px 16px 12px;}
.scroll-sheet-foot{display:flex; align-items:center; gap:12px; padding:10px 16px;
  border-top:1px solid var(--sc-line);}
.scroll-sheet-foot a{font-size:12px; color:var(--sc-accent); text-decoration:none; margin-left:auto;}
.scroll-sheet .form-select,.scroll-sheet .form-control{font-size:12px; padding:3px 8px;
  min-height:0; height:28px;}
.scroll-sheet .form-select{padding-right:24px;}
/* the moved per-panel look controls (Plot section): compact, full width */
.scroll-sheet .shiny-input-container{width:100%; margin-bottom:8px;}
.scroll-sheet .control-label,.scroll-sheet .form-label{font-size:12px; font-weight:600;
  margin-bottom:2px;}
.scroll-sheet .bslib-input-switch{margin-bottom:6px; font-size:13px;}
.scroll-sheet .irs{font-size:10px;}

/* tighter app bar on smaller screens (padding follows the layout's 28/20/16 steps) */
@media (max-width:1199.98px){.scroll-appbar{padding:10px 20px; column-gap:14px;}}
@media (max-width:900px){.scroll-appbar{padding:10px 16px;}}
@media (max-width:576px){.scroll-dataset{font-size:16px;} .scroll-cells{display:none;}}
/* Below 1400 the right control rail becomes an off-canvas DRAWER (positioning only;
   grid tracks are set per-range below, so they never fight on specificity). Toggled
   by .filters-open (the app-bar button). Including .controls-collapsed keeps a
   wide-mode hidden state from killing the drawer after a resize. */
@media (max-width:1399.98px){
  .scroll-layout.has-filters>.scroll-filters,
  .scroll-layout.controls-collapsed>.scroll-filters{
    position:fixed; top:var(--sc-appbar-h,70px); right:0; z-index:1050;
    display:flex; width:min(300px,90vw);
    /* explicit height (not top+bottom stretch) so overflow-y actually scrolls the
       content -- a tall filters list was overflowing off-screen and unreachable. */
    height:calc(100dvh - var(--sc-appbar-h,70px)); max-height:none; padding:16px 16px 24px;
    background:var(--sc-card); border-left:1px solid var(--sc-line);
    box-shadow:-14px 0 34px -18px rgba(17,24,38,.4);
    overflow-y:auto; -webkit-overflow-scrolling:touch; overscroll-behavior:contain;
    transform:translateX(105%); transition:transform .18s ease;}
  .scroll-layout.filters-open>.scroll-filters{transform:none;}
  /* thinner drawer -> theme sub-controls stack one per row (2-up is too wide here) */
  .scroll-filters .scroll-ctl-body .scroll-ctl-body{grid-template-columns:1fr;}
}
/* laptop band: drop the 300px filters track (it is a drawer now) */
@media (min-width:1200px) and (max-width:1399.98px){
  .scroll-layout.has-filters{grid-template-columns:var(--sc-rail-w) minmax(0,1fr); max-width:none;}
}
/* mid widths (portrait monitor / small landscape): a compact numeric rail frees
   ~136px for the plot; labels collapse to the number (hover shows the full label via
   the title attr). Applies to every layout -- filters are a drawer at this width. */
@media (min-width:900.02px) and (max-width:1199.98px){
  .scroll-layout,.scroll-layout.has-filters,.scroll-layout.controls-collapsed{
    grid-template-columns:64px minmax(0,1fr); gap:20px; padding:20px; max-width:none;}
  .scroll-rail-item{justify-content:center; padding:8px 0; gap:0;}
  .scroll-rail-item>span:last-child,.scroll-rail-item::before{display:none;}   /* keep .scroll-rail-num */
  .scroll-rail-foot{display:none;}
}
@media (max-width:900px){
  .scroll-layout,.scroll-layout.has-filters,.scroll-layout.controls-collapsed{
    grid-template-columns:1fr; gap:16px; padding:16px;}
  .scroll-rail{position:static; flex-direction:row; overflow-x:auto; top:auto;}
  .scroll-rail-item{cursor:default;}
  .scroll-rail-item::before{display:none;}     /* no drag on the horizontal strip */
  .scroll-rail-foot{display:none;}
}
"

.scroll_spy_js <- function() paste0(.scroll_tips_js(), "
// App-bar button labels: what the button does in its CURRENT state. Computed on
// hover/focus rather than set once, because the controls button changes role with
// the screen width (in-grid column vs drawer) and a resize would leave it stale.
function scrollTip(btn){
  var T=window.SCROLL_TIPS||{}, t=null, lay=null;
  if(btn.classList.contains('scroll-ctl-toggle')){
    var body=btn.closest('.scroll-appbar')&&btn.closest('.scroll-appbar').parentElement;
    lay=body?body.querySelector('.scroll-layout'):document.querySelector('.scroll-layout');
    var drawer=window.matchMedia('(max-width:1399.98px)').matches ||
               (lay&&lay.querySelector('.scroll-content.two-up'));
    t = drawer ? (lay&&lay.classList.contains('filters-open') ? T.drawer_close : T.drawer_open)
               : (lay&&lay.classList.contains('controls-collapsed') ? T.ctl_show : T.ctl_hide);
  } else if(btn.classList.contains('scroll-two-toggle')){
    t = btn.classList.contains('is-on') ? T.two_off : T.two_on;
  } else if(btn.classList.contains('scroll-info-toggle')){ t = T.info; }
  if(t){ btn.setAttribute('data-tip', t); btn.setAttribute('aria-label', t); }
}
['pointerover','focusin'].forEach(function(ev){
  document.addEventListener(ev, function(e){
    var b=e.target.closest&&e.target.closest('.scroll-appbar-actions [data-tip]');
    if(b) scrollTip(b);
  });
});
function scrollCloseDrawers(){
  document.querySelectorAll('.scroll-layout.filters-open').forEach(function(l){
    l.classList.remove('filters-open');
    var t=l.parentElement&&l.parentElement.querySelector('.scroll-ctl-toggle');
    if(t){t.classList.remove('is-open'); t.setAttribute('aria-expanded','false');}
  });
}
window.scrollToggleControls=function(btn){
  var body=btn.closest('.scroll-appbar').parentElement;
  var lay=body?body.querySelector('.scroll-layout'):document.querySelector('.scroll-layout');
  if(!lay) return;
  // filters live in a drawer below 1400, and also whenever two-up is on (any width)
  if(window.matchMedia('(max-width:1399.98px)').matches || lay.querySelector('.scroll-content.two-up')){
    scrollCloseSheets();                                 // the drawer and a sheet share the right edge
    var open=lay.classList.toggle('filters-open');
    btn.classList.toggle('is-open', open);
    btn.setAttribute('aria-expanded', String(open));
    scrollTip(btn);
    return;
  }
  var collapsed=lay.classList.toggle('controls-collapsed');  // wide: hide/show the in-grid rail
  btn.classList.toggle('is-collapsed', collapsed);
  btn.setAttribute('aria-expanded', String(!collapsed));
  scrollTip(btn);
};
window.scrollTogglePanelControls=function(btn){   // per-card controls collapse (narrow only)
  var card=btn.closest('.scroll-panel-card'); if(!card) return;
  var hidden=card.classList.toggle('scroll-ctl-hidden');
  btn.setAttribute('aria-expanded', String(!hidden));
  var ctl=card.querySelector('.scroll-panel.bslib-grid>.bslib-grid-item');
  if(ctl && window.jQuery) jQuery(ctl).trigger(hidden?'hidden':'shown');  // suspend/resume any outputs inside
};
// Per-plot Style sheets (R/style.R). One sheet open at a time; opening one closes
// the filters drawer. The FIRST open tells the server to insert the sheet's
// sections (seeded from the panel's current style). Sheets are non-modal and are
// NOT closed by an outside click -- the point is to click around while editing.
function scrollCloseSheets(){
  document.querySelectorAll('.scroll-sheet.is-open').forEach(function(sh){
    sh.classList.remove('is-open'); sh.setAttribute('inert','');
  });
  document.querySelectorAll('.scroll-style-btn.is-open').forEach(function(b){
    b.classList.remove('is-open'); b.setAttribute('aria-expanded','false');
  });
}
window.scrollCloseSheets=scrollCloseSheets;
window.scrollToggleStyle=function(btn){
  var sh=document.getElementById(btn.getAttribute('data-sheet')); if(!sh) return;
  var wasOpen=sh.classList.contains('is-open');
  scrollCloseSheets(); scrollCloseDrawers();
  if(wasOpen) return;
  sh.classList.add('is-open'); sh.removeAttribute('inert');
  btn.classList.add('is-open'); btn.setAttribute('aria-expanded','true');
  if(!sh.getAttribute('data-opened') && window.Shiny && Shiny.setInputValue){
    sh.setAttribute('data-opened','1');
    Shiny.setInputValue(btn.getAttribute('data-open'), Date.now(), {priority:'event'});
  }
  if(window.jQuery) jQuery(sh).trigger('shown');
  // keep the plot being styled on screen (an off-screen card stops rendering live)
  var card=btn.closest('.scroll-panel-card');
  if(card) card.scrollIntoView({block:'nearest', behavior:'smooth'});
};
document.addEventListener('keydown',function(e){
  if(e.key==='Escape'){ scrollCloseDrawers(); scrollCloseSheets(); } });
document.addEventListener('click',function(e){   // click outside the drawer (or its toggle) closes it
  if(e.target.closest('.scroll-filters')||e.target.closest('.scroll-ctl-toggle')||
     e.target.closest('.scroll-sheet')||e.target.closest('.scroll-style-btn')) return;
  scrollCloseDrawers();
});
/* Keep --sc-appbar-h in sync with the app bar's measured bottom, so sticky offsets
   (rail, filters, scroll-margin) and the drawer top stay correct when it wraps. The
   bottom already includes a multi-app sticky tab strip. Observes size only -> never
   changes a width, so it cannot feed the plot ResizeObserver loop. */
(function(){
  function setH(){ var bar=document.querySelector('.scroll-appbar'); if(!bar) return;
    document.documentElement.style.setProperty('--sc-appbar-h',
      Math.round(bar.getBoundingClientRect().bottom)+'px'); }
  function wire(){ var bar=document.querySelector('.scroll-appbar'); if(!bar) return;
    if(window.ResizeObserver) new ResizeObserver(setH).observe(bar); setH(); }
  window.addEventListener('resize', setH);
  if(document.readyState!=='loading') wire(); else document.addEventListener('DOMContentLoaded', wire);
})();
// Two panels per row: toggle .two-up on the (nearest) .scroll-content, remembered
// across reloads. auto-fit self-degrades to one column when there isn't room.
window.scrollToggleTwoUp=function(btn){
  var body=btn.closest('.scroll-appbar').parentElement;
  var c=body?body.querySelector('.scroll-content'):document.querySelector('.scroll-content');
  if(!c) return;
  var on=c.classList.toggle('two-up');
  btn.classList.toggle('is-on', on); btn.setAttribute('aria-pressed', String(on));
  scrollTip(btn);
  scrollCloseDrawers();                                 // filters drawer mode may change -> start closed
  try{ localStorage.setItem('scroll:two-up', on?'1':'0'); }catch(e){}
  scrollAdjustTables();                                 // cards changed width -> re-fit any DataTables
};
// A card's width flips (full <-> half) when two-up toggles, but a DataTable only
// lays out its column widths once; nudge every visible DataTable to recompute.
function scrollAdjustTables(){
  requestAnimationFrame(function(){
    if(!(window.jQuery && jQuery.fn && jQuery.fn.dataTable)) return;
    try{ jQuery.fn.dataTable.tables({visible:true, api:true}).columns.adjust(); }catch(e){}
  });
}
(function(){
  function apply(){ var v=null; try{ v=localStorage.getItem('scroll:two-up'); }catch(e){}
    if(v===null) return;                              // no stored pref -> keep the server (config) default
    var on = v==='1';
    document.querySelectorAll('.scroll-content').forEach(function(c){ c.classList.toggle('two-up', on); });
    document.querySelectorAll('.scroll-two-toggle').forEach(function(b){
      b.classList.toggle('is-on', on); b.setAttribute('aria-pressed', String(on));
      scrollTip(b); });
  }
  if(document.readyState!=='loading') apply(); else document.addEventListener('DOMContentLoaded',apply);
})();
(function(){
  function spy(){
    var items=document.querySelectorAll('.scroll-rail-item');
    var secs=document.querySelectorAll('.scroll-panel-card');
    if(!secs.length) return;
    var visible={};                                     // set of intersecting card ids
    var io=new IntersectionObserver(function(es){
      es.forEach(function(e){ if(e.isIntersecting) visible[e.target.id]=1; else delete visible[e.target.id]; });
      // highlight EVERY in-view panel -- so in two-up both side-by-side row-mates
      // light up, not just the first. Fall back to keeping the last active if the
      // centre band momentarily sees none (e.g. between cards).
      var any=false;
      items.forEach(function(it){ var on=!!visible[it.getAttribute('href').slice(1)];
        if(on) any=true; });
      if(any) items.forEach(function(it){
        it.classList.toggle('active', !!visible[it.getAttribute('href').slice(1)]); });
    },{rootMargin:'-45% 0px -45% 0px'});
    secs.forEach(function(s){io.observe(s);});
  }
  if(document.readyState!=='loading') spy(); else document.addEventListener('DOMContentLoaded',spy);
})();
// Drag-to-reorder panels via the section rail: drag a rail item, the panel cards
// mirror the new order, both are renumbered, and the order persists (keyed by the
// panel-id set, so it survives reloads and self-invalidates if the panels change).
// Moving the (Shiny-bound) card DOM does not re-render or rebind anything.
(function(){
  function idsOf(rail){ return [].slice.call(rail.querySelectorAll('.scroll-rail-item'))
    .map(function(a){ return a.getAttribute('href').slice(1); }); }
  function keyOf(rail){ return 'scroll-order:' + idsOf(rail).slice().sort().join(','); }
  function contentOf(rail){ var lay=rail.closest('.scroll-layout'); return lay&&lay.querySelector('.scroll-content'); }
  function renumber(rail){
    rail.querySelectorAll('.scroll-rail-item').forEach(function(a,i){
      var n=String(i+1).padStart(2,'0');
      var rn=a.querySelector('.scroll-rail-num'); if(rn) rn.textContent=n;
      var card=document.getElementById(a.getAttribute('href').slice(1));
      var cn=card&&card.querySelector('.scroll-num'); if(cn) cn.textContent=n;
    });
  }
  function mirror(rail){ var c=contentOf(rail); if(!c) return;
    idsOf(rail).forEach(function(id){ var card=document.getElementById(id); if(card) c.appendChild(card); });
    renumber(rail); }
  function persist(rail){ try{ localStorage.setItem(keyOf(rail), idsOf(rail).join(',')); }catch(e){} }
  function applyStored(){
    document.querySelectorAll('.scroll-rail').forEach(function(rail){
      var s; try{ s=localStorage.getItem(keyOf(rail)); }catch(e){ return; } if(!s) return;
      s.split(',').forEach(function(id){
        var a=rail.querySelector('.scroll-rail-item[href=\"#'+id+'\"]'); if(a) rail.appendChild(a); });
      mirror(rail);
    });
  }
  var dragEl=null, dragged=false;
  // FLIP: run `mutate` (a DOM reorder) and slide every rail item from its old
  // position to its new one, so the list visibly reflows as the dragged item is
  // pulled out and re-inserted (rather than a static before/after highlight line).
  function flip(rail, mutate){
    var items=[].slice.call(rail.querySelectorAll('.scroll-rail-item'));
    var y0={}; items.forEach(function(it){ y0[it.getAttribute('href')]=it.getBoundingClientRect().top; });
    mutate();
    items.forEach(function(it){
      var dy=y0[it.getAttribute('href')]-it.getBoundingClientRect().top; if(!dy) return;
      it.style.transition='none'; it.style.transform='translateY('+dy+'px)';
      requestAnimationFrame(function(){ it.style.transition='transform .16s ease'; it.style.transform=''; });
    });
  }
  function commit(){   // apply the (live-reordered) rail order to the cards + persist, then clean up
    var rail=dragEl && dragEl.closest('.scroll-rail');
    if(rail){ mirror(rail); persist(rail); }
    document.querySelectorAll('.scroll-rail-item').forEach(function(x){
      x.classList.remove('is-dragging'); x.style.transition=''; x.style.transform=''; });
    dragEl=null; setTimeout(function(){ dragged=false; }, 60);
  }
  document.addEventListener('dragstart', function(e){
    var it=e.target.closest('.scroll-rail-item'); if(!it) return;
    if(window.matchMedia('(max-width:900px)').matches){ e.preventDefault(); return; }  // horizontal strip
    dragEl=it; dragged=true;
    e.dataTransfer.effectAllowed='move';
    try{ e.dataTransfer.setData('text/plain', it.getAttribute('href')); e.dataTransfer.setDragImage(it,10,10); }catch(_){}
    requestAnimationFrame(function(){ if(dragEl) dragEl.classList.add('is-dragging'); });  // dim after the ghost is captured
  });
  document.addEventListener('dragover', function(e){
    if(!dragEl || e.target.closest('.scroll-rail')!==dragEl.closest('.scroll-rail')) return;
    e.preventDefault(); e.dataTransfer.dropEffect='move';   // ALWAYS allow a drop while over the rail
    var it=e.target.closest('.scroll-rail-item'); if(!it || it===dragEl) return;
    var rail=dragEl.closest('.scroll-rail'); var r=it.getBoundingClientRect();
    var ref=e.clientY > r.top + r.height/2 ? it.nextSibling : it;   // slot to drop into
    if(ref===dragEl || dragEl.nextSibling===ref) return;            // already there -> no-op (avoids FLIP thrash)
    flip(rail, function(){ rail.insertBefore(dragEl, ref); });      // live take-out & re-insert
  });
  document.addEventListener('drop', function(e){ if(dragEl){ e.preventDefault(); commit(); } });
  document.addEventListener('dragend', commit);   // fallback if drop didn't fire (dropped off a valid target)
  document.addEventListener('click', function(e){   // a click that ended a drag must not navigate
    if(dragged && e.target.closest('.scroll-rail-item')) e.preventDefault(); }, true);
  window.scrollResetOrder=function(el){
    var rail=el&&el.closest('.scroll-rail'); if(rail){ try{ localStorage.removeItem(keyOf(rail)); }catch(e){} }
    location.reload();
  };
  if(document.readyState!=='loading') applyStored(); else document.addEventListener('DOMContentLoaded', applyStored);
})();
")

# App-bar button labels, shared by the R markup (initial data-tip / aria-label) and
# the JS that recomputes them from live state (window.SCROLL_TIPS).
.SCROLL_TIPS <- list(
  info         = "Dataset details: cells, genes, assays and reductions",
  two_on       = "Show two panels side by side",
  two_off      = "Back to one panel per row",
  ctl_hide     = "Hide the filters column (more room for plots)",
  ctl_show     = "Show the filters column",
  drawer_open  = "Open the filters panel",
  drawer_close = "Close the filters panel",
  style        = "Style this plot: colours, theme, labels")

.scroll_tips_js <- function() {
  kv <- vapply(names(.SCROLL_TIPS), function(k)
    sprintf("%s:'%s'", k, gsub("'", "\\\\'", .SCROLL_TIPS[[k]])), "")
  sprintf("window.SCROLL_TIPS={%s};\n", paste(kv, collapse = ","))
}

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
