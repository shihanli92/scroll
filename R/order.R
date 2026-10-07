# --- Group label order ----------------------------------------------------------

# Categorical plots (DotPlot, Violin, Ridge, Heatmap, Proportions) put their groups
# in a natural order by default -- numbers by value, so cluster 2 comes before 10 --
# and the panel's Style sheet holds a drag-to-reorder list of the group labels for a
# manual order. The list is a small custom Shiny input (`scroll.order`, JS below);
# its value is "" until the user reorders (so an untouched list is not a saved
# style), then the full label order.

# Natural ("human") order of the distinct values of `x`: digit runs compare by
# value, case-insensitively (cluster_2 < cluster_10; "B cells" < "b2" < "T cells").
.scroll_natural_sort <- function(x) {
  x <- unique(as.character(x[!is.na(x)]))
  if (length(x) < 2L) return(x)
  num <- suppressWarnings(as.numeric(x))
  if (!anyNA(num)) return(x[order(num)])                # all numeric: by value
  key <- tolower(x)
  digits <- gregexpr("[0-9]+", key)
  runs <- regmatches(key, digits)
  regmatches(key, digits) <- lapply(runs, function(r) paste0(strrep("0", pmax(0, 20 - nchar(r))), r))
  x[order(key, x)]
}

# The levels of `lv` in display order: a manual order first (labels it names that
# are present), then the rest in natural order.
.scroll_level_order <- function(lv, manual = NULL) {
  lv <- .scroll_natural_sort(lv)
  manual <- as.character(unlist(manual))
  manual <- manual[nzchar(manual) & manual %in% lv]
  unique(c(manual, lv))
}

# A drag list's value as a character vector ("" -> character(0): no manual order).
.scroll_order_value <- function(x) {
  x <- as.character(unlist(x))
  x[!is.na(x) & nzchar(x)]
}

# Most groups a drag list will show; past this it says so instead (a list of
# hundreds of chips is no way to order anything).
.SCROLL_ORDER_MAX <- 80L

# The drag-to-reorder control. `id` is namespaced by the caller.
.scroll_order_input <- function(id, label = "Group order") {
  div(class = "form-group shiny-input-container scroll-order-wrap",
      tags$label(class = "control-label", `for` = id, label),
      div(id = id, class = "scroll-order", `data-max` = .SCROLL_ORDER_MAX,
          span(class = "scroll-order-empty", "Groups appear here.")),
      div(class = "scroll-order-foot",
          span(class = "scroll-order-hint", "Drag to reorder"),
          tags$a(href = "#", class = "scroll-order-reset", `data-for` = id, "Reset")))
}

# Keep a drag list's labels in step with the plot's groups. `levels_r` returns the
# current group labels (any order); the list shows them naturally sorted, and keeps
# a manual order across group changes (labels no longer present just drop out).
.scroll_bind_order <- function(session, id, levels_r) {
  observe({
    lv <- tryCatch(levels_r(), error = function(e) character(0))
    session$sendInputMessage(id, list(levels = as.list(.scroll_natural_sort(lv))))
  })
}

# The group labels of `cols` (one column, or the "a | b" interaction of several) on
# the active cells.
.scroll_group_levels <- function(cells, cols) {
  cols <- cols[cols %in% names(cells)]
  if (!length(cols)) return(character(0))
  .scroll_combo_levels(cells, cols)
}

.scroll_order_css <- function() "
.scroll-order{display:flex; flex-wrap:wrap; gap:4px; padding:6px; min-height:34px;
  border:1px solid var(--sc-line); border-radius:6px; background:var(--sc-card);}
.scroll-order-item{display:inline-flex; align-items:center; padding:2px 9px; font-size:12px;
  line-height:18px; border:1px solid var(--sc-line); border-radius:11px; background:var(--sc-ground);
  color:var(--sc-ink); cursor:grab; user-select:none; max-width:100%; overflow:hidden;
  text-overflow:ellipsis; white-space:nowrap;}
.scroll-order-item:focus{outline:2px solid var(--sc-accent); outline-offset:1px;}
.scroll-order-item.is-dragging{opacity:.35;}
.scroll-order-item.is-moved{border-color:var(--sc-accent); background:var(--sc-wash);}
.scroll-order-empty{font-size:11px; color:var(--sc-faint); padding:1px 2px;}
.scroll-order-foot{display:flex; justify-content:space-between; margin-top:3px;
  font-size:11px; color:var(--sc-faint);}
.scroll-order-reset{color:var(--sc-accent); text-decoration:none;}
"

# Input binding: value "" (default order) or the full label order once the user has
# moved something. Messages: {levels: [...]} from the server, {value: [...] | ""}
# from a saved style / Reset. Drag with the mouse, or focus a label and use the
# left / right arrow keys.
.scroll_order_js <- function() "
(function(){
  function want(el){ return Array.isArray(el.__want) ? el.__want : []; }
  function render(el){
    var lv = el.__levels || [], max = +el.getAttribute('data-max') || 80;
    var w = want(el).filter(function(v){ return lv.indexOf(v) >= 0; });
    var order = w.concat(lv.filter(function(v){ return w.indexOf(v) < 0; }));
    el.__custom = w.length > 0;
    el.innerHTML = '';
    if(!lv.length || lv.length > max){
      var s = document.createElement('span'); s.className = 'scroll-order-empty';
      s.textContent = lv.length ? (lv.length + ' groups: too many to order here.') : 'Groups appear here.';
      el.appendChild(s); el.__custom = false; return;
    }
    order.forEach(function(v){
      var c = document.createElement('span');
      c.className = 'scroll-order-item' + (w.indexOf(v) >= 0 ? ' is-moved' : '');
      c.textContent = v; c.title = v; c.setAttribute('data-value', v);
      c.draggable = true; c.tabIndex = 0;
      el.appendChild(c);
    });
  }
  function commit(el){
    el.__want = Array.prototype.map.call(el.querySelectorAll('.scroll-order-item'),
                                         function(c){ return c.getAttribute('data-value'); });
    render(el); $(el).trigger('change');
  }
  function wire(el){
    if(el.__wired) return; el.__wired = true;
    var drag = null;
    el.addEventListener('dragstart', function(e){
      var c = e.target.closest('.scroll-order-item'); if(!c) return;
      drag = c; c.classList.add('is-dragging');
      e.dataTransfer.effectAllowed = 'move';
      try { e.dataTransfer.setData('text/plain', c.getAttribute('data-value')); } catch(_){}
    });
    el.addEventListener('dragover', function(e){
      if(!drag) return; e.preventDefault();
      var c = e.target.closest('.scroll-order-item'); if(!c || c === drag) return;
      var r = c.getBoundingClientRect();
      var after = (e.clientY > r.bottom) || (e.clientY >= r.top && e.clientX > r.left + r.width / 2);
      el.insertBefore(drag, after ? c.nextSibling : c);
    });
    el.addEventListener('drop', function(e){ if(drag) e.preventDefault(); });
    el.addEventListener('dragend', function(){
      if(!drag) return; drag.classList.remove('is-dragging'); drag = null; commit(el);
    });
    el.addEventListener('keydown', function(e){
      var c = e.target.closest('.scroll-order-item'); if(!c) return;
      var v = c.getAttribute('data-value'), sib;
      if(e.key === 'ArrowLeft' || e.key === 'ArrowUp') sib = c.previousElementSibling;
      else if(e.key === 'ArrowRight' || e.key === 'ArrowDown') sib = c.nextElementSibling;
      else return;
      e.preventDefault(); if(!sib) return;
      el.insertBefore(c, (e.key === 'ArrowLeft' || e.key === 'ArrowUp') ? sib : sib.nextSibling);
      commit(el);
      var again = el.querySelector('.scroll-order-item[data-value=\"' + CSS.escape(v) + '\"]');
      if(again) again.focus();
    });
  }
  document.addEventListener('click', function(e){
    var a = e.target.closest && e.target.closest('.scroll-order-reset'); if(!a) return;
    e.preventDefault();
    var el = document.getElementById(a.getAttribute('data-for')); if(!el) return;
    el.__want = ''; render(el); $(el).trigger('change');
  });
  function reg(){
    if(window.__scrollOrderReg) return;
    if(!window.Shiny || !Shiny.InputBinding || !Shiny.inputBindings) return;
    window.__scrollOrderReg = true;
    var b = new Shiny.InputBinding();
    $.extend(b, {
      find: function(scope){ return $(scope).find('.scroll-order'); },
      initialize: function(el){ wire(el); },
      getValue: function(el){
        if(!el.__custom) return '';
        return Array.prototype.map.call(el.querySelectorAll('.scroll-order-item'),
                                        function(c){ return c.getAttribute('data-value'); });
      },
      setValue: function(el, v){ el.__want = v; render(el); },
      receiveMessage: function(el, m){
        if(m && m.hasOwnProperty('levels')) el.__levels = m.levels || [];
        if(m && m.hasOwnProperty('value')) el.__want = m.value;
        wire(el); render(el); $(el).trigger('change');
      },
      subscribe: function(el, cb){ $(el).on('change.scrollOrder', function(){ cb(); }); },
      unsubscribe: function(el){ $(el).off('.scrollOrder'); }
    });
    Shiny.inputBindings.register(b, 'scroll.order');
  }
  if(document.readyState !== 'loading') reg();
  document.addEventListener('DOMContentLoaded', reg);
})();
"
