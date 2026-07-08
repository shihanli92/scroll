// scroll: bridge closeread scroll position -> Shiny.
//
// closeread renders each narrative trigger as `<div class="trigger new-trigger"
// data-focus-on="...">` and drives them with scrollama. We can't rely on
// closeread's OJS internals from Shiny, so we run our own IntersectionObserver
// over the `.new-trigger` elements and, when one reaches the viewport centre,
// push its section id to Shiny as the global input `active_section`.
//
// Each authored trigger carries an inline `<span data-section="ID">` (see the
// scaffolded story.qmd); that ID names the config section this trigger should
// activate. This is the one place scroll position enters Shiny; the feature
// search box and toggles are ordinary Shiny inputs, so the drivers meet only
// inside the shared-visual reactive and never contend for the DOM.
(function () {
  function setActive(id) {
    if (window.Shiny && id) {
      Shiny.setInputValue("active_section", id, { priority: "event" });
    }
  }

  function sectionOf(triggerEl) {
    var span = triggerEl.querySelector("[data-section]");
    return span ? span.getAttribute("data-section") : null;
  }

  function wire(attempt) {
    var triggers = Array.prototype.slice.call(
      document.querySelectorAll(".cr-section .new-trigger")
    );
    if (!triggers.length) {
      // closeread may build/restructure triggers after load; retry briefly.
      if ((attempt || 0) < 20) setTimeout(function () { wire((attempt || 0) + 1); }, 150);
      return;
    }

    // Seed with the first trigger so the sticky has a section before any scroll.
    setActive(sectionOf(triggers[0]));

    var observer = new IntersectionObserver(
      function (entries) {
        entries.forEach(function (entry) {
          if (entry.isIntersecting) setActive(sectionOf(entry.target));
        });
      },
      // A thin band across the viewport centre: the trigger crossing it wins.
      { rootMargin: "-45% 0px -45% 0px", threshold: 0 }
    );
    triggers.forEach(function (t) { observer.observe(t); });
  }

  if (document.readyState !== "loading") wire();
  else document.addEventListener("DOMContentLoaded", wire);
})();
