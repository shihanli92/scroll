// scroll: bridge closeread scroll position -> Shiny.
//
// closeread marks the trigger paragraph in focus by toggling a class on the
// step element. We observe those steps and push the id of the centered one to
// Shiny as the global input `active_section`, which the single shared-visual
// reactive depends on. This is the one place scroll position enters Shiny; the
// feature search box is an ordinary Shiny input, so the two drivers meet only
// inside that reactive and never contend for the DOM.
(function () {
  function setActive(id) {
    if (window.Shiny && id) {
      Shiny.setInputValue("active_section", id, { priority: "event" });
    }
  }

  function triggerId(el) {
    // The authored id lives on the trigger span (e.g. {#overview}); fall back to
    // any id closeread stamps on the step element.
    var inner = el.querySelector("[id]");
    return (inner && inner.id) || el.id || null;
  }

  function wire() {
    // closeread trigger paragraphs carry the `.trigger` class within a
    // `.cr-section`. Observe whichever is nearest the viewport centre.
    var steps = document.querySelectorAll(".cr-section .trigger, .cr-section .step");
    if (!steps.length) return;

    var observer = new IntersectionObserver(
      function (entries) {
        entries.forEach(function (entry) {
          if (entry.isIntersecting) setActive(triggerId(entry.target));
        });
      },
      { rootMargin: "-45% 0px -45% 0px", threshold: 0 }
    );
    steps.forEach(function (s) { observer.observe(s); });
  }

  if (document.readyState !== "loading") wire();
  else document.addEventListener("DOMContentLoaded", wire);
})();
