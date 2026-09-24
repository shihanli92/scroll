#' @keywords internal
"_PACKAGE"

#' @importFrom ggplot2 .data
#' @importFrom methods is
#' @import shiny
NULL

# Soft-deprecation for exports being retired from the public API (0.2.19 -> 0.3.0).
# Warns only when called from OUTSIDE the package -- the built-in panels, scroll_de,
# and tests call these internally and must stay silent -- and once per session per
# function. `env` is the deprecated function's parent.frame() (its caller).
# `cast` is an arrow dplyr binding used inside dplyr::mutate() on an arrow query
utils::globalVariables("cast")

.scroll_dep_seen <- new.env(parent = emptyenv())
.scroll_soft_deprecate <- function(what, instead, env) {
  # internal = called from the namespace, or from the attached package env (where
  # testthat evaluates the package's own tests)
  if (environmentName(topenv(env)) %in% c("scroll", "package:scroll"))
    return(invisible(FALSE))
  if (isTRUE(.scroll_dep_seen[[what]])) return(invisible(FALSE))
  assign(what, TRUE, envir = .scroll_dep_seen)
  warning("`", what, "()` is deprecated as a public function and will no longer be ",
          "exported from scroll 0.3.0. ", instead, call. = FALSE)
  invisible(TRUE)
}
