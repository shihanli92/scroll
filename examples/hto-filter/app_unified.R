library(scroll)

# Cao et al. unified explorer, with an OPTIONAL password gate.
#
# Auth is enabled only when the SCROLL_PW environment variable is set (and
# shinymanager is installed); the password is read from the environment, so no
# secret lives in this file or the deploy bundle. Username defaults to "viewer"
# (override with SCROLL_USER). Leave SCROLL_PW unset for an open app (dev bypass).
#
#   Sys.setenv(SCROLL_USER = "viewer", SCROLL_PW = "…")   # then run the app
#
# The gate wraps only scroll's PUBLIC app object (its UI + server), so it works
# with any installed scroll version without reaching into internals.

app <- scroll_app(".")

# Load credentials from a local .Renviron beside this app.R if present (runApp()
# sets the working dir to the app folder, so R's startup .Renviron loading may miss
# it). Keeps SCROLL_USER / SCROLL_PW out of committed source — the file ships with
# the deploy bundle, not the repo.
if (file.exists(".Renviron")) readRenviron(".Renviron")

.pw <- Sys.getenv("SCROLL_PW")
if (nzchar(.pw) && requireNamespace("shinymanager", quietly = TRUE)) {
  .ui   <- environment(app$httpHandler)$ui
  .srv  <- app$serverFuncSource()
  .cred <- data.frame(user = Sys.getenv("SCROLL_USER", unset = "viewer"),
                      password = .pw, stringsAsFactors = FALSE)
  app <- shiny::shinyApp(
    ui = shinymanager::secure_app(.ui),
    server = function(input, output, session) {
      shinymanager::secure_server(
        check_credentials = shinymanager::check_credentials(.cred))
      .srv(input, output, session)
    })
}

app
