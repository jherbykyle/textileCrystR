#' Launch the textileCrystR interactive dashboard
#'
#' Opens a local Shiny web app for interactive Segal Crystallinity Index
#' analysis: upload an XRD file, pick a position preset, fine-tune the I200
#' and Iam positions with sliders or by clicking directly on the plot, and
#' download the plot or a summary report. This needs the `shiny` and
#' `rmarkdown` packages, which are not required for the rest of the package
#' and are only checked for when this function is actually called.
#'
#' @param launch Logical. If `TRUE` (default), starts the app with
#'   [shiny::runApp()]. If `FALSE`, returns the app object without starting
#'   it -- mainly useful for automated testing.
#' @param ... Further arguments passed to [shiny::runApp()] (for example
#'   `port` or `launch.browser`).
#' @param .app_dir Internal; the app directory to use. Not for end users --
#'   exposed only so tests can point at a directory that doesn't exist,
#'   without needing to reinstall the package to exercise that error path.
#' @return If `launch = TRUE`, called for its side effect (invisibly returns
#'   `NULL` after the app stops). If `launch = FALSE`, the `shiny.appobj`
#'   invisibly, without starting it.
#' @examples
#' \dontrun{
#' run_textile_app()
#' }
#' @export
run_textile_app <- function(launch = TRUE, ...,
                            .app_dir = system.file("shiny-app", package = "textileCrystR")) {
  needed <- c("shiny", "rmarkdown")
  have <- vapply(needed, requireNamespace, logical(1), quietly = TRUE)
  if (!all(have)) {
    stop(sprintf(paste0(
      "run_textile_app() needs the %s package(s), which %s not installed. Install %s with:\n",
      "  install.packages(c(%s))"),
      paste(needed[!have], collapse = " and "),
      if (sum(!have) > 1) "are" else "is",
      if (sum(!have) > 1) "them" else "it",
      paste(sprintf('"%s"', needed[!have]), collapse = ", ")), call. = FALSE)
  }
  if (!nzchar(.app_dir)) {
    stop("Could not find textileCrystR's Shiny app directory (inst/shiny-app). ",
         "Try reinstalling the package.", call. = FALSE)
  }
  app <- shiny::shinyAppDir(.app_dir)
  if (!isTRUE(launch)) return(invisible(app))
  shiny::runApp(app, ...)  # nocov: starts a real, blocking local server
}
