# Small, pure helpers behind the Shiny app's interactive controls. These hold
# the actual decisions (what position a preset or a plot click implies);
# app.R only wires them to updateSliderInput() calls. Kept here, rather than
# inline in inst/shiny-app/app.R, so they can be unit tested directly --
# Shiny's testServer() cannot verify update*Input() calls at all (its mock
# session's sendInputMessage() is a documented no-op), so any logic left
# inside an observeEvent() that only calls updateSliderInput() is untestable.

#' Clamp a preset's I200/Iam positions to a data range
#'
#' Used by the Shiny app to decide what to set the position sliders to when a
#' preset is chosen (or the file/sample changes): the preset's published
#' positions, clamped to the 2-theta range the uploaded file actually covers
#' so the sliders are never set outside the data.
#'
#' @param range Numeric vector `c(min, max)`, the data's 2-theta range.
#' @param preset A preset name from [segal_presets()], or an invalid/`NULL`
#'   name, in which case `"cellulose_I_segal"` is used.
#' @return A list with `i200` and `iam`, each clamped into `range`.
#' @keywords internal
.clamp_preset_positions <- function(range, preset) {
  if (!is.numeric(range) || length(range) != 2L || anyNA(range) || range[1] >= range[2]) {
    stop("`range` must be two increasing, non-missing numbers.", call. = FALSE)
  }
  p <- tryCatch(segal_presets(preset), error = function(e) segal_presets("cellulose_I_segal"))
  clamp <- function(x) min(max(x, range[1]), range[2])
  list(i200 = clamp(p$i200), iam = clamp(p$iam))
}

#' Which slider does a plot click move?
#'
#' @param click_target The app's `input$click_target` value.
#' @return `"iam"` if `click_target` is `"iam"`, otherwise `"i200"` (this is
#'   also the safe default if `click_target` is missing or unrecognized).
#' @keywords internal
.click_slider_id <- function(click_target) {
  if (isTRUE(click_target == "iam")) "iam" else "i200"
}
