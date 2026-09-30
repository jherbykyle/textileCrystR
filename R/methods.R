# S3 methods for textile_cryst, summary.textile_cryst, textile_cryst_sensitivity
# and textile_cryst_input objects.

.fnum <- function(x, digits = 2) formatC(x, format = "f", digits = digits, big.mark = ",")

#' @export
print.textile_cryst_input <- function(x, ...) {
  cat("<textile_cryst_input> validated XRD pattern\n")
  cat(sprintf("  Points   : %d used of %d supplied\n", x$n_used, x$n_input))
  cat(sprintf("  2-theta  : %s to %s deg (median step %s deg)\n",
              format(x$range[1]), format(x$range[2]), format(signif(x$median_step, 4))))
  if (length(x$messages)) cat(paste0("  Note     : ", x$messages, "\n"), sep = "")
  invisible(x)
}

#' Methods for `textile_cryst` objects
#'
#' * `print()` shows the Segal Crystallinity Index, the positions and
#'   intensities used, and any warnings.
#' * `summary()` adds the method, interpolation, preset, input checks and the
#'   interpretation caveat.
#' * `as.data.frame()` returns a one-row data frame suitable for binding
#'   several samples together.
#'
#' @param x,object A `textile_cryst` object from [segal_ci()].
#' @param row.names,optional Passed on to the data-frame constructor; see
#'   [as.data.frame()].
#' @param ... Unused.
#' @return `print()` returns `x` invisibly; `summary()` returns an object of
#'   class `summary.textile_cryst`; `as.data.frame()` a data frame with columns
#'   `sample`, `ci`, `i200_position`, `i200_intensity`, `iam_position`,
#'   `iam_intensity`, `interpolation`, `preset`, `n_used` and `n_warnings`.
#' @name textile_cryst-methods
#' @examples
#' tt <- seq(5, 40, by = 0.05)
#' y <- 300 + 4000 * exp(-((tt - 22.7) / 1.1)^2) + 1500 * exp(-((tt - 16) / 1.5)^2)
#' res <- segal_ci(tt, y, i200 = 22.6, iam = 18)
#' print(res)
#' summary(res)
#' as.data.frame(res)
NULL

#' @rdname textile_cryst-methods
#' @export
print.textile_cryst <- function(x, ...) {
  p <- x$points
  cat("<textile_cryst> Segal Crystallinity Index\n")
  if (!is.na(x$sample)) cat(sprintf("  Sample : %s\n", x$sample))
  cat(sprintf("  CI     : %s %%\n", .fnum(x$ci)))
  cat(sprintf("  I200   : %s at %s deg 2-theta\n", .fnum(p$intensity[1]), .fnum(p$position_used[1])))
  cat(sprintf("  Iam    : %s at %s deg 2-theta\n", .fnum(p$intensity[2]), .fnum(p$position_used[2])))
  cat(sprintf("  Interpolation: %s%s\n", x$method$interpolation,
              if (all(p$on_measured_point)) " (both positions on measured points)" else ""))
  if (length(x$warnings)) cat(paste0("  Warning: ", x$warnings, "\n"), sep = "")
  invisible(x)
}

#' @rdname textile_cryst-methods
#' @export
summary.textile_cryst <- function(object, ...) {
  structure(list(
    sample = object$sample, ci = object$ci, points = object$points,
    method = object$method, input = object$input, warnings = object$warnings,
    interpretation = paste0(
      "The Segal Crystallinity Index is an empirical, relative XRD index. It is meaningful for ",
      "comparing samples measured, prepared and processed the same way, not as an absolute ",
      "crystalline fraction. Report the positions used with the value.")
  ), class = "summary.textile_cryst")
}

#' @export
print.summary.textile_cryst <- function(x, ...) {
  cat("Segal Crystallinity Index (summary)\n")
  cat(strrep("-", 40), "\n", sep = "")
  if (!is.na(x$sample)) cat(sprintf("Sample        : %s\n", x$sample))
  cat(sprintf("Equation      : %s\n", x$method$equation))
  cat(sprintf("CI            : %s %%\n", .fnum(x$ci, 3)))
  cat("Points used:\n")
  p <- x$points
  for (i in seq_len(nrow(p))) {
    cat(sprintf("  %-4s requested %s, used %s deg; intensity %s (measured neighbours %s, %s)\n",
                p$point[i], .fnum(p$position_requested[i]), .fnum(p$position_used[i]),
                .fnum(p$intensity[i]), .fnum(p$lower_2theta[i], 3), .fnum(p$upper_2theta[i], 3)))
  }
  cat(sprintf("Interpolation : %s\n", x$method$interpolation))
  cat(sprintf("Preset        : %s\n", if (is.null(x$method$preset)) "none (explicit positions)"
                                      else x$method$preset))
  cat(sprintf("Input         : %d points used of %d (%s to %s deg)\n", x$input$n_used,
              x$input$n_input, format(x$input$range[1]), format(x$input$range[2])))
  if (length(x$input$messages)) cat(paste0("Input note    : ", x$input$messages, "\n"), sep = "")
  if (length(x$warnings)) cat(paste0("Warning       : ", x$warnings, "\n"), sep = "")
  cat("\n")
  cat(strwrap(x$interpretation, width = 78), sep = "\n")
  invisible(x)
}

#' @rdname textile_cryst-methods
#' @export
as.data.frame.textile_cryst <- function(x, row.names = NULL, optional = FALSE, ...) {
  p <- x$points
  data.frame(
    sample = x$sample, ci = x$ci,
    i200_position = p$position_used[1], i200_intensity = p$intensity[1],
    iam_position = p$position_used[2], iam_intensity = p$intensity[2],
    interpolation = x$method$interpolation,
    preset = if (is.null(x$method$preset)) NA_character_ else x$method$preset,
    n_used = x$input$n_used, n_warnings = length(x$warnings),
    stringsAsFactors = FALSE, row.names = row.names
  )
}

#' Methods for `textile_cryst_sensitivity` objects
#'
#' `print()` shows the sensitivity matrix and its summary; `as.data.frame()`
#' returns the long-form grid.
#'
#' @param x A `textile_cryst_sensitivity` object from [segal_sensitivity()].
#' @param row.names,optional See [as.data.frame()].
#' @param ... Unused.
#' @return `print()` returns `x` invisibly; `as.data.frame()` a data frame with
#'   one row per I200 / Iam combination.
#' @name textile_cryst_sensitivity-methods
#' @examples
#' tt <- seq(5, 40, by = 0.05)
#' y <- 300 + 4000 * exp(-((tt - 22.7) / 1.1)^2) + 1500 * exp(-((tt - 16) / 1.5)^2)
#' sens <- segal_sensitivity(tt, y)
#' print(sens)
#' head(as.data.frame(sens))
NULL

#' @rdname textile_cryst_sensitivity-methods
#' @export
print.textile_cryst_sensitivity <- function(x, ...) {
  cat("<textile_cryst_sensitivity> Segal Crystallinity Index over an I200 x Iam grid\n")
  if (!is.na(x$sample)) cat(sprintf("  Sample : %s\n", x$sample))
  cat("  CI (%), rows = Iam position, columns = I200 position (deg 2-theta):\n")
  m <- round(x$matrix, 2)
  print(m)
  s <- x$summary
  cat(sprintf("  Min %s, max %s, range %s percentage points\n", .fnum(s$ci_min), .fnum(s$ci_max),
              .fnum(s$ci_range)))
  cat(sprintf("  Central estimate %s %% (I200 %s, Iam %s deg)\n", .fnum(x$central$ci),
              .fnum(x$central$i200), .fnum(x$central$iam)))
  cat("  ", x$note, "\n", sep = "")
  if (length(x$warnings)) cat(paste0("  Warning: ", x$warnings, "\n"), sep = "")
  invisible(x)
}

#' @rdname textile_cryst_sensitivity-methods
#' @export
as.data.frame.textile_cryst_sensitivity <- function(x, row.names = NULL, optional = FALSE, ...) {
  out <- as.data.frame(x$grid, stringsAsFactors = FALSE)
  out <- cbind(sample = x$sample, out, stringsAsFactors = FALSE)
  if (!is.null(row.names)) rownames(out) <- row.names
  out
}
