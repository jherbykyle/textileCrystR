# S3 methods for the peak-area (deconvolution) result classes, and the
# plot_peak_fit() diagnostic plot.

.peak_colours <- function(n) {
  base <- c("#A5382C", "#1E7A5C", "#3B5BA5", "#B5722E", "#6E4C9E", "#2E8B8B")
  if (n <= length(base)) return(base[seq_len(n)])
  grDevices::colorRampPalette(base)(n)
}

#' Methods for `textile_amorphous_hump` objects
#'
#' @param x A `textile_amorphous_hump` from [fit_amorphous_hump()].
#' @param ... Unused.
#' @return `print()` and `plot()` return `x` invisibly (`plot()` after
#'   printing a `ggplot`).
#' @name textile_amorphous_hump-methods
#' @examples
#' tt <- seq(10, 40, by = 0.05)
#' y <- 200 + 600 * exp(-((tt - 20.5) / 3)^2)
#' hump <- fit_amorphous_hump(tt, y, window = c(15, 28))
#' hump
#' plot(hump)
NULL

#' @rdname textile_amorphous_hump-methods
#' @export
print.textile_amorphous_hump <- function(x, ...) {
  cat("<textile_amorphous_hump> fitted amorphous contribution\n")
  if (!is.na(x$sample)) cat(sprintf("  Sample   : %s\n", x$sample))
  cat(sprintf("  Model    : %s\n", x$model))
  cat(sprintf("  Center   : %s deg 2-theta\n", .fnum(x$center)))
  cat(sprintf("  Height   : %s\n", .fnum(x$height)))
  cat(sprintf("  FWHM     : %s deg\n", .fnum(x$fwhm_reported)))
  if (x$background$model != "none") {
    cat(sprintf("  Background: %s (intercept %s%s)\n", x$background$model, .fnum(x$background$intercept),
               if (x$background$model == "linear") sprintf(", slope %s", .fnum(x$background$slope, 4)) else ""))
  }
  cat(sprintf("  Area     : %s\n", .fnum(x$area)))
  cat(sprintf("  R-squared: %s   RMSE: %s\n", .fnum(x$r_squared, 4), .fnum(x$rmse)))
  if (!x$converged) cat("  Warning  : optimizer did not clearly converge\n")
  invisible(x)
}

#' @rdname textile_amorphous_hump-methods
#' @export
plot.textile_amorphous_hump <- function(x, ...) {
  d <- x$data
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$two_theta)) +
    ggplot2::geom_line(ggplot2::aes(y = .data$intensity), colour = "#243070", linewidth = 0.5) +
    ggplot2::geom_line(ggplot2::aes(y = .data$fitted), colour = "#D89B16", linewidth = 0.8,
                       linetype = "dashed") +
    ggplot2::labs(x = "2-theta (\u00b0)", y = "Intensity",
                  title = if (is.na(x$sample)) "Amorphous hump fit" else x$sample,
                  subtitle = sprintf("%s hump, R\u00b2 = %.3f", x$model, x$r_squared)) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.background = ggplot2::element_rect(fill = "white", colour = NA))
  print(p)
  invisible(x)
}

#' Methods for `textile_cryst_peakarea` objects
#'
#' * `print()` shows the crystallinity, fit quality and convergence.
#' * `summary()` adds the full peak table, metadata and any warnings.
#' * `as.data.frame()` returns a one-row overview (for binding several
#'   samples together); see [peak_fit_summary()] for the multi-row,
#'   per-peak table.
#' * `plot()` calls [plot_peak_fit()].
#'
#' @param x,object A `textile_cryst_peakarea` from [fit_crystalline_peaks()].
#' @param row.names,optional See [as.data.frame()].
#' @param ... For `plot()`, passed to [plot_peak_fit()]; otherwise unused.
#' @return `print()`/`plot()` return `x` invisibly; `summary()` returns a
#'   `summary.textile_cryst_peakarea` object; `as.data.frame()` a one-row
#'   data frame.
#' @name textile_cryst_peakarea-methods
#' @examples
#' tt <- seq(8, 40, by = 0.05)
#' y <- 200 + 500 * exp(-((tt - 20.5) / 3)^2) + 1800 * exp(-((tt - 22.6) / 0.7)^2)
#' fit <- fit_crystalline_peaks(tt, y, centers = 22.6, sample = "syn")
#' fit
#' summary(fit)
#' as.data.frame(fit)
NULL

#' @rdname textile_cryst_peakarea-methods
#' @export
print.textile_cryst_peakarea <- function(x, ...) {
  cat("<textile_cryst_peakarea> Peak-area (deconvolution) Crystallinity\n")
  if (!is.na(x$sample)) cat(sprintf("  Sample     : %s\n", x$sample))
  cat(sprintf("  CI         : %s %%\n", .fnum(x$ci)))
  cat(sprintf("  Model      : %s (%d crystalline peak(s) + amorphous hump)\n", x$model, nrow(x$peaks)))
  if (x$background$model != "none") {
    cat(sprintf("  Background : %s (intercept %s%s)\n", x$background$model, .fnum(x$background$intercept),
               if (x$background$model == "linear") sprintf(", slope %s", .fnum(x$background$slope, 4)) else ""))
  }
  cat(sprintf("  Fit quality: R-squared = %s, RMSE = %s\n", .fnum(x$r_squared, 4), .fnum(x$rmse)))
  if (!x$converged) cat(sprintf("  Warning    : %s\n", x$warnings[1]))
  invisible(x)
}

#' @rdname textile_cryst_peakarea-methods
#' @export
summary.textile_cryst_peakarea <- function(object, ...) {
  structure(list(
    sample = object$sample, ci = object$ci, peaks = object$peaks, amorphous = object$amorphous,
    background = object$background,
    r_squared = object$r_squared, rmse = object$rmse, model = object$model, window = object$window,
    metadata = object$metadata, warnings = object$warnings,
    interpretation = paste0(
      "Peak-area crystallinity depends on the chosen peak-shape model, the number and starting ",
      "positions of crystalline peaks, and the amorphous window -- it is not an absolute crystalline ",
      "fraction. Compare fits made the same way; report the model and positions used alongside the value.")
  ), class = "summary.textile_cryst_peakarea")
}

#' @export
print.summary.textile_cryst_peakarea <- function(x, ...) {
  cat("Peak-area (deconvolution) Crystallinity (summary)\n")
  cat(strrep("-", 45), "\n", sep = "")
  if (!is.na(x$sample)) cat(sprintf("Sample        : %s\n", x$sample))
  cat(sprintf("CI            : %s %%\n", .fnum(x$ci, 3)))
  cat(sprintf("Model         : %s\n", x$model))
  if (x$background$model != "none") {
    cat(sprintf("Background    : %s (intercept %s%s)\n", x$background$model, .fnum(x$background$intercept),
               if (x$background$model == "linear") sprintf(", slope %s per deg", .fnum(x$background$slope, 4))
               else ""))
  }
  cat(sprintf("Window        : %s to %s deg\n", format(x$window[1]), format(x$window[2])))
  cat(sprintf("Fit quality   : R-squared = %s, RMSE = %s\n", .fnum(x$r_squared, 4), .fnum(x$rmse)))
  cat("Peaks:\n")
  tab <- rbind(x$peaks[, c("label", "center", "height", "fwhm", "area", "pct_of_crystalline")],
              x$amorphous[, c("label", "center", "height", "fwhm", "area", "pct_of_crystalline")])
  for (i in seq_len(nrow(tab))) {
    cat(sprintf("  %-10s center %s  height %s  fwhm %s  area %s%s\n",
                tab$label[i], .fnum(tab$center[i]), .fnum(tab$height[i]), .fnum(tab$fwhm[i]),
                .fnum(tab$area[i]),
                if (is.na(tab$pct_of_crystalline[i])) "" else sprintf("  (%.1f%% of crystalline)",
                                                                      tab$pct_of_crystalline[i])))
  }
  cat(sprintf("Reference     : %s\n", x$metadata$reference))
  if (length(x$warnings)) cat(paste0("Warning       : ", x$warnings, "\n"), sep = "")
  cat("\n")
  cat(strwrap(x$interpretation, width = 78), sep = "\n")
  cat("\n")
  invisible(x)
}

#' @rdname textile_cryst_peakarea-methods
#' @export
as.data.frame.textile_cryst_peakarea <- function(x, row.names = NULL, optional = FALSE, ...) {
  data.frame(
    sample = x$sample, ci = x$ci, crystalline_area = x$crystalline_area,
    amorphous_area = x$amorphous_area, total_area = x$total_area, n_peaks = nrow(x$peaks),
    model = x$model, r_squared = x$r_squared, rmse = x$rmse, converged = x$converged,
    n_warnings = length(x$warnings), stringsAsFactors = FALSE, row.names = row.names
  )
}

#' @rdname textile_cryst_peakarea-methods
#' @export
plot.textile_cryst_peakarea <- function(x, ...) {
  plot_peak_fit(x, ...)
}

#' Plot a peak-area (deconvolution) fit
#'
#' Draws the raw diffraction pattern, each fitted crystalline peak, the
#' amorphous hump, the total fitted curve, and a residuals panel beneath.
#'
#' @param fit A `textile_cryst_peakarea` from [fit_crystalline_peaks()].
#' @param show_components Logical; draw the individual crystalline peaks and
#'   the amorphous hump as separate curves (default `TRUE`)?
#' @param show_residuals Logical; include the residuals panel (default
#'   `TRUE`)?
#' @param ... Unused.
#' @return A `patchwork`/`ggplot` object (when `show_residuals = TRUE`, a
#'   `patchwork` composite of two `ggplot` panels; each panel, and the
#'   composite itself, can still be modified with the usual `+` syntax).
#' @seealso [fit_crystalline_peaks()]
#' @examples
#' tt <- seq(8, 40, by = 0.05)
#' y <- 200 + 500 * exp(-((tt - 20.5) / 3)^2) + 1800 * exp(-((tt - 22.6) / 0.7)^2)
#' fit <- fit_crystalline_peaks(tt, y, centers = 22.6, sample = "syn")
#' plot_peak_fit(fit)
#' @export
plot_peak_fit <- function(fit, show_components = TRUE, show_residuals = TRUE, ...) {
  if (!inherits(fit, "textile_cryst_peakarea")) {
    stop("`fit` must be a textile_cryst_peakarea object from fit_crystalline_peaks().", call. = FALSE)
  }
  d <- fit$data
  cols <- .peak_colours(nrow(fit$peaks))

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$two_theta)) +
    ggplot2::geom_line(ggplot2::aes(y = .data$intensity), colour = "#5B6275", linewidth = 0.45)

  if (fit$background$model != "none") {
    p <- p + ggplot2::geom_line(ggplot2::aes(y = .data$background), colour = "#9AA0AE",
                                linewidth = 0.6, linetype = "dotdash")
  }

  if (isTRUE(show_components)) {
    nm <- .model_param_names(fit$model)
    comp_curves <- do.call(rbind, lapply(seq_len(nrow(fit$peaks)), function(i) {
      p_i <- stats::setNames(as.numeric(fit$peaks[i, nm]), nm)
      tibble::tibble(two_theta = d$two_theta, intensity = .eval_component(d$two_theta, p_i, fit$model),
                    label = fit$peaks$label[i])
    }))
    am_p <- stats::setNames(as.numeric(fit$amorphous[1, .model_param_names(fit$model)]),
                            .model_param_names(fit$model))
    am_curve <- tibble::tibble(two_theta = d$two_theta,
                               intensity = .eval_component(d$two_theta, am_p, fit$model),
                               label = "Amorphous")
    p <- p +
      ggplot2::geom_line(data = am_curve, ggplot2::aes(y = .data$intensity), colour = "#B0B4C2",
                         linewidth = 0.9, linetype = "dotted") +
      ggplot2::geom_line(data = comp_curves,
                         ggplot2::aes(y = .data$intensity, colour = .data$label), linewidth = 0.6,
                         linetype = "dashed") +
      ggplot2::scale_colour_manual(values = stats::setNames(cols, fit$peaks$label), name = "Peak")
  }

  p <- p +
    ggplot2::geom_line(ggplot2::aes(y = .data$fitted), colour = "#1E2757", linewidth = 0.8) +
    ggplot2::labs(x = NULL, y = "Intensity",
                  title = if (is.na(fit$sample)) "Peak-area fit" else fit$sample,
                  subtitle = sprintf("CI = %.2f%%  (%s, R\u00b2 = %.4f)", fit$ci, fit$model, fit$r_squared)) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "right", panel.grid.minor = ggplot2::element_blank(),
                   plot.background = ggplot2::element_rect(fill = "white", colour = NA))

  if (!isTRUE(show_residuals)) return(p)

  rp <- ggplot2::ggplot(fit$residuals, ggplot2::aes(x = .data$two_theta, y = .data$residual)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#B0B4C2", linewidth = 0.4) +
    ggplot2::geom_line(colour = "#A5382C", linewidth = 0.4) +
    ggplot2::labs(x = "2-theta (\u00b0)", y = "Residual") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.background = ggplot2::element_rect(fill = "white", colour = NA))

  patchwork::wrap_plots(p, rp, ncol = 1, heights = c(3, 1))
}
