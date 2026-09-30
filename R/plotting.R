#' Plot a Segal calculation
#'
#' Draws the XRD pattern with the I200 and Iam points, their labels and the
#' calculated Segal Crystallinity Index. The function returns an ordinary
#' [ggplot2::ggplot()] object, so it can be modified with the usual `+`
#' syntax (themes, scales, additional layers, `coord_cartesian()`, ...).
#'
#' @param x A `textile_cryst` object from [segal_ci()].
#' @param xlim Optional numeric vector of length two: the 2-theta range to
#'   display. The default shows the region from 10 degrees below the Iam
#'   position to 15 degrees above the I200 position, clipped to the measured
#'   range. The full pattern is always kept in the plot data.
#' @param show_ci Logical. Show the calculated index in the subtitle?
#' @param title Optional plot title; defaults to the sample name.
#' @param ... Unused; for compatibility with the [plot()] generic.
#' @return A `ggplot` object.
#' @examples
#' tt <- seq(5, 40, by = 0.05)
#' y <- 300 + 4000 * exp(-((tt - 22.7) / 1.1)^2) + 1500 * exp(-((tt - 16) / 1.5)^2)
#' res <- segal_ci(tt, y, i200 = 22.6, iam = 18, sample = "synthetic example")
#' p <- plot(res)
#' p + ggplot2::theme_bw()
#' @export
plot.textile_cryst <- function(x, xlim = NULL, show_ci = TRUE, title = NULL, ...) {
  d <- x$data
  pts <- x$points
  if (is.null(xlim)) {
    xlim <- c(max(min(d$two_theta), min(pts$position_used) - 10),
              min(max(d$two_theta), max(pts$position_used) + 15))
  }
  if (!is.numeric(xlim) || length(xlim) != 2L || anyNA(xlim) || xlim[1] >= xlim[2]) {
    stop("`xlim` must be two increasing numbers.", call. = FALSE)
  }
  vis <- d$intensity[d$two_theta >= xlim[1] & d$two_theta <= xlim[2]]
  if (!length(vis)) vis <- d$intensity
  ylim <- c(min(0, min(vis)), max(vis) * 1.08)
  span <- diff(ylim)

  pts$label <- sprintf("%s = %s at %.2f\u00b0", pts$point,
                       formatC(pts$intensity, format = "f", digits = 0, big.mark = ","),
                       pts$position_used)
  pts$label_x <- pts$position_used + c(0.6, 0)
  pts$label_y <- pts$intensity + c(0, -0.09 * span)
  pts$label_hjust <- c(0, 0.5)

  sub <- if (isTRUE(show_ci)) {
    sprintf("Segal Crystallinity Index = %.2f%%  (%s interpolation)", x$ci, x$method$interpolation)
  } else NULL
  ttl <- title %||% (if (is.na(x$sample)) "XRD pattern" else x$sample)

  ggplot2::ggplot(d, ggplot2::aes(x = .data$two_theta, y = .data$intensity)) +
    ggplot2::geom_line(colour = "#243070", linewidth = 0.5) +
    ggplot2::geom_segment(
      data = pts,
      ggplot2::aes(x = .data$position_used, xend = .data$position_used,
                   y = 0, yend = .data$intensity, colour = .data$point),
      linetype = "dashed", linewidth = 0.4, inherit.aes = FALSE) +
    ggplot2::geom_point(
      data = pts,
      ggplot2::aes(x = .data$position_used, y = .data$intensity, colour = .data$point),
      size = 3, inherit.aes = FALSE) +
    do.call(ggplot2::geom_label, c(
      list(data = pts,
           mapping = ggplot2::aes(x = .data$label_x, y = .data$label_y, label = .data$label,
                                  colour = .data$point, hjust = .data$label_hjust),
           size = 3.3, fill = "white", inherit.aes = FALSE, show.legend = FALSE),
      .label_border_zero()
    )) +
    ggplot2::scale_colour_manual(values = c(I200 = "#A5382C", Iam = "#D89B16"), guide = "none") +
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
    ggplot2::labs(x = "2-theta (\u00b0)", y = "Intensity", title = ttl, subtitle = sub) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "#5B6275", linewidth = 0.3),
      plot.title = ggplot2::element_text(face = "bold", colour = "#1E2757"),
      plot.subtitle = ggplot2::element_text(colour = "#5B6275"),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# geom_label()'s border-width parameter was renamed from `label.size` to
# `linewidth` in ggplot2 3.5.0 (label.size is soft-deprecated since then).
# Passing the wrong name for the installed version either triggers a
# deprecation warning (old name on new ggplot2) or is silently ignored,
# leaving the default nonzero border visible (new name on old ggplot2) --
# checked directly for both cases before relying on this, since neither
# failure mode raises an error that would otherwise catch it. Resolved once
# per plot call, based on the installed version, rather than hardcoding
# either name.
.label_border_zero <- function() {
  if (utils::packageVersion("ggplot2") >= "3.5.0") list(linewidth = 0) else list(label.size = 0)
}
