# The "one-shot" entry point: reads a file, computes the Segal index for
# every sample it contains, shows the plot, and optionally writes a plot
# file and/or an HTML report. It is a thin orchestration layer over
# import_xrd_file() + segal_ci(); it introduces no new science.

.batch_plot <- function(results, summary) {
  zoom_one <- function(r) {
    d <- r$data
    pos <- r$points$position_used
    xlim <- c(max(min(d$two_theta), min(pos) - 10), min(max(d$two_theta), max(pos) + 15))
    d <- d[d$two_theta >= xlim[1] & d$two_theta <= xlim[2], ]
    d$sample <- r$sample
    d
  }
  data <- do.call(rbind, lapply(results, zoom_one))
  pts <- do.call(rbind, lapply(results, function(r) {
    p <- r$points
    p$sample <- r$sample
    p
  }))
  ggplot2::ggplot(data, ggplot2::aes(x = .data$two_theta, y = .data$intensity)) +
    ggplot2::geom_line(colour = "#243070", linewidth = 0.4) +
    ggplot2::geom_point(
      data = pts,
      ggplot2::aes(x = .data$position_used, y = .data$intensity, colour = .data$point),
      size = 2, inherit.aes = FALSE) +
    ggplot2::scale_colour_manual(values = c(I200 = "#A5382C", Iam = "#D89B16"), name = NULL) +
    ggplot2::facet_wrap(~sample, scales = "free") +
    ggplot2::labs(x = "2-theta (\u00b0)", y = "Intensity",
                  title = "Segal Crystallinity Index by sample",
                  subtitle = "Composition/sample labels are metadata, not a crystallinity scale") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "top",
                   panel.grid.minor = ggplot2::element_blank(),
                   plot.background = ggplot2::element_rect(fill = "white", colour = NA))
}

#' Methods for `textile_cryst_batch` objects
#'
#' A `textile_cryst_batch`, returned by [analyze_textile_file()] when a file
#' contains more than one sample, bundles one [segal_ci()] result per sample.
#' `print()` shows a compact summary table; `plot()` draws every pattern in a
#' facet grid with its I200/Iam points marked; `as.data.frame()` returns the
#' tidy one-row-per-sample table.
#'
#' @param x A `textile_cryst_batch` from [analyze_textile_file()].
#' @param row.names,optional See [as.data.frame()].
#' @param ... Unused.
#' @return `print()` and `plot()` return `x` invisibly; `as.data.frame()`
#'   returns a data frame with one row per sample.
#' @name textile_cryst_batch-methods
NULL

#' @rdname textile_cryst_batch-methods
#' @export
print.textile_cryst_batch <- function(x, ...) {
  cat(sprintf("<textile_cryst_batch> %d sample(s), Segal Crystallinity Index\n", length(x$results)))
  print(as.data.frame(x)[, c("sample", "ci", "i200_position", "iam_position", "n_warnings")],
        row.names = FALSE)
  nw <- sum(vapply(x$results, function(r) length(r$warnings), integer(1)))
  if (nw > 0) {
    cat(sprintf("  %d sample(s) raised a warning; inspect individual results (x$results) for detail.\n",
                sum(vapply(x$results, function(r) length(r$warnings) > 0, logical(1)))))
  }
  invisible(x)
}

#' @rdname textile_cryst_batch-methods
#' @export
plot.textile_cryst_batch <- function(x, ...) {
  .batch_plot(x$results, x$summary)
}

#' @rdname textile_cryst_batch-methods
#' @export
as.data.frame.textile_cryst_batch <- function(x, row.names = NULL, optional = FALSE, ...) {
  out <- x$summary
  if (!is.null(row.names)) rownames(out) <- row.names
  out
}

#' Analyze an XRD file in one call
#'
#' Reads an XRD file, validates it, calculates the Segal Crystallinity Index
#' for every sample it contains, displays the plot, and optionally saves a
#' plot image and/or an HTML summary report -- all in a single command. This
#' is the quickest way to use the package; [import_xrd_file()], [segal_ci()]
#' and [segal_sensitivity()] remain available for anything more specific.
#'
#' @param path Path to a CSV/TSV/whitespace-delimited XRD file, or a file in
#'   the two-header-row multi-sample layout (see [import_xrd_file()]).
#' @param sample Optional: analyze only this sample name from a multi-sample
#'   file. If `NULL` (default) and the file has more than one sample, every
#'   sample is analyzed and a `textile_cryst_batch` is returned.
#' @param preset,i200,iam,interpolation Passed to [segal_ci()]. `preset`
#'   defaults to `"cellulose_I_segal"` when neither `i200` nor `iam` is given.
#' @param cols Passed to [import_xrd_file()] when the file needs its columns
#'   named explicitly.
#' @param plot Logical. Display the plot as a side effect? Default `TRUE`.
#' @param plot_file Optional path (e.g. ending in `.png` or `.pdf`) to save
#'   the plot to, via [ggplot2::ggsave()].
#' @param report_file Optional path ending in `.html` to write a self-contained
#'   summary report to, via [rmarkdown::render()] (requires the rmarkdown and
#'   knitr packages).
#' @param na.rm,duplicates Passed to [validate_crystallinity_input()] via
#'   [segal_ci()].
#' @return A `textile_cryst` (single sample) or `textile_cryst_batch`
#'   (multiple samples), returned invisibly if `report_file` or `plot_file`
#'   was written (so the console isn't cluttered when the point of the call
#'   was to produce a file), and visibly otherwise.
#' @examples
#' path <- system.file("extdata", "cotton.csv", package = "textileCrystR")
#' res <- analyze_textile_file(path, sample = "0%/1", plot = FALSE)
#' res
#'
#' \donttest{
#' batch <- analyze_textile_file(path, plot = FALSE)
#' batch
#' }
#' @export
analyze_textile_file <- function(path, sample = NULL, preset = NULL, i200 = NULL, iam = NULL,
                                 interpolation = c("linear", "nearest"), cols = NULL,
                                 plot = TRUE, plot_file = NULL, report_file = NULL,
                                 na.rm = FALSE, duplicates = c("error", "mean")) {
  if (is.null(preset) && is.null(i200) && is.null(iam)) preset <- "cellulose_I_segal"

  xrd <- tryCatch(
    import_xrd_file(path, cols = cols, sample = sample),
    error = function(e) {
      stop(sprintf(paste0(
        "analyze_textile_file() could not read '%s':\n%s"),
        path, conditionMessage(e)), call. = FALSE)
    }
  )

  samples <- unique(xrd$sample)
  if (!is.null(sample)) {
    if (!sample %in% samples) {
      stop(sprintf(paste0(
        "Sample '%s' was not found in '%s'. Samples in this file: %s."),
        sample, path, paste(samples, collapse = ", ")), call. = FALSE)
    }
    samples <- sample
  }

  run_one <- function(s) {
    d <- xrd[xrd$sample == s, ]
    tryCatch(
      segal_ci(d$two_theta, d$intensity, i200 = i200, iam = iam, preset = preset,
              interpolation = interpolation, sample = s, na.rm = na.rm, duplicates = duplicates),
      error = function(e) {
        stop(sprintf(paste0(
          "Could not calculate the Segal index for sample '%s' in '%s':\n%s\nIf this position ",
          "genuinely falls outside what was measured, check the file was exported over a wide ",
          "enough 2-theta range (it should cover at least the Iam position to the I200 position, ",
          "roughly 15-25 degrees for cellulose I)."), s, path, conditionMessage(e)), call. = FALSE)
      }
    )
  }

  wrote_file <- !is.null(plot_file) || !is.null(report_file)

  if (length(samples) == 1L) {
    res <- run_one(samples)
    if (isTRUE(plot) || !is.null(plot_file)) {
      p <- plot(res)
      if (isTRUE(plot)) print(p)
      if (!is.null(plot_file)) ggplot2::ggsave(plot_file, p, width = 7, height = 4.2, dpi = 150)
    }
    if (!is.null(report_file)) .render_report(list(result = res), report_file)
    return(if (wrote_file) invisible(res) else res)
  }

  results <- stats::setNames(lapply(samples, run_one), samples)
  summ <- do.call(rbind, lapply(results, as.data.frame))
  rownames(summ) <- NULL
  batch <- structure(list(results = results, summary = summ), class = "textile_cryst_batch")
  if (isTRUE(plot) || !is.null(plot_file)) {
    p <- plot(batch)
    if (isTRUE(plot)) print(p)
    if (!is.null(plot_file)) {
      ggplot2::ggsave(plot_file, p, width = 8, height = 5.5, dpi = 150)
    }
  }
  if (!is.null(report_file)) .render_report(list(result = batch), report_file)
  if (wrote_file) invisible(batch) else batch
}

# Shared by analyze_textile_file(report_file = ) and the Shiny app's
# "Download summary report" button, so both produce the same report.
.render_report <- function(params, output_file,
                           template = system.file("report", "report_template.Rmd",
                                                  package = "textileCrystR")) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop(paste0("Writing a report needs the rmarkdown package. Install it with:\n",
                '  install.packages("rmarkdown")'), call. = FALSE)
  }
  if (!nzchar(template)) {
    stop("Could not find textileCrystR's report template. Try reinstalling the package.",
         call. = FALSE)
  }
  output_file <- normalizePath(output_file, mustWork = FALSE)
  rmarkdown::render(template, output_file = output_file, params = params,
                    envir = new.env(parent = globalenv()), quiet = TRUE)
  invisible(output_file)
}
