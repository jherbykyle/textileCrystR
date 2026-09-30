# Input validation and file reading. Nothing in this file changes measured
# intensities: patterns are only re-ordered (sorting), reduced (dropping rows
# the user has explicitly agreed to drop) or averaged (duplicated 2-theta, only
# when the user asks for it), and every such action is recorded.

.POS_TOL <- 1e-6   # degrees; a requested position this close to a measured point
                   # is treated as that measured point (guards against floating
                   # point / single-precision storage of 2-theta values)
.MAX_GAP <- 0.5    # degrees; wider bracketing intervals trigger a warning

#' Validate an XRD pattern before a crystallinity calculation
#'
#' Checks that a 2-theta / intensity pair is usable for Segal peak-height
#' calculations and returns the cleaned pattern together with a record of
#' everything that was checked or changed. Measured intensities are never
#' altered: rows can only be re-ordered, dropped (if `na.rm = TRUE`) or, for
#' exactly repeated 2-theta values, averaged (if `duplicates = "mean"`).
#'
#' @section Rules:
#' * `two_theta` and `intensity` must be numeric vectors of equal length with
#'   at least two points.
#' * Missing or non-finite values are an error unless `na.rm = TRUE`, in which
#'   case the affected rows are dropped and counted.
#' * All 2-theta values must lie strictly between 0 and 180 degrees.
#' * Unsorted 2-theta values are sorted (ascending) and this is recorded.
#' * Repeated 2-theta values are an error unless `duplicates = "mean"`.
#' * Negative intensities are allowed (they can occur after background
#'   subtraction) but are counted and reported.
#'
#' @param two_theta Numeric vector of diffraction angles in degrees 2-theta.
#' @param intensity Numeric vector of intensities, same length as `two_theta`.
#' @param na.rm Logical. Drop rows with missing or non-finite values instead of
#'   stopping with an error? Default `FALSE`.
#' @param duplicates `"error"` (default) or `"mean"`. What to do when the same
#'   2-theta value occurs more than once.
#'
#' @return An object of class `textile_cryst_input`: a list with the cleaned
#'   `two_theta` and `intensity`, counts (`n_input`, `n_used`, `n_dropped`,
#'   `n_duplicates`, `n_negative`), `was_sorted`, the 2-theta `range`, the
#'   `median_step` and `max_step` between neighbouring points, and `messages`
#'   describing every action taken.
#' @examples
#' v <- validate_crystallinity_input(c(12, 10, 11, 13), c(5, 3, 4, 6))
#' v$two_theta
#' v$messages
#' @export
validate_crystallinity_input <- function(two_theta, intensity, na.rm = FALSE,
                                         duplicates = c("error", "mean")) {
  duplicates <- match.arg(duplicates)
  if (!is.logical(na.rm) || length(na.rm) != 1L || is.na(na.rm)) {
    stop("`na.rm` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(two_theta) || !is.numeric(intensity)) {
    stop("`two_theta` and `intensity` must both be numeric vectors.", call. = FALSE)
  }
  if (length(two_theta) != length(intensity)) {
    stop(sprintf("`two_theta` (length %d) and `intensity` (length %d) must have the same length.",
                 length(two_theta), length(intensity)), call. = FALSE)
  }
  tt <- as.numeric(two_theta)
  y <- as.numeric(intensity)
  n_input <- length(tt)
  msgs <- character()

  bad <- !is.finite(tt) | !is.finite(y)
  n_dropped <- sum(bad)
  if (n_dropped > 0L) {
    if (!na.rm) {
      stop(sprintf(paste0("%d missing or non-finite value(s) found in `two_theta` or `intensity`. ",
                          "Remove them or use na.rm = TRUE."), n_dropped), call. = FALSE)
    }
    tt <- tt[!bad]
    y <- y[!bad]
    msgs <- c(msgs, sprintf("Dropped %d row(s) with missing or non-finite values.", n_dropped))
  }
  if (length(tt) < 2L) {
    stop("At least two valid points are needed.", call. = FALSE)
  }
  if (any(tt <= 0 | tt >= 180)) {
    stop("All 2-theta values must lie strictly between 0 and 180 degrees.", call. = FALSE)
  }

  was_sorted <- is.unsorted(tt)
  if (was_sorted) {
    o <- order(tt)
    tt <- tt[o]
    y <- y[o]
    msgs <- c(msgs, "2-theta values were not in ascending order; the pattern was sorted.")
  }

  dup <- duplicated(tt)
  n_dup <- sum(dup)
  if (n_dup > 0L) {
    if (duplicates == "error") {
      stop(sprintf(paste0("%d duplicated 2-theta value(s) found. ",
                          "Use duplicates = \"mean\" to average their intensities."), n_dup),
           call. = FALSE)
    }
    u <- unique(tt)
    y <- as.numeric(tapply(y, match(tt, u), mean))
    tt <- u
    msgs <- c(msgs, sprintf("Averaged intensities at %d duplicated 2-theta value(s).", n_dup))
  }
  if (length(tt) < 2L) {
    stop("At least two distinct 2-theta values are needed.", call. = FALSE)
  }

  n_neg <- sum(y < 0)
  if (n_neg > 0L) {
    msgs <- c(msgs, sprintf(paste0("%d point(s) have negative intensity, which can result from ",
                                   "background subtraction."), n_neg))
  }
  steps <- diff(tt)
  structure(list(
    two_theta = tt, intensity = y,
    n_input = n_input, n_used = length(tt), n_dropped = n_dropped,
    n_duplicates = n_dup, n_negative = n_neg, was_sorted = was_sorted,
    range = range(tt), median_step = stats::median(steps), max_step = max(steps),
    messages = msgs
  ), class = "textile_cryst_input")
}

# Validate a single requested position (in degrees 2-theta).
.check_position <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 || x >= 180) {
    stop(sprintf("`%s` must be a single number of degrees 2-theta between 0 and 180.", name),
         call. = FALSE)
  }
  invisible(x)
}

#' Read a paired-column XRD file
#'
#' Reads a CSV file in which each sample occupies two adjacent columns
#' (2-theta, intensity), the first row holds the sample name above the first
#' column of each pair and the second row holds column units, as in the
#' `cotton.csv` file shipped with the package. Values are read exactly as
#' stored; no processing of any kind is applied.
#'
#' @param path Path to the CSV file.
#' @return A [tibble::tibble()] with columns `sample`, `two_theta` and
#'   `intensity`, samples in file order. Blank trailing cells of shorter
#'   series are dropped.
#' @examples
#' path <- system.file("extdata", "cotton.csv", package = "textileCrystR")
#' xrd <- read_xrd_pairs(path)
#' unique(xrd$sample)
#' @export
read_xrd_pairs <- function(path) {
  if (!is.character(path) || length(path) != 1L || !file.exists(path)) {
    stop("`path` must be the path of an existing file.", call. = FALSE)
  }
  # Read as text and strip a UTF-8 byte-order mark by bytes, so the result does
  # not depend on the session locale (unit text such as the degree sign is in
  # the skipped second row and is never interpreted).
  lines <- readLines(path, warn = FALSE)
  if (length(lines) < 3L) {
    stop("The file needs two header rows and at least one data row.", call. = FALSE)
  }
  lines[1] <- sub("^\xef\xbb\xbf", "", lines[1], useBytes = TRUE)
  hdr <- as.character(unlist(utils::read.csv(text = lines[1], header = FALSE,
                                             colClasses = "character",
                                             na.strings = character())[1, ]))
  raw <- utils::read.csv(text = lines[-(1:2)], header = FALSE, colClasses = "numeric",
                         na.strings = c("", "NA"))
  p <- ncol(raw)
  if (p < 2L || p %% 2L != 0L) {
    stop("The file must contain an even number of columns (2-theta / intensity pairs).",
         call. = FALSE)
  }
  starts <- seq(1L, p, by = 2L)
  samples <- trimws(hdr[starts])
  if (any(!nzchar(samples))) {
    stop("Every 2-theta / intensity pair needs a sample name in the first header row.",
         call. = FALSE)
  }
  parts <- lapply(seq_along(starts), function(k) {
    tt <- raw[[starts[k]]]
    y <- raw[[starts[k] + 1L]]
    keep <- !is.na(tt) & !is.na(y)
    tibble::tibble(sample = samples[k], two_theta = tt[keep], intensity = y[keep])
  })
  do.call(rbind, parts)
}
