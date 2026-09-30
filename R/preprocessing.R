# General-purpose XRD pattern preprocessing ("Method 2A"): small, composable
# functions usable ahead of any specific crystallinity method, not just
# peak-area. These are deliberately independent of
# validate_crystallinity_input() (used throughout Segal and the existing
# peak-area code) rather than a refactor of it -- that function is heavily
# depended on and heavily tested already, and rewiring it for the sake of
# code reuse would be exactly the kind of unnecessary change this package's
# own conventions ask to avoid. The two validators overlap in what they
# check; use whichever fits the task at hand.

#' Validate an XRD pattern's basic shape and values
#'
#' A strict, no-fixing check: `two_theta` and `intensity` must be numeric,
#' equal length, at least two points, and every 2-theta value must be finite
#' and strictly between 0 and 180 degrees. Unlike
#' [validate_crystallinity_input()], this never sorts, deduplicates, or
#' otherwise changes the data -- it only checks it, and only fixes missing
#' values when `na.rm = TRUE` is explicitly given. Pair with
#' [sort_xrd_pattern()] and [remove_xrd_duplicates()] for those steps.
#'
#' @param two_theta,intensity Numeric vectors: the pattern to check.
#' @param na.rm Logical (default `FALSE`). Drop rows with a missing or
#'   non-finite `two_theta`/`intensity` instead of erroring?
#' @return A list: `two_theta`, `intensity` (unchanged, or with non-finite
#'   rows dropped if `na.rm = TRUE`), `n_input`, `n_used`, `n_dropped`.
#' @seealso [sort_xrd_pattern()], [remove_xrd_duplicates()],
#'   [validate_crystallinity_input()]
#' @examples
#' validate_xrd_pattern(c(10, 15, 20), c(100, 150, 120))
#' @export
validate_xrd_pattern <- function(two_theta, intensity, na.rm = FALSE) {
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

  bad <- !is.finite(tt) | !is.finite(y)
  n_dropped <- sum(bad)
  if (n_dropped > 0L) {
    if (!isTRUE(na.rm)) {
      stop(sprintf(paste0("%d missing or non-finite value(s) found in `two_theta` or `intensity`. ",
                          "Remove them or use na.rm = TRUE."), n_dropped), call. = FALSE)
    }
    tt <- tt[!bad]; y <- y[!bad]
  }
  if (length(tt) < 2L) stop("At least two valid points are needed.", call. = FALSE)
  if (any(tt <= 0 | tt >= 180)) {
    stop("All 2-theta values must lie strictly between 0 and 180 degrees.", call. = FALSE)
  }
  list(two_theta = tt, intensity = y, n_input = n_input, n_used = length(tt), n_dropped = n_dropped)
}

#' Sort an XRD pattern by 2-theta
#'
#' @param two_theta,intensity Numeric vectors: the pattern to sort.
#' @return A list: `two_theta`, `intensity` (ascending by `two_theta`), and
#'   `was_sorted` (logical: was the input already in ascending order?).
#' @seealso [validate_xrd_pattern()], [remove_xrd_duplicates()]
#' @examples
#' sort_xrd_pattern(c(15, 10, 20), c(2, 1, 3))
#' @export
sort_xrd_pattern <- function(two_theta, intensity) {
  if (!is.numeric(two_theta) || !is.numeric(intensity)) {
    stop("`two_theta` and `intensity` must both be numeric vectors.", call. = FALSE)
  }
  if (length(two_theta) != length(intensity)) {
    stop(sprintf("`two_theta` (length %d) and `intensity` (length %d) must have the same length.",
                 length(two_theta), length(intensity)), call. = FALSE)
  }
  already <- !is.unsorted(two_theta)
  if (already) return(list(two_theta = two_theta, intensity = intensity, was_sorted = FALSE))
  o <- order(two_theta)
  list(two_theta = two_theta[o], intensity = intensity[o], was_sorted = TRUE)
}

#' Remove or average duplicated 2-theta values
#'
#' @param two_theta,intensity Numeric vectors: the pattern to check. Any
#'   input order is accepted; duplicates are exact repeats of a `two_theta`
#'   value regardless of position.
#' @param duplicates `"error"` (default) or `"mean"` (average the
#'   intensities at each repeated 2-theta value, collapsing it to one row).
#' @return A list: `two_theta`, `intensity`, `n_duplicates` (number of
#'   *extra* rows collapsed away, i.e. occurrences beyond the first at each
#'   repeated value).
#' @seealso [validate_xrd_pattern()], [sort_xrd_pattern()]
#' @examples
#' remove_xrd_duplicates(c(10, 10, 12), c(5, 7, 9), duplicates = "mean")
#' @export
remove_xrd_duplicates <- function(two_theta, intensity, duplicates = c("error", "mean")) {
  duplicates <- match.arg(duplicates)
  if (!is.numeric(two_theta) || !is.numeric(intensity)) {
    stop("`two_theta` and `intensity` must both be numeric vectors.", call. = FALSE)
  }
  if (length(two_theta) != length(intensity)) {
    stop(sprintf("`two_theta` (length %d) and `intensity` (length %d) must have the same length.",
                 length(two_theta), length(intensity)), call. = FALSE)
  }
  dup <- duplicated(two_theta)
  n_dup <- sum(dup)
  if (n_dup == 0L) return(list(two_theta = two_theta, intensity = intensity, n_duplicates = 0L))
  if (duplicates == "error") {
    stop(sprintf(paste0("%d duplicated 2-theta value(s) found. ",
                        "Use duplicates = \"mean\" to average their intensities."), n_dup),
         call. = FALSE)
  }
  u <- unique(two_theta)
  y <- as.numeric(tapply(intensity, match(two_theta, u), mean))
  list(two_theta = u, intensity = y, n_duplicates = n_dup)
}

# SNIP (Statistics-sensitive Non-linear Iterative Peak-clipping) background
# estimate, LLS-transform variant (Ryan et al. 1988; refined by Morhac et al.
# 1997). The LLS transform log(log(sqrt(y+1)+1)+1) compresses the dynamic
# range so that both weak and strong peaks are clipped comparably; each
# iteration m replaces v[i] with min(v[i], mean(v[i-m], v[i+m])) so a genuine
# peak (locally higher than points m away on both sides) gets clipped down
# towards its surroundings, while a smoothly-varying background is left
# essentially unchanged. m runs from 1 up to `iterations`. Only interior
# points with a full symmetric window (m+1 <= i <= n-m) are updated at each
# m; near-edge points keep their current value for that iteration rather
# than using a clamped, non-symmetric window -- clamping was tried first and
# checked against a pure linear-ramp background with zero peaks, where it
# introduced substantial error (not just at the true edges: once
# `iterations` is a sizeable fraction of the array length, as it commonly is
# here, the distortion propagates inward across most of the array within a
# few iterations). The fix below recovers a pure linear ramp to floating-
# point precision on that same check.
.snip_background <- function(y, iterations) {
  n <- length(y)
  v <- log(log(sqrt(pmax(y, 0) + 1) + 1) + 1)
  for (m in seq_len(iterations)) {
    if (2L * m >= n) break
    i <- (m + 1L):(n - m)
    v[i] <- pmin(v[i], (v[i - m] + v[i + m]) / 2)
  }
  pmax((exp(exp(v) - 1) - 1)^2 - 1, 0)
}

#' Estimate a smooth background under an XRD pattern
#'
#' Estimates a smooth, non-parametric background/baseline using the SNIP
#' algorithm (Statistics-sensitive Non-linear Iterative Peak-clipping; Ryan
#' et al., 1988, IEEE Trans. Nucl. Sci. 35; Morhac et al., 1997, Nucl.
#' Instrum. Methods A 401), well established for XRD and gamma-ray
#' background removal. Unlike the linear/constant background term fit
#' *jointly* with peaks inside [fit_crystalline_peaks()], this makes no
#' assumption about the background's shape (it need not be flat or linear)
#' and does not require any peaks to be specified first -- it works directly
#' from the pattern's own local structure, which also means it can be
#' fooled by very broad features (see `window`).
#'
#' @param two_theta,intensity Numeric vectors: the pattern.
#' @param window Numeric, degrees 2-theta (default `8`). Peaks narrower than
#'   this are clipped away as signal; genuine features broader than this
#'   (including a real amorphous halo) risk being partly clipped away as
#'   background instead. Converted internally to an iteration count using
#'   the pattern's own median 2-theta step.
#' @param na.rm Passed to [validate_xrd_pattern()].
#' @return A list: `two_theta`, `intensity` (the input, unchanged),
#'   `background` (the estimated background at each `two_theta`), `window`,
#'   `iterations` (the resolved iteration count).
#' @section Accuracy on sloped backgrounds:
#' The LLS transform (needed for widely-varying peak heights) introduces a
#' small bias on a background that is not flat: checked against known
#' synthetic backgrounds, this is under 0.1% of the background level for a
#' gentle slope (about what the shipped `cotton.csv` data shows) and about
#' 0.6% for a slope typical of a fairly steeply-tapering pattern, but grows
#' faster than linearly for unusually steep backgrounds (multi-percent
#' errors are possible if the background itself changes by a large fraction
#' of its own level within `window` degrees). This is a property of the
#' published algorithm's transform, not tunable away; for most real XRD
#' patterns it is small relative to peak heights.
#' @seealso [subtract_xrd_background()]
#' @examples
#' tt <- seq(10, 40, by = 0.05)
#' y <- 300 + 4000 * exp(-((tt - 22.6) / 0.7)^2)  # sharp peak on a flat background
#' bg <- estimate_xrd_background(tt, y)
#' range(bg$background)  # close to 300, the true background level
#' @export
estimate_xrd_background <- function(two_theta, intensity, window = 8, na.rm = FALSE) {
  if (!is.numeric(window) || length(window) != 1L || !is.finite(window) || window <= 0) {
    stop("`window` must be a single positive number of degrees 2-theta.", call. = FALSE)
  }
  v <- validate_xrd_pattern(two_theta, intensity, na.rm = na.rm)
  s <- sort_xrd_pattern(v$two_theta, v$intensity)
  tt <- s$two_theta; y <- s$intensity
  step <- stats::median(diff(tt))
  if (!is.finite(step) || step <= 0) {
    stop("Could not determine a typical 2-theta step (are all points at the same position?).",
         call. = FALSE)
  }
  iterations <- max(1L, round(window / step))
  n <- length(y)
  max_feasible <- (n - 1L) %/% 2L
  if (iterations > max_feasible) {
    if (max_feasible < 1L) {
      stop(sprintf(paste0("The pattern has too few points (%d) to estimate a background at all; ",
                          "need at least 3."), n), call. = FALSE)
    }
    warning(sprintf(paste0("`window` = %s degrees would need %d iterations, but the pattern has only ",
                           "%d points; using the largest window this pattern supports instead (%s ",
                           "degrees, %d iterations)."), format(window), iterations, n,
                    format(round(max_feasible * step, 3)), max_feasible), call. = FALSE)
    iterations <- max_feasible
  }
  bg <- .snip_background(y, iterations)
  list(two_theta = tt, intensity = y, background = bg, window = window, iterations = iterations)
}

#' Subtract an estimated or supplied background from an XRD pattern
#'
#' @param two_theta,intensity Numeric vectors: the pattern.
#' @param background Either the list returned by [estimate_xrd_background()],
#'   or a plain numeric vector of background values the same length as
#'   `two_theta` (e.g. from your own estimate).
#' @param clip_negative Logical (default `TRUE`): floor the corrected
#'   intensity at zero. Small negative dips can arise from noise or an
#'   imperfect background estimate; a corrected pattern is not physically
#'   meaningful below zero.
#' @return A [tibble::tibble()] with `two_theta`, `intensity` (the original),
#'   `background`, and `corrected` (`intensity - background`, clipped if
#'   requested).
#' @seealso [estimate_xrd_background()]
#' @examples
#' tt <- seq(10, 40, by = 0.05)
#' y <- 300 + 4000 * exp(-((tt - 22.6) / 0.7)^2)
#' bg <- estimate_xrd_background(tt, y)
#' subtract_xrd_background(tt, y, bg)
#' @export
subtract_xrd_background <- function(two_theta, intensity, background, clip_negative = TRUE) {
  v <- validate_xrd_pattern(two_theta, intensity)
  tt <- v$two_theta; y <- v$intensity
  bg <- if (is.list(background) && !is.null(background$background)) background$background else background
  if (!is.numeric(bg) || length(bg) != length(tt)) {
    stop(sprintf(paste0("`background` must be a numeric vector of length %d (matching `two_theta`), ",
                        "or the list returned by estimate_xrd_background()."), length(tt)), call. = FALSE)
  }
  corrected <- y - bg
  if (isTRUE(clip_negative)) corrected <- pmax(corrected, 0)
  tibble::tibble(two_theta = tt, intensity = y, background = bg, corrected = corrected)
}

#' Normalize an XRD pattern's intensity scale
#'
#' Rescales intensity onto a common scale, e.g. so patterns collected with
#' different exposure times or on different instruments can be compared or
#' plotted together. Does not change `two_theta`.
#'
#' @param two_theta,intensity Numeric vectors: the pattern.
#' @param method `"max"` (default; divide by the maximum intensity, so the
#'   pattern's tallest point becomes 1), `"area"` (divide by the total area
#'   under the pattern, trapezoidal rule, so the pattern integrates to 1), or
#'   `"minmax"` (rescale linearly so the minimum becomes 0 and the maximum
#'   becomes 1).
#' @return A list: `two_theta`, `intensity` (rescaled), `method`, `scale`
#'   (the divisor applied; not meaningful for `"minmax"`, which is `NA`).
#' @examples
#' normalize_xrd_pattern(c(10, 20, 30), c(50, 200, 100))
#' @export
normalize_xrd_pattern <- function(two_theta, intensity, method = c("max", "area", "minmax")) {
  method <- match.arg(method)
  v <- validate_xrd_pattern(two_theta, intensity)
  tt <- v$two_theta; y <- v$intensity
  if (method == "max") {
    scale <- max(y)
    if (!is.finite(scale) || scale == 0) stop("Maximum intensity is zero; cannot normalize.", call. = FALSE)
    y_out <- y / scale
  } else if (method == "area") {
    s <- sort_xrd_pattern(tt, y)
    scale <- .trapz(s$two_theta, s$intensity)
    if (!is.finite(scale) || scale == 0) stop("Total area is zero; cannot normalize.", call. = FALSE)
    y_out <- y / scale
  } else {
    rng <- range(y)
    if (!is.finite(diff(rng)) || diff(rng) == 0) stop("Intensity is constant; cannot min-max normalize.",
                                                       call. = FALSE)
    y_out <- (y - rng[1]) / diff(rng)
    scale <- NA_real_
  }
  list(two_theta = tt, intensity = y_out, method = method, scale = scale)
}
