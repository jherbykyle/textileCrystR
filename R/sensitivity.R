#' Sensitivity of the Segal index to the I200 and Iam positions
#'
#' Evaluates the Segal Crystallinity Index over a grid of I200 and Iam
#' positions and reports the full matrix together with its minimum, maximum,
#' range and a central estimate. The default grid is the validation grid
#' used in this package (I200 = 22.2, 22.4, 22.6, 22.8, 23.0 and Iam = 17.6,
#' 17.8, 18.0, 18.2, 18.4 degrees 2-theta), giving a 5 x 5 matrix.
#'
#' This is a **robustness analysis** that shows how much the empirical index
#' depends on where the two intensities are read. It is not a peak-selection
#' rule and does not identify a "best" position pair.
#'
#' @inheritParams segal_ci
#' @param i200,iam Numeric vectors of positions (degrees 2-theta) to evaluate.
#'   Every I200 position must be greater than every Iam position.
#' @param central Named numeric vector `c(i200 =, iam =)` giving the position
#'   pair used for the central estimate. Defaults to the median of each grid
#'   axis (22.6 and 18.0 degrees for the default grid).
#' @return An object of class `textile_cryst_sensitivity`, a list with
#' * `matrix`: numeric matrix of indices (rows = Iam positions, columns = I200
#'   positions),
#' * `grid`: the same values in long form with the intensities read,
#' * `summary`: one-row tibble with `ci_min`, `ci_max`, `ci_range`,
#'   `ci_central` and the positions at which the extremes occur,
#' * `central`, `interpolation`, `sample`, `input`, `warnings` and `note`.
#'
#' `print()` and `as.data.frame()` methods are provided.
#' @seealso [segal_ci()]
#' @examples
#' tt <- seq(5, 40, by = 0.05)
#' y <- 300 + 4000 * exp(-((tt - 22.7) / 1.1)^2) + 1500 * exp(-((tt - 16) / 1.5)^2)
#' sens <- segal_sensitivity(tt, y, sample = "synthetic example")
#' sens$matrix
#' sens$summary
#' @export
segal_sensitivity <- function(two_theta, intensity,
                              i200 = c(22.2, 22.4, 22.6, 22.8, 23.0),
                              iam = c(17.6, 17.8, 18.0, 18.2, 18.4),
                              central = c(i200 = stats::median(i200), iam = stats::median(iam)),
                              interpolation = c("linear", "nearest"), sample = NULL,
                              na.rm = FALSE, duplicates = c("error", "mean")) {
  interpolation <- .check_interpolation(interpolation)
  sample <- .check_sample(sample)
  for (v in list(list(i200, "i200"), list(iam, "iam"))) {
    if (!is.numeric(v[[1]]) || length(v[[1]]) < 1L) {
      stop(sprintf("`%s` must be a numeric vector of positions.", v[[2]]), call. = FALSE)
    }
    for (p in v[[1]]) .check_position(p, v[[2]])
    if (anyDuplicated(v[[1]])) {
      stop(sprintf("`%s` contains repeated positions.", v[[2]]), call. = FALSE)
    }
  }
  if (min(i200) <= max(iam)) {
    stop("Every I200 position must be greater than every Iam position.", call. = FALSE)
  }
  if (!is.numeric(central) || length(central) != 2L || is.null(names(central)) ||
      !setequal(names(central), c("i200", "iam"))) {
    stop("`central` must be a named numeric vector c(i200 = , iam = ).", call. = FALSE)
  }
  .check_position(central[["i200"]], "central['i200']")
  .check_position(central[["iam"]], "central['iam']")
  input <- validate_crystallinity_input(two_theta, intensity, na.rm = na.rm,
                                        duplicates = match.arg(duplicates))

  loc <- function(p) .locate(input$two_theta, input$intensity, p, interpolation, "position")
  top <- lapply(i200, loc)
  bot <- lapply(iam, loc)
  i200_int <- vapply(top, `[[`, numeric(1), "intensity")
  iam_int <- vapply(bot, `[[`, numeric(1), "intensity")
  if (any(i200_int <= 0)) {
    stop("An I200 intensity in the grid is not positive, so the Segal index is undefined.",
         call. = FALSE)
  }
  warn <- character()
  gaps <- c(vapply(top, function(z) if (z$on_measured_point) 0 else z$gap, numeric(1)),
            vapply(bot, function(z) if (z$on_measured_point) 0 else z$gap, numeric(1)))
  if (any(gaps > .MAX_GAP)) {
    warn <- c(warn, sprintf(paste0("Some grid positions were %s across a gap wider than %.1f ",
                                   "degrees between measured points."),
                            if (interpolation == "linear") "interpolated" else "read", .MAX_GAP))
  }
  if (any(iam_int < 0)) {
    warn <- c(warn, "Some Iam intensities in the grid are negative; check background handling.")
  }

  mat <- outer(iam_int, i200_int, function(a, b) (b - a) / b * 100)
  dimnames(mat) <- list(iam = format(iam), i200 = format(i200))
  grid <- tibble::tibble(
    i200_position = rep(i200, each = length(iam)),
    iam_position = rep(iam, times = length(i200)),
    i200_intensity = rep(i200_int, each = length(iam)),
    iam_intensity = rep(iam_int, times = length(i200)),
    ci = as.numeric(mat)
  )
  cen <- .segal_eval(input, central[["i200"]], central[["iam"]], interpolation)
  imin <- which.min(grid$ci)
  imax <- which.max(grid$ci)
  summ <- tibble::tibble(
    ci_min = grid$ci[imin], ci_max = grid$ci[imax],
    ci_range = grid$ci[imax] - grid$ci[imin], ci_central = cen$ci,
    i200_at_min = grid$i200_position[imin], iam_at_min = grid$iam_position[imin],
    i200_at_max = grid$i200_position[imax], iam_at_max = grid$iam_position[imax]
  )
  structure(list(
    matrix = mat, grid = grid, summary = summ,
    central = list(i200 = central[["i200"]], iam = central[["iam"]], ci = cen$ci),
    interpolation = interpolation, sample = sample,
    input = input[setdiff(names(input), c("two_theta", "intensity"))],
    warnings = c(warn, cen$warnings),
    note = paste0("Robustness analysis of the Segal Crystallinity Index to the choice of I200 ",
                  "and Iam positions; not a universal peak-selection rule.")
  ), class = "textile_cryst_sensitivity")
}
