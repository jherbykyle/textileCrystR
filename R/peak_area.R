# Peak-area (deconvolution) crystallinity: fit the crystalline peaks and the
# amorphous hump SIMULTANEOUSLY as one multi-component nonlinear model (this
# is the scientifically necessary approach -- the two overlap, so fitting
# either alone to the same data double-counts the other's contribution; see
# Salem et al. 2023 Sec. 2.1.2, "their sum will be adjusted to match the
# experimental pattern"). Crystallinity is then the crystalline peaks' area
# divided by the total (crystalline + amorphous) area.

# Literature-cited default crystalline peak positions for cellulose I,
# degrees 2-theta (Cu K-alpha), from the four-peak convention (012/102 often
# unresolved as a shoulder) described in Salem et al. (2023) Sec. 2.1.2 and
# used for cotton by Nam et al. (2016).
.CELLULOSE_I_PEAK_CENTERS <- c(`1-10` = 14.8, `110` = 16.5, `200` = 22.6, `004` = 34.5)
.CELLULOSE_I_PEAK_TOLERANCE <- 1.2  # deg; how far a fitted center may be labelled by its literature name

.label_peak <- function(center) {
  d <- abs(.CELLULOSE_I_PEAK_CENTERS - center)
  j <- which.min(d)
  if (d[j] <= .CELLULOSE_I_PEAK_TOLERANCE) names(.CELLULOSE_I_PEAK_CENTERS)[j] else NA_character_
}

.check_model <- function(model) match.arg(model, c("pseudo_voigt", "gaussian", "lorentzian", "voigt"))

.check_window <- function(window, name = "window") {
  if (!is.numeric(window) || length(window) != 2L || anyNA(window) || window[1] >= window[2]) {
    stop(sprintf("`%s` must be two increasing, non-missing numbers (degrees 2-theta).", name),
         call. = FALSE)
  }
  window
}

# Rough, data-driven starting FWHM for a peak near `center`: half-width at
# half the local maximum above the local minimum, searched within `radius`
# of `center`. Falls back to `default` when the local data doesn't support a
# clear estimate (too few points, or no local maximum). `clip` bounds the
# result to a plausible range for the kind of peak being estimated (tight
# for sharp crystalline peaks, wide for a broad amorphous hump).
.local_fwhm_guess <- function(x, y, center, radius = 2, default = 1.0, clip = c(0.3, 3)) {
  sel <- x >= center - radius & x <= center + radius
  if (sum(sel) < 5L) return(default)
  xx <- x[sel]; yy <- y[sel]
  i0 <- which.min(abs(xx - center))
  peak_y <- yy[i0]; base_y <- min(yy)
  if (peak_y <= base_y) return(default)
  half <- base_y + (peak_y - base_y) / 2
  n <- length(yy)
  li <- i0; while (li > 1L && yy[li] > half) li <- li - 1L
  ri <- i0; while (ri < n && yy[ri] > half) ri <- ri + 1L
  fwhm <- xx[ri] - xx[li]
  if (!is.finite(fwhm) || fwhm <= 0) return(default)
  min(max(fwhm, clip[1]), clip[2])
}

# Starting parameters and bounds for one component, as named numeric vectors
# in .model_param_names(model) order.
.component_start <- function(model, center, height, fwhm, center_tol, fwhm_max, fwhm_min = 0.05) {
  base <- c(center = center, height = height)
  extra <- switch(model,
    gaussian = ,
    lorentzian = c(fwhm = fwhm),
    pseudo_voigt = c(fwhm = fwhm, eta = 0.5),
    voigt = c(gaussian_fwhm = fwhm / 2, lorentzian_fwhm = fwhm / 2)
  )
  start <- c(base, extra)[.model_param_names(model)]
  lower <- start; upper <- start
  lower["center"] <- center - center_tol; upper["center"] <- center + center_tol
  lower["height"] <- 0; upper["height"] <- height * 5 + 1
  if (model %in% c("gaussian", "lorentzian", "pseudo_voigt")) {
    lower["fwhm"] <- fwhm_min; upper["fwhm"] <- fwhm_max
  } else {
    lower[c("gaussian_fwhm", "lorentzian_fwhm")] <- fwhm_min
    upper[c("gaussian_fwhm", "lorentzian_fwhm")] <- fwhm_max
  }
  if (model == "pseudo_voigt") { lower["eta"] <- 0; upper["eta"] <- 1 }
  list(start = start, lower = lower, upper = upper)
}

# Evaluate the sum of all components. par is a flat vector, k params per
# component, components in column order (matrix(par, nrow = k)).
.total_curve <- function(x, par, model, k) {
  m <- matrix(par, nrow = k)
  nm <- .model_param_names(model)
  out <- numeric(length(x))
  for (j in seq_len(ncol(m))) {
    p <- stats::setNames(m[, j], nm)
    out <- out + .eval_component(x, p, model)
  }
  out
}

# Evaluate the sum of all peak components, plus an optional background term
# appended at the end of `par` (n_bg = 0, 1 [constant] or 2 [linear,
# intercept + slope about x_mid] extra trailing values). par is a flat
# vector, k params per component, components in column order for the
# non-background part (matrix(par[1:(length(par)-n_bg)], nrow = k)).
.total_curve <- function(x, par, model, k, n_bg = 0, x_mid = 0) {
  n_comp_par <- length(par) - n_bg
  m <- matrix(par[seq_len(n_comp_par)], nrow = k)
  nm <- .model_param_names(model)
  out <- numeric(length(x))
  for (j in seq_len(ncol(m))) {
    p <- stats::setNames(m[, j], nm)
    out <- out + .eval_component(x, p, model)
  }
  if (n_bg >= 1L) out <- out + par[n_comp_par + 1L]
  if (n_bg >= 2L) out <- out + par[n_comp_par + 2L] * (x - x_mid)
  out
}

.fit_joint <- function(x, y, starts, lowers, uppers, model, max_iter,
                       bg_start = numeric(0), bg_lower = numeric(0), bg_upper = numeric(0),
                       x_mid = 0) {
  k <- length(.model_param_names(model))
  n_bg <- length(bg_start)
  par0 <- c(unlist(starts), bg_start); lower <- c(unlist(lowers), bg_lower)
  upper <- c(unlist(uppers), bg_upper)
  resid_fn <- function(par) y - .total_curve(x, par, model, k, n_bg, x_mid)
  fit <- minpack.lm::nls.lm(par = par0, lower = lower, upper = upper, fn = resid_fn,
                            control = minpack.lm::nls.lm.control(maxiter = max_iter, maxfev = max_iter * 200))
  n_comp_par <- length(par0) - n_bg
  list(par = matrix(fit$par[seq_len(n_comp_par)], nrow = k, dimnames = list(.model_param_names(model), NULL)),
       bg = if (n_bg > 0) fit$par[(n_comp_par + 1L):(n_comp_par + n_bg)] else numeric(0),
       info = fit$info, message = fit$message, niter = fit$niter,
       converged = fit$info %in% 1:4)
}

#' Fit a broad amorphous hump to an XRD pattern
#'
#' Fits a single, broad peak (using the same peak-shape models as
#' [fit_crystalline_peaks()]) to a specified 2-theta window, as a standalone
#' fit. Used on its own -- for example on a ball-milled or otherwise fully
#' amorphous reference pattern -- or as a source of starting values for the
#' amorphous component of [fit_crystalline_peaks()]'s joint fit (where its
#' parameters are then refined further, jointly with the crystalline peaks).
#'
#' @param two_theta,intensity Numeric vectors: the diffraction pattern.
#' @param window Numeric `c(min, max)`, degrees 2-theta: the region the hump
#'   is fit within. Default `c(18, 25)`, the amorphous-dominated region
#'   commonly used for cellulose (its maximum lies near 20-21 degrees;
#'   Salem et al. 2023 Sec. 2.1.2).
#' @param model One of `"pseudo_voigt"` (default), `"gaussian"`,
#'   `"lorentzian"`, `"voigt"`.
#' @param center Optional starting center (degrees 2-theta); default the
#'   position of the maximum intensity within `window`.
#' @param background One of `"linear"` (default), `"constant"` or `"none"`;
#'   see [fit_crystalline_peaks()] for why this exists. Excluded from `area`.
#' @param sample Optional sample name stored in the result.
#' @param na.rm,duplicates Passed to [validate_crystallinity_input()].
#' @param max_iter Maximum optimizer iterations (default 200).
#'
#' @return An object of class `textile_amorphous_hump`: a list with the
#'   fitted `center`, `height`, width parameter(s) (`fwhm`, or
#'   `gaussian_fwhm`/`lorentzian_fwhm` for `model = "voigt"`), `eta` (for
#'   `pseudo_voigt`), `area` (background excluded), `background` (`model`,
#'   `intercept`, `slope`, `x_mid`), `model`, `window`, `r_squared`, `rmse`,
#'   `converged`, `sample` and the windowed `data` used (with a `background`
#'   column).
#' @seealso [fit_crystalline_peaks()]
#' @examples
#' tt <- seq(10, 40, by = 0.05)
#' y <- 200 + 600 * exp(-((tt - 20.5) / 3)^2)
#' hump <- fit_amorphous_hump(tt, y, window = c(15, 28))
#' hump$center
#' hump$r_squared
#' @export
fit_amorphous_hump <- function(two_theta, intensity, window = c(18, 25),
                               model = c("pseudo_voigt", "gaussian", "lorentzian", "voigt"),
                               center = NULL, background = c("linear", "constant", "none"),
                               sample = NULL, na.rm = FALSE,
                               duplicates = c("error", "mean"), max_iter = 200) {
  model <- .check_model(model)
  window <- .check_window(window)
  sample <- .check_sample(sample)
  background <- match.arg(background)
  input <- validate_crystallinity_input(two_theta, intensity, na.rm = na.rm,
                                        duplicates = match.arg(duplicates))
  sel <- input$two_theta >= window[1] & input$two_theta <= window[2]
  if (sum(sel) < 5L) {
    stop(sprintf(paste0("Only %d measured point(s) fall inside window = c(%s, %s); need at least 5 ",
                        "to fit a peak. Check the window against the pattern's 2-theta range (%s to %s)."),
                sum(sel), format(window[1]), format(window[2]), format(input$range[1]), format(input$range[2])),
         call. = FALSE)
  }
  x <- input$two_theta[sel]; y <- input$intensity[sel]

  if (is.null(center)) center <- x[which.max(y)]
  .check_position(center, "center")
  x_mid <- mean(window)
  n_bg <- c(none = 0L, constant = 1L, linear = 2L)[[background]]
  bg0_start <- max(as.numeric(stats::quantile(y, 0.05)), 0)
  height0 <- max(max(y) - bg0_start, 1)
  fwhm0 <- .local_fwhm_guess(x, y, center, radius = diff(window) / 2, default = diff(window) * 0.6,
                             clip = c(1, diff(window) * 1.5))
  bg_start <- if (n_bg >= 1L) c(bg0_start, if (n_bg == 2L) 0) else numeric(0)
  bg_lower <- if (n_bg >= 1L) c(0, if (n_bg == 2L) -Inf) else numeric(0)
  bg_upper <- if (n_bg >= 1L) c(max(y), if (n_bg == 2L) Inf) else numeric(0)

  sp <- .component_start(model, center, height0, fwhm0, center_tol = diff(window) / 2,
                         fwhm_max = diff(window) * 4)
  fit <- .fit_joint(x, y, list(sp$start), list(sp$lower), list(sp$upper), model, max_iter,
                    bg_start = bg_start, bg_lower = bg_lower, bg_upper = bg_upper, x_mid = x_mid)
  p <- stats::setNames(fit$par[, 1], .model_param_names(model))
  peak_curve <- .eval_component(x, p, model)
  background_curve <- if (n_bg == 0L) rep(0, length(x))
                       else if (n_bg == 1L) rep(fit$bg[1], length(x))
                       else fit$bg[1] + fit$bg[2] * (x - x_mid)
  fitted_curve <- peak_curve + background_curve

  out <- as.list(p)
  out$area <- .component_area(p, model)
  out$fwhm_reported <- .component_fwhm(p, model)
  out$background <- list(model = background, x_mid = x_mid,
                         intercept = if (n_bg >= 1L) unname(fit$bg[1]) else NA_real_,
                         slope = if (n_bg == 2L) unname(fit$bg[2]) else NA_real_)
  out$model <- model
  out$window <- window
  out$r_squared <- .r_squared(y, fitted_curve)
  out$rmse <- .rmse(y, fitted_curve)
  out$converged <- fit$converged
  out$optimizer_message <- fit$message
  out$sample <- sample
  out$data <- tibble::tibble(two_theta = x, intensity = y, fitted = fitted_curve,
                             background = background_curve)
  structure(out, class = "textile_amorphous_hump")
}

#' Fit crystalline peaks and the amorphous hump jointly (peak-area method)
#'
#' Fits several sharp crystalline peaks and one broad amorphous hump
#' **simultaneously** to an XRD pattern -- the scientifically necessary
#' approach, since the two overlap and fitting either alone to the same data
#' would double-count the other's contribution (Salem et al., 2023, Sec.
#' 2.1.2). Crystallinity is the crystalline peaks' combined area divided by
#' the total (crystalline + amorphous) area, following the peak-area /
#' deconvolution approach reviewed there and in Nam et al. (2016).
#'
#' @param two_theta,intensity Numeric vectors: the diffraction pattern.
#' @param centers Optional numeric vector of starting positions (degrees
#'   2-theta) for the crystalline peaks. Default the literature four-peak
#'   cellulose I convention, `c(14.8, 16.5, 22.6, 34.5)` -- the (1-10), (110),
#'   (200) and (004) reflections (Salem et al., 2023; Nam et al., 2016); the
#'   (012/102) shoulder near 20.5 degrees is often unresolved and is not
#'   included by default. Pass your own vector for other polymorphs or fibres.
#' @param model One of `"pseudo_voigt"` (default), `"gaussian"`,
#'   `"lorentzian"`, `"voigt"`. All crystalline peaks and the amorphous hump
#'   use the same model.
#' @param window Numeric `c(min, max)`, degrees 2-theta, the region fit as a
#'   whole. Default `c(10, 40)`, clipped to the data's own range.
#' @param amorphous Either a numeric `c(min, max)` window (default
#'   `c(18, 25)`) used to fit a preliminary amorphous hump via
#'   [fit_amorphous_hump()] for starting values, or an existing
#'   `textile_amorphous_hump` object (e.g. fit on a separate amorphous
#'   reference pattern) whose parameters are used as the starting point
#'   instead. Either way, the amorphous component is then refined jointly
#'   with the crystalline peaks, not held fixed.
#' @param sample Optional sample name stored in the result.
#' @param background One of `"linear"` (default), `"constant"` or `"none"`.
#'   Real XRD patterns typically carry a broad, roughly flat-to-sloping
#'   background (air scatter, detector dark counts, ...) well outside any
#'   peak that a "crystalline peaks + one amorphous hump" model has no way to
#'   represent; left unmodelled, peaks are pulled towards implausible widths
#'   trying to absorb it (this was observed on the package's own experimental
#'   validation data before adding this option). `"linear"` fits an
#'   intercept and slope, `"constant"` an intercept only, `"none"` fits no
#'   background term at all (matching the "amorphous hump absorbs everything
#'   else" assumption used in some published deconvolutions). The background
#'   is **never** counted as crystalline or amorphous area; it is excluded
#'   from both before `ci` is calculated.
#' @param na.rm,duplicates Passed to [validate_crystallinity_input()].
#' @param max_iter Maximum optimizer iterations for the joint fit (default
#'   500).
#'
#' @return An object of class `textile_cryst_peakarea`, a list with:
#' * `ci`: crystallinity (percent) = crystalline area / total area * 100,
#'   with the background (if any) excluded from both,
#' * `crystalline_area`, `amorphous_area`, `total_area`,
#' * `background`: a list with `model` (`"none"`/`"constant"`/`"linear"`),
#'   `intercept`, `slope` (`NA` if not fit), and `x_mid` (the slope is about
#'   this 2-theta value, for numerical stability),
#' * `peaks`: a tibble, one row per crystalline peak (`label`, `center`,
#'   `height`, and the model's width parameter(s), `area`,
#'   `pct_of_crystalline`),
#' * `amorphous`: the same kind of row for the amorphous hump,
#' * `r_squared`, `rmse`: fit quality over the full `window` (background
#'   included, since that is what was actually fit to the data),
#' * `residuals`: a tibble (`two_theta`, `observed`, `fitted`, `background`,
#'   `residual`),
#' * `converged`, `optimizer_message`,
#' * `data`: the windowed pattern with the total fitted curve,
#' * `metadata`: `method`, `equation`, `assumptions`, `preprocessing`,
#'   `parameters`, `reference`, `package_version`,
#' * `sample`, `model`, `window`.
#'
#' It has `print()`, `summary()`, `plot()` and `as.data.frame()` methods.
#' @seealso [fit_amorphous_hump()], [peak_area_crystallinity()],
#'   [peak_fit_summary()], [plot_peak_fit()]
#' @examples
#' tt <- seq(8, 40, by = 0.05)
#' y <- 200 + 500 * exp(-((tt - 20.5) / 3)^2) +
#'   1800 * exp(-((tt - 22.6) / 0.7)^2) + 600 * exp(-((tt - 16.3) / 0.8)^2) +
#'   400 * exp(-((tt - 14.8) / 0.7)^2) + 300 * exp(-((tt - 34.5) / 1)^2)
#' fit <- fit_crystalline_peaks(tt, y, sample = "synthetic example")
#' fit
#' @export
fit_crystalline_peaks <- function(two_theta, intensity, centers = NULL,
                                  model = c("pseudo_voigt", "gaussian", "lorentzian", "voigt"),
                                  window = c(10, 40), amorphous = c(18, 25), sample = NULL,
                                  background = c("linear", "constant", "none"),
                                  na.rm = FALSE, duplicates = c("error", "mean"), max_iter = 500) {
  model <- .check_model(model)
  window <- .check_window(window)
  sample <- .check_sample(sample)
  background <- match.arg(background)
  duplicates <- match.arg(duplicates)
  if (is.null(centers)) centers <- unname(.CELLULOSE_I_PEAK_CENTERS[c("1-10", "110", "200", "004")])
  if (!is.numeric(centers) || length(centers) < 1L || anyNA(centers)) {
    stop("`centers` must be a numeric vector of at least one 2-theta position.", call. = FALSE)
  }

  input <- validate_crystallinity_input(two_theta, intensity, na.rm = na.rm, duplicates = duplicates)
  window <- c(max(window[1], input$range[1]), min(window[2], input$range[2]))
  window <- .check_window(window)
  sel <- input$two_theta >= window[1] & input$two_theta <= window[2]
  n_par_per <- length(.model_param_names(model))
  if (sum(sel) < (length(centers) + 1L) * n_par_per) {
    stop(sprintf(paste0("Only %d measured point(s) fall inside window = c(%s, %s), too few to fit %d ",
                        "crystalline peak(s) plus the amorphous hump. Widen `window` or reduce the ",
                        "number of `centers`."), sum(sel), format(window[1]), format(window[2]),
                length(centers)), call. = FALSE)
  }
  x <- input$two_theta[sel]; y <- input$intensity[sel]
  interp_y <- function(pos) stats::approx(x, y, xout = pos, rule = 2)$y

  # Preliminary amorphous estimate (used only as a starting point; refined
  # jointly below).
  if (inherits(amorphous, "textile_amorphous_hump")) {
    am_hump <- amorphous
    if (am_hump$model != model) {
      stop(sprintf(paste0("`amorphous` was fit with model = \"%s\" but this call uses model = \"%s\". ",
                          "Fit it again with fit_amorphous_hump(..., model = \"%s\") first."),
                  am_hump$model, model, model), call. = FALSE)
    }
  } else {
    am_window <- .check_window(amorphous, "amorphous")
    am_hump <- tryCatch(
      fit_amorphous_hump(x, y, window = am_window, model = model, background = background,
                        max_iter = max_iter),
      error = function(e) {
        stop(sprintf(paste0("Could not fit a preliminary amorphous hump in window = c(%s, %s): %s\n",
                            "Pass a wider `amorphous` window, or a `textile_amorphous_hump` object ",
                            "fit separately."), format(am_window[1]), format(am_window[2]),
                    conditionMessage(e)), call. = FALSE)
      }
    )
  }

  center_tol <- max(diff(range(centers)) / max(length(centers) - 1L, 1L) / 2, 0.6)
  crystalline_starts <- lapply(centers, function(cc) {
    local_fwhm0 <- .local_fwhm_guess(x, y, cc)
    .component_start(model, cc, max(interp_y(cc), 1), fwhm = local_fwhm0, center_tol = center_tol,
                     fwhm_max = 3, fwhm_min = 0.05)
  })
  am_p <- am_hump[.model_param_names(model)]
  am_fwhm0 <- if (model == "voigt") am_p$gaussian_fwhm + am_p$lorentzian_fwhm else am_p$fwhm
  # The amorphous hump is broad by physical assumption -- its FWHM floor is
  # kept well above the crystalline peaks' ceiling so the two components
  # cannot trade places during optimization (a real failure mode: without
  # this separation, a wide "crystalline" peak and a narrow "amorphous" one
  # can swap roles and converge to a self-consistent but wrong local optimum;
  # caught here by validating against known-parameter synthetic data before
  # this fix, not merely by inspecting the code).
  am_start <- .component_start(model, am_p$center, am_p$height, fwhm = am_fwhm0,
                               center_tol = min(diff(window) / 6, 2.5), fwhm_max = diff(window) * 1.2,
                               fwhm_min = 4)
  am_start$start <- unlist(am_p)[.model_param_names(model)]

  x_mid <- mean(window)
  n_bg <- c(none = 0L, constant = 1L, linear = 2L)[[background]]
  bg0_start <- max(as.numeric(stats::quantile(y, 0.05)), 0)
  bg_start <- if (n_bg >= 1L) c(bg0_start, if (n_bg == 2L) 0) else numeric(0)
  bg_lower <- if (n_bg >= 1L) c(0, if (n_bg == 2L) -Inf) else numeric(0)
  bg_upper <- if (n_bg >= 1L) c(max(y), if (n_bg == 2L) Inf) else numeric(0)

  starts <- c(crystalline_starts, list(am_start))
  fit <- .fit_joint(x, y, lapply(starts, `[[`, "start"), lapply(starts, `[[`, "lower"),
                    lapply(starts, `[[`, "upper"), model, max_iter,
                    bg_start = bg_start, bg_lower = bg_lower, bg_upper = bg_upper, x_mid = x_mid)
  if (!fit$converged) {
    warning_msg <- sprintf(paste0("The optimizer did not clearly converge (message: \"%s\"). Results ",
                                  "may be unreliable; consider different starting `centers`, a wider ",
                                  "`window`, or more `max_iter`."), fit$message)
  } else {
    warning_msg <- NULL
  }

  nm <- .model_param_names(model)
  n_cryst <- length(centers)
  comp_list <- lapply(seq_len(ncol(fit$par)), function(j) stats::setNames(fit$par[, j], nm))
  cryst_comps <- comp_list[seq_len(n_cryst)]
  am_comp <- comp_list[[n_cryst + 1L]]

  cryst_areas <- vapply(cryst_comps, .component_area, numeric(1), model = model)
  am_area <- .component_area(am_comp, model)
  crystalline_area <- sum(cryst_areas)
  total_area <- crystalline_area + am_area
  ci <- crystalline_area / total_area * 100

  peak_row <- function(p, area, label) {
    tibble::tibble(
      label = label, center = p[["center"]], height = p[["height"]],
      fwhm = .component_fwhm(p, model),
      gaussian_fwhm = if (model == "voigt") p[["gaussian_fwhm"]] else NA_real_,
      lorentzian_fwhm = if (model == "voigt") p[["lorentzian_fwhm"]] else NA_real_,
      eta = if (model == "pseudo_voigt") p[["eta"]] else NA_real_,
      area = area
    )
  }
  peaks <- do.call(rbind, Map(function(p, a) peak_row(p, a, .label_peak(p[["center"]]) %||%
                                                           sprintf("Peak %d", which(vapply(cryst_comps, identical, logical(1), p)))),
                              cryst_comps, cryst_areas))
  peaks$label[is.na(peaks$label)] <- sprintf("Peak %d", seq_along(peaks$label))[is.na(peaks$label)]
  peaks$pct_of_crystalline <- peaks$area / crystalline_area * 100
  amorphous_row <- peak_row(am_comp, am_area, "Amorphous")
  amorphous_row$pct_of_crystalline <- NA_real_

  total_fitted_no_bg <- .total_curve(x, unlist(comp_list), model, n_par_per)
  background_curve <- if (n_bg == 0L) rep(0, length(x))
                       else if (n_bg == 1L) rep(fit$bg[1], length(x))
                       else fit$bg[1] + fit$bg[2] * (x - x_mid)
  total_fitted <- total_fitted_no_bg + background_curve
  residuals <- tibble::tibble(two_theta = x, observed = y, fitted = total_fitted,
                              background = background_curve, residual = y - total_fitted)
  background_info <- list(model = background, x_mid = x_mid,
                          intercept = if (n_bg >= 1L) unname(fit$bg[1]) else NA_real_,
                          slope = if (n_bg == 2L) unname(fit$bg[2]) else NA_real_)

  equation <- paste0("CI = Area(crystalline peaks) / [Area(crystalline peaks) + Area(amorphous hump)] * 100, ",
                     "all peaks fit jointly as ", model, " profiles",
                     if (background != "none") paste0(" with a ", background, " background term (excluded",
                                                       " from crystalline and amorphous area)") else "")
  metadata <- list(
    method = "Peak area / deconvolution",
    equation = equation,
    assumptions = c(
      "The pattern is fully described by the sum of the specified crystalline peaks, one amorphous hump,",
      if (background == "none") "and nothing else (no background term is fit)."
      else sprintf("and a %s background term.", background),
      "The background (if fit) is excluded from both the crystalline and amorphous area.",
      "All peak components share the same peak-shape model.",
      "Peak centers are constrained near their starting positions, not free to relabel to a different reflection."
    ),
    preprocessing = sprintf("No smoothing is applied before fitting; background handling is \"%s\" (see `background` field).",
                            background),
    parameters = list(model = model, window = window, centers = centers, n_crystalline_peaks = n_cryst,
                      amorphous_window = if (inherits(amorphous, "textile_amorphous_hump")) NA else amorphous,
                      background = background, max_iter = max_iter),
    reference = paste0(
      "Salem, K. S. et al. (2023) Chem. Soc. Rev., Sec. 2.1.2 (peak fitting / deconvolution method); ",
      "Nam, S. et al. (2016) Carbohydr. Polym. 135, 1-9 (cellulose I peak positions for cotton)."),
    package_version = as.character(utils::packageVersion("textileCrystR"))
  )

  structure(list(
    ci = ci, crystalline_area = crystalline_area, amorphous_area = am_area, total_area = total_area,
    background = background_info,
    peaks = peaks, amorphous = amorphous_row,
    r_squared = .r_squared(y, total_fitted), rmse = .rmse(y, total_fitted),
    residuals = residuals, converged = fit$converged, optimizer_message = fit$message,
    data = tibble::tibble(two_theta = x, intensity = y, fitted = total_fitted, background = background_curve),
    metadata = metadata, sample = sample, model = model, window = window,
    warnings = if (is.null(warning_msg)) character() else warning_msg
  ), class = "textile_cryst_peakarea")
}

#' Crystallinity from a peak-area (deconvolution) fit
#'
#' A convenience wrapper: fits the crystalline peaks and amorphous hump
#' jointly with [fit_crystalline_peaks()] and returns the same
#' `textile_cryst_peakarea` result. Provided as a directly discoverable,
#' crystallinity-named entry point alongside [segal_ci()]; the two functions
#' do the same work; use whichever name you find clearer.
#'
#' @inheritParams fit_crystalline_peaks
#' @return As [fit_crystalline_peaks()]: a `textile_cryst_peakarea` object.
#' @seealso [fit_crystalline_peaks()]
#' @examples
#' tt <- seq(8, 40, by = 0.05)
#' y <- 200 + 500 * exp(-((tt - 20.5) / 3)^2) + 1800 * exp(-((tt - 22.6) / 0.7)^2)
#' peak_area_crystallinity(tt, y, centers = 22.6)$ci
#' @export
peak_area_crystallinity <- function(two_theta, intensity, centers = NULL,
                                    model = c("pseudo_voigt", "gaussian", "lorentzian", "voigt"),
                                    window = c(10, 40), amorphous = c(18, 25), sample = NULL,
                                    background = c("linear", "constant", "none"),
                                    na.rm = FALSE, duplicates = c("error", "mean"), max_iter = 500) {
  fit_crystalline_peaks(two_theta, intensity, centers = centers, model = model, window = window,
                        amorphous = amorphous, sample = sample, background = background, na.rm = na.rm,
                        duplicates = duplicates, max_iter = max_iter)
}

#' Tidy per-component summary of a peak-area fit
#'
#' @param fit A `textile_cryst_peakarea` from [fit_crystalline_peaks()].
#' @return A tibble, one row per crystalline peak plus one row for the
#'   amorphous hump: `label`, `center`, `height`, `fwhm` (composite FWHM for
#'   Voigt), `area`, `pct_of_crystalline` (`NA` for the amorphous row).
#' @seealso [fit_crystalline_peaks()]
#' @examples
#' tt <- seq(8, 40, by = 0.05)
#' y <- 200 + 500 * exp(-((tt - 20.5) / 3)^2) + 1800 * exp(-((tt - 22.6) / 0.7)^2)
#' fit <- fit_crystalline_peaks(tt, y, centers = 22.6)
#' peak_fit_summary(fit)
#' @export
peak_fit_summary <- function(fit) {
  if (!inherits(fit, "textile_cryst_peakarea")) {
    stop("`fit` must be a textile_cryst_peakarea object from fit_crystalline_peaks().", call. = FALSE)
  }
  cols <- c("label", "center", "height", "fwhm", "area", "pct_of_crystalline")
  rbind(fit$peaks[, cols], fit$amorphous[, cols])
}
