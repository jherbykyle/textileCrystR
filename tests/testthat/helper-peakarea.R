# Synthetic ground-truth pattern generators for peak-area tests. Each returns
# list(two_theta, intensity, true_ci, true_crystalline_area, true_amorphous_area)
# for a *clean* pattern (no injected background) built from the package's own
# internal peak-shape functions, so recovery can be checked against exactly
# known values.

.synthetic_peakarea <- function(model = c("gaussian", "lorentzian", "pseudo_voigt", "voigt"),
                                noise_sd = 3, seed = 1, bg0 = 0, bg1 = 0) {
  model <- match.arg(model)
  tt <- seq(8, 40, by = 0.02)
  centers <- c(14.8, 16.5, 22.6, 34.5)
  heights <- c(400, 650, 1800, 320)
  fwhm <- c(0.9, 0.85, 0.75, 1.1)
  am_center <- 20.3; am_height <- 480; am_fwhm <- 7.5

  if (model %in% c("gaussian", "lorentzian")) {
    fn <- if (model == "gaussian") textileCrystR:::.gaussian_peak else textileCrystR:::.lorentzian_peak
    area_fn <- if (model == "gaussian") textileCrystR:::.gaussian_area else textileCrystR:::.lorentzian_area
    y_cryst <- Reduce(`+`, Map(function(c, h, f) fn(tt, c, h, f), centers, heights, fwhm))
    y_am <- fn(tt, am_center, am_height, am_fwhm)
    true_ca <- sum(mapply(area_fn, heights, fwhm))
    true_aa <- area_fn(am_height, am_fwhm)
  } else if (model == "pseudo_voigt") {
    eta <- c(0.3, 0.5, 0.4, 0.6); am_eta <- 0.5
    fn <- textileCrystR:::.pseudo_voigt_peak; area_fn <- textileCrystR:::.pseudo_voigt_area
    y_cryst <- Reduce(`+`, Map(function(c, h, f, e) fn(tt, c, h, f, e), centers, heights, fwhm, eta))
    y_am <- fn(tt, am_center, am_height, am_fwhm, am_eta)
    true_ca <- sum(mapply(area_fn, heights, fwhm, eta))
    true_aa <- area_fn(am_height, am_fwhm, am_eta)
  } else {
    gfwhm <- c(0.7, 0.6, 0.55, 0.9); lfwhm <- c(0.4, 0.5, 0.3, 0.6)
    am_gfwhm <- 6; am_lfwhm <- 3
    fn <- textileCrystR:::.voigt_peak; area_fn <- textileCrystR:::.voigt_area
    y_cryst <- Reduce(`+`, Map(function(c, h, gf, lf) fn(tt, c, h, gf, lf), centers, heights, gfwhm, lfwhm))
    y_am <- fn(tt, am_center, am_height, am_gfwhm, am_lfwhm)
    true_ca <- sum(mapply(area_fn, heights, gfwhm, lfwhm))
    true_aa <- area_fn(am_height, am_gfwhm, am_lfwhm)
  }
  set.seed(seed)
  y_bg <- bg0 + bg1 * (tt - mean(range(tt)))
  y <- y_cryst + y_am + y_bg + rnorm(length(tt), sd = noise_sd)
  list(two_theta = tt, intensity = y, centers = centers,
       true_ci = true_ca / (true_ca + true_aa) * 100,
       true_crystalline_area = true_ca, true_amorphous_area = true_aa)
}
