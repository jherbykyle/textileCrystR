test_that("peak_area_crystallinity is a thin wrapper equivalent to fit_crystalline_peaks", {
  s <- .synthetic_peakarea("gaussian", seed = 1)
  a <- fit_crystalline_peaks(s$two_theta, s$intensity, centers = s$centers, model = "gaussian",
                             amorphous = c(17, 24), sample = "x")
  b <- peak_area_crystallinity(s$two_theta, s$intensity, centers = s$centers, model = "gaussian",
                               amorphous = c(17, 24), sample = "x")
  expect_s3_class(b, "textile_cryst_peakarea")
  expect_equal(a$ci, b$ci)
  expect_equal(a$peaks, b$peaks)
})

test_that("peak_fit_summary returns one row per crystalline peak plus the amorphous hump", {
  s <- .synthetic_peakarea("gaussian", seed = 1)
  fit <- fit_crystalline_peaks(s$two_theta, s$intensity, centers = s$centers, model = "gaussian",
                               amorphous = c(17, 24))
  tab <- peak_fit_summary(fit)
  expect_equal(nrow(tab), 5L)  # 4 crystalline + 1 amorphous
  expect_true(all(c("label", "center", "height", "fwhm", "area", "pct_of_crystalline") %in% names(tab)))
  expect_equal(tab$label[5], "Amorphous")
  expect_true(is.na(tab$pct_of_crystalline[5]))
  expect_equal(sum(tab$pct_of_crystalline[1:4]), 100, tolerance = 1e-6)
})

test_that("peak_fit_summary requires a textile_cryst_peakarea object", {
  expect_error(peak_fit_summary(1), "textile_cryst_peakarea")
  expect_error(peak_fit_summary(list(a = 1)), "textile_cryst_peakarea")
})

# ---- input validation / error handling ----

test_that("fit_crystalline_peaks validates model, window and centers", {
  s <- .synthetic_peakarea("gaussian", seed = 1)
  f <- function(...) fit_crystalline_peaks(s$two_theta, s$intensity, ...)
  expect_error(f(model = "spline"), "should be one of")
  expect_error(f(window = c(20, 10)), "increasing")
  expect_error(f(window = c(20, 20)), "increasing")
  expect_error(f(centers = character()), "numeric vector")
  expect_error(f(centers = NA_real_), "numeric vector")
  expect_error(f(background = "quadratic"), "should be one of")
})

test_that("fit_crystalline_peaks errors clearly when too few points fall in the window", {
  tt <- seq(10, 40, by = 5)
  y <- rep(500, length(tt))
  expect_error(fit_crystalline_peaks(tt, y, centers = c(14.8, 16.5, 22.6, 34.5)), "too few")
})

test_that("fit_crystalline_peaks errors clearly when the preliminary amorphous fit fails", {
  s <- .synthetic_peakarea("gaussian", seed = 1)
  expect_error(
    fit_crystalline_peaks(s$two_theta, s$intensity, centers = s$centers, amorphous = c(50, 60)),
    "amorphous"
  )
})

test_that("fit_crystalline_peaks accepts a pre-fit textile_amorphous_hump and enforces matching models", {
  s <- .synthetic_peakarea("gaussian", seed = 1)
  hump <- fit_amorphous_hump(s$two_theta, s$intensity, window = c(17, 24), model = "gaussian")
  fit <- fit_crystalline_peaks(s$two_theta, s$intensity, centers = s$centers, model = "gaussian",
                               amorphous = hump)
  expect_s3_class(fit, "textile_cryst_peakarea")
  expect_equal(fit$ci, s$true_ci, tolerance = 0.5)

  hump_lorentzian <- fit_amorphous_hump(s$two_theta, s$intensity, window = c(17, 24), model = "lorentzian")
  expect_error(
    fit_crystalline_peaks(s$two_theta, s$intensity, centers = s$centers, model = "gaussian",
                          amorphous = hump_lorentzian),
    "lorentzian.*gaussian"
  )
})

test_that("fit_crystalline_peaks passes na.rm/duplicates through to input validation", {
  s <- .synthetic_peakarea("gaussian", seed = 1)
  tt <- s$two_theta; y <- s$intensity
  tt_dup <- c(tt, tt[1]); y_dup <- c(y, y[1] + 1)
  expect_error(fit_crystalline_peaks(tt_dup, y_dup, centers = s$centers), "duplicated")
  fit <- fit_crystalline_peaks(tt_dup, y_dup, centers = s$centers, duplicates = "mean")
  expect_s3_class(fit, "textile_cryst_peakarea")

  y_na <- y; y_na[1] <- NA
  expect_error(fit_crystalline_peaks(tt, y_na, centers = s$centers), "missing")
  fit2 <- fit_crystalline_peaks(tt, y_na, centers = s$centers, na.rm = TRUE)
  expect_s3_class(fit2, "textile_cryst_peakarea")
})

test_that("an unconverged fit is flagged with a warning message, not silently trusted", {
  # Deliberately implausible centers (nowhere near real peaks) with a very
  # low iteration cap, to reliably trigger non-convergence.
  s <- .synthetic_peakarea("gaussian", seed = 1)
  fit <- suppressWarnings(fit_crystalline_peaks(s$two_theta, s$intensity, centers = c(11, 12, 13, 38),
                                                amorphous = c(17, 24), max_iter = 1))
  if (!fit$converged) {
    expect_length(fit$warnings, 1L)
    expect_match(fit$warnings, "did not clearly converge")
  }
})

# ---- S3 methods for textile_cryst_peakarea ----

test_that("print/summary/as.data.frame/plot methods work for textile_cryst_peakarea", {
  s <- .synthetic_peakarea("gaussian", seed = 1)
  fit <- fit_crystalline_peaks(s$two_theta, s$intensity, centers = s$centers, model = "gaussian",
                               amorphous = c(17, 24), sample = "syn")

  out <- capture.output(print(fit))
  expect_true(any(grepl("textile_cryst_peakarea", out)))
  expect_true(any(grepl("Sample     : syn", out)))
  expect_true(any(grepl("CI", out)))
  ret <- NULL
  capture.output(ret <- print(fit))
  expect_identical(ret, fit)

  s3 <- summary(fit)
  expect_s3_class(s3, "summary.textile_cryst_peakarea")
  expect_equal(s3$ci, fit$ci)
  sout <- capture.output(print(s3))
  expect_true(any(grepl("Peak-area", sout)))
  expect_true(any(grepl("Reference", sout)))
  expect_true(any(grepl("empirical, relative fraction|not an absolute crystalline", paste(sout, collapse = " "))))

  df <- as.data.frame(fit)
  expect_equal(nrow(df), 1L)
  expect_true(all(c("sample", "ci", "crystalline_area", "amorphous_area", "total_area", "n_peaks",
                    "model", "r_squared", "rmse", "converged", "n_warnings") %in% names(df)))
  expect_equal(df$ci, fit$ci)
  expect_equal(df$n_peaks, 4L)

  p <- plot(fit)
  expect_true(inherits(p, "patchwork") || inherits(p, "ggplot"))
})

test_that("plot_peak_fit requires a textile_cryst_peakarea object and honours show_* options", {
  expect_error(plot_peak_fit(1), "textile_cryst_peakarea")
  s <- .synthetic_peakarea("gaussian", seed = 1)
  fit <- fit_crystalline_peaks(s$two_theta, s$intensity, centers = s$centers, amorphous = c(17, 24))

  p_full <- plot_peak_fit(fit)
  expect_s3_class(p_full, "patchwork")

  p_nores <- plot_peak_fit(fit, show_residuals = FALSE)
  expect_s3_class(p_nores, "ggplot")
  expect_false(inherits(p_nores, "patchwork"))

  p_nocomp <- plot_peak_fit(fit, show_components = FALSE, show_residuals = FALSE)
  expect_s3_class(p_nocomp, "ggplot")
})

test_that("non-convergence is surfaced in print() for both amorphous hump and peak-area fits", {
  tt <- seq(10, 40, by = 0.02)
  y <- 300 + 600 * exp(-4 * log(2) * (tt - 20.5)^2 / 6^2)
  hump <- suppressWarnings(fit_amorphous_hump(tt, y, window = c(12, 30), center = 25, max_iter = 1))
  if (!hump$converged) {
    out <- capture.output(print(hump))
    expect_true(any(grepl("optimizer did not clearly converge", out)))
  }

  s <- .synthetic_peakarea("gaussian", seed = 1)
  fit <- suppressWarnings(fit_crystalline_peaks(s$two_theta, s$intensity, centers = c(11, 12, 13, 38),
                                                amorphous = c(17, 24), max_iter = 1))
  if (!fit$converged) {
    out <- capture.output(print(fit))
    expect_true(any(grepl("Warning", out)))
    sout <- capture.output(print(summary(fit)))
    expect_true(any(grepl("Warning", sout)))
  }
})

test_that("plot_peak_fit's colour palette extends beyond its 6 named defaults", {
  # Seven crystalline peaks packed into a narrow range, forcing the
  # colorRampPalette() fallback in .peak_colours() (the base palette has
  # exactly 6 colours) rather than the direct lookup.
  tt <- seq(8, 40, by = 0.02)
  centers <- seq(12, 17, length.out = 7)
  y <- 300 + Reduce(`+`, lapply(centers, function(c) 200 * exp(-4 * log(2) * (tt - c)^2 / 0.25^2))) +
    250 * exp(-4 * log(2) * (tt - 20.5)^2 / 6^2)
  fit <- fit_crystalline_peaks(tt, y, centers = centers, amorphous = c(18, 25))
  expect_equal(nrow(fit$peaks), 7L)
  p <- plot_peak_fit(fit)
  expect_s3_class(p, "patchwork")
})

test_that(".local_fwhm_guess falls back to `default` when data can't support an estimate", {
  expect_equal(.local_fwhm_guess(1:3, c(1, 2, 3), 2, default = 1.23), 1.23)  # too few points
  x <- seq(10, 30, by = 0.5)
  expect_equal(.local_fwhm_guess(x, rep(5, length(x)), 20, default = 2.5), 2.5)  # flat: peak_y == base_y
})

test_that(".model_param_names errors on an unreachable/unknown model", {
  expect_error(.model_param_names("not_a_model"), "Unreachable")
})

test_that("metadata fields are present and populated as documented", {
  s <- .synthetic_peakarea("gaussian", seed = 1)
  fit <- fit_crystalline_peaks(s$two_theta, s$intensity, centers = s$centers, amorphous = c(17, 24))
  m <- fit$metadata
  expect_true(all(c("method", "equation", "assumptions", "preprocessing", "parameters", "reference",
                    "package_version") %in% names(m)))
  expect_match(m$method, "Peak area")
  expect_true(length(m$assumptions) >= 3L)
  expect_match(m$reference, "Salem")
  expect_match(m$reference, "Nam")
  expect_equal(m$package_version, as.character(utils::packageVersion("textileCrystR")))
  expect_equal(m$parameters$model, fit$model)
})
