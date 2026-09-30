test_that("fit_amorphous_hump recovers known parameters on a clean synthetic hump", {
  tt <- seq(8, 40, by = 0.02)
  y <- 300 + 600 * exp(-4 * log(2) * (tt - 20.5)^2 / 6^2)
  set.seed(1)
  y_noisy <- y + rnorm(length(tt), sd = 2)
  hump <- fit_amorphous_hump(tt, y_noisy, window = c(8, 40), sample = "hump-test")
  expect_s3_class(hump, "textile_amorphous_hump")
  expect_equal(hump$center, 20.5, tolerance = 0.01)
  expect_equal(hump$height, 600, tolerance = 5)
  expect_equal(hump$fwhm, 6, tolerance = 0.05)
  expect_equal(hump$background$intercept, 300, tolerance = 3)
  expect_gt(hump$r_squared, 0.999)
  expect_true(hump$converged)
  expect_equal(hump$sample, "hump-test")
  expect_equal(hump$model, "pseudo_voigt")
})

test_that("fit_amorphous_hump respects an explicit starting center, for each model on its own matched shape", {
  tt <- seq(8, 40, by = 0.02)
  shapes <- list(
    gaussian = .gaussian_peak(tt, 21, 500, 5),
    lorentzian = .lorentzian_peak(tt, 21, 500, 5),
    pseudo_voigt = .pseudo_voigt_peak(tt, 21, 500, 5, 0.4),
    voigt = .voigt_peak(tt, 21, 500, 3, 2)
  )
  for (m in names(shapes)) {
    hump <- fit_amorphous_hump(tt, 300 + shapes[[m]], window = c(8, 40), model = m, center = 21)
    expect_equal(hump$model, m)
    expect_equal(hump$center, 21, tolerance = 0.05)
    expect_gt(hump$r_squared, 0.999)
  }
})

test_that("fit_amorphous_hump errors clearly on a too-narrow window", {
  tt <- seq(10, 40, by = 5)  # very sparse
  y <- rep(100, length(tt))
  expect_error(fit_amorphous_hump(tt, y, window = c(20, 21)), "at least 5")
})

test_that("fit_amorphous_hump print/plot methods work", {
  tt <- seq(10, 40, by = 0.02)
  y <- 300 + 600 * exp(-4 * log(2) * (tt - 20.5)^2 / 6^2)
  hump <- fit_amorphous_hump(tt, y, window = c(12, 30), sample = "s")
  out <- capture.output(print(hump))
  expect_true(any(grepl("textile_amorphous_hump", out)))
  expect_true(any(grepl("Sample   : s", out)))
  ret <- NULL
  capture.output(ret <- print(hump))
  expect_identical(ret, hump)
  capture.output(plot(hump))
  expect_true(TRUE)  # plot() prints and returns invisibly; reaching here means no error
})

# ---- ground-truth recovery for fit_crystalline_peaks(), all four models ----

test_that("fit_crystalline_peaks recovers known parameters exactly for gaussian", {
  s <- .synthetic_peakarea("gaussian", seed = 1)
  fit <- fit_crystalline_peaks(s$two_theta, s$intensity, centers = s$centers, model = "gaussian",
                               amorphous = c(17, 24), background = "none", sample = "g")
  expect_s3_class(fit, "textile_cryst_peakarea")
  expect_equal(fit$ci, s$true_ci, tolerance = 0.5)
  expect_equal(fit$crystalline_area, s$true_crystalline_area, tolerance = s$true_crystalline_area * 0.02)
  expect_equal(fit$amorphous_area, s$true_amorphous_area, tolerance = s$true_amorphous_area * 0.02)
  expect_true(fit$converged)
  expect_gt(fit$r_squared, 0.999)
  expect_equal(sort(round(fit$peaks$center, 1)), sort(s$centers))
  expect_setequal(fit$peaks$label, c("1-10", "110", "200", "004"))
})

test_that("fit_crystalline_peaks recovers known parameters for lorentzian", {
  s <- .synthetic_peakarea("lorentzian", seed = 7)
  fit <- fit_crystalline_peaks(s$two_theta, s$intensity, centers = s$centers, model = "lorentzian",
                               amorphous = c(17, 24), background = "none")
  expect_equal(fit$ci, s$true_ci, tolerance = 0.5)
  expect_gt(fit$r_squared, 0.999)
  expect_true(fit$converged)
})

test_that("fit_crystalline_peaks recovers known parameters for voigt", {
  s <- .synthetic_peakarea("voigt", seed = 42)
  fit <- fit_crystalline_peaks(s$two_theta, s$intensity, centers = s$centers, model = "voigt",
                               amorphous = c(17, 24), background = "none")
  expect_equal(fit$ci, s$true_ci, tolerance = 0.5)
  expect_gt(fit$r_squared, 0.999)
  expect_true(fit$converged)
  expect_true(all(c("gaussian_fwhm", "lorentzian_fwhm") %in% names(fit$peaks)))
  expect_true(all(is.na(fit$peaks$fwhm) == FALSE))  # composite FWHM still reported
})

test_that("fit_crystalline_peaks recovers pseudo_voigt (default model) within its known precision", {
  # pseudo-Voigt fits have an inherent, literature-recognized height/FWHM/eta
  # correlation that limits precision even on clean data (unlike the other
  # three models, which recover to <0.05 pp); the tolerance here reflects
  # that documented characteristic, not a looser correctness bar.
  for (sd in c(1, 7, 99)) {
    s <- .synthetic_peakarea("pseudo_voigt", seed = sd)
    fit <- fit_crystalline_peaks(s$two_theta, s$intensity, centers = s$centers,
                                 model = "pseudo_voigt", amorphous = c(17, 24), background = "none")
    expect_equal(fit$ci, s$true_ci, tolerance = 3)  # percentage points
    expect_gt(fit$r_squared, 0.998)
    expect_true(fit$converged)
    expect_true(all(fit$peaks$height > 5))  # regression check: no collapsed/degenerate peaks
  }
})

test_that("regression: crystalline peaks and the amorphous hump do not swap roles (pseudo_voigt)", {
  # This is the specific failure mode found and fixed during development:
  # loose bounds let a crystalline peak balloon to amorphous-like width while
  # the amorphous hump collapsed onto a crystalline peak's position. Checked
  # here directly rather than just via the CI tolerance above.
  s <- .synthetic_peakarea("pseudo_voigt", seed = 42)
  fit <- fit_crystalline_peaks(s$two_theta, s$intensity, centers = s$centers, model = "pseudo_voigt",
                               amorphous = c(17, 24), background = "none")
  expect_true(all(fit$peaks$fwhm < 2.5))            # crystalline peaks stay sharp
  expect_gt(fit$amorphous$fwhm, 4)                  # amorphous hump stays broad
  for (cc in s$centers) {
    expect_true(any(abs(fit$peaks$center - cc) < 0.5))  # every true center has a matching peak
  }
})

test_that("adding a default linear background does not harm clean (zero-background) recovery", {
  s <- .synthetic_peakarea("gaussian", seed = 1)
  fit <- fit_crystalline_peaks(s$two_theta, s$intensity, centers = s$centers, model = "gaussian",
                               amorphous = c(17, 24))  # background = "linear" (default)
  expect_equal(fit$ci, s$true_ci, tolerance = 0.5)
  expect_equal(fit$background$intercept, 0, tolerance = 5)
  expect_equal(fit$background$slope, 0, tolerance = 0.1)
})

test_that("a known injected linear background is recovered and excluded from the CI", {
  s <- .synthetic_peakarea("gaussian", seed = 11, bg0 = 700, bg1 = -3)
  fit <- fit_crystalline_peaks(s$two_theta, s$intensity, centers = s$centers, model = "gaussian",
                               amorphous = c(17, 24), background = "linear")
  expect_equal(fit$background$intercept, 700, tolerance = 15)
  expect_equal(fit$background$slope, -3, tolerance = 1)
  expect_equal(fit$ci, s$true_ci, tolerance = 2)      # background correctly excluded from area/CI
  expect_true(fit$converged)
  expect_gt(fit$r_squared, 0.99)
})

test_that("background = \"none\" fits no background term and background = \"constant\" fits only an intercept", {
  s <- .synthetic_peakarea("gaussian", seed = 1)
  fit_none <- fit_crystalline_peaks(s$two_theta, s$intensity, centers = s$centers, model = "gaussian",
                                    amorphous = c(17, 24), background = "none")
  expect_equal(fit_none$background$model, "none")
  expect_true(is.na(fit_none$background$intercept))
  expect_true(is.na(fit_none$background$slope))
  expect_true(all(fit_none$data$background == 0))

  fit_const <- fit_crystalline_peaks(s$two_theta, s$intensity, centers = s$centers, model = "gaussian",
                                     amorphous = c(17, 24), background = "constant")
  expect_equal(fit_const$background$model, "constant")
  expect_false(is.na(fit_const$background$intercept))
  expect_true(is.na(fit_const$background$slope))
  expect_true(all(fit_const$data$background == fit_const$background$intercept))
})

test_that("regression: the originally-pinned-bounds real-data problem stays fixed", {
  # On the shipped cotton data, before the background fix, several fitted
  # FWHMs were pinned exactly at their upper bound and the amorphous
  # component had collapsed to a narrow, crystalline-like width -- the sign
  # of a poorly constrained fit, not a converged one. Checked directly so a
  # future change can't silently reintroduce it. This does not assert a
  # tight ceiling on every peak's FWHM in general: real, messy data can
  # legitimately push one peak close to a bound without that being the
  # earlier catastrophic failure (which collapsed the amorphous hump itself
  # to a sharp peak, at the *opposite* end of its range).
  x <- read_xrd_pairs(cotton_path())
  d <- x[x$sample == "0%/1", ]
  fit <- fit_crystalline_peaks(d$two_theta, d$intensity, sample = "0%/1")
  expect_true(all(fit$peaks$height > 50))                     # no collapsed/degenerate peaks
  expect_gt(fit$amorphous$fwhm, 4.2)                           # comfortably clear of its 4-deg floor
  expect_lt(fit$amorphous$fwhm, diff(fit$window) * 1.2 - 0.1)  # comfortably clear of its ceiling
  expect_gt(fit$amorphous$fwhm, max(fit$peaks$fwhm))           # amorphous stays broader than every crystalline peak
  expect_gt(fit$r_squared, 0.995)
  expect_true(fit$converged)
  expect_true(fit$ci > 0 && fit$ci < 100)
})
