test_that(".faddeeva matches the exact w(iy) = exp(y^2)*erfc(y) relation for real y", {
  for (y in c(0.001, 0.05, 0.3, 1, 2, 5, 12, 20)) {
    got <- Re(.faddeeva(complex(real = 0, imaginary = y)))
    exact <- exp(y^2) * 2 * pnorm(-y * sqrt(2))
    expect_equal(got, exact, tolerance = 1e-5)
  }
})

test_that(".faddeeva(0) = 1", {
  expect_equal(Re(.faddeeva(complex(real = 0, imaginary = 0))), 1, tolerance = 1e-6)
})

test_that("peak shapes equal height at their own center", {
  expect_equal(.gaussian_peak(5, 5, 123, 2), 123)
  expect_equal(.lorentzian_peak(5, 5, 123, 2), 123)
  expect_equal(.pseudo_voigt_peak(5, 5, 123, 2, 0.4), 123)
  expect_equal(.voigt_peak(5, 5, 123, 2, 1), 123, tolerance = 1e-9)
})

test_that("pseudo-Voigt is the eta-weighted mixture of Gaussian and Lorentzian", {
  x <- seq(-4, 4, by = 0.1)
  g <- .gaussian_peak(x, 0, 100, 1.5)
  l <- .lorentzian_peak(x, 0, 100, 1.5)
  for (eta in c(0, 0.3, 0.7, 1)) {
    expect_equal(.pseudo_voigt_peak(x, 0, 100, 1.5, eta), eta * l + (1 - eta) * g)
  }
})

test_that("Voigt reduces to Gaussian as the Lorentzian width shrinks to zero", {
  x <- seq(-5, 5, by = 0.05)
  g <- .gaussian_peak(x, 0, 100, 2)
  v <- .voigt_peak(x, 0, 100, 2, 1e-6)
  expect_lt(max(abs(g - v)), 1e-2)  # absolute; height scale is 100
})

test_that("Voigt reduces to Lorentzian as the Gaussian width shrinks to zero", {
  x <- seq(-5, 5, by = 0.05)
  l <- .lorentzian_peak(x, 0, 100, 2)
  v <- .voigt_peak(x, 0, 100, 1e-4, 2)
  expect_lt(max(abs(l - v)), 1e-2)
})

test_that("analytic area formulas match wide-window numerical integration", {
  x <- seq(-3000, 3000, by = 0.02)
  expect_equal(.trapz(x, .gaussian_peak(x, 0, 100, 1.5)), .gaussian_area(100, 1.5), tolerance = 1e-4)
  expect_equal(.trapz(x, .lorentzian_peak(x, 0, 100, 1.5)), .lorentzian_area(100, 1.5),
               tolerance = 5e-4)
  expect_equal(.trapz(x, .pseudo_voigt_peak(x, 0, 100, 1.5, 0.4)),
               .pseudo_voigt_area(100, 1.5, 0.4), tolerance = 5e-4)
  expect_equal(.trapz(x, .voigt_peak(x, 0, 100, 1.5, 0.8)), .voigt_area(100, 1.5, 0.8),
               tolerance = 5e-4)
})

test_that(".component_fwhm reports the stored FWHM directly except for Voigt", {
  p <- c(center = 20, height = 100, fwhm = 1.3)
  expect_equal(.component_fwhm(p, "gaussian"), 1.3)
  pv <- c(center = 20, height = 100, fwhm = 1.3, eta = 0.5)
  expect_equal(.component_fwhm(pv, "pseudo_voigt"), 1.3)
  pvo <- c(center = 20, height = 100, gaussian_fwhm = 0.8, lorentzian_fwhm = 0.5)
  fwhm_v <- .component_fwhm(pvo, "voigt")
  expect_equal(fwhm_v, 0.5346 * 0.5 + sqrt(0.2166 * 0.5^2 + 0.8^2))
  expect_gt(fwhm_v, max(0.8, 0.5))  # composite Voigt FWHM exceeds either component's own width
})

test_that(".trapz handles edge cases", {
  expect_equal(.trapz(1, 5), 0)
  expect_equal(.trapz(numeric(0), numeric(0)), 0)
  expect_equal(.trapz(c(0, 1, 2), c(0, 1, 0)), 1)  # triangle, base 2, height 1 -> area 1
})

test_that(".r_squared and .rmse behave sensibly", {
  y <- c(1, 2, 3, 4)
  expect_equal(.r_squared(y, y), 1)
  expect_equal(.rmse(y, y), 0)
  expect_equal(.r_squared(y, rep(mean(y), 4)), 0)
  expect_true(is.na(.r_squared(rep(5, 4), rep(5, 4))))  # zero total variance
  expect_equal(.rmse(c(0, 0), c(3, 4)), 3.5355339, tolerance = 1e-6)
})
