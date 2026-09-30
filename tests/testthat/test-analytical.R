# Validation component 1: analytical ground truth. These are software-unit tests of
# the Segal arithmetic on patterns whose answer is known exactly.

test_that("Segal arithmetic gives the exact analytical results", {
  cases <- data.frame(i200 = c(1000, 850, 1250, 500, 1000, 640),
                      iam = c(250, 340, 100, 500, 0, 160),
                      expected = c(75, 60, 92, 0, 100, 75))
  for (k in seq_len(nrow(cases))) {
    d <- mini_pattern(cases$i200[k], cases$iam[k])
    res <- segal_ci(d$two_theta, d$intensity, i200 = 22.6, iam = 18)
    expect_equal(res$ci, cases$expected[k], tolerance = 1e-12)
  }
})

test_that("the result is unchanged when intensities are scaled", {
  d <- mini_pattern(1000, 250)
  base <- segal_ci(d$two_theta, d$intensity, i200 = 22.6, iam = 18)$ci
  for (f in c(0.001, 7.3, 1e6)) {
    expect_equal(segal_ci(d$two_theta, d$intensity * f, i200 = 22.6, iam = 18)$ci, base,
                 tolerance = 1e-12)
  }
})

test_that("a constant offset changes the index exactly as the equation predicts", {
  d <- mini_pattern(1000, 250)
  res <- segal_ci(d$two_theta, d$intensity + 100, i200 = 22.6, iam = 18)
  expect_equal(res$ci, (1100 - 350) / 1100 * 100, tolerance = 1e-12)
})

test_that("indices outside 0-100 % follow the equation and are flagged", {
  d <- mini_pattern(1000, 1200)
  expect_warning(res <- segal_ci(d$two_theta, d$intensity, i200 = 22.6, iam = 18), NA)
  expect_equal(res$ci, -20)
  expect_true(any(grepl("negative", res$warnings)))
  d <- mini_pattern(1000, -50)
  res <- segal_ci(d$two_theta, d$intensity, i200 = 22.6, iam = 18)
  expect_equal(res$ci, 105)
  expect_true(any(grepl("negative", res$warnings)))
})

test_that("a non-positive I200 is an error", {
  d <- mini_pattern(0, 10)
  expect_error(segal_ci(d$two_theta, d$intensity, i200 = 22.6, iam = 18), "not positive")
  d <- mini_pattern(-5, 10)
  expect_error(segal_ci(d$two_theta, d$intensity, i200 = 22.6, iam = 18), "not positive")
})

test_that("interpolation is exact for linear data and respects the method", {
  tt <- 5:40
  y <- 2 * tt + 10                       # linear, so linear interpolation is exact
  res <- segal_ci(tt, y, i200 = 22.6, iam = 18.25)
  expect_equal(res$points$intensity, c(2 * 22.6 + 10, 2 * 18.25 + 10), tolerance = 1e-12)
  expect_equal(res$ci, (55.2 - 46.5) / 55.2 * 100, tolerance = 1e-12)
  expect_false(any(res$points$on_measured_point))
  expect_equal(res$points$lower_2theta, c(22, 18))
  expect_equal(res$points$upper_2theta, c(23, 19))

  near <- segal_ci(tt, y, i200 = 22.6, iam = 18.25, interpolation = "nearest")
  expect_equal(near$points$position_used, c(23, 18))
  expect_equal(near$points$intensity, c(56, 46))
  expect_equal(near$points$position_requested, c(22.6, 18.25))
})

test_that("a position on a measured point uses the measured intensity exactly", {
  d <- mini_pattern(1000, 250)
  p <- segal_points(d$two_theta, d$intensity, i200 = 22.6 + 5e-7, iam = 18)
  expect_equal(p$i200_intensity, 1000)
  expect_true(p$i200_on_measured_point && p$iam_on_measured_point)
})
