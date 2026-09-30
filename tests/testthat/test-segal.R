test_that("segal_ci returns a documented textile_cryst object", {
  d <- synthetic_xrd()
  res <- segal_ci(d$two_theta, d$intensity, i200 = 22.6, iam = 18, sample = "syn")
  expect_s3_class(res, "textile_cryst")
  expect_named(res, c("ci", "points", "method", "sample", "data", "input", "warnings"))
  expect_true(is.numeric(res$ci) && length(res$ci) == 1L)
  expect_equal(res$sample, "syn")
  expect_equal(res$points$point, c("I200", "Iam"))
  expect_equal(res$points$position_requested, c(22.6, 18))
  expect_equal(res$method$interpolation, "linear")
  expect_equal(res$method$name, "Segal peak height")
  expect_match(res$method$equation, "I200")
  expect_equal(nrow(res$data), nrow(d))
  expect_equal(res$data$intensity, d$intensity)          # measured data retained unchanged
  expect_equal(res$input$n_used, nrow(d))
  expect_equal(res$ci, (res$points$intensity[1] - res$points$intensity[2]) /
                 res$points$intensity[1] * 100)
})

test_that("unnamed samples are allowed and sample is checked", {
  d <- synthetic_xrd()
  res <- segal_ci(d$two_theta, d$intensity, i200 = 22.6, iam = 18)
  expect_true(is.na(res$sample))
  expect_error(segal_ci(d$two_theta, d$intensity, 22.6, 18, sample = 1), "sample")
  expect_error(segal_ci(d$two_theta, d$intensity, 22.6, 18, sample = c("a", "b")), "sample")
})

test_that("positions must be supplied and valid", {
  d <- synthetic_xrd()
  f <- function(...) segal_ci(d$two_theta, d$intensity, ...)
  expect_error(f(), "preset")
  expect_error(f(i200 = 22.6), "preset")
  expect_error(f(iam = 18), "preset")
  expect_error(f(i200 = NA_real_, iam = 18), "i200")
  expect_error(f(i200 = "a", iam = 18), "i200")
  expect_error(f(i200 = c(22, 23), iam = 18), "i200")
  expect_error(f(i200 = -1, iam = 18), "i200")
  expect_error(f(i200 = 22.6, iam = 200), "iam")
  expect_error(f(i200 = 22.6, iam = Inf), "iam")
  expect_error(f(i200 = 18, iam = 18), "different")
  expect_error(f(i200 = 50, iam = 18), "outside the measured")
  expect_error(f(i200 = 22.6, iam = 2), "outside the measured")
  expect_error(f(i200 = 22.6, iam = 18, interpolation = "cubic"), "should be one of")
})

test_that("swapped positions give a warning", {
  d <- synthetic_xrd()
  res <- segal_ci(d$two_theta, d$intensity, i200 = 18, iam = 22.6)
  expect_true(any(grepl("lower angle", res$warnings)))
})

test_that("interpolating across a wide gap warns", {
  tt <- c(5, 10, 15, 20, 25, 30, 40)
  y <- c(100, 200, 300, 500, 900, 400, 100)
  res <- segal_ci(tt, y, i200 = 22.6, iam = 12.5)
  expect_true(any(grepl("gap", res$warnings)))
  res <- segal_ci(tt, y, i200 = 22.6, iam = 12.5, interpolation = "nearest")
  expect_true(any(grepl("gap", res$warnings)))
})

test_that("presets supply positions and explicit positions override them", {
  d <- synthetic_xrd()
  a <- segal_ci(d$two_theta, d$intensity, preset = "cellulose_I_segal")
  b <- segal_ci(d$two_theta, d$intensity, i200 = 22.7, iam = 18)
  expect_equal(a$ci, b$ci)
  expect_equal(a$method$preset, "cellulose_I_segal")
  expect_false(a$method$preset_overridden)
  expect_match(a$method$preset_reference, "Segal")
  expect_length(a$warnings, 0L)

  o <- segal_ci(d$two_theta, d$intensity, preset = "cellulose_I_segal", i200 = 22.6)
  expect_true(o$method$preset_overridden)
  expect_equal(o$points$position_requested, c(22.6, 18))
  expect_true(any(grepl("differ from preset", o$warnings)))

  same <- segal_ci(d$two_theta, d$intensity, preset = "cellulose_I_segal", i200 = 22.7)
  expect_false(same$method$preset_overridden)
  half <- segal_ci(d$two_theta, d$intensity, preset = "cellulose_I_segal", iam = 17.9)
  expect_true(half$method$preset_overridden)

  expect_error(segal_ci(d$two_theta, d$intensity, preset = "nope"), "Unknown preset")
})

test_that("segal_points returns positions, intensities, method and input information", {
  d <- synthetic_xrd()
  p <- segal_points(d$two_theta, d$intensity, i200 = 22.62, iam = 18.013, sample = "s")
  expect_s3_class(p, "tbl_df")
  expect_equal(nrow(p), 1L)
  expect_true(all(c("i200_position", "i200_intensity", "iam_position", "iam_intensity",
                    "interpolation", "n_input", "n_used", "two_theta_min", "two_theta_max",
                    "median_step") %in% names(p)))
  expect_equal(p$interpolation, "linear")
  expect_equal(p$sample, "s")
  expect_false(p$iam_on_measured_point)
  expect_equal(p$iam_lower, 18)
  expect_equal(p$iam_upper, 18.05)
  expect_equal(p$i200_intensity, stats::approx(d$two_theta, d$intensity, 22.62)$y)
  pn <- segal_points(d$two_theta, d$intensity, preset = "cotton_validation",
                     interpolation = "nearest")
  expect_equal(pn$preset, "cotton_validation")
  expect_equal(pn$i200_position, 22.6, tolerance = 1e-6)
  # consistent with segal_ci
  res <- segal_ci(d$two_theta, d$intensity, i200 = 22.62, iam = 18.013)
  expect_equal(res$points$intensity, c(p$i200_intensity, p$iam_intensity))
})

test_that("unsorted, duplicated and missing data are handled as documented", {
  d <- mini_pattern(1000, 250)
  o <- c(3, 1, 6, 2, 5, 4)
  res <- segal_ci(d$two_theta[o], d$intensity[o], i200 = 22.6, iam = 18)
  expect_equal(res$ci, 75)
  expect_true(res$input$was_sorted)
  expect_true(any(grepl("sorted", res$input$messages)))
  expect_equal(res$data$two_theta, d$two_theta)

  dd <- rbind(d, data.frame(two_theta = 22.6, intensity = 1000))
  expect_error(segal_ci(dd$two_theta, dd$intensity, i200 = 22.6, iam = 18), "duplicated")
  res <- segal_ci(dd$two_theta, dd$intensity, i200 = 22.6, iam = 18, duplicates = "mean")
  expect_equal(res$ci, 75)
  dd$intensity[nrow(dd)] <- 800
  res <- segal_ci(dd$two_theta, dd$intensity, i200 = 22.6, iam = 18, duplicates = "mean")
  expect_equal(res$points$intensity[1], 900)
  expect_equal(res$input$n_duplicates, 1L)

  dm <- d
  dm$intensity[2] <- NA
  expect_error(segal_ci(dm$two_theta, dm$intensity, i200 = 22.6, iam = 18), "missing")
  res <- segal_ci(dm$two_theta, dm$intensity, i200 = 22.6, iam = 18, na.rm = TRUE)
  expect_equal(res$ci, 75)
  expect_equal(res$input$n_dropped, 1L)
  expect_equal(res$input$n_used, 5L)
})

test_that("an interpolation neighbour that was dropped is reflected in the neighbours", {
  d <- mini_pattern(1000, 250)
  d$intensity[3] <- NA                       # the point at 18 degrees is missing
  res <- segal_ci(d$two_theta, d$intensity, i200 = 22.6, iam = 18, na.rm = TRUE)
  expect_false(res$points$on_measured_point[2])
  expect_equal(res$points$lower_2theta[2], 10)
  expect_equal(res$points$upper_2theta[2], 22.6)
})
