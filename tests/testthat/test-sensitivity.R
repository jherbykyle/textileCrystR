test_that("the default grid gives the complete 5 x 5 matrix with summary statistics", {
  d <- synthetic_xrd()
  s <- segal_sensitivity(d$two_theta, d$intensity, sample = "syn")
  expect_s3_class(s, "textile_cryst_sensitivity")
  expect_equal(dim(s$matrix), c(5L, 5L))
  expect_equal(trimws(rownames(s$matrix)), c("17.6", "17.8", "18.0", "18.2", "18.4"))
  expect_equal(trimws(colnames(s$matrix)), c("22.2", "22.4", "22.6", "22.8", "23.0"))
  expect_equal(nrow(s$grid), 25L)
  expect_true(all(is.finite(s$matrix)))

  # every cell equals the single-point calculation
  i200 <- c(22.2, 22.4, 22.6, 22.8, 23.0)
  iam <- c(17.6, 17.8, 18.0, 18.2, 18.4)
  for (i in 1:5) for (j in 1:5) {
    one <- segal_ci(d$two_theta, d$intensity, i200 = i200[j], iam = iam[i])$ci
    expect_equal(unname(s$matrix[i, j]), one, tolerance = 1e-10)
  }
  m <- s$matrix
  expect_equal(s$summary$ci_min, min(m))
  expect_equal(s$summary$ci_max, max(m))
  expect_equal(s$summary$ci_range, max(m) - min(m))
  expect_equal(s$summary$ci_central, unname(m["18.0", "22.6"]))
  expect_equal(s$central$ci, s$summary$ci_central)
  expect_equal(c(s$central$i200, s$central$iam), c(22.6, 18))
  imin <- unname(which(m == min(m), arr.ind = TRUE)[1, ])
  imax <- unname(which(m == max(m), arr.ind = TRUE)[1, ])
  expect_equal(c(s$summary$i200_at_min, s$summary$iam_at_min),
               c(i200[imin[2]], iam[imin[1]]))
  expect_equal(c(s$summary$i200_at_max, s$summary$iam_at_max),
               c(i200[imax[2]], iam[imax[1]]))
  expect_match(s$note, "not a universal peak-selection rule")
})

test_that("custom grids, central estimate and interpolation are supported", {
  d <- synthetic_xrd()
  s <- segal_sensitivity(d$two_theta, d$intensity, i200 = c(22.5, 23), iam = c(17.5, 18, 18.5),
                         central = c(iam = 18.25, i200 = 22.7), interpolation = "nearest")
  expect_equal(dim(s$matrix), c(3L, 2L))
  expect_equal(s$interpolation, "nearest")
  expect_equal(c(s$central$i200, s$central$iam), c(22.7, 18.25))
  one <- segal_ci(d$two_theta, d$intensity, i200 = 22.7, iam = 18.25, interpolation = "nearest")
  expect_equal(s$central$ci, one$ci)
  single <- segal_sensitivity(d$two_theta, d$intensity, i200 = 22.6, iam = 18)
  expect_equal(dim(single$matrix), c(1L, 1L))
  expect_equal(single$summary$ci_range, 0)
})

test_that("invalid grids and central values are rejected", {
  d <- synthetic_xrd()
  f <- function(...) segal_sensitivity(d$two_theta, d$intensity, ...)
  expect_error(f(i200 = "a"), "numeric vector")
  expect_error(f(iam = numeric(0)), "numeric vector")
  expect_error(f(i200 = c(22.2, 22.2)), "repeated")
  expect_error(f(iam = c(18, 18)), "repeated")
  expect_error(f(i200 = c(22.2, NA)), "i200")
  expect_error(f(iam = 250), "iam")
  expect_error(f(i200 = c(17, 22.6)), "greater than every Iam")
  expect_error(f(central = c(22.6, 18)), "named numeric")
  expect_error(f(central = c(a = 22.6, b = 18)), "named numeric")
  expect_error(f(central = c(i200 = 22.6)), "named numeric")
  expect_error(f(central = "x"), "named numeric")
  expect_error(f(central = c(i200 = 22.6, iam = 400)), "iam")
  expect_error(f(central = c(i200 = 400, iam = 18)), "i200")
  expect_error(f(i200 = 60), "outside the measured")
  expect_error(f(sample = 5), "sample")
  expect_error(f(interpolation = "spline"), "should be one of")
})

test_that("non-positive I200 stops and problematic patterns warn", {
  tt <- seq(5, 40, by = 0.05)
  y <- rep(100, length(tt))
  y[abs(tt - 22.6) < 1e-6] <- 0
  expect_error(segal_sensitivity(tt, y), "not positive")

  y <- 300 + 4000 * exp(-((tt - 22.7) / 1.1)^2)
  y[abs(tt - 18) < 1e-6] <- -10
  s <- segal_sensitivity(tt, y)
  expect_true(any(grepl("negative", s$warnings)))

  sparse <- seq(5, 40, by = 2.5)
  ys <- 300 + 4000 * exp(-((sparse - 22.7) / 1.1)^2)
  s <- segal_sensitivity(sparse, ys)
  expect_true(any(grepl("gap", s$warnings)))
  s <- segal_sensitivity(sparse, ys, interpolation = "nearest")
  expect_true(any(grepl("gap", s$warnings)))
})

test_that("print() and as.data.frame() work for sensitivity results", {
  d <- synthetic_xrd()
  s <- segal_sensitivity(d$two_theta, d$intensity, sample = "syn")
  out <- capture.output(print(s))
  expect_true(any(grepl("Sample : syn", out)))
  expect_true(any(grepl("Min ", out)))
  expect_true(any(grepl("Central estimate", out)))
  expect_true(any(grepl("not a universal peak-selection rule", out)))
  expect_identical(print(s), s)
  df <- as.data.frame(s)
  expect_equal(nrow(df), 25L)
  expect_true(all(c("sample", "i200_position", "iam_position", "ci") %in% names(df)))
  expect_equal(unique(df$sample), "syn")
  expect_equal(rownames(as.data.frame(s, row.names = as.character(1:25)))[25], "25")

  unnamed <- segal_sensitivity(d$two_theta, d$intensity)
  expect_false(any(grepl("Sample :", capture.output(print(unnamed)))))
  sp <- seq(5, 40, by = 2.5)
  w <- segal_sensitivity(sp, 300 + 4000 * exp(-((sp - 22.7) / 1.1)^2))
  expect_true(any(grepl("Warning", capture.output(print(w)))))
})
