# Tests for Method 2A: reusable XRD pattern preprocessing. These functions
# are deliberately independent of validate_crystallinity_input() (see
# R/preprocessing.R); test-analytical.R, test-validation.R etc. continue to
# cover that existing function on its own.

# ---- validate_xrd_pattern() ----

test_that("validate_xrd_pattern accepts a clean pattern and reports counts", {
  v <- validate_xrd_pattern(c(10, 15, 20), c(1, 2, 3))
  expect_equal(v$two_theta, c(10, 15, 20))
  expect_equal(v$n_input, 3L)
  expect_equal(v$n_used, 3L)
  expect_equal(v$n_dropped, 0L)
})

test_that("validate_xrd_pattern does not sort or deduplicate", {
  v <- validate_xrd_pattern(c(15, 10, 15), c(1, 2, 3))
  expect_equal(v$two_theta, c(15, 10, 15))  # unchanged order, duplicate kept
})

test_that("validate_xrd_pattern enforces type, length, count and range checks", {
  expect_error(validate_xrd_pattern("a", 1), "numeric vectors")
  expect_error(validate_xrd_pattern(1:3, 1:2), "same length")
  expect_error(validate_xrd_pattern(10, 1), "At least two")
  expect_error(validate_xrd_pattern(c(0, 10), c(1, 2)), "between 0 and 180")
  expect_error(validate_xrd_pattern(c(10, 180), c(1, 2)), "between 0 and 180")
  expect_error(validate_xrd_pattern(c(10, NA), c(1, 2)), "missing or non-finite")
  expect_error(validate_xrd_pattern(c(10, 20), c(1, 2), na.rm = "yes"), "TRUE or FALSE")
})

test_that("validate_xrd_pattern drops non-finite rows only when na.rm = TRUE", {
  v <- validate_xrd_pattern(c(10, NA, 20, Inf), c(1, 2, 3, 4), na.rm = TRUE)
  expect_equal(v$two_theta, c(10, 20))
  expect_equal(v$n_input, 4L)
  expect_equal(v$n_dropped, 2L)
})

# ---- sort_xrd_pattern() ----

test_that("sort_xrd_pattern sorts ascending and reports was_sorted correctly", {
  s <- sort_xrd_pattern(c(15, 10, 20), c(2, 1, 3))
  expect_equal(s$two_theta, c(10, 15, 20))
  expect_equal(s$intensity, c(1, 2, 3))
  expect_true(s$was_sorted)

  already <- sort_xrd_pattern(c(10, 15, 20), c(1, 2, 3))
  expect_false(already$was_sorted)
  expect_equal(already$two_theta, c(10, 15, 20))
})

test_that("sort_xrd_pattern validates type and length", {
  expect_error(sort_xrd_pattern("a", 1), "numeric vectors")
  expect_error(sort_xrd_pattern(1:3, 1:2), "same length")
})

# ---- remove_xrd_duplicates() ----

test_that("remove_xrd_duplicates averages repeated 2-theta values regardless of input order", {
  d <- remove_xrd_duplicates(c(10, 12, 10), c(4, 9, 6), duplicates = "mean")
  expect_equal(d$two_theta, c(10, 12))
  expect_equal(d$intensity, c(5, 9))  # mean(4,6) = 5
  expect_equal(d$n_duplicates, 1L)
})

test_that("remove_xrd_duplicates errors by default and passes clean data through unchanged", {
  expect_error(remove_xrd_duplicates(c(10, 12, 10), c(4, 9, 6)), "duplicated")
  clean <- remove_xrd_duplicates(c(10, 12, 14), c(4, 9, 6))
  expect_equal(clean$n_duplicates, 0L)
  expect_equal(clean$two_theta, c(10, 12, 14))
})

test_that("remove_xrd_duplicates validates type and length", {
  expect_error(remove_xrd_duplicates("a", 1), "numeric vectors")
  expect_error(remove_xrd_duplicates(1:3, 1:2), "same length")
})

# ---- estimate_xrd_background() / subtract_xrd_background(): ground truth ----

test_that("estimate_xrd_background recovers a flat background under a sharp peak almost exactly", {
  tt <- seq(10, 40, by = 0.02)
  y <- 300 + 4000 * exp(-((tt - 22.6) / 0.7)^2)
  bg <- estimate_xrd_background(tt, y, window = 8)
  expect_lt(max(abs(bg$background - 300)), 1e-6)
})

test_that("estimate_xrd_background recovers a sloped background under multiple sharp peaks", {
  tt <- seq(10, 40, by = 0.02)
  true_bg <- 700 + 5 * (tt - 25)
  y <- true_bg + 1800 * exp(-((tt - 22.6) / 0.7)^2) + 600 * exp(-((tt - 16.3) / 0.8)^2) +
    400 * exp(-((tt - 14.8) / 0.7)^2) + 300 * exp(-((tt - 34.5) / 1)^2)
  bg <- estimate_xrd_background(tt, y, window = 8)
  # The LLS transform (needed for widely-varying peak heights) has a small,
  # known bias on sloped backgrounds -- see the regression test below and
  # ?estimate_xrd_background. For a slope typical of real data (matching
  # this test and the shipped cotton.csv), that bias is under 1% of the
  # background level (~700-850 here); this checks well clear of that, not
  # an unbiased/exact recovery.
  expect_lt(max(abs(bg$background - true_bg)), 10)
  expect_lt(mean(abs(bg$background - true_bg)), 5)
})

test_that("regression: SNIP background estimate is not badly biased by a sloped background alone", {
  # Found and fixed during development: clamping the window index at the
  # array edges (rather than only updating points with a full symmetric
  # window available) polluted a large fraction of the array once
  # `iterations` was a sizeable fraction of the pattern length -- which is
  # the common case here, not a rare edge condition. Checked directly on a
  # pure linear ramp with *zero* peaks, where any bias can only come from
  # the algorithm itself, not from genuine peak content. Before the fix,
  # this specific case had max error ~150 (about 21% of the background
  # level); the LLS transform itself contributes a separate, much smaller,
  # documented bias on sloped backgrounds (~4-5 here) that remains after the
  # fix, so the bound below is set well clear of the *old bug's* magnitude
  # rather than at zero.
  tt <- seq(10, 40, by = 0.02)
  true_bg <- 700 + 5 * (tt - 25)
  bg <- estimate_xrd_background(tt, true_bg, window = 8)
  expect_lt(max(abs(bg$background - true_bg)), 20)
})

test_that("estimate_xrd_background warns and falls back when window exceeds what the pattern supports", {
  tt <- seq(15, 25, by = 0.5)
  y <- 300 + 2000 * exp(-((tt - 20) / 1)^2)
  expect_warning(bg <- estimate_xrd_background(tt, y, window = 8), "using the largest window")
  expect_lt(bg$iterations, round(8 / 0.5))
})

test_that("estimate_xrd_background errors clearly on a too-short pattern", {
  expect_error(estimate_xrd_background(c(10, 11), c(1, 2)), "too few points")
})

test_that("estimate_xrd_background errors clearly when the 2-theta step can't be determined", {
  # heavily duplicated 2-theta values can degenerate the median step to zero
  expect_error(estimate_xrd_background(c(10, 10, 10, 10, 20), c(1, 2, 3, 4, 5)),
              "Could not determine a typical 2-theta step")
})

test_that(".snip_background's own safety break is reachable directly (belt-and-suspenders check)", {
  # estimate_xrd_background() always pre-caps `iterations` so this branch is
  # not reached via the public function; checked directly since it is not
  # otherwise exercised, and future refactoring could change that guarantee.
  y <- c(10, 50, 12, 8, 60, 9, 11)
  expect_no_error(out <- .snip_background(y, iterations = 10L))  # far more than length(y) supports
  expect_length(out, length(y))
})

test_that("estimate_xrd_background validates `window` and sorts unsorted input first", {
  expect_error(estimate_xrd_background(c(10, 20), c(1, 2), window = -1), "positive number")
  expect_error(estimate_xrd_background(c(10, 20), c(1, 2), window = "a"), "positive number")

  tt <- seq(10, 40, by = 0.02)
  y <- 300 + 4000 * exp(-((tt - 22.6) / 0.7)^2)
  o <- sample(length(tt))
  bg_shuffled <- estimate_xrd_background(tt[o], y[o], window = 8)
  expect_false(is.unsorted(bg_shuffled$two_theta))
  expect_lt(max(abs(bg_shuffled$background - 300)), 1e-6)
})

test_that("subtract_xrd_background accepts either the estimate_xrd_background() list or a plain vector", {
  tt <- seq(10, 40, by = 0.02)
  y <- 300 + 4000 * exp(-((tt - 22.6) / 0.7)^2)
  bg <- estimate_xrd_background(tt, y, window = 8)

  d1 <- subtract_xrd_background(tt, y, bg)
  d2 <- subtract_xrd_background(tt, y, bg$background)
  expect_equal(d1$corrected, d2$corrected)
  expect_equal(d1$corrected, pmax(y - bg$background, 0))
})

test_that("subtract_xrd_background's clip_negative controls flooring at zero", {
  tt <- c(10, 15, 20)
  y <- c(50, 40, 60)
  bg <- c(45, 45, 45)  # 15 -> negative before clipping
  clipped <- subtract_xrd_background(tt, y, bg, clip_negative = TRUE)
  unclipped <- subtract_xrd_background(tt, y, bg, clip_negative = FALSE)
  expect_equal(clipped$corrected, c(5, 0, 15))
  expect_equal(unclipped$corrected, c(5, -5, 15))
})

test_that("subtract_xrd_background validates the background vector's length", {
  expect_error(subtract_xrd_background(c(10, 20, 30), c(1, 2, 3), c(1, 2)), "length 3")
})

# ---- normalize_xrd_pattern() ----

test_that("normalize_xrd_pattern: max method scales the tallest point to 1", {
  n <- normalize_xrd_pattern(c(10, 20, 30), c(50, 200, 100), method = "max")
  expect_equal(max(n$intensity), 1)
  expect_equal(n$intensity, c(50, 200, 100) / 200)
  expect_equal(n$scale, 200)
})

test_that("normalize_xrd_pattern: area method integrates to 1", {
  tt <- seq(1, 11, by = 0.01)
  y <- rep(5, length(tt))  # constant 5 over a 10-degree span -> area 50
  n <- normalize_xrd_pattern(tt, y, method = "area")
  expect_equal(.trapz(n$two_theta, n$intensity), 1, tolerance = 1e-6)
  expect_equal(n$scale, 50, tolerance = 1e-6)
})

test_that("normalize_xrd_pattern: minmax method rescales to [0, 1]", {
  n <- normalize_xrd_pattern(c(10, 20, 30), c(50, 200, 100), method = "minmax")
  expect_equal(min(n$intensity), 0)
  expect_equal(max(n$intensity), 1)
  expect_true(is.na(n$scale))
})

test_that("normalize_xrd_pattern errors clearly on degenerate (constant/zero) intensity", {
  expect_error(normalize_xrd_pattern(c(10, 20), c(0, 0), method = "max"), "zero")
  expect_error(normalize_xrd_pattern(c(10, 20), c(5, 5), method = "minmax"), "constant")
  expect_error(normalize_xrd_pattern(c(10, 20), c(0, 0), method = "area"), "area is zero")
})

test_that("normalize_xrd_pattern rejects an unknown method", {
  expect_error(normalize_xrd_pattern(c(10, 20), c(1, 2), method = "bogus"), "should be one of")
})

# ---- real cotton data smoke test ----

test_that("background estimation and subtraction give plausible results on real cotton data", {
  x <- read_xrd_pairs(cotton_path())
  d <- x[x$sample == "0%/1", ]
  bg <- estimate_xrd_background(d$two_theta, d$intensity, window = 8)
  expect_true(all(bg$background >= 0))
  expect_true(all(bg$background <= bg$intensity + 1e-6))  # background never exceeds the raw pattern
  expect_true(max(bg$background) < max(bg$intensity))     # doesn't eat into the tallest peak

  sub <- subtract_xrd_background(d$two_theta, d$intensity, bg)
  expect_true(all(sub$corrected >= 0))
  expect_gt(max(sub$corrected), 0.9 * max(d$intensity))   # the crystalline peak survives, largely intact
})
