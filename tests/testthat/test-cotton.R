# Validation components 2 and 4: independent implementation check and the seven
# experimental cotton patterns.

cotton_ci <- function(i200 = 22.6, iam = 18) {
  x <- read_xrd_pairs(cotton_path())
  do.call(rbind, lapply(unique(x$sample), function(s) {
    d <- x[x$sample == s, ]
    as.data.frame(segal_ci(d$two_theta, d$intensity, i200 = i200, iam = iam, sample = s))
  }))
}

test_that("the shipped cotton.csv is the unmodified experimental file", {
  expect_equal(unname(tools::md5sum(cotton_path())), "cf33ee9c15cd55f1638c093ee28b804a")
})

test_that("all seven cotton patterns load and validate", {
  x <- read_xrd_pairs(cotton_path())
  expect_setequal(unique(x$sample), names(cotton_expected_2dp))
  expect_equal(nrow(x), 7L * 8701L)
  for (s in unique(x$sample)) {
    d <- x[x$sample == s, ]
    v <- validate_crystallinity_input(d$two_theta, d$intensity)
    expect_equal(v$n_used, 8701L)
    expect_equal(v$range, c(3, 90))
    expect_false(v$was_sorted)
    expect_length(v$messages, 0L)
    expect_equal(v$median_step, 0.01, tolerance = 1e-4)
  }
})

test_that("Segal CI for all seven patterns is calculated from the raw data", {
  res <- cotton_ci()
  expect_equal(nrow(res), 7L)
  expect_true(all(res$interpolation == "linear"))
  expect_equal(res$i200_position, rep(22.6, 7))
  expect_equal(res$iam_position, rep(18, 7))
  # approximately the values in the project brief (2 d.p.)
  got <- setNames(round(res$ci, 2), res$sample)
  expect_equal(got[names(cotton_expected_2dp)], cotton_expected_2dp)
})

test_that("independent (spreadsheet) reference agrees with the package result", {
  ref <- utils::read.csv(system.file("extdata", "reference_cotton_segal.csv",
                                     package = "textileCrystR", mustWork = TRUE))
  res <- cotton_ci()
  m <- merge(res, ref, by = "sample")
  expect_equal(nrow(m), 7L)
  diff <- abs(m$ci - m$ci_reference)
  expect_lt(max(diff), 1e-8)                              # percentage points
  expect_equal(m$i200_intensity.x, m$i200_intensity.y, tolerance = 1e-9)
  expect_equal(m$iam_intensity.x, m$iam_intensity.y, tolerance = 1e-9)
})

test_that("the 5 x 5 sensitivity grid agrees with the spreadsheet for three samples", {
  ref <- utils::read.csv(system.file("extdata", "reference_cotton_sensitivity.csv",
                                     package = "textileCrystR", mustWork = TRUE))
  expect_equal(nrow(ref), 75L)
  x <- read_xrd_pairs(cotton_path())
  for (s in unique(ref$sample)) {
    d <- x[x$sample == s, ]
    sens <- segal_sensitivity(d$two_theta, d$intensity, sample = s)
    r <- ref[ref$sample == s, ]
    m <- merge(as.data.frame(sens), r, by = c("i200_position", "iam_position"))
    expect_equal(nrow(m), 25L)
    expect_lt(max(abs(m$ci - m$ci_reference)), 1e-8)
  }
})

test_that("interpolation and nearest-point reading agree on the measured grid", {
  a <- cotton_ci()
  x <- read_xrd_pairs(cotton_path())
  d <- x[x$sample == "50%/1", ]
  n <- segal_ci(d$two_theta, d$intensity, i200 = 22.6, iam = 18, interpolation = "nearest")
  expect_equal(n$ci, a$ci[a$sample == "50%/1"])
  # an off-grid position is interpolated between neighbouring measured points
  off <- segal_ci(d$two_theta, d$intensity, i200 = 22.603, iam = 18)
  expect_false(off$points$on_measured_point[1])
  expect_gt(off$points$intensity[1], min(d$intensity[d$two_theta > 22.6 & d$two_theta < 22.62]) - 1)
})

test_that("gravimetric composition is not used by the calculation", {
  x <- read_xrd_pairs(cotton_path())
  d <- x[x$sample == "40%/1", ]
  a <- segal_ci(d$two_theta, d$intensity, i200 = 22.6, iam = 18, sample = "40%/1")$ci
  b <- segal_ci(d$two_theta, d$intensity, i200 = 22.6, iam = 18, sample = "anything")$ci
  expect_identical(a, b)
})
