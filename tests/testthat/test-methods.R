make_res <- function(...) {
  d <- synthetic_xrd()
  segal_ci(d$two_theta, d$intensity, i200 = 22.6, iam = 18, sample = "syn", ...)
}

test_that("print() reports the Segal Crystallinity Index and never calls it absolute", {
  res <- make_res()
  out <- capture.output(print(res))
  expect_true(any(grepl("Segal Crystallinity Index", out)))
  expect_true(any(grepl("CI", out)))
  expect_true(any(grepl("I200", out)))
  expect_true(any(grepl("Iam", out)))
  expect_true(any(grepl("Sample : syn", out)))
  expect_false(any(grepl("absolute crystallinity", out, ignore.case = TRUE)))
  expect_identical(print(res), res)                       # returns invisibly
  expect_true(any(grepl("both positions on measured points", out)))
  d <- synthetic_xrd()
  off <- segal_ci(d$two_theta, d$intensity, i200 = 22.62, iam = 18)
  out2 <- capture.output(print(off))
  expect_false(any(grepl("both positions", out2)))
  expect_false(any(grepl("Sample :", out2)))
  bad <- segal_ci(d$two_theta, d$intensity, i200 = 18, iam = 22.6)
  expect_true(any(grepl("Warning", capture.output(print(bad)))))
})

test_that("summary() gives method, points, input checks and the interpretation caveat", {
  res <- make_res()
  s <- summary(res)
  expect_s3_class(s, "summary.textile_cryst")
  expect_equal(s$ci, res$ci)
  out <- capture.output(print(s))
  expect_true(any(grepl("Segal Crystallinity Index \\(summary\\)", out)))
  expect_true(any(grepl("CI = \\(I200 - Iam\\) / I200 \\* 100", out)))
  expect_true(any(grepl("Interpolation : linear", out)))
  expect_true(any(grepl("none \\(explicit positions\\)", out)))
  expect_true(any(grepl("empirical, relative", out)))
  expect_false(any(grepl("absolute crystallinity", out, ignore.case = TRUE)))

  d <- synthetic_xrd()
  o <- c(2:nrow(d), 1)
  pres <- segal_ci(d$two_theta[o], d$intensity[o], preset = "cellulose_I_segal",
                   i200 = 22.61)
  out <- capture.output(print(summary(pres)))
  expect_true(any(grepl("cellulose_I_segal", out)))
  expect_true(any(grepl("Input note", out)))
  expect_true(any(grepl("Warning", out)))
})

test_that("as.data.frame() returns one tidy row", {
  res <- make_res()
  df <- as.data.frame(res)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 1L)
  expect_named(df, c("sample", "ci", "i200_position", "i200_intensity", "iam_position",
                     "iam_intensity", "interpolation", "preset", "n_used", "n_warnings"))
  expect_equal(df$ci, res$ci)
  expect_true(is.na(df$preset))
  expect_equal(df$n_warnings, 0L)
  d <- synthetic_xrd()
  pres <- segal_ci(d$two_theta, d$intensity, preset = "cotton_validation")
  expect_equal(as.data.frame(pres)$preset, "cotton_validation")
  two <- rbind(df, as.data.frame(pres))
  expect_equal(nrow(two), 2L)
  expect_equal(rownames(as.data.frame(res, row.names = "a")), "a")
})

test_that("plot() returns a modifiable ggplot object showing the points and the index", {
  res <- make_res()
  p <- plot(res)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "syn")
  expect_match(p$labels$subtitle, "Segal Crystallinity Index = ")
  expect_match(p$labels$subtitle, sprintf("%.2f%%", res$ci), fixed = TRUE)
  expect_gte(length(p$layers), 4L)
  labs_layer <- p$layers[[4]]$data
  expect_equal(labs_layer$point, c("I200", "Iam"))
  expect_match(labs_layer$label[1], "I200 = ")
  expect_match(labs_layer$label[2], "Iam = ")

  # user modifications keep working and build
  p2 <- p + ggplot2::theme_bw() + ggplot2::labs(title = "changed") +
    ggplot2::geom_hline(yintercept = 1000)
  expect_s3_class(p2, "ggplot")
  expect_equal(p2$labels$title, "changed")
  expect_no_error(ggplot2::ggplot_build(p2))
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("plot() options: xlim, show_ci, title, and validation of xlim", {
  res <- make_res()
  p <- plot(res, xlim = c(10, 30), show_ci = FALSE, title = "custom")
  expect_null(p$labels$subtitle)
  expect_equal(p$labels$title, "custom")
  expect_equal(p$coordinates$limits$x, c(10, 30))
  expect_equal(nrow(p$data), nrow(res$data))               # full pattern is kept
  expect_error(plot(res, xlim = c(30, 10)), "xlim")
  expect_error(plot(res, xlim = 10), "xlim")
  expect_error(plot(res, xlim = c(NA, 10)), "xlim")
  expect_no_error(ggplot2::ggplot_build(p))

  d <- synthetic_xrd()
  unnamed <- segal_ci(d$two_theta, d$intensity, i200 = 22.6, iam = 18)
  expect_equal(plot(unnamed)$labels$title, "XRD pattern")
  auto <- plot(unnamed)$coordinates$limits$x
  expect_equal(auto, c(8, 37.6))
  # window outside the data falls back to all intensities for the y range
  expect_no_error(ggplot2::ggplot_build(plot(unnamed, xlim = c(100, 120))))
})

test_that("plot() works with nearest interpolation and negative intensities", {
  d <- synthetic_xrd()
  res <- segal_ci(d$two_theta, d$intensity - 400, i200 = 22.62, iam = 18.01,
                  interpolation = "nearest")
  expect_no_error(ggplot2::ggplot_build(plot(res)))
})

test_that("regression: plot() does not trigger geom_label's ggplot2 3.5.0 label.size deprecation warning", {
  # geom_label()'s border-width argument was renamed from label.size to
  # linewidth in ggplot2 3.5.0 (label.size soft-deprecated since); found on a
  # real user's newer ggplot2 install, where every plot() call warned. The
  # wrong argument name for the installed version either warns (old name on
  # new ggplot2) or is silently ignored, leaving a visible border (new name
  # on old ggplot2) -- so both directions are checked here via mocking,
  # rather than assuming which branch is "current" in this environment.
  res <- make_res()
  expect_no_warning(plot(res))

  testthat::with_mocked_bindings(
    packageVersion = function(pkg) as.package_version("3.5.1"),
    .package = "utils",
    expect_equal(.label_border_zero(), list(linewidth = 0))
  )
  testthat::with_mocked_bindings(
    packageVersion = function(pkg) as.package_version("3.4.0"),
    .package = "utils",
    expect_equal(.label_border_zero(), list(label.size = 0))
  )
})
