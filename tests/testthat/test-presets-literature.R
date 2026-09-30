# Validation component 3: literature context. The Nam et al. (2016) and Salem et al.
# (2023) papers define the equation and the reference positions; they are NOT used
# as targets for the experimental values.

test_that("presets store every required field and are documented", {
  p <- segal_presets()
  expect_s3_class(p, "tbl_df")
  expect_true(all(c("preset", "material", "polymorph", "i200", "iam", "i200_range",
                    "iam_range", "radiation", "reference", "notes", "limitations") %in% names(p)))
  expect_true(all(c("cellulose_I_segal", "cellulose_II_extended", "cotton_validation") %in% p$preset))
  for (col in c("material", "i200_range", "iam_range", "reference", "notes", "limitations")) {
    expect_true(all(nzchar(p[[col]]) & !is.na(p[[col]])), info = col)
  }
  expect_true(all(is.finite(p$i200) & is.finite(p$iam) & p$i200 > p$iam))
})

test_that("preset positions match the supplied literature", {
  # Nam et al. (2016): I200 at 22.7 deg and Iam at 18 deg for cellulose I; (020) at
  # 21.7 deg and Iam at 16 deg for cellulose II. Salem et al. (2023): (200) maximum
  # between 22 and 23 deg, minimum near 18 deg.
  one <- segal_presets("cellulose_I_segal")
  expect_equal(c(one$i200, one$iam), c(22.7, 18.0))
  expect_match(one$i200_range, "22 and 23")
  expect_match(one$iam_range, "18\\.6")
  expect_match(one$reference, "Segal")
  expect_match(one$reference, "Nam")
  expect_match(one$reference, "Salem")
  two <- segal_presets("cellulose_II_extended")
  expect_equal(c(two$i200, two$iam), c(21.7, 16.0))
  expect_match(two$limitations, "partially mercerized")
  val <- segal_presets("cotton_validation")
  expect_equal(c(val$i200, val$iam), c(22.6, 18.0))
  expect_match(val$reference, "Fig. 8")
  expect_match(val$limitations, "not the crystalline fraction")
  expect_error(segal_presets("unknown"), "Unknown preset")
  expect_error(segal_presets(3), "Unknown preset")
})

test_that("the implemented equation is the Segal equation of both references", {
  # Salem et al. (2023) Eq. 1: CrI = 100 x (I200 - IAM) / I200
  # Nam et al. (2016) Eq. 2:   CI = (It - Ia) / It x 100
  tt <- c(10, 18, 22.7, 30)
  y <- c(50, 420, 1350, 80)
  res <- segal_ci(tt, y, preset = "cellulose_I_segal")
  expect_equal(res$ci, 100 * (1350 - 420) / 1350)
  expect_equal(res$ci, (1350 - 420) / 1350 * 100)
})

test_that("the reference-position exercise is independent of the experimental values", {
  # Nam et al. report Segal CIs of 84.7 % (control) and 68.8 % (mercerized) for their own
  # simulated/experimental patterns. Those numbers are context, and the package's
  # cotton results are deliberately not compared against them here.
  x <- read_xrd_pairs(cotton_path())
  d <- x[x$sample == "0%/1", ]
  ci <- segal_ci(d$two_theta, d$intensity, preset = "cotton_validation")$ci
  expect_false(isTRUE(all.equal(ci, 84.7, tolerance = 0.01)))
})
