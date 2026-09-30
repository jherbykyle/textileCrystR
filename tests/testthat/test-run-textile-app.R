# Pure-function tests and run_textile_app() itself. The dashboard's reactive
# logic (file upload, sample pooling, method dispatch, single/batch modes)
# is tested in test-dashboard.R.

test_that(".clamp_preset_positions clamps a preset's positions to the data range", {
  # preset's published positions already inside the range: unchanged
  v <- .clamp_preset_positions(c(5, 40), "cellulose_I_segal")
  expect_equal(v, list(i200 = 22.7, iam = 18.0))

  # range narrower than the preset: both positions clamped to the nearest edge
  v2 <- .clamp_preset_positions(c(15, 20), "cellulose_I_segal")
  expect_equal(v2$i200, 20)
  expect_equal(v2$iam, 18)

  v3 <- .clamp_preset_positions(c(19, 25), "cellulose_I_segal")
  expect_equal(v3$i200, 22.7)
  expect_equal(v3$iam, 19)

  # cellulose_II_extended preset (16 / 21.7) clamped the other direction
  v4 <- .clamp_preset_positions(c(17, 20), "cellulose_II_extended")
  expect_equal(v4$i200, 20)
  expect_equal(v4$iam, 17)
})

test_that(".clamp_preset_positions falls back to cellulose_I_segal for an unknown preset", {
  v <- .clamp_preset_positions(c(5, 40), "not_a_real_preset")
  expect_equal(v, list(i200 = 22.7, iam = 18.0))
  v2 <- .clamp_preset_positions(c(5, 40), NULL)
  expect_equal(v2, list(i200 = 22.7, iam = 18.0))
})

test_that(".clamp_preset_positions validates its range argument", {
  expect_error(.clamp_preset_positions(c(5, 5), "cellulose_I_segal"), "two increasing")
  expect_error(.clamp_preset_positions(c(10, 5), "cellulose_I_segal"), "two increasing")
  expect_error(.clamp_preset_positions(c(NA, 5), "cellulose_I_segal"), "two increasing")
  expect_error(.clamp_preset_positions(5, "cellulose_I_segal"), "two increasing")
})

test_that(".click_slider_id picks the slider named by click_target, defaulting to i200", {
  expect_equal(.click_slider_id("iam"), "iam")
  expect_equal(.click_slider_id("i200"), "i200")
  expect_equal(.click_slider_id("anything_else"), "i200")
  expect_equal(.click_slider_id(NULL), "i200")
  expect_equal(.click_slider_id(character(0)), "i200")
})

test_that("run_textile_app(launch = FALSE) returns a shiny app object without starting a server", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("rmarkdown")
  app <- run_textile_app(launch = FALSE)
  expect_s3_class(app, "shiny.appobj")
})

test_that("the shipped app directory contains app.R", {
  app_dir <- system.file("shiny-app", package = "textileCrystR")
  skip_if(!nzchar(app_dir), "package not installed (running via load_all); app_dir lookup differs")
  expect_true(file.exists(file.path(app_dir, "app.R")))
})

test_that("run_textile_app() names the missing package(s) and how to install them", {
  testthat::local_mocked_bindings(
    requireNamespace = function(package, ...) package != "shiny",
    .package = "base"
  )
  err <- tryCatch(run_textile_app(launch = FALSE), error = function(e) e)
  expect_match(conditionMessage(err), "shiny")
  expect_match(conditionMessage(err), 'install.packages\\(c\\("shiny"\\)\\)')
})

test_that("run_textile_app() reports both packages when both are missing", {
  testthat::local_mocked_bindings(
    requireNamespace = function(package, ...) FALSE,
    .package = "base"
  )
  err <- tryCatch(run_textile_app(launch = FALSE), error = function(e) e)
  expect_match(conditionMessage(err), "shiny.*and.*rmarkdown|rmarkdown.*and.*shiny")
})

test_that("run_textile_app() reports a missing app directory clearly", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("rmarkdown")
  err <- tryCatch(run_textile_app(launch = FALSE, .app_dir = ""), error = function(e) e)
  expect_match(conditionMessage(err), "Could not find textileCrystR's Shiny app directory")
})
