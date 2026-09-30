# Tests for the redesigned dashboard: multi-file upload, multi-sample pooling,
# a method selector (Segal / peak area), single-sample interactive mode, and
# batch mode. See test-run-textile-app.R for run_textile_app() itself and the
# pure helper functions (.clamp_preset_positions, .click_slider_id); this file
# covers the app's own reactive logic.

.upload_df <- function(names, paths) {
  data.frame(name = names, size = 1, type = "text/csv", datapath = paths, stringsAsFactors = FALSE)
}

test_that("uploading the shipped cotton.csv pools all 7 samples with clean (no file-prefix) labels", {
  skip_if_not_installed("shiny")
  app_dir <- system.file("shiny-app", package = "textileCrystR")
  skip_if(!nzchar(app_dir), "package not installed (running via load_all)")

  shiny::testServer(app_dir, {
    session$setInputs(files = .upload_df("cotton.csv", cotton_path()))
    expect_false(multi_file())
    sc <- sample_choices()
    expect_setequal(names(sc), c("0%/1", "20%/1", "40%/1", "50%/1", "60%/1", "80%/1", "100%/1"))
    expect_false(any(grepl("::", names(sc))))       # single file: no compound-key prefix in labels
    expect_true(all(grepl("cotton\\.csv ::", unname(sc))))  # but the underlying keys are still unique
  })
})

test_that("uploading two files makes multi_file() TRUE and labels carry the file name", {
  skip_if_not_installed("shiny")
  app_dir <- system.file("shiny-app", package = "textileCrystR")
  skip_if(!nzchar(app_dir), "package not installed (running via load_all)")

  f2 <- tempfile(fileext = ".csv")
  tt <- seq(10, 40, by = 0.1)
  writeLines(c("Angle,Counts", paste(tt, round(1000 + 50 * sin(tt), 1), sep = ",")), f2)
  on.exit(unlink(f2), add = TRUE)

  shiny::testServer(app_dir, {
    session$setInputs(files = .upload_df(c("cotton.csv", "extra.csv"), c(cotton_path(), f2)))
    expect_true(multi_file())
    sc <- sample_choices()
    expect_true(all(grepl("::", names(sc))))
    expect_equal(length(sc), 8L)  # 7 cotton samples + 1 from the extra file
  })
})

test_that("a file that fails to import surfaces its own friendly message via shiny::validate()", {
  skip_if_not_installed("shiny")
  app_dir <- system.file("shiny-app", package = "textileCrystR")
  skip_if(!nzchar(app_dir), "package not installed (running via load_all)")

  bad <- tempfile(fileext = ".csv")
  writeLines(c("no structure here", "still none", "nope"), bad)
  on.exit(unlink(bad), add = TRUE)

  shiny::testServer(app_dir, {
    session$setInputs(files = .upload_df("bad.csv", bad))
    err <- tryCatch(all_data(), error = function(e) e, condition = function(c) c)
    expect_s3_class(err, "shiny.silent.error")
    expect_match(conditionMessage(err), "bad.csv")
    expect_match(conditionMessage(err), "Could not work out how")
  })
})

test_that("Segal single-sample mode computes the same result as segal_ci() directly", {
  skip_if_not_installed("shiny")
  app_dir <- system.file("shiny-app", package = "textileCrystR")
  skip_if(!nzchar(app_dir), "package not installed (running via load_all)")

  d0 <- read_xrd_pairs(cotton_path())
  d0 <- d0[d0$sample == "0%/1", ]

  shiny::testServer(app_dir, {
    session$setInputs(files = .upload_df("cotton.csv", cotton_path()))
    key0 <- unname(sample_choices()["0%/1"])
    session$setInputs(view_mode = "single", sample_single = key0, method = "segal",
                      preset = "cotton_validation", i200 = 22.6, iam = 18)
    res <- result()
    expect_s3_class(res, "textile_cryst")
    expect_equal(res$ci, segal_ci(d0$two_theta, d0$intensity, i200 = 22.6, iam = 18)$ci)
    expect_equal(res$sample, "0%/1")  # the clean label, not the compound key
  })
})

test_that("peak-area single-sample mode computes the same result as fit_crystalline_peaks() directly", {
  skip_if_not_installed("shiny")
  app_dir <- system.file("shiny-app", package = "textileCrystR")
  skip_if(!nzchar(app_dir), "package not installed (running via load_all)")

  d0 <- read_xrd_pairs(cotton_path())
  d0 <- d0[d0$sample == "0%/1", ]

  shiny::testServer(app_dir, {
    session$setInputs(files = .upload_df("cotton.csv", cotton_path()))
    key0 <- unname(sample_choices()["0%/1"])
    session$setInputs(view_mode = "single", sample_single = key0, method = "peakarea",
                      model = "pseudo_voigt", centers = "", am_win_lo = 18, am_win_hi = 25,
                      background = "linear")
    session$setInputs(run_fit = 1)
    res <- result()
    expect_s3_class(res, "textile_cryst_peakarea")
    expect_equal(res$ci, fit_crystalline_peaks(d0$two_theta, d0$intensity, sample = "0%/1")$ci)
  })
})

test_that("custom peak centers are parsed and used, and malformed centers give a friendly error", {
  skip_if_not_installed("shiny")
  app_dir <- system.file("shiny-app", package = "textileCrystR")
  skip_if(!nzchar(app_dir), "package not installed (running via load_all)")

  shiny::testServer(app_dir, {
    session$setInputs(files = .upload_df("cotton.csv", cotton_path()))
    key0 <- unname(sample_choices()["0%/1"])
    session$setInputs(view_mode = "single", sample_single = key0, method = "peakarea",
                      model = "gaussian", centers = "14.8, 16.5, 22.6, 34.5",
                      am_win_lo = 18, am_win_hi = 25, background = "none")
    session$setInputs(run_fit = 1)
    expect_equal(nrow(result()$peaks), 4L)

    session$setInputs(centers = "abc, def")
    session$setInputs(run_fit = 2)
    err <- tryCatch(result(), error = function(e) e, condition = function(c) c)
    expect_s3_class(err, "shiny.silent.error")
    expect_match(conditionMessage(err), "comma-separated list of numbers")
  })
})

test_that("batch Segal mode uses the preset directly and matches analyze_textile_file()'s batch results", {
  skip_if_not_installed("shiny")
  app_dir <- system.file("shiny-app", package = "textileCrystR")
  skip_if(!nzchar(app_dir), "package not installed (running via load_all)")

  expected <- analyze_textile_file(cotton_path(), preset = "cellulose_I_segal", plot = FALSE)

  shiny::testServer(app_dir, {
    session$setInputs(files = .upload_df("cotton.csv", cotton_path()))
    # deliberately do NOT set i200/iam -- batch mode must not depend on them (see regression test below)
    session$setInputs(view_mode = "batch", samples_batch = character(0), method = "segal",
                      preset = "cellulose_I_segal")
    session$setInputs(run_batch = 1)
    b <- batch_result()
    expect_equal(b$method, "segal")
    expect_equal(nrow(b$summary), 7L)
    expect_length(b$errors, 0L)
    got <- stats::setNames(b$summary$ci, b$summary$sample)
    exp <- stats::setNames(expected$summary$ci, expected$summary$sample)
    expect_equal(got[names(exp)], exp)
  })
})

test_that("regression: batch Segal does not depend on the single-sample i200/iam sliders", {
  # This is the specific bug found during development: the batch loop used
  # to read input$i200/input$iam (tied to whichever one sample was last
  # viewed in single mode, and never initialized if the user goes straight
  # to batch mode), so every batch run failed with "Give both i200 and iam".
  # Fixed by using the preset directly, recomputed per sample.
  skip_if_not_installed("shiny")
  app_dir <- system.file("shiny-app", package = "textileCrystR")
  skip_if(!nzchar(app_dir), "package not installed (running via load_all)")

  shiny::testServer(app_dir, {
    session$setInputs(files = .upload_df("cotton.csv", cotton_path()))
    session$setInputs(view_mode = "batch", samples_batch = character(0), method = "segal",
                      preset = "cotton_validation")
    # input$i200 and input$iam are unset (NULL) at this point -- batch must still work
    session$setInputs(run_batch = 1)
    b <- batch_result()
    expect_equal(nrow(b$summary), 7L)
    expect_length(b$errors, 0L)
  })
})

test_that("batch mode respects a specific (non-empty) sample selection", {
  skip_if_not_installed("shiny")
  app_dir <- system.file("shiny-app", package = "textileCrystR")
  skip_if(!nzchar(app_dir), "package not installed (running via load_all)")

  shiny::testServer(app_dir, {
    session$setInputs(files = .upload_df("cotton.csv", cotton_path()))
    chosen <- unname(sample_choices()[c("0%/1", "100%/1")])
    session$setInputs(view_mode = "batch", samples_batch = chosen, method = "peakarea",
                      model = "pseudo_voigt", centers = "", am_win_lo = 18, am_win_hi = 25,
                      background = "linear")
    session$setInputs(run_batch = 1)
    b <- batch_result()
    expect_equal(nrow(b$summary), 2L)
    expect_setequal(b$summary$sample, c("0%/1", "100%/1"))
  })
})

test_that("batch mode reports partial failures without discarding the successful results", {
  skip_if_not_installed("shiny")
  app_dir <- system.file("shiny-app", package = "textileCrystR")
  skip_if(!nzchar(app_dir), "package not installed (running via load_all)")

  narrow <- tempfile(fileext = ".csv")
  tt <- seq(19, 21, by = 0.05)  # too narrow to cover the Segal preset positions
  writeLines(c("Angle,Counts", paste(tt, round(1000 + 50 * sin(tt), 1), sep = ",")), narrow)
  on.exit(unlink(narrow), add = TRUE)

  shiny::testServer(app_dir, {
    session$setInputs(files = .upload_df(c("cotton.csv", "narrow.csv"), c(cotton_path(), narrow)))
    session$setInputs(view_mode = "batch", samples_batch = character(0), method = "segal",
                      preset = "cellulose_I_segal")
    session$setInputs(run_batch = 1)
    b <- batch_result()
    expect_equal(nrow(b$summary), 7L)       # all 7 cotton samples still succeed
    expect_length(b$errors, 1L)             # the narrow-range file's one sample fails
    expect_match(b$errors, "narrow.csv")
    expect_match(b$errors, "outside the measured")
  })
})

test_that("batch mode errors clearly when every selected sample fails", {
  skip_if_not_installed("shiny")
  app_dir <- system.file("shiny-app", package = "textileCrystR")
  skip_if(!nzchar(app_dir), "package not installed (running via load_all)")

  narrow <- tempfile(fileext = ".csv")
  tt <- seq(19, 21, by = 0.05)
  writeLines(c("Angle,Counts", paste(tt, round(1000 + 50 * sin(tt), 1), sep = ",")), narrow)
  on.exit(unlink(narrow), add = TRUE)

  shiny::testServer(app_dir, {
    session$setInputs(files = .upload_df("narrow.csv", narrow))
    session$setInputs(view_mode = "batch", samples_batch = character(0), method = "segal",
                      preset = "cellulose_I_segal")
    session$setInputs(run_batch = 1)
    err <- tryCatch(batch_result(), error = function(e) e, condition = function(c) c)
    expect_s3_class(err, "shiny.silent.error")
    expect_match(conditionMessage(err), "Every selected sample failed")
  })
})

test_that("downloads (plot and report) work for both methods from the app's own reactive state", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("rmarkdown")
  app_dir <- system.file("shiny-app", package = "textileCrystR")
  skip_if(!nzchar(app_dir), "package not installed (running via load_all)")

  shiny::testServer(app_dir, {
    session$setInputs(files = .upload_df("cotton.csv", cotton_path()))
    key0 <- unname(sample_choices()["0%/1"])
    session$setInputs(view_mode = "single", sample_single = key0, method = "segal",
                      preset = "cotton_validation", i200 = 22.6, iam = 18)
    res <- result()
    rep1 <- tempfile(fileext = ".html")
    textileCrystR:::.render_report(list(result = res), rep1)
    expect_true(file.exists(rep1))
    plot1 <- tempfile(fileext = ".png")
    ggplot2::ggsave(plot1, plot(res), width = 7, height = 4.2, dpi = 100)
    expect_gt(file.info(plot1)$size, 0)

    session$setInputs(method = "peakarea", model = "pseudo_voigt", centers = "",
                      am_win_lo = 18, am_win_hi = 25, background = "linear")
    session$setInputs(run_fit = 1)
    res2 <- result()
    expect_s3_class(res2, "textile_cryst_peakarea")
    rep2 <- tempfile(fileext = ".html")
    textileCrystR:::.render_report(list(result = res2), rep2)
    expect_true(file.exists(rep2))
    plot2 <- tempfile(fileext = ".png")
    ggplot2::ggsave(plot2, plot(res2), width = 8, height = 6, dpi = 100)
    expect_gt(file.info(plot2)$size, 0)
  })
})

test_that("batch CSV export contains the expected columns and values", {
  skip_if_not_installed("shiny")
  app_dir <- system.file("shiny-app", package = "textileCrystR")
  skip_if(!nzchar(app_dir), "package not installed (running via load_all)")

  shiny::testServer(app_dir, {
    session$setInputs(files = .upload_df("cotton.csv", cotton_path()))
    session$setInputs(view_mode = "batch", samples_batch = character(0), method = "segal",
                      preset = "cotton_validation")
    session$setInputs(run_batch = 1)
    b <- batch_result()
    f <- tempfile(fileext = ".csv")
    utils::write.csv(b$summary, f, row.names = FALSE)
    back <- utils::read.csv(f)
    expect_setequal(back$sample, c("0%/1", "20%/1", "40%/1", "50%/1", "60%/1", "80%/1", "100%/1"))
    expect_equal(round(back$ci[back$sample == "0%/1"], 2), 71.19)
  })
})

# ---- background subtraction ----

test_that("background subtraction is off by default: current_data() equals the raw data", {
  skip_if_not_installed("shiny")
  app_dir <- system.file("shiny-app", package = "textileCrystR")
  skip_if(!nzchar(app_dir), "package not installed (running via load_all)")

  d0 <- read_xrd_pairs(cotton_path()); d0 <- d0[d0$sample == "0%/1", ]

  shiny::testServer(app_dir, {
    session$setInputs(files = .upload_df("cotton.csv", cotton_path()))
    key0 <- unname(sample_choices()["0%/1"])
    session$setInputs(view_mode = "single", sample_single = key0, method = "segal",
                      preset = "cotton_validation", i200 = 22.6, iam = 18)
    expect_identical(current_data()$intensity, current_data_raw()$intensity)
    expect_equal(result()$ci, segal_ci(d0$two_theta, d0$intensity, i200 = 22.6, iam = 18)$ci)
  })
})

test_that("enabling background subtraction changes current_data() to match a direct computation", {
  skip_if_not_installed("shiny")
  app_dir <- system.file("shiny-app", package = "textileCrystR")
  skip_if(!nzchar(app_dir), "package not installed (running via load_all)")

  d0 <- read_xrd_pairs(cotton_path()); d0 <- d0[d0$sample == "0%/1", ]
  bg <- estimate_xrd_background(d0$two_theta, d0$intensity, window = 8)
  corrected <- subtract_xrd_background(d0$two_theta, d0$intensity, bg, clip_negative = TRUE)

  shiny::testServer(app_dir, {
    session$setInputs(files = .upload_df("cotton.csv", cotton_path()))
    key0 <- unname(sample_choices()["0%/1"])
    session$setInputs(view_mode = "single", sample_single = key0, method = "segal",
                      preset = "cotton_validation", i200 = 22.6, iam = 18,
                      subtract_background = TRUE, bg_window = 8)
    expect_equal(current_data()$intensity, corrected$corrected)

    expected <- segal_ci(corrected$two_theta, corrected$corrected, i200 = 22.6, iam = 18)
    expect_equal(result()$ci, expected$ci)
    expect_false(isTRUE(all.equal(result()$ci,
                                  segal_ci(d0$two_theta, d0$intensity, i200 = 22.6, iam = 18)$ci)))
  })
})

test_that("background subtraction also applies to peak-area single-sample results", {
  skip_if_not_installed("shiny")
  app_dir <- system.file("shiny-app", package = "textileCrystR")
  skip_if(!nzchar(app_dir), "package not installed (running via load_all)")

  d0 <- read_xrd_pairs(cotton_path()); d0 <- d0[d0$sample == "0%/1", ]
  bg <- estimate_xrd_background(d0$two_theta, d0$intensity, window = 8)
  corrected <- subtract_xrd_background(d0$two_theta, d0$intensity, bg, clip_negative = TRUE)

  shiny::testServer(app_dir, {
    session$setInputs(files = .upload_df("cotton.csv", cotton_path()))
    key0 <- unname(sample_choices()["0%/1"])
    session$setInputs(view_mode = "single", sample_single = key0, method = "peakarea",
                      model = "pseudo_voigt", centers = "", am_win_lo = 18, am_win_hi = 25,
                      background = "linear", subtract_background = TRUE, bg_window = 8)
    session$setInputs(run_fit = 1)
    res <- result()
    expect_s3_class(res, "textile_cryst_peakarea")
    expected <- fit_crystalline_peaks(corrected$two_theta, corrected$corrected, sample = "0%/1")
    expect_equal(res$ci, expected$ci)
  })
})

test_that("a window too wide for the pattern is captured as a warning note, not an error", {
  skip_if_not_installed("shiny")
  app_dir <- system.file("shiny-app", package = "textileCrystR")
  skip_if(!nzchar(app_dir), "package not installed (running via load_all)")

  short <- tempfile(fileext = ".csv")
  tt <- seq(15, 25, by = 0.5)
  writeLines(c("Angle,Counts", paste(tt, round(300 + 2000 * exp(-((tt - 20) / 1)^2), 1), sep = ",")), short)
  on.exit(unlink(short), add = TRUE)

  shiny::testServer(app_dir, {
    session$setInputs(files = .upload_df("short.csv", short))
    session$setInputs(view_mode = "single", sample_single = unname(sample_choices())[1],
                      method = "segal", i200 = 20, iam = 17,
                      subtract_background = TRUE, bg_window = 8)
    bg <- current_background()
    expect_false(inherits(bg, "error"))
    expect_match(attr(bg, "warning"), "using the largest window")
    expect_no_error(current_data())  # still produces a usable (capped-window) result, not an error
  })
})

test_that("a background estimate that genuinely fails surfaces a clear error, not a crash", {
  skip_if_not_installed("shiny")
  app_dir <- system.file("shiny-app", package = "textileCrystR")
  skip_if(!nzchar(app_dir), "package not installed (running via load_all)")

  degenerate <- tempfile(fileext = ".csv")
  writeLines(c("Angle,Counts", "10,5", "10,6"), degenerate)
  on.exit(unlink(degenerate), add = TRUE)

  shiny::testServer(app_dir, {
    session$setInputs(files = .upload_df("degenerate.csv", degenerate))
    session$setInputs(view_mode = "single", sample_single = unname(sample_choices())[1],
                      method = "segal", subtract_background = TRUE, bg_window = 8)
    err <- tryCatch(current_data(), error = function(e) e, condition = function(c) c)
    expect_s3_class(err, "shiny.silent.error")
    expect_match(conditionMessage(err), "Could not determine a typical 2-theta step")
  })
})

test_that("batch mode applies background subtraction consistently across all selected samples", {
  skip_if_not_installed("shiny")
  app_dir <- system.file("shiny-app", package = "textileCrystR")
  skip_if(!nzchar(app_dir), "package not installed (running via load_all)")

  d0 <- read_xrd_pairs(cotton_path()); d0 <- d0[d0$sample == "0%/1", ]
  bg <- estimate_xrd_background(d0$two_theta, d0$intensity, window = 8)
  corrected <- subtract_xrd_background(d0$two_theta, d0$intensity, bg, clip_negative = TRUE)
  expected <- segal_ci(corrected$two_theta, corrected$corrected, preset = "cotton_validation")

  shiny::testServer(app_dir, {
    session$setInputs(files = .upload_df("cotton.csv", cotton_path()))
    session$setInputs(view_mode = "batch", samples_batch = character(0), method = "segal",
                      preset = "cotton_validation", subtract_background = TRUE, bg_window = 8)
    session$setInputs(run_batch = 1)
    b <- batch_result()
    expect_equal(nrow(b$summary), 7L)
    expect_length(b$errors, 0L)
    expect_equal(b$summary$ci[b$summary$sample == "0%/1"], expected$ci)
  })
})

test_that("regression: batch mode reports a per-sample background failure without discarding the rest", {
  skip_if_not_installed("shiny")
  app_dir <- system.file("shiny-app", package = "textileCrystR")
  skip_if(!nzchar(app_dir), "package not installed (running via load_all)")

  degenerate <- tempfile(fileext = ".csv")
  writeLines(c("Angle,Counts", "10,5", "10,6"), degenerate)
  on.exit(unlink(degenerate), add = TRUE)

  shiny::testServer(app_dir, {
    session$setInputs(files = .upload_df(c("cotton.csv", "degenerate.csv"), c(cotton_path(), degenerate)))
    session$setInputs(view_mode = "batch", samples_batch = character(0), method = "segal",
                      preset = "cotton_validation", subtract_background = TRUE, bg_window = 8)
    session$setInputs(run_batch = 1)
    b <- batch_result()
    expect_equal(nrow(b$summary), 7L)   # the 7 good cotton.csv samples still succeed
    expect_length(b$errors, 1L)
    expect_match(b$errors, "degenerate.csv")
  })
})

test_that("toggling background subtraction off after batch analysis reverts to the original values", {
  skip_if_not_installed("shiny")
  app_dir <- system.file("shiny-app", package = "textileCrystR")
  skip_if(!nzchar(app_dir), "package not installed (running via load_all)")

  shiny::testServer(app_dir, {
    session$setInputs(files = .upload_df("cotton.csv", cotton_path()))
    session$setInputs(view_mode = "batch", samples_batch = character(0), method = "segal",
                      preset = "cotton_validation", subtract_background = TRUE, bg_window = 8)
    session$setInputs(run_batch = 1)
    with_bg <- batch_result()$summary$ci[batch_result()$summary$sample == "0%/1"]

    session$setInputs(subtract_background = FALSE)
    session$setInputs(run_batch = 2)
    without_bg <- batch_result()$summary$ci[batch_result()$summary$sample == "0%/1"]

    expect_equal(round(without_bg, 2), 71.19)
    expect_false(isTRUE(all.equal(with_bg, without_bg)))
  })
})
