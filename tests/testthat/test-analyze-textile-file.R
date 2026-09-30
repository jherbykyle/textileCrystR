test_that("a single-sample selection returns an ordinary textile_cryst, computed correctly", {
  res <- analyze_textile_file(cotton_path(), sample = "0%/1", preset = "cotton_validation",
                              plot = FALSE)
  expect_s3_class(res, "textile_cryst")
  expect_equal(round(res$ci, 2), unname(cotton_expected_2dp["0%/1"]))
  expect_equal(res$sample, "0%/1")
})

test_that("defaulting preset to cellulose_I_segal happens only when i200/iam are both unset", {
  a <- analyze_textile_file(cotton_path(), sample = "0%/1", plot = FALSE)
  expect_equal(a$method$preset, "cellulose_I_segal")
  b <- analyze_textile_file(cotton_path(), sample = "0%/1", i200 = 22.6, iam = 18, plot = FALSE)
  expect_null(b$method$preset)
  cc <- analyze_textile_file(cotton_path(), sample = "0%/1", preset = "cotton_validation",
                            plot = FALSE)
  expect_equal(cc$method$preset, "cotton_validation")
})

test_that("no sample given on a multi-sample file returns a textile_cryst_batch for all samples", {
  batch <- analyze_textile_file(cotton_path(), preset = "cotton_validation", plot = FALSE)
  expect_s3_class(batch, "textile_cryst_batch")
  expect_length(batch$results, 7L)
  expect_setequal(names(batch$results), names(cotton_expected_2dp))
  got <- setNames(round(batch$summary$ci, 2), batch$summary$sample)
  expect_equal(got[names(cotton_expected_2dp)], cotton_expected_2dp)
})

test_that("an explicit sample= on a single-pattern file still works and is not treated as a batch", {
  p <- simple_pattern()
  f <- .write_fixture("single_named.csv", paste(p$two_theta, p$intensity, sep = ","))
  res <- analyze_textile_file(f, sample = "single_named", i200 = 22.6, iam = 18, plot = FALSE)
  expect_s3_class(res, "textile_cryst")
})

test_that("plot = FALSE suppresses the side-effect plot without affecting the result", {
  res1 <- analyze_textile_file(cotton_path(), sample = "0%/1", preset = "cotton_validation",
                               plot = FALSE)
  res2 <- analyze_textile_file(cotton_path(), sample = "0%/1", preset = "cotton_validation",
                               plot = TRUE)
  expect_equal(res1$ci, res2$ci)
})

test_that("plot_file writes a plot image for a single sample and returns the result invisibly", {
  f <- tempfile(fileext = ".png")
  on.exit(unlink(f), add = TRUE)
  expect_false(file.exists(f))
  val <- withVisible(analyze_textile_file(cotton_path(), sample = "0%/1",
                                          preset = "cotton_validation", plot = FALSE,
                                          plot_file = f))
  expect_true(file.exists(f))
  expect_gt(file.info(f)$size, 0)
  expect_false(val$visible)
  expect_s3_class(val$value, "textile_cryst")
})

test_that("plot_file writes a combined plot for a batch result", {
  f <- tempfile(fileext = ".png")
  on.exit(unlink(f), add = TRUE)
  val <- withVisible(analyze_textile_file(cotton_path(), preset = "cotton_validation",
                                          plot = FALSE, plot_file = f))
  expect_true(file.exists(f))
  expect_false(val$visible)
  expect_s3_class(val$value, "textile_cryst_batch")
})

test_that("report_file writes a self-contained HTML report for single and batch results", {
  skip_if_not_installed("rmarkdown")
  f1 <- tempfile(fileext = ".html")
  on.exit(unlink(f1), add = TRUE)
  analyze_textile_file(cotton_path(), sample = "0%/1", preset = "cotton_validation",
                       plot = FALSE, report_file = f1)
  expect_true(file.exists(f1))
  html1 <- paste(readLines(f1, warn = FALSE), collapse = "\n")
  expect_match(html1, "Segal Crystallinity Index")
  expect_match(html1, "0%/1", fixed = TRUE)

  f2 <- tempfile(fileext = ".html")
  on.exit(unlink(f2), add = TRUE)
  analyze_textile_file(cotton_path(), preset = "cotton_validation", plot = FALSE,
                       report_file = f2)
  expect_true(file.exists(f2))
  html2 <- paste(readLines(f2, warn = FALSE), collapse = "\n")
  expect_match(html2, "7 sample")
})

test_that("an unrecognized sample name lists the samples that do exist", {
  err <- tryCatch(analyze_textile_file(cotton_path(), sample = "999%/1", plot = FALSE),
                  error = function(e) e)
  expect_match(conditionMessage(err), "999%/1", fixed = TRUE)
  expect_match(conditionMessage(err), "80%/1", fixed = TRUE)
})

test_that("a read failure is re-raised with analyze_textile_file() context and the original message", {
  f <- .write_fixture("garbage2.txt", c("no structure", "still none", "nope"))
  err <- tryCatch(analyze_textile_file(f, plot = FALSE), error = function(e) e)
  expect_match(conditionMessage(err), "analyze_textile_file\\(\\) could not read")
  expect_match(conditionMessage(err), "Could not work out how")
})

test_that("a position outside the measured range is re-raised with the sample/file in context", {
  err <- tryCatch(
    analyze_textile_file(cotton_path(), sample = "0%/1", i200 = 150, iam = 18, plot = FALSE),
    error = function(e) e
  )
  expect_match(conditionMessage(err), "sample '0%/1'", fixed = TRUE)
  expect_match(conditionMessage(err), "outside the measured")
  expect_match(conditionMessage(err), "wide enough 2-theta range")
})

test_that("textile_cryst_batch print/plot/as.data.frame methods behave as documented", {
  batch <- analyze_textile_file(cotton_path(), preset = "cotton_validation", plot = FALSE)

  out <- capture.output(print(batch))
  expect_true(any(grepl("7 sample", out)))
  expect_true(any(grepl("0%/1", out)))
  expect_identical(print(batch), batch)

  df <- as.data.frame(batch)
  expect_equal(nrow(df), 7L)
  expect_true(all(c("sample", "ci", "i200_position", "iam_position") %in% names(df)))
  expect_equal(df, batch$summary)

  p <- plot(batch)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_match(p$labels$title, "Segal Crystallinity Index by sample")
})

test_that("plot = TRUE prints the batch plot as a side effect without changing the result", {
  res_quiet <- analyze_textile_file(cotton_path(), preset = "cotton_validation", plot = FALSE)
  res_shown <- analyze_textile_file(cotton_path(), preset = "cotton_validation", plot = TRUE)
  expect_equal(res_quiet$summary$ci, res_shown$summary$ci)
})

test_that("as.data.frame.textile_cryst_batch respects an explicit row.names", {
  batch <- analyze_textile_file(cotton_path(), preset = "cotton_validation", plot = FALSE)
  df <- as.data.frame(batch, row.names = LETTERS[seq_len(nrow(batch$summary))])
  expect_equal(rownames(df), LETTERS[seq_len(nrow(batch$summary))])
})

test_that(".render_report reports a missing rmarkdown installation clearly", {
  testthat::local_mocked_bindings(
    requireNamespace = function(package, ...) package != "rmarkdown",
    .package = "base"
  )
  res <- analyze_textile_file(cotton_path(), sample = "0%/1", preset = "cotton_validation",
                              plot = FALSE)
  err <- tryCatch(.render_report(list(result = res), tempfile(fileext = ".html")),
                  error = function(e) e)
  expect_match(conditionMessage(err), "needs the rmarkdown package")
  expect_match(conditionMessage(err), 'install.packages\\("rmarkdown"\\)')
})

test_that(".render_report reports a missing template file clearly", {
  res <- analyze_textile_file(cotton_path(), sample = "0%/1", preset = "cotton_validation",
                              plot = FALSE)
  err <- tryCatch(.render_report(list(result = res), tempfile(fileext = ".html"), template = ""),
                  error = function(e) e)
  expect_match(conditionMessage(err), "Could not find textileCrystR's report template")
})

test_that("a sample with a warning is reflected in the batch print output", {
  tt <- c(5, 10, 18, 22.6, 30, 40)
  y <- c(100, 200, 250, 1000, 300, 100)
  y2 <- c(100, 200, 1200, 1000, 300, 100)  # Iam > I200 -> warning, negative CI
  results <- list(
    ok = segal_ci(tt, y, i200 = 22.6, iam = 18, sample = "ok"),
    warned = segal_ci(tt, y2, i200 = 22.6, iam = 18, sample = "warned")
  )
  summ <- do.call(rbind, lapply(results, as.data.frame))
  rownames(summ) <- NULL
  batch <- structure(list(results = results, summary = summ), class = "textile_cryst_batch")
  expect_gt(length(results$warned$warnings), 0L)
  out <- capture.output(print(batch))
  expect_true(any(grepl("1 sample\\(s\\) raised a warning", out)))
})
