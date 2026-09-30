test_that("import_xrd_file detects the paired multi-sample layout (cotton.csv)", {
  path <- cotton_path()
  xrd <- import_xrd_file(path)
  expect_s3_class(xrd, "tbl_df")
  expect_named(xrd, c("sample", "two_theta", "intensity"))
  expect_setequal(unique(xrd$sample), names(cotton_expected_2dp))
  expect_equal(nrow(xrd), 7L * 8701L)
  # must agree exactly with the dedicated paired-layout reader
  expect_equal(xrd, read_xrd_pairs(path))
})

test_that("a single pattern with a recognizable header is read correctly, in any column order", {
  p <- simple_pattern()
  f <- .write_fixture("named_std.csv",
                      c("2-theta,Intensity", paste(p$two_theta, p$intensity, sep = ",")))
  d <- import_xrd_file(f)
  expect_equal(d$two_theta, p$two_theta)
  expect_equal(d$intensity, p$intensity)
  expect_equal(unique(d$sample), "named_std")

  f2 <- .write_fixture("named_swapped.csv",
                       c("Counts,Angle", paste(p$intensity, p$two_theta, sep = ",")))
  d2 <- import_xrd_file(f2)
  expect_equal(d2$two_theta, p$two_theta)
  expect_equal(d2$intensity, p$intensity)
})

test_that("a header-less single pattern falls back to column order 1, 2", {
  p <- simple_pattern()
  f <- .write_fixture("noheader.csv", paste(p$two_theta, p$intensity, sep = ","))
  d <- import_xrd_file(f)
  expect_equal(d$two_theta, p$two_theta)
  expect_equal(d$intensity, p$intensity)
})

test_that("tab and whitespace delimited files are auto-detected", {
  p <- simple_pattern()
  f_tab <- .write_fixture("tabbed.txt", paste(p$two_theta, p$intensity, sep = "\t"))
  d_tab <- import_xrd_file(f_tab)
  expect_equal(d_tab$two_theta, p$two_theta)

  f_ws <- .write_fixture("whitespace.xy",
                         paste(sprintf("%.3f", p$two_theta), sprintf("%.2f", p$intensity)))
  d_ws <- import_xrd_file(f_ws)
  expect_equal(d_ws$two_theta, p$two_theta, tolerance = 1e-6)
  expect_equal(d_ws$intensity, p$intensity, tolerance = 1e-6)

  f_semi <- .write_fixture("semicolon.csv", paste(p$two_theta, p$intensity, sep = ";"))
  d_semi <- import_xrd_file(f_semi)
  expect_equal(d_semi$two_theta, p$two_theta)
})

test_that("explicit cols overrides both header-name and positional guessing", {
  p <- simple_pattern()
  f <- .write_fixture("explicit.csv",
                      c("Angle,Counts", paste(p$two_theta, p$intensity, sep = ",")))
  d1 <- import_xrd_file(f, cols = c("Angle", "Counts"))
  expect_equal(d1$two_theta, p$two_theta)
  d2 <- import_xrd_file(f, cols = c(1, 2))
  expect_equal(d2$two_theta, p$two_theta)
})

test_that("a custom `sample` name overrides the file-name default", {
  p <- simple_pattern()
  f <- .write_fixture("anyname.csv", paste(p$two_theta, p$intensity, sep = ","))
  d <- import_xrd_file(f, sample = "my sample")
  expect_equal(unique(d$sample), "my sample")
})

test_that("a nonexistent path gives a plain-English error naming the path", {
  expect_error(import_xrd_file("/no/such/path/here.csv"), "does not exist")
  expect_error(import_xrd_file("/no/such/path/here.csv"), "/no/such/path/here.csv", fixed = TRUE)
})

test_that("a directory instead of a file is refused clearly", {
  expect_error(import_xrd_file(.fixture_dir()), "folder, not a file")
})

test_that("`path` must be a single string", {
  expect_error(import_xrd_file(1), "single file path")
  expect_error(import_xrd_file(c("a", "b")), "single file path")
})

test_that("a file with fewer than two non-blank lines is refused", {
  f <- .write_fixture("one_line.csv", "only,one,line")
  expect_error(import_xrd_file(f), "fewer than two non-blank lines")
})

test_that("a file with no detectable delimiter is refused with a helpful preview", {
  f <- .write_fixture("garbage.txt", c("nothing structured here", "still nothing", "nope"))
  err <- tryCatch(import_xrd_file(f), error = function(e) e)
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "Could not work out how")
  expect_match(conditionMessage(err), "nothing structured here", fixed = TRUE)
})

test_that("an unrecognized header plus an implausible 2-theta column is refused", {
  f <- .write_fixture("mislabeled.csv",
                      c("Sample,Reading", "501,3200", "502,3305", "503,3100", "504,2950"))
  err <- tryCatch(import_xrd_file(f), error = function(e) e)
  expect_match(conditionMessage(err), "does not look like a 2-theta axis")
  expect_match(conditionMessage(err), "did not clearly name a 2-theta")
  expect_match(conditionMessage(err), "cols = c\\(two_theta_column, intensity_column\\)")
})

test_that("a header-less file with an implausible 2-theta column is refused", {
  f <- .write_fixture("noheader_bad.csv", paste(5000:5010, 400 + (1:11), sep = ","))
  err <- tryCatch(import_xrd_file(f), error = function(e) e)
  expect_match(conditionMessage(err), "does not look like a 2-theta axis")
  expect_match(conditionMessage(err), "no header row")
})

test_that("explicit cols with an unknown column name is refused, listing the real header", {
  f <- .write_fixture("named2.csv", c("Angle,Counts", "10,100", "11,110"))
  err <- tryCatch(import_xrd_file(f, cols = c("TwoTheta", "Counts")), error = function(e) e)
  expect_match(conditionMessage(err), "'TwoTheta'.*not found", perl = TRUE)
  expect_match(conditionMessage(err), "Angle, Counts", fixed = TRUE)
})

test_that("explicit cols naming a column that does not exist in a header-less file errors clearly", {
  f <- .write_fixture("noheader2.csv", "10,100\n11,110")
  expect_error(import_xrd_file(f, cols = c("Angle", "Counts")), "no header row to match against")
})

test_that("`cols` must have exactly two elements", {
  f <- .write_fixture("simple3.csv", "10,100\n11,110")
  expect_error(import_xrd_file(f, cols = 1), "exactly two columns")
  expect_error(import_xrd_file(f, cols = c(1, 2, 3)), "exactly two columns")
})

test_that("an out-of-range explicit column number is refused", {
  f <- .write_fixture("simple4.csv", "10,100\n11,110")
  expect_error(import_xrd_file(f, cols = c(1, 5)), "only has 2 column")
})

test_that("all-non-numeric data in the chosen columns is refused", {
  f <- .write_fixture("allna.csv", c("Angle,Counts", "a,b", "c,d"))
  err <- tryCatch(import_xrd_file(f), error = function(e) e)
  expect_match(conditionMessage(err), "contains no numeric values")
})

test_that("a broken multi-sample file names the underlying read_xrd_pairs() problem", {
  # even number of columns, looks paired, but the sample-name row is blank
  # for one pair -> read_xrd_pairs() itself should reject it
  f <- .write_fixture("broken_paired.csv",
                      c(",,B/1,", "u,u,u,u", "10,1,10,7", "11,2,11,8"))
  err <- tryCatch(import_xrd_file(f), error = function(e) e)
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "multi-sample file")
  expect_match(conditionMessage(err), "single XRD pattern")
})
