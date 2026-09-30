test_that("validate_crystallinity_input returns a cleaned pattern and a record", {
  v <- validate_crystallinity_input(c(12, 10, 11, 13), c(5, 3, 4, 6))
  expect_s3_class(v, "textile_cryst_input")
  expect_equal(v$two_theta, 10:13)
  expect_equal(v$intensity, c(3, 4, 5, 6))
  expect_true(v$was_sorted)
  expect_equal(v$n_input, 4L)
  expect_equal(v$n_used, 4L)
  expect_equal(v$range, c(10, 13))
  expect_equal(v$median_step, 1)
  expect_equal(v$max_step, 1)
  expect_output(print(v), "validated XRD pattern")
  expect_output(print(v), "sorted")
  clean <- validate_crystallinity_input(1:5, c(1, 2, 3, 2, 1))
  expect_length(clean$messages, 0L)
  expect_output(print(clean), "5 used of 5")
})

test_that("bad input types and lengths are rejected", {
  expect_error(validate_crystallinity_input("a", 1), "numeric")
  expect_error(validate_crystallinity_input(1:3, c("a", "b", "c")), "numeric")
  expect_error(validate_crystallinity_input(1:3, 1:2), "same length")
  expect_error(validate_crystallinity_input(1, 1), "two valid points")
  expect_error(validate_crystallinity_input(numeric(0), numeric(0)), "two valid points")
  expect_error(validate_crystallinity_input(1:3, 1:3, na.rm = NA), "na.rm")
  expect_error(validate_crystallinity_input(1:3, 1:3, duplicates = "drop"), "should be one of")
})

test_that("missing and non-finite values are rejected or dropped on request", {
  tt <- c(10, 11, NA, 13, 14)
  y <- c(1, 2, 3, Inf, 5)
  expect_error(validate_crystallinity_input(tt, y), "2 missing or non-finite")
  v <- validate_crystallinity_input(tt, y, na.rm = TRUE)
  expect_equal(v$n_dropped, 2L)
  expect_equal(v$two_theta, c(10, 11, 14))
  expect_error(validate_crystallinity_input(c(NA, NA, 1), c(1, 2, NA), na.rm = TRUE),
               "two valid points")
})

test_that("2-theta must lie between 0 and 180 degrees", {
  expect_error(validate_crystallinity_input(c(0, 10, 20), 1:3), "between 0 and 180")
  expect_error(validate_crystallinity_input(c(10, 20, 180), 1:3), "between 0 and 180")
  expect_error(validate_crystallinity_input(c(-5, 10, 20), 1:3), "between 0 and 180")
})

test_that("duplicates are an error or are averaged, and negatives are reported", {
  tt <- c(10, 11, 11, 12)
  y <- c(1, 2, 4, -3)
  expect_error(validate_crystallinity_input(tt, y), "1 duplicated")
  v <- validate_crystallinity_input(tt, y, duplicates = "mean")
  expect_equal(v$two_theta, c(10, 11, 12))
  expect_equal(v$intensity, c(1, 3, -3))
  expect_equal(v$n_negative, 1L)
  expect_true(any(grepl("negative intensity", v$messages)))
  expect_true(any(grepl("Averaged", v$messages)))
  expect_error(validate_crystallinity_input(c(5, 5, 5), c(1, 2, 3), duplicates = "mean"),
               "two distinct")
})

test_that("read_xrd_pairs reads paired columns exactly and handles short series", {
  f <- tempfile(fileext = ".csv")
  bom <- as.raw(c(0xef, 0xbb, 0xbf))
  body <- charToRaw(paste0(
    "A/1,,B/1,\r\n",
    "2theta,counts,2theta,counts\r\n",
    "10,1.5,10,7\r\n",
    "11,2.5,11,8\r\n",
    "12,3.5,,\r\n"))
  writeBin(c(bom, body), f)
  x <- read_xrd_pairs(f)
  expect_s3_class(x, "tbl_df")
  expect_named(x, c("sample", "two_theta", "intensity"))
  expect_equal(x$sample, c("A/1", "A/1", "A/1", "B/1", "B/1"))
  expect_equal(x$two_theta, c(10, 11, 12, 10, 11))
  expect_equal(x$intensity, c(1.5, 2.5, 3.5, 7, 8))
})

test_that("read_xrd_pairs rejects malformed files", {
  expect_error(read_xrd_pairs("no_such_file.csv"), "existing file")
  expect_error(read_xrd_pairs(1), "existing file")
  f <- tempfile(fileext = ".csv")
  writeLines(c("A,", "u,u"), f)
  expect_error(read_xrd_pairs(f), "two header rows")
  writeLines(c("A,,B", "u,u,u", "1,2,3", "4,5,6"), f)
  expect_error(read_xrd_pairs(f), "even number")
  writeLines(c(",,B,", "u,u,u,u", "1,2,3,4", "4,5,6,7"), f)
  expect_error(read_xrd_pairs(f), "sample name")
})
