test_that(".file_hint suggests similarly-typed sibling files when the folder exists", {
  d <- file.path(.fixture_dir(), "hint_test")
  dir.create(d, showWarnings = FALSE)
  writeLines("x", file.path(d, "pattern_a.csv"))
  writeLines("x", file.path(d, "pattern_b.txt"))
  writeLines("x", file.path(d, "readme.md"))  # not a recognized XRD extension
  hint <- .file_hint(file.path(d, "does_not_exist.csv"))
  expect_match(hint, "pattern_a.csv", fixed = TRUE)
  expect_match(hint, "pattern_b.txt", fixed = TRUE)
  expect_false(grepl("readme.md", hint, fixed = TRUE))
})

test_that(".file_hint returns an empty string when the folder has no similarly-typed files", {
  d <- file.path(.fixture_dir(), "hint_empty")
  dir.create(d, showWarnings = FALSE)
  writeLines("x", file.path(d, "unrelated.md"))
  expect_equal(.file_hint(file.path(d, "missing.csv")), "")
})

test_that(".file_hint returns an empty string when the folder itself does not exist", {
  expect_equal(.file_hint("/no/such/folder/at/all/missing.csv"), "")
})

test_that(".sniff_delim falls through to the whitespace candidate when no fixed delimiter fits", {
  lines <- c("10.0 100.0", "11.0 110.0", "12.0 120.0", "13.0 130.0")
  info <- .sniff_delim(lines)
  expect_equal(info$fixed, FALSE)
  expect_equal(info$delim, "\\s+")
  expect_equal(info$n_cols, 2L)
})

test_that(".sniff_delim returns NULL when nothing is consistent, even as whitespace", {
  lines <- c("no structure here", "still none", "nope")
  expect_null(.sniff_delim(lines))
})

test_that(".count_fields_fixed returns 0 on a genuinely unparseable line, and .sniff_delim skips it", {
  # An unterminated quote makes read.csv() error out for this line specifically.
  bad_line <- 'a,"unterminated'
  expect_equal(.count_fields_fixed(bad_line, ","), 0L)
  # .sniff_delim should still find comma from the other, well-formed lines.
  lines <- c(bad_line, "1,2", "3,4", "5,6", "7,8")
  info <- .sniff_delim(lines)
  expect_equal(info$delim, ",")
  expect_equal(info$n_cols, 2L)
})

test_that(".sniff_delim falls through to whitespace when every fixed delimiter fails on every line", {
  # An unescaped quote breaks read.csv() the same way regardless of which
  # separator is tried, so comma, tab and semicolon all fail on every line
  # here (not just one) -- this is the "every candidate entirely unusable"
  # case, distinct from the single-bad-line case above.
  lines <- sprintf('%d.0 "%s', 10:13, letters[1:4])
  info <- .sniff_delim(lines)
  expect_equal(info$fixed, FALSE)
  expect_equal(info$delim, "\\s+")
})

test_that(".looks_like_paired_layout requires at least 3 lines to consider the layout at all", {
  expect_false(.looks_like_paired_layout(c("A,,B,", "1,2,3,4"), ",", TRUE))
  expect_false(.looks_like_paired_layout(character(0), ",", TRUE))
})

test_that(".looks_like_paired_layout recognizes the real layout and rejects look-alikes", {
  paired <- c("A/1,,B/1,", "u,u,u,u", "1,2,3,4", "5,6,7,8")
  expect_true(.looks_like_paired_layout(paired, ",", TRUE))

  # ordinary 2-column file: header + data, no alternating blanks
  ordinary <- c("Angle,Counts", "10,100", "11,110", "12,120")
  expect_false(.looks_like_paired_layout(ordinary, ",", TRUE))

  # even column count but not alternating blank (both columns named)
  not_alternating <- c("A,B,C,D", "1,2,3,4", "5,6,7,8")
  expect_false(.looks_like_paired_layout(not_alternating, ",", TRUE))
})

test_that(".split_line preserves trailing empty fields for fixed delimiters", {
  expect_equal(.split_line("a,,b,", ",", TRUE), c("a", "", "b", ""))
  expect_equal(.split_line("10  20", "\\s+", FALSE), c("10", "20"))
})
