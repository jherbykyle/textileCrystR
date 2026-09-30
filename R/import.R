# A forgiving file reader for single- or multi-sample XRD text exports. This
# is deliberately more permissive than read_xrd_pairs(): it sniffs delimiter,
# header and column roles, and only ever falls back to read_xrd_pairs() for
# the specific two-header-row layout that function expects. It never guesses
# at ambiguous numbers silently -- when it cannot be sure, it stops with a
# plain-English explanation of what it saw and how to fix it.

.ANGLE_NAME_RX <- "(?i)^\\s*(2\\s*[-_.]?\\s*theta|two[-_. ]?theta|angle|2th(eta)?|theta)\\b"
.INTENSITY_NAME_RX <- "(?i)^\\s*(intensity|counts?|cps|y)\\b"

.file_hint <- function(path) {
  dir <- dirname(path)
  if (!dir.exists(dir)) return("")
  siblings <- list.files(dir, full.names = FALSE)
  siblings <- siblings[grepl("\\.(csv|txt|xy|dat|tsv)$", siblings, ignore.case = TRUE)]
  if (!length(siblings)) return("")
  sprintf(" Files of a similar type in that folder: %s.",
          paste(utils::head(siblings, 8), collapse = ", "))
}

# Try candidate delimiters and keep the one that yields the most columns
# consistently across the first few lines (ties broken by the order given,
# which favours comma, the most common plain-text export).
.count_fields_fixed <- function(line, delim) {
  tryCatch(
    ncol(utils::read.csv(text = line, header = FALSE, sep = delim, colClasses = "character",
                         na.strings = character())),
    error = function(e) 0L
  )
}

# Try candidate delimiters and keep the one that yields the most columns
# consistently across the first few lines (ties broken by the order given,
# which favours comma, the most common plain-text export). Field counts use
# read.csv-style parsing (via .count_fields_fixed / .split_line), not raw
# strsplit(), because strsplit() silently drops a trailing empty field after
# a final delimiter (e.g. "a,,b," has 3 fields to strsplit but 4 to read.csv),
# which would otherwise undercount files with a trailing separator.
.sniff_delim <- function(lines) {
  sample_lines <- lines[nzchar(trimws(lines))]
  sample_lines <- utils::head(sample_lines, 20L)
  candidates <- list(comma = ",", tab = "\t", semicolon = ";")
  best <- NULL
  for (nm in names(candidates)) {
    counts <- vapply(sample_lines, .count_fields_fixed, integer(1), delim = candidates[[nm]])
    counts <- counts[counts > 0L]
    if (!length(counts)) next
    tab <- table(counts)
    n_cols <- as.integer(names(tab)[which.max(tab)])
    consistency <- max(tab) / length(counts)
    if (n_cols >= 2L && consistency >= 0.8 &&
        (is.null(best) || n_cols > best$n_cols)) {
      best <- list(delim = candidates[[nm]], n_cols = n_cols, fixed = TRUE)
    }
  }
  if (!is.null(best)) return(best)
  # Whitespace-separated as a last resort (common for .xy / .dat exports).
  counts <- lengths(strsplit(trimws(sample_lines), "\\s+"))
  counts <- counts[counts > 0L]
  if (length(counts)) {
    tab <- table(counts)
    n_cols <- as.integer(names(tab)[which.max(tab)])
    consistency <- max(tab) / length(counts)
    if (n_cols >= 2L && consistency >= 0.8) {
      return(list(delim = "\\s+", n_cols = n_cols, fixed = FALSE))
    }
  }
  NULL
}

.split_line <- function(line, delim, fixed) {
  if (fixed) {
    row <- utils::read.csv(text = line, header = FALSE, sep = delim, colClasses = "character",
                           na.strings = character())
    trimws(as.character(unlist(row[1, ])))
  } else {
    trimws(strsplit(trimws(line), delim, fixed = fixed)[[1]])
  }
}

# Does this look like the two-header-row, sample-per-pair layout that
# read_xrd_pairs() expects (as in the shipped cotton.csv)? Its defining
# feature is alternating blank columns -- "Name1,,Name2,,Name3,..." -- not
# merely "non-numeric text in both of the first two rows", which an ordinary
# single-pattern file with a bad header or a garbage data row could also
# match. Requiring at least two sample-pairs (>= 4 columns) also sidesteps
# the genuinely ambiguous single-pair case, which is far more likely to be
# an ordinary two-column file.
.looks_like_paired_layout <- function(lines, delim, fixed) {
  if (length(lines) < 3L) return(FALSE)
  f1 <- .split_line(lines[1], delim, fixed)
  f2 <- .split_line(lines[2], delim, fixed)
  if (length(f1) < 4L || length(f1) %% 2L != 0L) return(FALSE)
  n1 <- suppressWarnings(as.numeric(f1))
  n2 <- suppressWarnings(as.numeric(f2))
  name_like <- mean(is.na(n1[nzchar(f1)])) > 0.5 && mean(is.na(n2[nzchar(f2)])) > 0.5
  even_idx <- seq(2L, length(f1), by = 2L)
  alternating_blank <- mean(!nzchar(f1[even_idx])) > 0.5
  name_like && alternating_blank
}

# Decide whether the first line is a header (non-numeric fields) and, if so,
# try to name the two_theta and intensity columns from it.
.guess_columns <- function(lines, delim, fixed, cols) {
  first_fields <- .split_line(lines[1], delim, fixed)
  first_numeric <- suppressWarnings(as.numeric(first_fields))
  has_header <- any(is.na(first_numeric) & nzchar(first_fields))
  header <- if (has_header) first_fields else NULL
  skip <- if (has_header) 1L else 0L

  if (!is.null(cols)) {
    if (length(cols) != 2L) {
      stop("`cols` must name exactly two columns: c(two_theta_column, intensity_column).",
           call. = FALSE)
    }
    resolve_one <- function(col, what) {
      if (is.character(col)) {
        if (is.null(header)) {
          stop(sprintf(paste0("`cols` names a column ('%s') by name for %s, but the file has no ",
                              "header row to match against. Use column numbers instead, e.g. ",
                              "cols = c(1, 2)."), col, what), call. = FALSE)
        }
        j <- match(col, header)
        if (is.na(j)) {
          stop(sprintf(paste0("Column '%s' (requested for %s) was not found in the file's header. ",
                              "Header columns found: %s."),
                       col, what, paste(header, collapse = ", ")), call. = FALSE)
        }
        j
      } else {
        as.integer(col)
      }
    }
    return(list(two_theta_col = resolve_one(cols[1], "two_theta"),
                intensity_col = resolve_one(cols[2], "intensity"),
                skip = skip, header = header, guessed = FALSE))
  }

  if (!is.null(header)) {
    tcol <- which(grepl(.ANGLE_NAME_RX, header, perl = TRUE))[1]
    icol <- which(grepl(.INTENSITY_NAME_RX, header, perl = TRUE))[1]
    if (!is.na(tcol) && !is.na(icol) && tcol != icol) {
      return(list(two_theta_col = tcol, intensity_col = icol, skip = skip, header = header,
                 guessed = FALSE))
    }
  }
  list(two_theta_col = 1L, intensity_col = 2L, skip = skip, header = header, guessed = TRUE)
}

#' Read an XRD text file, auto-detecting its layout
#'
#' A forgiving reader for XRD exports as CSV, TSV or whitespace-delimited
#' plain text: it detects the delimiter, whether the first row is a header,
#' and (from a header, if present) which column is 2-theta and which is
#' intensity. Files in the shipped two-header-row layout (one sample name per
#' pair of columns, as in `cotton.csv`) are detected and handed to
#' [read_xrd_pairs()] unchanged. Nothing about the numbers themselves is
#' altered; only column and row selection is inferred.
#'
#' @param path Path to a `.csv`, `.txt`, `.tsv`, `.dat` or `.xy` file.
#' @param cols Optional length-2 vector naming the 2-theta and intensity
#'   columns explicitly, as `c(two_theta, intensity)`, either both column
#'   numbers or (if the file has a header) both column names. Use this when
#'   auto-detection picks the wrong columns or the file has no clear header.
#' @param sample Optional sample name to use for a single-pattern file
#'   (default: the file name without its extension). Ignored for files in the
#'   multi-sample layout, which supply their own sample names.
#' @return A [tibble::tibble()] with columns `sample`, `two_theta` and
#'   `intensity`, one or more samples in file order.
#' @examples
#' path <- system.file("extdata", "cotton.csv", package = "textileCrystR")
#' xrd <- import_xrd_file(path)
#' unique(xrd$sample)
#' @export
import_xrd_file <- function(path, cols = NULL, sample = NULL) {
  if (!is.character(path) || length(path) != 1L) {
    stop("`path` must be a single file path (as a character string).", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop(sprintf("The file '%s' does not exist.%s", path, .file_hint(path)), call. = FALSE)
  }
  if (file.info(path)$isdir) {
    stop(sprintf("'%s' is a folder, not a file. Point `path` at the XRD file itself.", path),
         call. = FALSE)
  }
  lines <- readLines(path, warn = FALSE)
  lines[1] <- sub("^\xef\xbb\xbf", "", lines[1], useBytes = TRUE)
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) < 2L) {
    stop(sprintf(paste0("'%s' has fewer than two non-blank lines, so there is no data to read. ",
                        "Check that the file was saved/exported correctly."), path), call. = FALSE)
  }

  delim_info <- .sniff_delim(lines)
  if (is.null(delim_info)) {
    preview <- paste(utils::head(lines, 3), collapse = "\n  ")
    stop(sprintf(paste0(
      "Could not work out how '%s' is separated into columns (tried comma, tab, semicolon and ",
      "whitespace, but none gave a consistent number of columns across the first lines). The ",
      "first lines of the file look like:\n  %s\nIf this is really tabular data, try re-saving it ",
      "as a plain comma-separated (.csv) file."), path, preview), call. = FALSE)
  }

  if (is.null(cols) && .looks_like_paired_layout(lines, delim_info$delim, delim_info$fixed)) {
    out <- tryCatch(
      read_xrd_pairs(path),
      error = function(e) {
        stop(sprintf(paste0(
          "'%s' looks like a multi-sample file (one sample name per pair of columns), but it could ",
          "not be read that way: %s\nIf this is really a single XRD pattern, save it as a plain ",
          "two-column (2-theta, intensity) CSV instead."), path, conditionMessage(e)), call. = FALSE)
      }
    )
    return(out)
  }

  col_info <- .guess_columns(lines, delim_info$delim, delim_info$fixed, cols)
  raw <- utils::read.table(text = paste(lines, collapse = "\n"), sep = if (delim_info$fixed) delim_info$delim else "",
                           header = FALSE, skip = col_info$skip, fill = TRUE,
                           stringsAsFactors = FALSE, strip.white = TRUE)
  n_cols_needed <- max(col_info$two_theta_col, col_info$intensity_col)
  if (ncol(raw) < n_cols_needed) {
    stop(sprintf(paste0(
      "'%s' only has %d column(s) after splitting on '%s', but column %d was requested. Pass ",
      "`cols = c(two_theta_column, intensity_column)` to say explicitly which columns to use."),
      path, ncol(raw), if (delim_info$fixed) delim_info$delim else "whitespace", n_cols_needed),
      call. = FALSE)
  }
  two_theta <- suppressWarnings(as.numeric(raw[[col_info$two_theta_col]]))
  intensity <- suppressWarnings(as.numeric(raw[[col_info$intensity_col]]))

  if (col_info$guessed) {
    # No confident column choice was made (no cols= given, and either no
    # header or a header whose names didn't match known angle/intensity
    # terms): sanity-check the first-two-columns fallback before committing
    # to it.
    tt_ok <- mean(is.finite(two_theta) & two_theta > 0 & two_theta < 180, na.rm = TRUE) > 0.9
    if (!tt_ok) {
      preview <- utils::head(raw, 3)
      header_note <- if (is.null(col_info$header)) {
        "the file has no header row, so"
      } else {
        sprintf("the header (%s) did not clearly name a 2-theta/angle column, so",
                paste(col_info$header, collapse = ", "))
      }
      stop(sprintf(paste0(
        "In '%s', %s column 1 was assumed to be 2-theta and column 2 intensity -- but column 1 ",
        "does not look like a 2-theta axis (values should mostly lie between 0 and 180 degrees). ",
        "The first rows read as:\n%s\nPass `cols = c(two_theta_column, intensity_column)` (by ",
        "column number or, if there's a header, by name) to tell import_xrd_file() which columns ",
        "to use."),
        path, header_note, paste(utils::capture.output(print(preview)), collapse = "\n")),
        call. = FALSE)
    }
  }
  if (all(is.na(two_theta)) || all(is.na(intensity))) {
    stop(sprintf(paste0(
      "Column %d (2-theta) or column %d (intensity) in '%s' contains no numeric values at all. ",
      "Check that `cols` (if given) points at the right columns, and that the file does not have ",
      "extra header or footer lines that need skipping."),
      col_info$two_theta_col, col_info$intensity_col, path), call. = FALSE)
  }

  keep <- !is.na(two_theta) & !is.na(intensity)
  sample <- sample %||% tools::file_path_sans_ext(basename(path))
  tibble::tibble(sample = sample, two_theta = two_theta[keep], intensity = intensity[keep])
}
