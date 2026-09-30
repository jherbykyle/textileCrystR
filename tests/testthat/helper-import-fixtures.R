# Small on-disk fixtures for import_xrd_file() tests. Written fresh into a
# per-session temp dir so tests never depend on files left by earlier runs.

.fixture_dir <- function() {
  d <- file.path(tempdir(), "textileCrystR-import-fixtures")
  dir.create(d, showWarnings = FALSE)
  d
}

.write_fixture <- function(name, lines) {
  path <- file.path(.fixture_dir(), name)
  writeLines(lines, path)
  path
}

simple_pattern <- function() {
  tt <- seq(5, 40, by = 0.1)
  y <- 100 + 800 * exp(-((tt - 22.6) / 1.1)^2) + 250 * exp(-((tt - 18) / 1.5)^2)
  list(two_theta = tt, intensity = round(y, 2))
}
