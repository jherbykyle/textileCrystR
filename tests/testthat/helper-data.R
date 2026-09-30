cotton_path <- function() {
  system.file("extdata", "cotton.csv", package = "textileCrystR", mustWork = TRUE)
}

# Small pattern with exactly prescribed intensities at 18 and 22.6 degrees.
mini_pattern <- function(i200 = 1000, iam = 250) {
  data.frame(two_theta = c(5, 10, 18, 22.6, 30, 40),
             intensity = c(100, 200, iam, i200, 300, 100))
}

# Smooth synthetic pattern (cellulose-I-like) on a regular grid.
synthetic_xrd <- function(step = 0.05) {
  tt <- seq(5, 40, by = step)
  y <- 300 + 4000 * exp(-((tt - 22.7) / 1.1)^2) + 1500 * exp(-((tt - 16) / 1.5)^2)
  data.frame(two_theta = tt, intensity = y)
}

# Expected Segal indices (2 d.p.) quoted in the project brief for the cotton series.
cotton_expected_2dp <- c("0%/1" = 71.19, "20%/1" = 67.65, "40%/1" = 69.36, "50%/1" = 69.24,
                         "60%/1" = 65.65, "80%/1" = 62.92, "100%/1" = 56.96)
