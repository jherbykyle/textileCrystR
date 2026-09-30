# Small, dependency-free numeric helpers. These are simple enough (a
# trapezoidal rule, R-squared, RMSE) that writing them directly is clearer
# and lighter than adding a dependency for them.

# Trapezoidal integral of y over x (x need not be evenly spaced; must be
# sorted ascending, as validate_crystallinity_input() already guarantees for
# anything this package passes in).
.trapz <- function(x, y) {
  n <- length(x)
  if (n < 2L) return(0)
  sum((x[-1] - x[-n]) * (y[-1] + y[-n]) / 2)
}

# Coefficient of determination between observed and fitted values.
.r_squared <- function(observed, fitted) {
  ss_res <- sum((observed - fitted)^2)
  ss_tot <- sum((observed - mean(observed))^2)
  if (ss_tot <= 0) return(NA_real_)
  1 - ss_res / ss_tot
}

.rmse <- function(observed, fitted) {
  sqrt(mean((observed - fitted)^2))
}
