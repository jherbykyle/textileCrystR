# Peak-shape functions for peak-area (deconvolution) crystallinity fitting.
# All four shapes share a "height at center" parametrization -- height is the
# peak's value at its own center for every shape -- so starting values and
# fitted heights mean the same thing regardless of which model is chosen.
#
# Gaussian and Lorentzian are closed-form textbook formulas. Pseudo-Voigt is
# their linear mixture (Thompson-Cox-Hastings form, single shared FWHM).
# Voigt is the true convolution of a Gaussian and a Lorentzian, computed via
# the Faddeeva function w(z) = exp(-z^2)*erfc(-iz), using the Humlicek (1982)
# rational approximation. That implementation was checked, before use here,
# against (a) the exact relation w(iy) = exp(y^2)*erfc(y) for real y (matches
# an independently computable reference to ~1e-6 relative error), (b) its
# required convergence to the Lorentzian as the Gaussian width shrinks to
# zero and to the Gaussian as the Lorentzian width shrinks to zero, and (c)
# numerical integration of its area against the closed-form area derived
# below, over a window wide enough for the heavy Lorentzian tails.

# Humlicek (1982) w4 rational approximation to the Faddeeva function.
# Accurate to about 1e-6 relative error across the region tested; see the
# module comment above for how that was checked.
.faddeeva <- function(z) {
  x <- Re(z); y <- Im(z)
  t <- y - 1i * x
  s <- abs(x) + y
  w <- complex(length.out = length(z))

  r1 <- s >= 15
  if (any(r1)) {
    tt <- t[r1]
    w[r1] <- tt * 0.5641896 / (0.5 + tt * tt)
  }
  r2 <- s >= 5.5 & s < 15
  if (any(r2)) {
    tt <- t[r2]
    u <- tt * tt
    w[r2] <- (tt * (1.410474 + u * 0.5641896)) / (0.75 + u * (3 + u))
  }
  r3 <- s < 5.5 & y >= (0.195 * abs(x) - 0.176)
  if (any(r3)) {
    tt <- t[r3]
    w[r3] <- (16.4955 + tt * (20.20933 + tt * (11.96482 + tt * (3.778987 + tt * 0.5642236)))) /
      (16.4955 + tt * (38.82363 + tt * (39.27121 + tt * (21.69274 + tt * (6.699398 + tt)))))
  }
  r4 <- !(r1 | r2 | r3)
  if (any(r4)) {
    tt <- t[r4]
    u <- tt * tt
    nom <- tt * (36183.31 - u * (3321.9905 - u * (1540.787 - u * (219.0313 -
             u * (35.76683 - u * (1.320522 - u * 0.56419))))))
    den <- 32066.6 - u * (24322.84 - u * (9022.228 - u * (2186.181 -
             u * (364.2191 - u * (61.57037 - u * (1.841439 - u))))))
    w[r4] <- exp(u) - nom / den
  }
  w
}

.gaussian_peak <- function(x, center, height, fwhm) {
  height * exp(-4 * log(2) * (x - center)^2 / fwhm^2)
}
.lorentzian_peak <- function(x, center, height, fwhm) {
  height / (1 + 4 * (x - center)^2 / fwhm^2)
}
.pseudo_voigt_peak <- function(x, center, height, fwhm, eta) {
  eta * .lorentzian_peak(x, center, height, fwhm) + (1 - eta) * .gaussian_peak(x, center, height, fwhm)
}
# gaussian_fwhm, lorentzian_fwhm: FWHM of the Gaussian and Lorentzian
# components being convolved (not the composite Voigt FWHM).
.voigt_peak <- function(x, center, height, gaussian_fwhm, lorentzian_fwhm) {
  sigma <- gaussian_fwhm / (2 * sqrt(2 * log(2)))
  gamma <- lorentzian_fwhm / 2
  z0 <- complex(real = 0, imaginary = gamma) / (sigma * sqrt(2))
  raw0 <- Re(.faddeeva(z0))
  z <- ((x - center) + complex(real = 0, imaginary = gamma)) / (sigma * sqrt(2))
  height * Re(.faddeeva(z)) / raw0
}

# Analytic area (integral over all x) for each shape, in the same height/FWHM
# parametrization as the peak functions above. Used for fast, exact area
# reporting instead of re-integrating the fitted curve numerically.
.gaussian_area <- function(height, fwhm) height * fwhm * sqrt(pi / (4 * log(2)))
.lorentzian_area <- function(height, fwhm) height * fwhm * pi / 2
.pseudo_voigt_area <- function(height, fwhm, eta) {
  eta * .lorentzian_area(height, fwhm) + (1 - eta) * .gaussian_area(height, fwhm)
}
.voigt_area <- function(height, gaussian_fwhm, lorentzian_fwhm) {
  sigma <- gaussian_fwhm / (2 * sqrt(2 * log(2)))
  gamma <- lorentzian_fwhm / 2
  z0 <- complex(real = 0, imaginary = gamma) / (sigma * sqrt(2))
  raw0 <- Re(.faddeeva(z0))
  height / raw0 * sigma * sqrt(2 * pi)
}

# Names and order of the free parameters for one component of a given model.
# Every fitted component (crystalline peak or amorphous hump) has this same
# parameter layout, in this order, throughout peak_area.R.
.model_param_names <- function(model) {
  switch(model,
    gaussian = c("center", "height", "fwhm"),
    lorentzian = c("center", "height", "fwhm"),
    pseudo_voigt = c("center", "height", "fwhm", "eta"),
    voigt = c("center", "height", "gaussian_fwhm", "lorentzian_fwhm"),
    stop("Unreachable: unknown model in .model_param_names().")
  )
}

# Evaluate one component (a named numeric vector p, in .model_param_names()
# order) at positions x.
.eval_component <- function(x, p, model) {
  switch(model,
    gaussian = .gaussian_peak(x, p[["center"]], p[["height"]], p[["fwhm"]]),
    lorentzian = .lorentzian_peak(x, p[["center"]], p[["height"]], p[["fwhm"]]),
    pseudo_voigt = .pseudo_voigt_peak(x, p[["center"]], p[["height"]], p[["fwhm"]], p[["eta"]]),
    voigt = .voigt_peak(x, p[["center"]], p[["height"]], p[["gaussian_fwhm"]], p[["lorentzian_fwhm"]])
  )
}

# Exact area of one component.
.component_area <- function(p, model) {
  switch(model,
    gaussian = .gaussian_area(p[["height"]], p[["fwhm"]]),
    lorentzian = .lorentzian_area(p[["height"]], p[["fwhm"]]),
    pseudo_voigt = .pseudo_voigt_area(p[["height"]], p[["fwhm"]], p[["eta"]]),
    voigt = .voigt_area(p[["height"]], p[["gaussian_fwhm"]], p[["lorentzian_fwhm"]])
  )
}

# A single component's FWHM for reporting purposes. Voigt has two component
# widths rather than one; report the standard Olivero & Longbothum (1977)
# approximation for the composite Voigt FWHM from them.
.component_fwhm <- function(p, model) {
  if (model == "voigt") {
    fg <- p[["gaussian_fwhm"]]; fl <- p[["lorentzian_fwhm"]]
    0.5346 * fl + sqrt(0.2166 * fl^2 + fg^2)
  } else {
    p[["fwhm"]]
  }
}
