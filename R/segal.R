#' textileCrystR: Segal XRD crystallinity index for textile fibres and cellulose
#'
#' `textileCrystR` implements one method, the **Segal peak-height X-ray
#' diffraction crystallinity index**, in a transparent and reproducible way:
#' positions and intensities used are always retained, interpolation is explicit,
#' parameter choices can be stress-tested with [segal_sensitivity()], and every
#' result is labelled the *Segal Crystallinity Index*, never an absolute
#' crystallinity.
#'
#' @section Main functions:
#' * [segal_ci()] calculates the index and returns a `textile_cryst` object with
#'   `print()`, `summary()`, `plot()` and `as.data.frame()` methods.
#' * [segal_points()] returns only the intensities read at the two positions.
#' * [segal_presets()] lists literature-supported position pairs.
#' * [segal_sensitivity()] evaluates a grid of I200 and Iam positions.
#' * [validate_crystallinity_input()] checks a pattern before use.
#' * [read_xrd_pairs()] reads paired-column XRD files such as the shipped
#'   `cotton.csv`.
#'
#' @section Scientific limitations:
#' The Segal index is an **empirical, relative** XRD index. It was introduced to
#' compare cotton cellulose samples treated in different ways on the same
#' instrument (Segal et al., 1959) and is best used for relative changes within
#' one consistently measured sample set (Salem et al., 2023). It is not an
#' absolute crystalline mass fraction. In particular:
#'
#' * **Peak overlap.** The method treats all intensity at the amorphous
#'   position as amorphous scatter, but the broad neighbouring crystalline
#'   peaks of small cellulose crystallites also contribute there. Simulated
#'   patterns of perfect, cotton-sized crystals give a Segal index well below
#'   100 % (about 90 % in Salem et al., 2023), and real samples never reach 0 %.
#' * **Crystallite size / FWHM.** Broader peaks (smaller crystallites) raise the
#'   intensity at the amorphous position and change the index even when the
#'   amorphous fraction is identical (Nam et al., 2016).
#' * **Cellulose polymorph.** The index depends on the polymorph. For cellulose
#'   II the amorphous position used in the literature (16 degrees) misses most
#'   of the amorphous scattering, which is centred near the (110) reflection, so
#'   the amorphous fraction is underestimated (Nam et al., 2016). Values for
#'   different polymorphs are not directly comparable.
#' * **Preferred orientation.** Non-random crystallite orientation changes the
#'   I200 : Iam ratio at identical crystallinity; the (012)/(102) reflections
#'   near 20.3-20.6 degrees, close to the amorphous position, are suppressed or
#'   exposed depending on sample form and geometry (Nam et al., 2016; Salem et
#'   al., 2023).
#' * **Background handling.** Segal originally recommended no background
#'   correction; Salem et al. (2023) prefer subtracting a blank run, and note
#'   that subtraction routines that pull minima to zero are unsuited to
#'   cellulose. Indices calculated with different background treatment are not
#'   comparable. `textileCrystR` applies **no** background correction: whatever
#'   is in the intensity vector is used.
#' * **Choice of I200 and Iam positions.** The result depends on the positions
#'   chosen (about 22-23 degrees for the (200) maximum; about 18 degrees, with
#'   minima reported closer to 18.6 degrees in modern work). Use
#'   [segal_sensitivity()] to quantify this dependence and report the positions
#'   with every value. Radiation is assumed to be Cu K-alpha unless positions
#'   are converted.
#' * **Sample preparation and measurement.** Sample form (pressed powder,
#'   nonwoven, fabric), sample spinning and diffractometer geometry can shift the
#'   index by large amounts (Nam et al., 2016; Salem et al., 2023).
#' * **Partially mercerized cotton.** In blends of cellulose I and II, the
#'   (1-10) and (110) reflections of cellulose I remain near the Segal amorphous
#'   position of cellulose II, so the extended cellulose II Segal position is not
#'   appropriate for samples that are not completely mercerized (Nam et al.,
#'   2016). A single fixed position pair applied across a blend series
#'   measures a peak-height ratio, not the crystalline fraction of either
#'   polymorph.
#'
#' @section Scope:
#' Segal peak height and peak area / deconvolution (see
#' [fit_crystalline_peaks()]) are both implemented. Amorphous subtraction and
#' Rietveld-result integration, two further XRD approaches reviewed by Salem
#' et al. (2023), are intentionally not included -- the package's scope is
#' kept to these two well-validated methods rather than growing indefinitely.
#'
#' @references
#' Segal, L., Creely, J. J., Martin, A. E., Conrad, C. M. (1959). An empirical
#' method for estimating the degree of crystallinity of native cellulose using
#' the X-ray diffractometer. *Textile Research Journal* 29, 786-794.
#'
#' Nam, S., French, A. D., Condon, B. D., Concha, M. (2016). Segal
#' crystallinity index revisited by the simulation of X-ray diffraction patterns
#' of cotton cellulose I-beta and cellulose II. *Carbohydrate Polymers* 135,
#' 1-9. \doi{10.1016/j.carbpol.2015.08.035}
#'
#' Salem, K. S., Kasera, N. K., Rahman, M. A., Jameel, H., Habibi, Y., Eichhorn,
#' S. J., French, A. D., Pal, L., Lucia, L. A. (2023). Comparison and assessment
#' of methods for cellulose crystallinity determination. *Chemical Society
#' Reviews*. \doi{10.1039/d2cs00569g}
#'
#' @importFrom ggplot2 .data
#' @keywords internal
"_PACKAGE"

# ---- presets --------------------------------------------------------------------

#' Literature-supported Segal position presets
#'
#' A small table of I200 / Iam position pairs, each tied to a source in the
#' supplied literature. Positions are **not** universal rules: every preset
#' records its material, position or range, reference, notes and limitations,
#' and results should always be reported together with the positions used.
#'
#' @section Presets:
#' * `cellulose_I_segal`: cellulose I (native cotton), I200 at 22.7 degrees and
#'   Iam at 18.0 degrees. This is the Segal (1959) definition as stated by Nam
#'   et al. (2016) and Salem et al. (2023).
#' * `cellulose_II_extended`: the extension of the Segal method to cellulose II
#'   described by Nam et al. (2016), with the (020) reflection at 21.7 degrees
#'   and Iam at 16.0 degrees. For this preset the `i200` argument of
#'   [segal_ci()] holds the total intensity of the (020) reflection.
#' * `cotton_validation`: 22.6 and 18.0 degrees, the convention used for the
#'   experimental cotton validation in this package. 22.6 degrees is the (200)
#'   position of the simulated perfect cellulose I-beta crystal in Nam et al.
#'   (2016, Fig. 8) and lies inside the 22-23 degree range of Segal's method; it
#'   is a study convention, not a published Segal position.
#'
#' @param preset Optional preset name. If `NULL` (default) all presets are
#'   returned.
#' @return A [tibble::tibble()] with columns `preset`, `material`, `polymorph`,
#'   `i200`, `iam` (degrees 2-theta), `i200_range`, `iam_range`, `radiation`,
#'   `reference`, `notes` and `limitations`.
#' @examples
#' segal_presets()
#' segal_presets("cellulose_I_segal")$reference
#' @export
segal_presets <- function(preset = NULL) {
  tab <- tibble::tibble(
    preset = c("cellulose_I_segal", "cellulose_II_extended", "cotton_validation"),
    material = c("Native cellulose I (cotton)",
                 "Mercerized / regenerated cellulose II",
                 "Cotton, raw and mercerized blends (package validation convention)"),
    polymorph = c("cellulose I", "cellulose II", "cellulose I positions"),
    i200 = c(22.7, 21.7, 22.6),
    iam = c(18.0, 16.0, 18.0),
    i200_range = c(
      "Maximum of the (200) peak between 22 and 23 degrees 2-theta",
      "Point value of the (020) reflection; no range given",
      "22.2-23.0 degrees (robustness grid, not a published range)"),
    iam_range = c(
      "About 18 degrees; the minimum is closer to 18.6 degrees in modern work",
      "Point value; no range given",
      "17.6-18.4 degrees (robustness grid, not a published range)"),
    radiation = "Cu K-alpha",
    reference = c(
      paste0("Segal et al. (1959) Text. Res. J. 29, 786-794, as defined in Nam et al. (2016) ",
             "Carbohydr. Polym. 135, 1-9 (Sec. 1 and Sec. 2.4, Eq. 2) and Salem et al. (2023) ",
             "Chem. Soc. Rev. (Sec. 2.1.1, Eq. 1)"),
      paste0("Nam et al. (2016) Carbohydr. Polym. 135, 1-9 (Sec. 1 and Sec. 2.4), ",
             "extension of the Segal method to cellulose II"),
      paste0("Package validation convention; 22.6 degrees is the simulated (200) position ",
             "in Nam et al. (2016) Fig. 8; 18 degrees follows Segal et al. (1959) as ",
             "described by Nam et al. (2016) and Salem et al. (2023)")),
    notes = c(
      paste0("Original definition for native cotton cellulose. Salem et al. (2023) stress ",
             "that the amorphous value is the trough between the (110) and (200) peaks, not ",
             "the (110) peak height."),
      paste0("Uses 21.7 degrees for the (020) reflection and 16 degrees for Iam. Values are ",
             "not comparable with cellulose I Segal indices."),
      paste0("Same positions applied to every pattern of a blend series so the indices are ",
             "comparable with each other.")),
    limitations = c(
      paste0("Empirical relative index; assumes all intensity at the amorphous position is ",
             "amorphous and ignores peak overlap, crystallite size, preferred orientation ",
             "and background treatment."),
      paste0("Underestimates the amorphous fraction because the amorphous scatter of ",
             "cellulose II is centred near the (110) reflection; Nam et al. (2016) judge the ",
             "extended method inappropriate for partially mercerized samples."),
      paste0("Cellulose I positions applied to samples that contain cellulose II: the result ",
             "is a peak-height ratio at fixed positions, not the crystalline fraction of ",
             "either polymorph. Positions are a convention, not a peak-selection rule."))
  )
  if (is.null(preset)) return(tab)
  if (!is.character(preset) || length(preset) != 1L || !preset %in% tab$preset) {
    stop(sprintf("Unknown preset '%s'. Available presets: %s.", format(preset),
                 paste(tab$preset, collapse = ", ")), call. = FALSE)
  }
  tab[tab$preset == preset, , drop = FALSE]
}

# ---- internal helpers -------------------------------------------------------------

# Combine explicit positions with an optional preset. Explicit arguments win.
.resolve_positions <- function(i200, iam, preset) {
  info <- NULL
  overridden <- FALSE
  if (!is.null(preset)) {
    info <- segal_presets(preset)
    if (is.null(i200)) i200 <- info$i200 else overridden <- overridden || !isTRUE(all.equal(i200, info$i200))
    if (is.null(iam)) iam <- info$iam else overridden <- overridden || !isTRUE(all.equal(iam, info$iam))
  }
  if (is.null(i200) || is.null(iam)) {
    stop(paste0("Give both `i200` and `iam` (degrees 2-theta), or a `preset` ",
                "from segal_presets()."), call. = FALSE)
  }
  .check_position(i200, "i200")
  .check_position(iam, "iam")
  if (isTRUE(all.equal(i200, iam))) {
    stop("`i200` and `iam` must be different positions.", call. = FALSE)
  }
  list(i200 = i200, iam = iam, preset = info, overridden = overridden)
}

# Intensity at one requested position of a validated pattern.
.locate <- function(tt, y, pos, interpolation, label = "position") {
  n <- length(tt)
  if (pos < tt[1] - .POS_TOL || pos > tt[n] + .POS_TOL) {
    stop(sprintf(paste0("%s = %s is outside the measured 2-theta range (%s to %s); ",
                        "extrapolation is not performed."),
                 label, format(pos), format(tt[1]), format(tt[n])), call. = FALSE)
  }
  d <- abs(tt - pos)
  j <- which.min(d)
  on_grid <- d[j] <= .POS_TOL
  if (on_grid) {
    lower <- upper <- tt[j]
    value <- y[j]
    used <- if (interpolation == "nearest") tt[j] else pos
  } else {
    k <- min(findInterval(pos, tt), n - 1L)
    lower <- tt[k]
    upper <- tt[k + 1L]
    if (interpolation == "linear") {
      value <- y[k] + (y[k + 1L] - y[k]) * (pos - lower) / (upper - lower)
      used <- pos
    } else {
      jn <- if (pos - lower <= upper - pos) k else k + 1L
      value <- y[jn]
      used <- tt[jn]
    }
  }
  list(position_requested = pos, position_used = used, intensity = value,
       lower = lower, upper = upper, on_measured_point = on_grid,
       gap = upper - lower)
}

# Core calculation on an already validated pattern.
.segal_eval <- function(input, i200, iam, interpolation, preset_info = NULL) {
  top <- .locate(input$two_theta, input$intensity, i200, interpolation, "i200")
  bot <- .locate(input$two_theta, input$intensity, iam, interpolation, "iam")
  if (top$intensity <= 0) {
    stop(sprintf(paste0("The I200 intensity is not positive (%s at %s degrees), so the Segal ",
                        "index is undefined."), format(signif(top$intensity, 4)),
                 format(top$position_used)), call. = FALSE)
  }
  ci <- (top$intensity - bot$intensity) / top$intensity * 100

  warn <- character()
  if (i200 < iam) {
    warn <- c(warn, "`i200` is at a lower angle than `iam`; check that the positions are not swapped.")
  }
  if (bot$intensity < 0) {
    warn <- c(warn, sprintf(paste0("The Iam intensity is negative (%s), so the index exceeds ",
                                   "100 %%. Check background handling."),
                            format(signif(bot$intensity, 4))))
  }
  if (bot$intensity > top$intensity) {
    warn <- c(warn, paste0("Iam exceeds I200, so the index is negative. Check the positions ",
                           "and the sample."))
  }
  for (p in list(list(top, "I200"), list(bot, "Iam"))) {
    if (!p[[1]]$on_measured_point && p[[1]]$gap > .MAX_GAP) {
      warn <- c(warn, sprintf(paste0("%s was %s across a %.2f degree gap between measured points; ",
                                     "the value may be unreliable."),
                              p[[2]], if (interpolation == "linear") "interpolated" else "read",
                              p[[1]]$gap))
    }
  }
  points <- tibble::tibble(
    point = c("I200", "Iam"),
    position_requested = c(top$position_requested, bot$position_requested),
    position_used = c(top$position_used, bot$position_used),
    intensity = c(top$intensity, bot$intensity),
    lower_2theta = c(top$lower, bot$lower),
    upper_2theta = c(top$upper, bot$upper),
    on_measured_point = c(top$on_measured_point, bot$on_measured_point)
  )
  list(ci = ci, points = points, warnings = warn)
}

.check_interpolation <- function(interpolation) {
  match.arg(interpolation, c("linear", "nearest"))
}

.check_sample <- function(sample) {
  if (is.null(sample)) return(NA_character_)
  if (!is.character(sample) || length(sample) != 1L || is.na(sample)) {
    stop("`sample` must be a single character string.", call. = FALSE)
  }
  sample
}

# ---- user-facing functions ----------------------------------------------------------

#' Segal Crystallinity Index
#'
#' Calculates the Segal peak-height crystallinity index
#' \deqn{CI = \frac{I_{200} - I_{am}}{I_{200}} \times 100}{CI = (I200 - Iam) / I200 * 100}
#' from a measured X-ray diffraction pattern. `I200` is the intensity of the
#' (200) reflection of cellulose I and `Iam` the intensity at the amorphous
#' position, both read **at the positions you give**. If a position was not
#' measured exactly, the intensity is interpolated (never extrapolated).
#'
#' The result is the *Segal Crystallinity Index*: an empirical, relative index,
#' **not** an absolute crystallinity. See the section *Limitations* and
#' [textileCrystR-package].
#'
#' @section Limitations:
#' The Segal index depends on peak overlap, crystallite size (FWHM), cellulose
#' polymorph, preferred orientation, background handling, the positions chosen
#' for I200 and Iam, sample preparation, and is not appropriate for partially
#' mercerized cotton when the cellulose II positions are used. No background
#' correction is applied by this function. Details and references are given in
#' [textileCrystR-package].
#'
#' @param two_theta Numeric vector of diffraction angles (degrees 2-theta).
#' @param intensity Numeric vector of intensities, same length as `two_theta`.
#' @param i200 Position (degrees 2-theta) at which I200 is read. May be omitted
#'   when `preset` is given.
#' @param iam Position (degrees 2-theta) at which Iam is read. May be omitted
#'   when `preset` is given.
#' @param preset Optional name of a preset from [segal_presets()]. Explicit
#'   `i200` / `iam` values override the preset's positions.
#' @param interpolation `"linear"` (default) reads the intensity by linear
#'   interpolation between the two neighbouring measured points; `"nearest"`
#'   uses the closest measured point. A requested position that coincides with a
#'   measured point (within 1e-6 degrees) uses that measured intensity exactly.
#' @param sample Optional sample name stored in the result.
#' @param na.rm,duplicates Passed to [validate_crystallinity_input()].
#'
#' @return An object of class `textile_cryst`, a list with
#' * `ci`: the Segal Crystallinity Index (percent),
#' * `points`: a tibble with the requested and used positions, intensities and
#'   neighbouring measured 2-theta values for `I200` and `Iam`,
#' * `method`: name, equation, interpolation, preset and radiation,
#' * `sample`, `data` (the validated pattern used), `input` (validation record)
#'   and `warnings`.
#'
#' It has `print()`, `summary()`, `plot()` and `as.data.frame()` methods.
#' @seealso [segal_points()], [segal_sensitivity()], [segal_presets()]
#' @examples
#' tt <- seq(5, 40, by = 0.05)
#' y <- 300 + 4000 * exp(-((tt - 22.7) / 1.1)^2) + 1500 * exp(-((tt - 16) / 1.5)^2)
#' res <- segal_ci(tt, y, i200 = 22.6, iam = 18.0, sample = "synthetic example")
#' res
#' summary(res)
#' as.data.frame(res)
#' @export
segal_ci <- function(two_theta, intensity, i200 = NULL, iam = NULL, preset = NULL,
                     interpolation = c("linear", "nearest"), sample = NULL,
                     na.rm = FALSE, duplicates = c("error", "mean")) {
  interpolation <- .check_interpolation(interpolation)
  sample <- .check_sample(sample)
  pos <- .resolve_positions(i200, iam, preset)
  input <- validate_crystallinity_input(two_theta, intensity, na.rm = na.rm,
                                        duplicates = match.arg(duplicates))
  ev <- .segal_eval(input, pos$i200, pos$iam, interpolation)
  if (pos$overridden) {
    ev$warnings <- c(ev$warnings, sprintf(
      "Explicit positions differ from preset '%s'; the positions actually used are reported.", preset))
  }
  structure(list(
    ci = ev$ci,
    points = ev$points,
    method = list(name = "Segal peak height",
                  equation = "CI = (I200 - Iam) / I200 * 100",
                  interpolation = interpolation,
                  preset = preset,
                  preset_overridden = pos$overridden,
                  preset_reference = if (is.null(pos$preset)) NA_character_ else pos$preset$reference,
                  radiation = "Cu K-alpha positions assumed"),
    sample = sample,
    data = tibble::tibble(two_theta = input$two_theta, intensity = input$intensity),
    input = input[setdiff(names(input), c("two_theta", "intensity"))],
    warnings = ev$warnings
  ), class = "textile_cryst")
}

#' Intensities read at the Segal positions
#'
#' Returns the four numbers that define a Segal calculation (position and
#' intensity for I200 and Iam) together with the interpolation method and basic
#' facts about the input pattern, without computing the index. Requested
#' positions need not exist in the data.
#'
#' @inheritParams segal_ci
#' @return A one-row [tibble::tibble()] with `sample`, `i200_position`,
#'   `i200_intensity`, `iam_position`, `iam_intensity`, `interpolation`,
#'   `preset`, whether each position `..._on_measured_point`, the measured
#'   neighbours `..._lower` / `..._upper`, and input information (`n_input`,
#'   `n_used`, `two_theta_min`, `two_theta_max`, `median_step`).
#' @examples
#' tt <- seq(5, 40, by = 0.05)
#' y <- 300 + 4000 * exp(-((tt - 22.7) / 1.1)^2)
#' segal_points(tt, y, i200 = 22.62, iam = 18.01)
#' @export
segal_points <- function(two_theta, intensity, i200 = NULL, iam = NULL, preset = NULL,
                         interpolation = c("linear", "nearest"), sample = NULL,
                         na.rm = FALSE, duplicates = c("error", "mean")) {
  interpolation <- .check_interpolation(interpolation)
  sample <- .check_sample(sample)
  pos <- .resolve_positions(i200, iam, preset)
  input <- validate_crystallinity_input(two_theta, intensity, na.rm = na.rm,
                                        duplicates = match.arg(duplicates))
  top <- .locate(input$two_theta, input$intensity, pos$i200, interpolation, "i200")
  bot <- .locate(input$two_theta, input$intensity, pos$iam, interpolation, "iam")
  tibble::tibble(
    sample = sample,
    i200_position = top$position_used, i200_intensity = top$intensity,
    iam_position = bot$position_used, iam_intensity = bot$intensity,
    interpolation = interpolation,
    preset = if (is.null(preset)) NA_character_ else preset,
    i200_on_measured_point = top$on_measured_point,
    iam_on_measured_point = bot$on_measured_point,
    i200_lower = top$lower, i200_upper = top$upper,
    iam_lower = bot$lower, iam_upper = bot$upper,
    n_input = input$n_input, n_used = input$n_used,
    two_theta_min = input$range[1], two_theta_max = input$range[2],
    median_step = input$median_step
  )
}
