#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------------
# Experimental cotton validation of textileCrystR 0.1.0
#
#   Rscript scripts/validate_cotton.R [output_dir]
#
# Run from the package root (or with the package installed). The raw file
# inst/extdata/cotton.csv is only read, never written or modified. Nothing is
# hard-coded: every Segal index is calculated from the raw patterns.
#
# Steps
#   1. load cotton.csv                        5. composition-vs-CI plot
#   2. validate the data                      6. sensitivity analysis (5 x 5 grid)
#   3. process all seven patterns             7. comparison with the independent
#   4. tidy results table                        (spreadsheet) reference calculation
# Gravimetric composition is a label of the sample series. It is NOT a ground-truth
# crystallinity and is never used in the calculation.
# ---------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1L) args[[1]] else file.path("scripts", "validation_output")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

if (requireNamespace("textileCrystR", quietly = TRUE)) {
  suppressPackageStartupMessages(library(textileCrystR))
  extdata <- function(f) system.file("extdata", f, package = "textileCrystR", mustWork = TRUE)
} else {
  pkgload::load_all(".", quiet = TRUE)
  extdata <- function(f) file.path("inst", "extdata", f)
}

I200 <- 22.6    # degrees 2-theta, package validation convention (see segal_presets())
IAM <- 18.0
TOL <- 1e-6     # agreement threshold with the reference, percentage points
RAW_MD5 <- "cf33ee9c15cd55f1638c093ee28b804a"

say <- function(...) cat(sprintf(...), "\n", sep = "")
status <- character()

# ---- 1. load ------------------------------------------------------------------------
say("== 1. Load cotton.csv")
path <- extdata("cotton.csv")
md5 <- unname(tools::md5sum(path))
say("   file: %s\n   md5 : %s (%s)", path, md5,
    if (md5 == RAW_MD5) "matches the supplied raw file" else "DIFFERS from the supplied raw file")
status["raw file unchanged"] <- if (md5 == RAW_MD5) "PASS" else "FAIL"
xrd <- read_xrd_pairs(path)
samples <- unique(xrd$sample)
composition <- as.numeric(sub("%.*$", "", samples))
samples <- samples[order(composition)]
say("   %d patterns, %d rows in total", length(samples), nrow(xrd))

# ---- 2. validate --------------------------------------------------------------------
say("== 2. Validate the data")
checks <- do.call(rbind, lapply(samples, function(s) {
  d <- xrd[xrd$sample == s, ]
  v <- validate_crystallinity_input(d$two_theta, d$intensity)
  data.frame(sample = s, n_used = v$n_used, two_theta_min = v$range[1],
             two_theta_max = v$range[2], median_step = v$median_step,
             n_negative = v$n_negative, actions = paste(v$messages, collapse = " "))
}))
print(checks[, c("sample", "n_used", "two_theta_min", "two_theta_max", "median_step")],
      row.names = FALSE)
utils::write.csv(checks, file.path(outdir, "cotton_input_checks.csv"), row.names = FALSE)
status["input validation (7 patterns)"] <- if (nrow(checks) == 7L && all(checks$actions == "")) "PASS" else "FAIL"

# ---- 3-4. Segal CI and tidy table ---------------------------------------------------
say("== 3-4. Segal Crystallinity Index at I200 = %.1f, Iam = %.1f deg", I200, IAM)
fits <- lapply(samples, function(s) {
  d <- xrd[xrd$sample == s, ]
  segal_ci(d$two_theta, d$intensity, i200 = I200, iam = IAM, sample = s)
})
names(fits) <- samples
res <- do.call(rbind, lapply(fits, as.data.frame))
rownames(res) <- NULL
res$composition_label_pct <- as.numeric(sub("%.*$", "", res$sample))
res$ci_change_vs_first_pp <- res$ci - res$ci[1]
res <- res[, c("sample", "composition_label_pct", "ci", "ci_change_vs_first_pp",
               "i200_position", "i200_intensity", "iam_position", "iam_intensity",
               "interpolation", "n_used", "n_warnings")]
print(transform(res, ci = round(ci, 2), ci_change_vs_first_pp = round(ci_change_vs_first_pp, 2),
                i200_intensity = round(i200_intensity, 1), iam_intensity = round(iam_intensity, 1))
      [, c("sample", "ci", "ci_change_vs_first_pp", "i200_intensity", "iam_intensity")],
      row.names = FALSE)
utils::write.csv(res, file.path(outdir, "cotton_segal_results.csv"), row.names = FALSE)
status["no warnings on the seven patterns"] <- if (all(res$n_warnings == 0L)) "PASS" else "FAIL"

# ---- 5. composition vs CI plot ------------------------------------------------------
say("== 5. Composition-vs-CI plot")
p <- ggplot2::ggplot(res, ggplot2::aes(x = composition_label_pct, y = ci)) +
  ggplot2::geom_line(colour = "#A9B4D0", linewidth = 0.6) +
  ggplot2::geom_point(colour = "#243070", size = 3) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", ci)), vjust = -1, size = 3.2,
                     colour = "#1E2757") +
  ggplot2::scale_x_continuous(breaks = res$composition_label_pct) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.12))) +
  ggplot2::labs(
    x = "Composition label of the sample series (%)",
    y = "Segal Crystallinity Index (%)",
    title = "Segal Crystallinity Index of the seven cotton patterns",
    subtitle = sprintf("I200 = %.1f deg, Iam = %.1f deg (linear interpolation). Composition is a sample label, not a crystallinity scale.",
                       I200, IAM)) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 plot.subtitle = ggplot2::element_text(colour = "#5B6275", size = 9),
                 plot.title = ggplot2::element_text(face = "bold", colour = "#1E2757"),
                 plot.background = ggplot2::element_rect(fill = "white", colour = NA))
ggplot2::ggsave(file.path(outdir, "cotton_composition_vs_ci.png"), p, width = 7.5, height = 4.6,
                dpi = 150)

# ---- 6. sensitivity -----------------------------------------------------------------
say("== 6. Sensitivity analysis (I200 22.2-23.0, Iam 17.6-18.4)")
sens <- lapply(samples, function(s) {
  d <- xrd[xrd$sample == s, ]
  segal_sensitivity(d$two_theta, d$intensity, sample = s)
})
names(sens) <- samples
sens_summary <- do.call(rbind, lapply(sens, function(z) cbind(sample = z$sample, z$summary)))
rownames(sens_summary) <- NULL
sens_grid <- do.call(rbind, lapply(sens, as.data.frame))
rownames(sens_grid) <- NULL
print(data.frame(sample = sens_summary$sample, min = round(sens_summary$ci_min, 2),
                 max = round(sens_summary$ci_max, 2), range = round(sens_summary$ci_range, 2),
                 central = round(sens_summary$ci_central, 2)), row.names = FALSE)
utils::write.csv(sens_summary, file.path(outdir, "cotton_sensitivity_summary.csv"), row.names = FALSE)
utils::write.csv(sens_grid, file.path(outdir, "cotton_sensitivity_grid.csv"), row.names = FALSE)

# Does the ordering of the seven patterns survive the choice of positions?
wide <- reshape(sens_grid[, c("sample", "i200_position", "iam_position", "ci")],
                idvar = c("i200_position", "iam_position"), timevar = "sample",
                direction = "wide")
cols <- paste0("ci.", samples)
central_ci <- res$ci
rho <- apply(wide[, cols], 1, function(v) stats::cor(v, central_ci, method = "spearman"))
first_last <- wide[[cols[length(cols)]]] < wide[[cols[1]]]
say("   Rank agreement with the central (22.6/18.0) ordering across 25 cells: Spearman rho %.2f-%.2f",
    min(rho), max(rho))
say("   %s (%s) is lower than %s (%s) in %d of 25 grid cells", samples[length(samples)],
    "highest label", samples[1], "lowest label", sum(first_last))
rank_stability <- data.frame(i200_position = wide$i200_position, iam_position = wide$iam_position,
                             spearman_vs_central = rho, last_lower_than_first = first_last)
utils::write.csv(rank_stability, file.path(outdir, "cotton_sensitivity_rank_stability.csv"),
                 row.names = FALSE)

pr <- ggplot2::ggplot(sens_summary, ggplot2::aes(y = factor(sample, levels = rev(samples)))) +
  ggplot2::geom_segment(ggplot2::aes(x = ci_min, xend = ci_max, yend = factor(sample, levels = rev(samples))),
                        colour = "#A9B4D0", linewidth = 3, lineend = "round") +
  ggplot2::geom_point(ggplot2::aes(x = ci_central), colour = "#A5382C", size = 3) +
  ggplot2::labs(x = "Segal Crystallinity Index (%)", y = NULL,
                title = "Sensitivity to the I200 / Iam positions",
                subtitle = "Bar: minimum to maximum over the 5 x 5 grid; point: central estimate (22.6 / 18.0 deg)") +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 plot.subtitle = ggplot2::element_text(colour = "#5B6275", size = 9),
                 plot.title = ggplot2::element_text(face = "bold", colour = "#1E2757"),
                 plot.background = ggplot2::element_rect(fill = "white", colour = NA))
ggplot2::ggsave(file.path(outdir, "cotton_sensitivity_range.png"), pr, width = 7.5, height = 4.2,
                dpi = 150)

# Diagnostic only (not a selection rule): where is the (200) maximum in the 22-23 deg window
# that Segal's peak-height definition refers to, relative to the fixed 22.6 deg position?
apex <- do.call(rbind, lapply(samples, function(s) {
  d <- xrd[xrd$sample == s & xrd$two_theta >= 22 & xrd$two_theta <= 23, ]
  j <- which.max(d$intensity)
  data.frame(sample = s, apex_2theta = d$two_theta[j], apex_intensity = d$intensity[j])
}))
apex$i200_at_fixed_position <- res$i200_intensity
apex$fixed_over_apex <- apex$i200_at_fixed_position / apex$apex_intensity
utils::write.csv(apex, file.path(outdir, "cotton_apex_diagnostic.csv"), row.names = FALSE)
say("   Diagnostic: (200) maximum within 22-23 deg lies at %.2f-%.2f deg 2-theta (fixed I200 at %.1f deg)",
    min(apex$apex_2theta), max(apex$apex_2theta), I200)

# ---- 7. independent reference comparison -----------------------------------------------
say("== 7. Comparison with the independent spreadsheet calculation")
ref <- utils::read.csv(extdata("reference_cotton_segal.csv"))
cmp <- merge(data.frame(sample = res$sample, ci_R = res$ci), ref[, c("sample", "ci_reference")],
             by = "sample")
cmp <- cmp[match(samples, cmp$sample), ]
cmp$abs_difference <- abs(cmp$ci_R - cmp$ci_reference)
print(data.frame(sample = cmp$sample, CI_R = round(cmp$ci_R, 6),
                 CI_reference = round(cmp$ci_reference, 6),
                 abs_difference = signif(cmp$abs_difference, 3)), row.names = FALSE)
utils::write.csv(cmp, file.path(outdir, "cotton_reference_comparison.csv"), row.names = FALSE)
status["independent reference (7 CI values)"] <-
  if (nrow(cmp) == 7L && max(cmp$abs_difference) < TOL) "PASS" else "FAIL"

ref_s <- utils::read.csv(extdata("reference_cotton_sensitivity.csv"))
cmp_s <- merge(sens_grid, ref_s, by = c("sample", "i200_position", "iam_position"))
cmp_s$abs_difference <- abs(cmp_s$ci - cmp_s$ci_reference)
say("   Sensitivity grids (3 samples x 25 cells): max abs difference %.2e percentage points",
    max(cmp_s$abs_difference))
utils::write.csv(cmp_s, file.path(outdir, "cotton_sensitivity_reference_comparison.csv"),
                 row.names = FALSE)
status["independent reference (75 sensitivity cells)"] <-
  if (nrow(cmp_s) == 75L && max(cmp_s$abs_difference) < TOL) "PASS" else "FAIL"

# Agreement with the values quoted in the brief (2 d.p.); informational cross-check.
brief <- c("0%/1" = 71.19, "20%/1" = 67.65, "40%/1" = 69.36, "50%/1" = 69.24,
           "60%/1" = 65.65, "80%/1" = 62.92, "100%/1" = 56.96)
status["approx. values from the brief (2 d.p.)"] <-
  if (all(round(setNames(res$ci, res$sample)[names(brief)], 2) == brief)) "PASS" else "FAIL"

# ---- summary ---------------------------------------------------------------------------
say("== Summary")
for (k in names(status)) say("   [%s] %s", status[[k]], k)
say("Outputs written to %s", normalizePath(outdir))
if (any(status != "PASS")) quit(status = 1L)
