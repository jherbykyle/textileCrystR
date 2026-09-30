# Extract the independent (spreadsheet) Segal calculation from the validation
# workbook and store it as plain CSV in inst/extdata/.
#
# The workbook computes Segal CI with Excel formulas of the form
#   =(INDEX(Raw_XRD!col, row_I200) - INDEX(Raw_XRD!col, row_Iam)) / INDEX(...) * 100
# i.e. by reading measured rows directly, with no R code and no interpolation.
# Only the values Excel calculated (cached results) are copied; nothing here
# recomputes the Segal index. Run from the package root:
#   Rscript data-raw/make_reference_from_workbook.R

wb <- "data-raw/textileCrystR_cotton_validation.xlsx"
stopifnot(file.exists(wb))

# --- main validation table (Validation!A7:F14) ---------------------------------
v <- readxl::read_excel(wb, sheet = "Validation", range = "A7:F14")
names(v) <- c("sample", "i200_position", "i200_intensity", "iam_position",
              "iam_intensity", "ci_reference")
v$composition_label_pct <- as.numeric(sub("%.*$", "", v$sample))
v <- v[, c("sample", "composition_label_pct", "i200_position", "i200_intensity",
           "iam_position", "iam_intensity", "ci_reference")]
utils::write.csv(v, "inst/extdata/reference_cotton_segal.csv", row.names = FALSE)

# --- sensitivity grids (Sensitivity!A7:R12), three samples ----------------------
s <- readxl::read_excel(wb, sheet = "Sensitivity", range = "A7:R12", col_names = FALSE)
blocks <- list("0%/1" = 1:6, "50%/1" = 7:12, "100%/1" = 13:18)
out <- lapply(names(blocks), function(nm) {
  b <- as.data.frame(s[, blocks[[nm]]])
  i200 <- as.numeric(b[1, -1])
  iam <- as.numeric(b[-1, 1])
  vals <- as.matrix(b[-1, -1])
  data.frame(sample = nm,
             i200_position = rep(i200, each = length(iam)),
             iam_position = rep(iam, times = length(i200)),
             ci_reference = as.numeric(vals))
})
out <- do.call(rbind, out)
utils::write.csv(out, "inst/extdata/reference_cotton_sensitivity.csv", row.names = FALSE)
message("Wrote ", nrow(v), " main rows and ", nrow(out), " sensitivity rows.")
