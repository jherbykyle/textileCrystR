# textileCrystR 0.1.0

Initial release.

* **Segal peak-height index**: `segal_ci()`, `segal_points()`,
  `segal_presets()` and `segal_sensitivity()`, with explicit interpolation,
  literature-supported position presets, and a robustness analysis over a
  grid of I200/Iam positions. Validated against an independent spreadsheet
  calculation and a shipped experimental cotton data set.
* **Peak-area / deconvolution method**: `fit_crystalline_peaks()`,
  `fit_amorphous_hump()`, `peak_area_crystallinity()`, `peak_fit_summary()`
  and `plot_peak_fit()`. Crystalline peaks and the amorphous hump are fit
  jointly (Gaussian, Lorentzian, pseudo-Voigt or Voigt profiles), with an
  optional linear or constant background term that is excluded from the
  crystallinity calculation. Checked against synthetic patterns with known
  parameters for every peak model.
* **General XRD preprocessing**: `validate_xrd_pattern()`,
  `sort_xrd_pattern()`, `remove_xrd_duplicates()`,
  `estimate_xrd_background()` (SNIP algorithm), `subtract_xrd_background()`
  and `normalize_xrd_pattern()`. The accuracy of the background estimate on
  sloped backgrounds is quantified in `?estimate_xrd_background`.
* **File input**: `import_xrd_file()` reads CSV/TSV/whitespace-delimited
  exports with automatic delimiter, header and column detection, and
  `read_xrd_pairs()` reads the paired multi-sample layout.
* **One-shot analysis**: `analyze_textile_file()` reads a file, validates it,
  computes the Segal index, plots it, and can write a plot image and a
  self-contained HTML report; multiple samples give a `textile_cryst_batch`.
* **Interactive dashboard**: `run_textile_app()` (needs the `shiny` and
  `rmarkdown` packages) supports both methods, multiple uploaded files and
  samples, single-sample tuning or batch analysis, optional background
  subtraction with a preview, and plot, report and CSV downloads.

Every method is treated as an empirical, relative XRD index and results
from different methods are never averaged together; see `?textileCrystR`
for the scientific limitations and references.
