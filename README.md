# textileCrystR

A small, reproducible R package for X-ray diffraction (XRD) crystallinity
determination of textile fibres and cellulose materials. It implements
**distinct, independent published methods** rather than claiming a single
"true" crystallinity, and never averages them into one number.

## Supported methods

* **Segal peak height** -- `segal_ci()` and friends. The classic, simplest
  empirical index:

  ```
  CI = (I200 - Iam) / I200 * 100
  ```

* **Peak area / deconvolution** -- `fit_crystalline_peaks()` and friends.
  Jointly fits several sharp crystalline peaks and one broad amorphous hump
  (Gaussian, Lorentzian, pseudo-Voigt or Voigt profiles), and reports
  crystallinity as the crystalline peak area over the total area.

The package's scope is intentionally kept to these two well-validated
methods (see `?textileCrystR`).

## Install

```r
# from the package source tarball
install.packages("textileCrystR_0.1.0.tar.gz", repos = NULL, type = "source")

# optional, only needed for the interactive dashboard:
install.packages(c("shiny", "rmarkdown"))
```

## Use it in one line (Segal)

```r
library(textileCrystR)

path <- system.file("extdata", "cotton.csv", package = "textileCrystR")
res <- analyze_textile_file(path, sample = "0%/1")
res
```

`analyze_textile_file()` reads the file (auto-detecting delimiter, header
and columns via `import_xrd_file()`, forgivingly enough for real messy
exports), validates it, calculates the Segal index, shows the plot, and can
write a plot image and/or a self-contained HTML summary report in the same
call:

```r
analyze_textile_file(path, sample = "0%/1",
                     plot_file = "segal_plot.png", report_file = "segal_report.html")
```

Leave `sample` unset on a file with several patterns (like the shipped
`cotton.csv`) to analyze all of them at once and get a `textile_cryst_batch`
with its own `print()`/`plot()`/`as.data.frame()` methods.

## Peak-area / deconvolution

```r
xrd <- read_xrd_pairs(path)
d <- xrd[xrd$sample == "0%/1", ]

fit <- fit_crystalline_peaks(d$two_theta, d$intensity, sample = "0%/1")
fit
plot_peak_fit(fit)          # raw pattern + peaks + amorphous hump + residuals
peak_fit_summary(fit)       # tidy one-row-per-peak table
```

Fits four crystalline reflections (the cellulose I (1-10), (110), (200) and
(004) peaks, by default) plus the amorphous hump **simultaneously** -- the
scientifically necessary approach, since they overlap. Includes optional
linear background handling for real (non-baseline-corrected) patterns.

## Or go interactive

```r
run_textile_app()
```

Launches a local Shiny dashboard covering **both methods**: upload one or
more files at once (every sample from every file is pooled into a single
selector), pick Segal or peak area / deconvolution from a method selector,
and either fine-tune one sample live (sliders and click-to-set for Segal; a
"Run fit" button for peak area) or switch to batch mode to analyze any
number of the pooled samples together, producing a results table, a
comparison chart, and a downloadable CSV. Partial batch failures (e.g. a
sample whose range doesn't cover the chosen positions) are reported
individually without losing the ones that succeeded.

## The building blocks

For more control than the one-shot wrappers give you:

* `segal_ci()` / `segal_points()` / `segal_presets()` / `segal_sensitivity()`
  -- the Segal method and its robustness analysis.
* `fit_crystalline_peaks()` / `fit_amorphous_hump()` /
  `peak_area_crystallinity()` / `peak_fit_summary()` / `plot_peak_fit()` --
  the peak-area / deconvolution method.
* `validate_crystallinity_input()` / `import_xrd_file()` / `read_xrd_pairs()`
  -- input checking and file reading, usable on their own.

## Learn more

* `vignette("getting-started")` -- a practical, start-to-finish tour of the
  package as a piece of software, covering both methods and the dashboard.
* `vignette("segal-cotton-validation")` -- the Segal method itself, validated
  end to end against the shipped experimental cotton data set, including an
  independent (spreadsheet) reference calculation.

## Scientific scope

Every method here is **empirical and relative**, not an absolute
crystallinity — see the *Limitations* section of `?textileCrystR` for the
full list of caveats (peak overlap, crystallite size, cellulose polymorph,
preferred orientation, background handling, position/model choice, sample
preparation, and partially mercerized cotton) and their literature sources
(Segal et al., 1959; Nam et al., 2016; Salem et al., 2023).
