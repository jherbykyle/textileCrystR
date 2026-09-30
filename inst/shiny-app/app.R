# textileCrystR interactive dashboard.
#
# Launched via textileCrystR::run_textile_app(); not meant to be sourced
# directly by users. Internal package helpers (.render_report,
# .clamp_preset_positions, .click_slider_id) are reached with ::: since this
# script runs outside the package namespace.

library(shiny)
library(textileCrystR)

# need()'s second argument is evaluated eagerly (it is an ordinary function
# argument, not lazy), so `need(!inherits(x, "error"), conditionMessage(x))`
# calls conditionMessage() on x even when x is NOT an error -- which errors
# out immediately, since conditionMessage() only accepts condition objects.
# This helper only touches conditionMessage() on the error path.
.stop_if_error <- function(x) {
  msg <- if (inherits(x, "error")) conditionMessage(x) else NULL
  validate(need(!inherits(x, "error"), msg))
}

# Parse a comma/whitespace-separated list of numbers (crystalline peak
# centers). Returns NULL (meaning "use the package default") on blank input,
# or a plain-English error via validate() on unparseable input.
.parse_centers <- function(text) {
  text <- trimws(text)
  if (!nzchar(text)) return(NULL)
  parts <- strsplit(text, "[,\\s]+", perl = TRUE)[[1]]
  parts <- parts[nzchar(parts)]
  nums <- suppressWarnings(as.numeric(parts))
  validate(need(!anyNA(nums) && length(nums) >= 1L,
               paste0("Peak centers should be a comma-separated list of numbers, e.g. ",
                      "\"14.8, 16.5, 22.6, 34.5\". Got: \"", text, "\".")))
  nums
}

# A compact, method-agnostic key identifying one sample within one uploaded
# file, so the same sample name in two different files never collides.
.sample_key <- function(file, sample) paste(file, sample, sep = " :: ")

.safe_filename <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

# ---- UI ---------------------------------------------------------------------------

ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { background-color: #FAFAFB; }
    .well { background-color: #FFFFFF; border-color: #E3E6EF; }
    h2.app-title { color: #1E2757; font-weight: 700; margin-bottom: 0; }
    .app-subtitle { color: #5B6275; margin-top: 2px; margin-bottom: 18px; }
    .btn-download { width: 100%; margin-bottom: 8px; }
    .method-note { color: #5B6275; font-size: 90%; }
  "))),
  titlePanel(title = tagList(
    h2("textileCrystR", class = "app-title"),
    div("XRD crystallinity -- interactive dashboard (Segal and peak-area methods)",
        class = "app-subtitle")
  ), windowTitle = "textileCrystR"),

  sidebarLayout(
    sidebarPanel(
      width = 4,
      fileInput("files", "XRD file(s) (CSV / TXT / TSV)", multiple = TRUE,
                accept = c(".csv", ".txt", ".tsv", ".dat", ".xy")),
      uiOutput("sample_ui"),

      tags$hr(),
      checkboxInput("subtract_background", "Subtract estimated background before analysis",
                   value = FALSE),
      conditionalPanel(
        "input.subtract_background",
        numericInput("bg_window", "Background window (\u00b0 2-theta)", value = 8, min = 1, step = 1),
        helpText("Estimates a smooth, non-parametric background (SNIP algorithm) and subtracts it ",
                 "before running the selected method. Too narrow a window follows real peaks too ",
                 "closely (subtracting signal); too wide risks clipping into a genuine amorphous ",
                 "halo. Preview below in single-sample mode before trusting it.")
      ),

      tags$hr(),
      radioButtons("method", "Method",
                   choices = c("Segal peak height" = "segal",
                               "Peak area / deconvolution" = "peakarea"),
                   selected = "segal"),

      conditionalPanel(
        "input.method == 'segal'",
        selectInput("preset", "Position preset",
                    choices = stats::setNames(segal_presets()$preset, segal_presets()$material)),
        helpText("Choosing a preset sets the sliders below to its published positions, ",
                 "clamped to what the sample actually covers. Drag the sliders, or click ",
                 "directly on the peak/trough in the plot (single-sample mode), to fine-tune."),
        sliderInput("i200", "I200 position (crystalline peak, \u00b0 2-theta)",
                    min = 0, max = 1, value = 0.5, step = 0.05),
        sliderInput("iam", "Iam position (amorphous minimum, \u00b0 2-theta)",
                    min = 0, max = 1, value = 0.5, step = 0.05),
        radioButtons("click_target", "Clicking the plot sets:",
                     choices = c("I200 (crystalline peak)" = "i200",
                                 "Iam (amorphous minimum)" = "iam"),
                     selected = "i200", inline = TRUE),
        actionButton("reset_preset", "Reset sliders to preset", icon = icon("rotate-left"))
      ),

      conditionalPanel(
        "input.method == 'peakarea'",
        selectInput("model", "Peak-shape model",
                    choices = c("Pseudo-Voigt (default)" = "pseudo_voigt", "Gaussian" = "gaussian",
                                "Lorentzian" = "lorentzian", "Voigt" = "voigt")),
        textInput("centers", "Crystalline peak centers (\u00b0 2-theta, comma-separated)",
                  value = "", placeholder = "default: 14.8, 16.5, 22.6, 34.5"),
        fluidRow(
          column(6, numericInput("am_win_lo", "Amorphous window: from", value = 18, step = 0.5)),
          column(6, numericInput("am_win_hi", "to", value = 25, step = 0.5))
        ),
        selectInput("background", "Background",
                    choices = c("Linear (default)" = "linear", "Constant" = "constant", "None" = "none")),
        helpText("Peak fitting is a nonlinear optimization and takes a moment -- press \"Run fit\" ",
                 "after changing settings (single-sample mode only; batch mode always uses the ",
                 "current settings when you run it)."),
        actionButton("run_fit", "Run fit", icon = icon("play"), class = "btn-primary")
      ),

      tags$hr(),
      conditionalPanel(
        "input.view_mode == 'single' || typeof input.view_mode === 'undefined'",
        downloadButton("dl_plot", "Download plot (PNG)", class = "btn-download"),
        downloadButton("dl_report", "Download summary report (HTML)", class = "btn-download")
      )
    ),

    mainPanel(
      width = 8,
      uiOutput("main_ui")
    )
  )
)

# ---- Server -----------------------------------------------------------------------

server <- function(input, output, session) {

  # ---- data: read every uploaded file, pool all samples together ----------------

  all_data <- reactive({
    req(input$files)
    parts <- lapply(seq_len(nrow(input$files)), function(i) {
      out <- tryCatch(import_xrd_file(input$files$datapath[i]), error = function(e) e)
      list(name = input$files$name[i], data = out)
    })
    errs <- Filter(function(p) inherits(p$data, "error"), parts)
    if (length(errs) > 0) {
      msgs <- vapply(errs, function(p) sprintf("%s: %s", p$name, conditionMessage(p$data)), character(1))
      validate(need(FALSE, paste(msgs, collapse = "\n\n")))
    }
    do.call(rbind, lapply(parts, function(p) {
      d <- p$data
      d$file <- p$name
      d$key <- .sample_key(p$name, d$sample)
      d
    }))
  })

  multi_file <- reactive(length(unique(all_data()$file)) > 1L)

  sample_choices <- reactive({
    d <- all_data()
    keys <- unique(d$key)
    labels <- if (isTRUE(multi_file())) keys else vapply(keys, function(k) d$sample[d$key == k][1], character(1))
    stats::setNames(keys, labels)
  })

  output$sample_ui <- renderUI({
    req(sample_choices())
    tagList(
      radioButtons("view_mode", "View",
                   choices = c("Single sample (interactive)" = "single",
                               "Batch (compare samples)" = "batch"),
                   selected = "single", inline = TRUE),
      conditionalPanel(
        "input.view_mode == 'single'",
        selectInput("sample_single", "Sample", choices = sample_choices())
      ),
      conditionalPanel(
        "input.view_mode == 'batch'",
        selectizeInput("samples_batch", "Samples", choices = sample_choices(), multiple = TRUE,
                       options = list(placeholder = "Select samples, or leave empty for all")),
        helpText("Batch mode uses the current preset (Segal) or model/window/background settings ",
                 "(peak area) applied fresh to every selected sample -- not the fine-tuned sliders ",
                 "from single-sample mode."),
        actionButton("run_batch", "Analyze selected samples", icon = icon("play"), class = "btn-primary")
      )
    )
  })

  view_mode <- reactive(if (is.null(input$view_mode)) "single" else input$view_mode)

  current_data_raw <- reactive({
    req(input$sample_single)
    d <- all_data()
    d[d$key == input$sample_single, ]
  })
  current_label <- reactive({
    req(input$sample_single)
    names(sample_choices())[sample_choices() == input$sample_single][1]
  })

  # Background estimate for the currently selected single sample. NULL when
  # the toggle is off. A warning (e.g. the requested window being capped for
  # a short pattern) is captured for display without interrupting the
  # estimate itself -- unlike a plain tryCatch(warning = ...), which would
  # abort the computation and return only the warning condition.
  current_background <- reactive({
    if (!isTRUE(input$subtract_background)) return(NULL)
    d <- current_data_raw()
    req(nrow(d) > 0)
    warn_msg <- NULL
    out <- tryCatch(
      withCallingHandlers(
        estimate_xrd_background(d$two_theta, d$intensity, window = input$bg_window),
        warning = function(w) { warn_msg <<- conditionMessage(w); invokeRestart("muffleWarning") }
      ),
      error = function(e) e
    )
    if (!inherits(out, "error")) attr(out, "warning") <- warn_msg
    out
  })

  current_data <- reactive({
    raw <- current_data_raw()
    bg <- current_background()
    if (is.null(bg)) return(raw)
    .stop_if_error(bg)
    corrected <- subtract_xrd_background(raw$two_theta, raw$intensity, bg, clip_negative = TRUE)
    data.frame(sample = raw$sample[1], file = raw$file[1], key = raw$key[1],
              two_theta = corrected$two_theta, intensity = corrected$corrected,
              stringsAsFactors = FALSE)
  })

  batch_keys <- reactive({
    chosen <- input$samples_batch
    if (is.null(chosen) || length(chosen) == 0) unname(sample_choices()) else chosen
  })

  # ---- Segal: live sliders, reset-to-preset -------------------------------------

  observeEvent(list(input$preset, current_data_raw(), input$reset_preset), {
    d <- current_data_raw()
    req(nrow(d) > 0)
    rng <- range(d$two_theta)
    req(diff(rng) > 0)  # a degenerate (zero-width) range can't be clamped into; leave sliders as-is
    vals <- textileCrystR:::.clamp_preset_positions(rng, input$preset)
    updateSliderInput(session, "i200", min = round(rng[1], 2), max = round(rng[2], 2),
                      value = round(vals$i200, 2), step = 0.05)
    updateSliderInput(session, "iam", min = round(rng[1], 2), max = round(rng[2], 2),
                      value = round(vals$iam, 2), step = 0.05)
  }, ignoreInit = FALSE)

  observeEvent(input$plot_click, {
    req(input$plot_click$x, input$method == "segal")
    target <- textileCrystR:::.click_slider_id(input$click_target)
    updateSliderInput(session, target, value = round(input$plot_click$x, 2))
  })

  # ---- single-sample result: dispatches on method; plot()/summary() are S3-generic ----

  segal_result <- reactive({
    d <- current_data()
    req(nrow(d) > 0, input$i200, input$iam)
    tryCatch(segal_ci(d$two_theta, d$intensity, i200 = input$i200, iam = input$iam,
                      sample = current_label()),
            error = function(e) e)
  })

  peakarea_result <- eventReactive(input$run_fit, {
    d <- current_data()
    req(nrow(d) > 0)
    centers <- .parse_centers(input$centers)
    tryCatch(
      fit_crystalline_peaks(d$two_theta, d$intensity, centers = centers, model = input$model,
                            amorphous = c(input$am_win_lo, input$am_win_hi),
                            background = input$background, sample = current_label()),
      error = function(e) e
    )
  })

  result <- reactive({
    if (identical(input$method, "peakarea")) {
      req(input$run_fit > 0)
      peakarea_result()
    } else {
      segal_result()
    }
  })

  # ---- batch result: loops the same per-sample function used above -------------

  batch_result <- eventReactive(input$run_batch, {
    d <- all_data()
    keys <- batch_keys()
    validate(need(length(keys) >= 1L, "Select at least one sample."))
    method <- isolate(input$method)
    centers <- if (method == "peakarea") isolate(.parse_centers(input$centers)) else NULL
    preset <- isolate(input$preset)
    subtract_bg <- isolate(isTRUE(input$subtract_background))
    bg_window <- isolate(input$bg_window)
    rows <- list(); errs <- character()
    for (k in keys) {
      dk <- d[d$key == k, ]
      label <- names(sample_choices())[sample_choices() == k][1]
      res <- tryCatch({
        if (subtract_bg) {
          bg_k <- estimate_xrd_background(dk$two_theta, dk$intensity, window = bg_window)
          corrected_k <- subtract_xrd_background(dk$two_theta, dk$intensity, bg_k, clip_negative = TRUE)
          dk <- data.frame(two_theta = corrected_k$two_theta, intensity = corrected_k$corrected)
        }
        if (method == "peakarea") {
          fit_crystalline_peaks(dk$two_theta, dk$intensity, centers = centers,
                                model = isolate(input$model),
                                amorphous = c(isolate(input$am_win_lo), isolate(input$am_win_hi)),
                                background = isolate(input$background), sample = label)
        } else {
          # Uses the preset directly (recomputed fresh for each sample's own
          # measured range), not the single-sample sliders: those are tied
          # to whichever one sample was last viewed there and may not even
          # be initialized if the user goes straight to batch mode.
          segal_ci(dk$two_theta, dk$intensity, preset = preset, sample = label)
        }
      }, error = function(e) e)
      if (inherits(res, "error")) {
        errs <- c(errs, sprintf("%s: %s", label, conditionMessage(res)))
      } else {
        rows[[label]] <- res
      }
    }
    validate(need(length(rows) >= 1L, paste0("Every selected sample failed:\n\n",
                                             paste(errs, collapse = "\n\n"))))
    list(method = method, results = rows, errors = errs,
        summary = do.call(rbind, lapply(rows, as.data.frame)))
  })

  # ---- outputs: single-sample panel ---------------------------------------------

  output$main_ui <- renderUI({
    if (identical(view_mode(), "batch")) {
      tagList(
        h4("Batch results"),
        tableOutput("batch_table"),
        uiOutput("batch_errors_ui"),
        tags$hr(),
        plotOutput("batch_plot", height = "420px"),
        tags$hr(),
        downloadButton("dl_batch_csv", "Download batch results (CSV)", class = "btn-download"),
        tags$p(tags$em(paste0(
          "Every XRD crystallinity method here is empirical and relative -- see ?textileCrystR ",
          "for the full list of caveats.")), class = "method-note")
      )
    } else {
      tagList(
        uiOutput("background_preview_ui"),
        plotOutput("xrd_plot", click = "plot_click", height = "440px"),
        tags$hr(),
        h4("Summary"),
        verbatimTextOutput("summary_text"),
        tags$p(tags$em(paste0(
          "Every XRD crystallinity method here is empirical and relative, not an absolute ",
          "crystallinity -- see ?textileCrystR for the full list of caveats.")), class = "method-note")
      )
    }
  })

  output$xrd_plot <- renderPlot({
    res <- result()
    .stop_if_error(res)
    plot(res)
  })

  output$background_preview_ui <- renderUI({
    req(isTRUE(input$subtract_background))
    tagList(
      h4("Background preview"),
      plotOutput("background_preview_plot", height = "300px"),
      uiOutput("background_warning_ui"),
      tags$hr()
    )
  })

  output$background_preview_plot <- renderPlot({
    req(isTRUE(input$subtract_background))
    raw <- current_data_raw()
    bg <- current_background()
    .stop_if_error(bg)
    d <- subtract_xrd_background(raw$two_theta, raw$intensity, bg, clip_negative = TRUE)
    ggplot2::ggplot(d, ggplot2::aes(x = .data$two_theta)) +
      ggplot2::geom_line(ggplot2::aes(y = .data$intensity, colour = "Raw"), linewidth = 0.5) +
      ggplot2::geom_line(ggplot2::aes(y = .data$background, colour = "Background"),
                         linetype = "dashed", linewidth = 0.8) +
      ggplot2::geom_line(ggplot2::aes(y = .data$corrected, colour = "Corrected"), linewidth = 0.6) +
      ggplot2::scale_colour_manual(name = NULL,
                                  values = c(Raw = "#5B6275", Background = "#D89B16",
                                            Corrected = "#A5382C"),
                                  breaks = c("Raw", "Background", "Corrected")) +
      ggplot2::labs(title = "Estimated background (SNIP)", x = "2-theta (\u00b0)", y = "Intensity") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(),
                     plot.background = ggplot2::element_rect(fill = "white", colour = NA))
  })

  output$background_warning_ui <- renderUI({
    bg <- current_background()
    if (is.null(bg) || inherits(bg, "error")) return(NULL)
    w <- attr(bg, "warning")
    if (is.null(w)) return(NULL)
    tags$p(tags$strong("Note: "), w, class = "method-note")
  })

  output$summary_text <- renderPrint({
    res <- result()
    .stop_if_error(res)
    summary(res)
  })

  # ---- outputs: batch panel ------------------------------------------------------

  output$batch_table <- renderTable({
    b <- batch_result()
    tab <- b$summary
    if (b$method == "peakarea") {
      data.frame(Sample = tab$sample, `CI (%)` = round(tab$ci, 2), Model = tab$model,
                `R-squared` = round(tab$r_squared, 4), `N peaks` = tab$n_peaks,
                Warnings = tab$n_warnings, check.names = FALSE)
    } else {
      data.frame(Sample = tab$sample, `CI (%)` = round(tab$ci, 2),
                `I200 (deg)` = round(tab$i200_position, 2), `Iam (deg)` = round(tab$iam_position, 2),
                Warnings = tab$n_warnings, check.names = FALSE)
    }
  })

  output$batch_errors_ui <- renderUI({
    b <- batch_result()
    if (length(b$errors) == 0) return(NULL)
    tagList(tags$p(tags$strong(sprintf("%d sample(s) could not be analyzed:", length(b$errors)))),
           tags$pre(paste(b$errors, collapse = "\n\n")))
  })

  output$batch_plot <- renderPlot({
    b <- batch_result()
    tab <- b$summary
    tab$sample <- factor(tab$sample, levels = tab$sample)
    ggplot2::ggplot(tab, ggplot2::aes(x = .data$sample, y = .data$ci)) +
      ggplot2::geom_col(fill = "#243070") +
      ggplot2::labs(x = NULL, y = "Crystallinity Index (%)",
                    title = sprintf("%s -- %d sample(s)",
                                    if (b$method == "peakarea") "Peak area / deconvolution" else "Segal",
                                    nrow(tab))) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
                     panel.grid.minor = ggplot2::element_blank(),
                     plot.background = ggplot2::element_rect(fill = "white", colour = NA))
  })

  output$dl_batch_csv <- downloadHandler(
    filename = function() sprintf("textileCrystR_batch_%s.csv", batch_result()$method),
    content = function(file) utils::write.csv(batch_result()$summary, file, row.names = FALSE)
  )

  # ---- downloads: single-sample plot/report --------------------------------------

  output$dl_plot <- downloadHandler(
    filename = function() sprintf("%s_%s_plot.png", .safe_filename(current_label()), input$method),
    content = function(file) {
      res <- result()
      .stop_if_error(res)
      ggplot2::ggsave(file, plot(res), width = 8, height = if (input$method == "peakarea") 6 else 4.2,
                      dpi = 150)
    }
  )

  output$dl_report <- downloadHandler(
    filename = function() sprintf("%s_%s_report.html", .safe_filename(current_label()), input$method),
    content = function(file) {
      res <- result()
      .stop_if_error(res)
      textileCrystR:::.render_report(list(result = res), file)
    }
  )
}

shinyApp(ui, server)
