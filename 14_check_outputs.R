# ==============================================================================
# 14_check_outputs.R
# ------------------------------------------------------------------------------
# Verifies that a run actually produced everything it should, and says precisely
# what is missing if not.
#
# This exists because a long run can fail part-way - a network write error, a
# reboot - and leave an output folder that looks populated but is not complete.
# Comparing what is on disk against the expected manifest is the only reliable
# way to know whether results can be used.
#
# Run on its own at any time:  source("14_check_outputs.R")
# run_all.R calls it automatically at the end.
#
# Outputs: logs/output_check.csv   (one row per expected file, with status)
# ==============================================================================

expected_outputs <- function() {

  e <- list()

  add <- function(step, subdir, files) {
    if (!length(files)) return(invisible(NULL))
    e[[length(e) + 1]] <<- data.frame(step = step, subdir = subdir,
                                      file = files, stringsAsFactors = FALSE)
  }

  add("02_descriptives", "csv",
      c("descriptives_by_region.csv", "descriptives_national.csv",
        "sample_size.csv", "exposure_days_per_year_by_region.csv",
        "fhs_coverage_by_region_year.csv",
        "tableA1_summary_stats_long.csv"))
  add("02_descriptives", "tables", "tableA1_summary_stats.csv")
  add("02_descriptives", "figures", "fig1_fhs_coverage_trend.png")

  # Figure 1 itself.  The maps need geobr and sf, which 02_descriptives.R
  # installs on demand; if that fails the script still writes the coverage data
  # and skips the drawing, so only the csv outputs are treated as required and a
  # run is not reported incomplete for want of a mapping package.
  add("02_descriptives", "csv",
      c("fig1_coverage_by_municipality.csv",
        "fig1_coverage_weighted_by_year.csv"))

  # 03 writes three files per region x outcome, plus the combined files
  grid <- expand.grid(reg = REGIONS, out = OUTCOMES, stringsAsFactors = FALSE)
  add("03_main_dlm", "csv",
      c(sprintf("cumulative_nointeraction_%s_%s.csv", grid$reg, grid$out),
        sprintf("predicted_by_coverage_%s_%s.csv",    grid$reg, grid$out),
        "cumulative_nointeraction_ALL.csv", "predicted_by_coverage_ALL.csv"))
  add("03_main_dlm", "models", sprintf("fit_%s_%s.rds", grid$reg, grid$out))
  add("03_main_dlm", "figures", "fig2_temperature_by_coverage_rate_total.png")

  add("04_non_interacted", "figures",
      "appendixB_temperature_mortality_no_interaction.png")

  add("05_table1", "tables",
      c("table1_fhs90_vs_fhs0_very_hot_days.csv", "table1_formatted.csv"))

  add("06_gdp_split", "csv", "appendixD_gdp_split_ALL.csv")
  add("06_gdp_split", "figures",
      sprintf("appendixD_gdp_split_%s.png",
              gsub("[^A-Za-z]", "", unname(REGION_LABELS[SENS_REGIONS]))))

  add("07_deaths_averted", "csv",
      c("averted_deaths_grand_total.csv", "averted_deaths_by_year.csv",
        "averted_deaths_by_bin.csv", "averted_deaths_by_age.csv",
        "averted_deaths_decomposition_check.csv"))

  # Appendix F: Southeast and South only, one figure per region, panels by age
  # group - the structure of the submitted appendix.  SENS_REGIONS and
  # SENS_OUTCOMES in 00_config.R set the scope.
  sens_reg <- gsub("[^A-Za-z]", "", unname(REGION_LABELS[SENS_REGIONS]))

  add("08_placebo_future_change", "csv",
      c("appendixF_placebo_future_change.csv",
        "appendixF_placebo_future_change_diagnostics.csv"))
  add("08_placebo_future_change", "figures",
      sprintf("appendixF_placebo_future_change_%s.png", sens_reg))

  add("09_lag_windows",     "csv", "appendixF_lag_windows.csv")
  add("09_lag_windows",     "figures",
      sprintf("appendixF_lag_windows_%s.png", sens_reg))

  add("10_coverage_bands",  "csv",
      c("appendixF_coverage_bands.csv",
        "appendixF_coverage_band_composition.csv"))
  add("10_coverage_bands",  "figures",
      sprintf("appendixF_coverage_bands_%s.png", sens_reg))

  add("11_exclude_large",   "csv", "appendixF_exclude_large_muni.csv")
  add("11_exclude_large",   "figures",
      sprintf("appendixF_exclude_large_muni_%s.png", sens_reg))

  add("12_percentile_bins", "csv",
      c("appendixF_percentile_bins.csv", "appendixF_percentile_cutpoints.csv"))
  add("12_percentile_bins", "figures",
      sprintf("appendixF_percentile_bins_%s.png", sens_reg))
  add("13_session_info",    "logs", "session_info.txt")

  do.call(rbind, e)
}

check_outputs <- function(verbose = TRUE) {

  man <- expected_outputs()
  man$path   <- file.path(OUT_DIR, man$subdir, man$file)
  man$exists <- file.exists(man$path)
  man$bytes  <- ifelse(man$exists, file.size(man$path), NA_real_)
  man$ok     <- man$exists & !is.na(man$bytes) & man$bytes > 0

  by_step <- do.call(rbind, lapply(split(man, man$step), function(d)
    data.frame(step = d$step[1], expected = nrow(d), present = sum(d$ok),
               missing = sum(!d$ok), stringsAsFactors = FALSE)))
  by_step <- by_step[order(by_step$step), ]

  if (verbose) {
    cat("\n================ OUTPUT CHECK ================\n")
    print(by_step, row.names = FALSE)
    if (any(!man$ok)) {
      cat("\nMissing or empty (", sum(!man$ok), " files):\n", sep = "")
      miss <- man[!man$ok, c("step", "subdir", "file")]
      print(utils::head(miss, 40), row.names = FALSE)
      if (nrow(miss) > 40) cat("... and", nrow(miss) - 40, "more\n")
      cat("\nRe-run run_all.R: with SKIP_EXISTING = TRUE, completed work is\n",
          "reused and only the missing pieces are computed.\n", sep = "")
    } else {
      cat("\nComplete: every expected output is present and non-empty.\n")
    }
  }

  write.csv(man[, c("step", "subdir", "file", "exists", "bytes", "ok")],
            file.path(dir_out("logs"), "output_check.csv"), row.names = FALSE)

  invisible(list(manifest = man, by_step = by_step, complete = all(man$ok)))
}

if (!exists(".sourced_by_run_all")) invisible(check_outputs())
