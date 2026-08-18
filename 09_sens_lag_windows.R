# ==============================================================================
# 09_sens_lag_windows.R
# ------------------------------------------------------------------------------
# Appendix F(ii), Figures F3 (Southeast) and F4 (South): alternative
# distributed-lag windows.  The main analysis cumulates temperature effects over
# lags 0-30; here the same model is estimated cumulating over 0-7 and 0-14 days.
#
# A shorter window captures less mortality displacement (harvesting): if a hot
# day mainly brings forward deaths that would have occurred within a fortnight,
# the 7-day cumulative effect is larger than the 30-day one.  Comparing the
# three windows is therefore informative about displacement, not only about
# robustness.
#
# STRUCTURE: one figure per region, columns = lag window, rows = age group, to
# match the submitted appendix.  The 30-day column is NOT refitted here - it is
# the main model, so it is read back from the CSVs 03_main_dlm.R wrote.
#
# Scope comes from SENS_REGIONS and SENS_OUTCOMES in 00_config.R.
#
# Inputs : panel (via get_panel()); csv/predicted_by_coverage_<reg>_<out>.csv
# Outputs: csv/appendixF_lag_windows.csv
#          figures/appendixF_lag_windows_<region>.png
# ==============================================================================

panel <- get_panel()

results <- list()

for (k in LAG_MAX_SENS) {
  for (reg in SENS_REGIONS) {

    ref_bin  <- REF_BIN[[reg]]
    lag_vars <- lag_vars_for(k)

    # Which outcomes still need fitting for this window and region?
    todo <- Filter(function(out) {
      f <- file.path(dir_out("models"),
                     sprintf("sens_lag%02d_%s_%s.rds", k, reg, out))
      if (SKIP_EXISTING && file.exists(f) && file.size(f) > 0) {
        results[[paste(k, reg, out)]] <<- readRDS(f)
        FALSE
      } else TRUE
    }, SENS_OUTCOMES)

    if (!length(todo)) {
      say("=== ", k, "-day window: ", REGION_LABELS[[reg]],
          " - all outcomes already done")
      next
    }

    say("=== ", k, "-day window: ", REGION_LABELS[[reg]], " ===")
    dt_reg <- prepare_region(panel, reg, k_max = k)

    for (out in todo) {
      say("  outcome: ", out)
      t0  <- Sys.time()
      fit <- cache_vcov(fit_dlm(dt_reg, out, lag_vars, ref_bin, interact = TRUE))

      res <- purrr::map_dfr(FHS_LEVELS, function(cv)
        purrr::map_dfr(BIN_LEVELS, ~ pred_at_coverage(fit, lag_vars, .x,
                                                      ref_bin, cv))) %>%
        mutate(region = REGION_LABELS[[reg]], outcome = out, lag_window = k,
               ref_bin = ref_bin, n_obs = model_nobs(fit))

      saveRDS(res, file.path(dir_out("models"),
                             sprintf("sens_lag%02d_%s_%s.rds", k, reg, out)))
      results[[paste(k, reg, out)]] <- res
      say(sprintf("    done in %.1f minutes",
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
      rm(fit); gc()
    }
    rm(dt_reg); gc()
  }
}

# ------------------------------------------------------------------------------
# The 30-day window is the main model.  Read it back rather than refitting, so
# the third column of Figures F3 and F4 is exactly the estimate the paper
# reports elsewhere.
# ------------------------------------------------------------------------------
main30 <- list()
for (reg in SENS_REGIONS) for (out in SENS_OUTCOMES) {
  f <- file.path(dir_out("csv"),
                 sprintf("predicted_by_coverage_%s_%s.csv", reg, out))
  if (!file.exists(f)) {
    warning("Missing ", f, " - run 03_main_dlm.R first; the 30-day column of ",
            "the lag-window figure will be empty.", call. = FALSE)
    next
  }
  main30[[paste(reg, out)]] <- readr::read_csv(f, show_col_types = FALSE) %>%
    mutate(lag_window = LAG_MAX, outcome = out) %>%
    select(any_of(c("temp_bin", "coverage", "estimate", "se", "lwr", "upr",
                    "region", "outcome", "lag_window", "ref_bin")))
}

res <- bind_rows(bind_rows(results), bind_rows(main30)) %>%
  mutate(temp_bin = factor(as.character(temp_bin), levels = BIN_LEVELS),
         coverage_label = factor(paste0("FHS = ", coverage),
                                 levels = paste0("FHS = ", FHS_LEVELS)),
         region = factor(region, levels = unname(REGION_LABELS)),
         age_group = factor(age_label(outcome),
                            levels = age_label(SENS_OUTCOMES)),
         window = factor(paste0(lag_window, "-day cumulative effect"),
                         levels = paste0(sort(c(LAG_MAX_SENS, LAG_MAX)),
                                         "-day cumulative effect")))

write_out(res, "csv", "appendixF_lag_windows.csv")

for (r in levels(droplevels(res$region))) {
  df <- res %>% filter(region == r)
  p <- ggplot(df, aes(temp_bin, estimate, colour = coverage_label)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(position = position_dodge(0.5), size = 2) +
    geom_errorbar(aes(ymin = lwr, ymax = upr),
                  position = position_dodge(0.5), width = 0.2) +
    scale_colour_manual(values = FHS_COLOURS) +
    facet_grid(age_group ~ window, scales = "free_y") +
    theme_minimal() +
    labs(x = "Temperature bin (ºC)",
         y = "Cumulative impact on mortality rate (deaths per 100 000)",
         colour = "FHS coverage",
         caption = paste0("The ", LAG_MAX, "-day column is the main model. ",
                          "Shorter windows capture less mortality ",
                          "displacement, so a larger short-window effect is ",
                          "evidence of harvesting rather than of instability."))
  save_fig(p, sprintf("appendixF_lag_windows_%s.png", gsub("[^A-Za-z]", "", r)),
           width = 11, height = 8)
}

say("09_sens_lag_windows.R done")
