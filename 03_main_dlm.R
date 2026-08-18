# ==============================================================================
# 03_main_dlm.R
# ------------------------------------------------------------------------------
# Main analysis: the distributed-lag model of daily mortality on temperature
# bins, with and without the FHS coverage interaction, estimated separately for
# each macro-region and each outcome.
#
# This script produces:
#   * Figure 2  - cumulative temperature effect by bin at FHS coverage of
#                 0%, 50% and 90% (all-cause, adults 20+), one panel per region
#   * Appendix B - the same relationship from the model WITHOUT the interaction
#   * Appendix C - cause- and age-specific versions of the interacted results
#   * the fitted models, saved for reuse by 05_table1.R and 07_deaths_averted.R
#
# Predictions are evaluated at the three coverage levels shown in Figure 2 -
# 0%, 50% and 90% - so every series in the figure comes straight from the model.
# The 50% interval is sqrt(L'VL) at loading 50, which is NOT what interpolating
# between the 0% and 90% intervals would give: the standard error is a
# square-root quadratic in coverage, not a linear function of it.  Estimating
# each level directly is what makes the figure's intervals correct, and it means
# no post-processing step stands between the model and the plot.
#
# Inputs : panel (via get_panel())
# Outputs: csv/cumulative_nointeraction_<region>_<outcome>.csv
#          csv/predicted_by_coverage_<region>_<outcome>.csv
#          models/fit_<region>_<outcome>.rds
#          figures/fig2_temperature_by_coverage_<outcome>.png
# ==============================================================================

panel <- get_panel()

OUTCOMES_MAIN <- OUTCOMES        # all outcomes; Figure 2 uses rate_total

all_pred <- list()
all_cum  <- list()

for (reg in REGIONS) {

  ref_bin  <- REF_BIN[[reg]]
  lag_vars <- lag_vars_for(LAG_MAX)

  say("=== Region ", REGION_LABELS[[reg]], " (reference bin ", ref_bin, ") ===")
  dt_reg <- prepare_region(panel, reg, k_max = LAG_MAX)

  for (out in OUTCOMES_MAIN) {

    # ---- Resume support -------------------------------------------------------
    # This script is the long pole: a full run is many hours, and an interrupted
    # one should not repeat work already on disk.  A region-outcome is skipped
    # only when all three of its outputs exist and are non-empty; the CSVs are
    # read back so the combined files and Figure 2 are still complete.
    outs <- main_dlm_outputs(reg, out)
    if (SKIP_EXISTING && outputs_exist(outs)) {
      say("  outcome: ", out, " - already done, reusing files on disk")
      all_cum[[paste(reg, out)]]  <- readr::read_csv(outs[1], show_col_types = FALSE)
      all_pred[[paste(reg, out)]] <- readr::read_csv(outs[2], show_col_types = FALSE)
      next
    }

    say("  outcome: ", out)
    t_out <- Sys.time()

    # ---- Model without the interaction (Appendix B) --------------------------
    # This is the temperature-mortality relationship pooled across all observed
    # coverage levels, and is what the literature usually reports.
    fit_bin <- cache_vcov(fit_dlm(dt_reg, out, lag_vars, ref_bin, interact = FALSE))

    cum <- purrr::map_dfr(BIN_LEVELS, ~ cum_effect(fit_bin, lag_vars, .x, ref_bin)) %>%
      mutate(region = REGION_LABELS[[reg]], outcome = out, ref_bin = ref_bin)

    write_out(cum, "csv", sprintf("cumulative_nointeraction_%s_%s.csv", reg, out))
    all_cum[[paste(reg, out)]] <- cum

    # ---- Interacted model ----------------------------------------------------
    fit_int <- cache_vcov(fit_dlm(dt_reg, out, lag_vars, ref_bin, interact = TRUE))

    # Saved slim: coefficients and the cluster-robust covariance only.  That is
    # everything 05_table1.R needs, and it keeps each file at a megabyte or two
    # instead of the hundreds of megabytes a full fixest object would occupy.
    saveRDS(list(fit = slim_fit(fit_int), ref_bin = ref_bin, lag_vars = lag_vars,
                 region = reg, outcome = out,
                 n_obs = model_nobs(fit_int)),
            file.path(dir_out("models"), sprintf("fit_%s_%s.rds", reg, out)))

    pred <- purrr::map_dfr(FHS_LEVELS, function(cv)
      purrr::map_dfr(BIN_LEVELS, ~ pred_at_coverage(fit_int, lag_vars, .x,
                                                    ref_bin, cv))) %>%
      mutate(region = REGION_LABELS[[reg]], outcome = out, ref_bin = ref_bin,
             coverage_label = factor(paste0("FHS = ", coverage),
                                     levels = paste0("FHS = ", FHS_LEVELS)),
             temp_bin = factor(temp_bin, levels = BIN_LEVELS))

    write_out(pred, "csv", sprintf("predicted_by_coverage_%s_%s.csv", reg, out))
    all_pred[[paste(reg, out)]] <- pred

    say(sprintf("  %s finished in %.1f minutes", out,
                as.numeric(difftime(Sys.time(), t_out, units = "mins"))))
  }

  rm(dt_reg); gc()
}

pred_all <- bind_rows(all_pred) %>%
  mutate(temp_bin = factor(as.character(temp_bin), levels = BIN_LEVELS),
         coverage_label = factor(paste0("FHS = ", coverage),
                                 levels = paste0("FHS = ", FHS_LEVELS)))
cum_all  <- bind_rows(all_cum)
write_out(pred_all, "csv", "predicted_by_coverage_ALL.csv")
write_out(cum_all,  "csv", "cumulative_nointeraction_ALL.csv")

# ------------------------------------------------------------------------------
# Figure 2 - one panel per region, all-cause mortality among adults 20+
# ------------------------------------------------------------------------------
fig2_df <- pred_all %>%
  filter(outcome == "rate_total") %>%
  mutate(region = factor(region, levels = unname(REGION_LABELS)))

p2 <- ggplot(fig2_df, aes(temp_bin, estimate, colour = coverage_label)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(position = position_dodge(0.5), size = 2.2) +
  geom_errorbar(aes(ymin = lwr, ymax = upr),
                position = position_dodge(0.5), width = 0.2) +
  scale_colour_manual(values = FHS_COLOURS) +
  facet_wrap(~ region, scales = "free_y") +
  theme_minimal() +
  labs(x = "Temperature bin (ºC)",
       y = paste0(LAG_MAX, "-day cumulative impact on mortality rate ",
                  "(deaths per 100 000)"),
       colour = "FHS coverage",
       caption = paste0("Reference bin is the region-specific minimum-mortality ",
                        "bin: 20-25 ºC in the Southeast and South, 25-30 ºC ",
                        "elsewhere. Vertical lines are 95% confidence intervals ",
                        "from municipality-clustered standard errors."))

save_fig(p2, "fig2_temperature_by_coverage_rate_total.png", width = 11, height = 7)

# Cause- and age-specific versions (Appendix C)
for (out in setdiff(OUTCOMES_MAIN, "rate_total")) {
  df <- pred_all %>% filter(outcome == out) %>%
    mutate(region = factor(region, levels = unname(REGION_LABELS)))
  if (!nrow(df)) next
  save_fig(p2 %+% df, sprintf("appendixC_temperature_by_coverage_%s.png", out),
           width = 11, height = 7)
}

say("03_main_dlm.R done")
