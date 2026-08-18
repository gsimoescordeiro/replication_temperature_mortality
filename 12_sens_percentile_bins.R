# ==============================================================================
# 12_sens_percentile_bins.R
# ------------------------------------------------------------------------------
# Appendix F(v), Figures F9 (Southeast) and F10 (South): region-specific
# percentile temperature bins instead of absolute ones.  A day above 30 ºC is
# routine in the North and exceptional in the South, so absolute bins mean
# different things in different places.  Percentile bins (<p10, p10-p25,
# p25-p75, p75-p90, >p90 of each region's own distribution) hold the LOCAL
# rarity of the exposure constant instead of the temperature itself.
#
# INTERPRETATION NOTE: the two binnings are not measuring the same exposure.
# The Southeast's p90 is about 26.7 ºC and the South's about 25.6 ºC, both well
# below the 30 ºC threshold used in the main analysis, so the top percentile bin
# contains many days that the absolute analysis treats as unremarkable.  A
# weaker gradient here is expected and is not by itself evidence against the
# main result.  The cut-points are written out so this can be stated.
#
# CHANGE FROM THE ORIGINAL SCRIPT: coefficients are now selected by exact name.
# The original selected them with grep() on the bin label, and the labels were
# "ext_cold", "cold", "mild", "hot", "ext_hot" - so the pattern "cold" also
# matched "ext_cold" and "hot" also matched "ext_hot".  The reported cold and
# hot effects were therefore the sum of two bins, with a standard error for that
# wrong combination.  Exact matching is in bin_coef_index() in 01_helpers.R and
# is used by every script, so the same mistake cannot recur.
#
# STRUCTURE: one figure per region, panels = age group, matching the submitted
# appendix.  Scope comes from SENS_REGIONS and SENS_OUTCOMES in 00_config.R.
#
# Inputs : panel (via get_panel())
# Outputs: csv/appendixF_percentile_bins.csv
#          csv/appendixF_percentile_cutpoints.csv
#          figures/appendixF_percentile_bins_<region>.png
# ==============================================================================

panel <- get_panel()

results <- list(); cutpoints <- list()

for (reg in SENS_REGIONS) {

  lag_vars <- lag_vars_for(LAG_MAX)

  todo <- Filter(function(out) {
    f <- file.path(dir_out("models"), sprintf("sens_pct_%s_%s.rds", reg, out))
    if (SKIP_EXISTING && file.exists(f) && file.size(f) > 0) {
      st <- readRDS(f)
      results[[paste(reg, out)]] <<- st$res
      cutpoints[[reg]] <<- st$cut
      FALSE
    } else TRUE
  }, SENS_OUTCOMES)

  if (!length(todo)) {
    say("=== Percentile bins: ", REGION_LABELS[[reg]],
        " - all outcomes already done")
    next
  }

  say("=== Percentile bins: ", REGION_LABELS[[reg]], " ===")
  dt_reg <- prepare_region(panel, reg, k_max = LAG_MAX, bins = "percentile")

  cut <- tibble(
    region = REGION_LABELS[[reg]],
    percentile = c("p10", "p25", "p75", "p90"),
    temperature_c = as.numeric(quantile(dt_reg$DAT,
                                        c(.10, .25, .75, .90), na.rm = TRUE)))
  cutpoints[[reg]] <- cut

  for (out in todo) {
    say("  outcome: ", out)
    t0  <- Sys.time()
    fit <- cache_vcov(fit_dlm(dt_reg, out, lag_vars, PCT_REF_BIN, interact = TRUE))

    res <- purrr::map_dfr(FHS_LEVELS, function(cv)
      purrr::map_dfr(PCT_BIN_LEVELS, ~ pred_at_coverage(fit, lag_vars, .x,
                                                        PCT_REF_BIN, cv))) %>%
      mutate(region = REGION_LABELS[[reg]], outcome = out,
             ref_bin = PCT_REF_BIN, n_obs = model_nobs(fit))

    saveRDS(list(res = res, cut = cut),
            file.path(dir_out("models"), sprintf("sens_pct_%s_%s.rds", reg, out)))
    results[[paste(reg, out)]] <- res
    say(sprintf("    done in %.1f minutes",
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    rm(fit); gc()
  }
  rm(dt_reg); gc()
}

if (!length(results) || !nrow(bind_rows(results)))
  stop("No results were produced. Check the messages above: every model for 12_sens_percentile_bins.R was skipped or failed.", call. = FALSE)

res <- bind_rows(results) %>%
  mutate(temp_bin = factor(as.character(temp_bin), levels = PCT_BIN_LEVELS),
         coverage_label = factor(paste0("FHS = ", coverage),
                                 levels = paste0("FHS = ", FHS_LEVELS)),
         region = factor(region, levels = unname(REGION_LABELS)),
         age_group = factor(age_label(outcome),
                            levels = age_label(SENS_OUTCOMES)))

write_out(res, "csv", "appendixF_percentile_bins.csv")
write_out(bind_rows(cutpoints), "csv", "appendixF_percentile_cutpoints.csv")

cuts <- bind_rows(cutpoints)

for (r in levels(droplevels(res$region))) {
  df  <- res %>% filter(region == r)
  p90 <- cuts %>% filter(region == r, percentile == "p90") %>% pull(temperature_c)
  p <- ggplot(df, aes(temp_bin, estimate, colour = coverage_label)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(position = position_dodge(0.5), size = 2) +
    geom_errorbar(aes(ymin = lwr, ymax = upr),
                  position = position_dodge(0.5), width = 0.2) +
    scale_colour_manual(values = FHS_COLOURS) +
    facet_wrap(~ age_group, scales = "free_y", ncol = 1) +
    theme_minimal() +
    labs(x = "Temperature percentile bin (region-specific)",
         y = paste0(LAG_MAX, "-day cumulative impact (deaths per 100 000)"),
         colour = "FHS coverage",
         caption = paste0("Reference bin: ", PCT_REF_BIN, ". Cut-points are ",
                          "computed on this region's own distribution of daily ",
                          "mean temperature; its 90th percentile is ",
                          if (length(p90)) sprintf("%.1f ºC", p90) else "reported",
                          ", so the top bin is not the same exposure as the ",
                          ">30 ºC bin in the main analysis."))
  save_fig(p, sprintf("appendixF_percentile_bins_%s.png",
                      gsub("[^A-Za-z]", "", r)), width = 9, height = 9)
}

say("12_sens_percentile_bins.R done")
