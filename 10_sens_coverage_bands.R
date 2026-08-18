# ==============================================================================
# 10_sens_coverage_bands.R
# ------------------------------------------------------------------------------
# Appendix F(iii), Figures F5 (Southeast) and F6 (South): split-sample
# estimation by FHS coverage band, as a check on the linearity that the
# interacted model imposes.  The interaction assumes the effect of coverage is
# linear in coverage; estimating the temperature-mortality relationship
# separately within bands lets the shape differ between them.
#
# The bands are FIXED CUT-POINTS at 25% and 75% (set in FHS_BANDS), not
# empirical terciles: FHS coverage piles up at 0 and at 100, so equal-sized
# groups would have cut-points that differ by region and would be hard to
# interpret.  Describe them as coverage bands, not terciles - the submitted
# appendix caption says "terciles" and needs changing.
#
# The split is applied at MUNICIPALITY-YEAR level, so a municipality that
# expanded coverage contributes its early years to the low band and its later
# years to the high band.  That is intended - the comparison is between
# municipality-years at different coverage levels - but it means the bands also
# differ in calendar composition, which should be said when the results are
# discussed.  appendixF_coverage_band_composition.csv reports the mean year and
# the sample share of each band for exactly that purpose.
#
# Models here are estimated WITHOUT the coverage interaction: within a band,
# coverage barely varies, so the interaction is not identified.
#
# STRUCTURE: one figure per region, panels = age group, matching the submitted
# appendix.  Scope comes from SENS_REGIONS and SENS_OUTCOMES in 00_config.R.
#
# Inputs : panel (via get_panel())
# Outputs: csv/appendixF_coverage_bands.csv
#          csv/appendixF_coverage_band_composition.csv
#          figures/appendixF_coverage_bands_<region>.png
# ==============================================================================

panel <- get_panel()

band_of <- function(x) {
  dplyr::case_when(
    x <  FHS_BANDS$low[2]  ~ sprintf("FHS < %g%%", FHS_BANDS$low[2]),
    x <= FHS_BANDS$mid[2]  ~ sprintf("%g%% <= FHS <= %g%%",
                                     FHS_BANDS$mid[1], FHS_BANDS$mid[2]),
    x >  FHS_BANDS$high[1] ~ sprintf("FHS > %g%%", FHS_BANDS$high[1])
  )
}

BAND_LEVELS <- c(sprintf("FHS < %g%%", FHS_BANDS$low[2]),
                 sprintf("%g%% <= FHS <= %g%%", FHS_BANDS$mid[1], FHS_BANDS$mid[2]),
                 sprintf("FHS > %g%%", FHS_BANDS$high[1]))

results <- list(); shares <- list()

for (reg in SENS_REGIONS) {

  ref_bin  <- REF_BIN[[reg]]
  lag_vars <- lag_vars_for(LAG_MAX)

  todo <- Filter(function(out) {
    f <- file.path(dir_out("models"), sprintf("sens_band_%s_%s.rds", reg, out))
    if (SKIP_EXISTING && file.exists(f) && file.size(f) > 0) {
      st <- readRDS(f)
      results[[paste(reg, out)]] <<- st$res
      shares[[reg]] <<- st$share
      FALSE
    } else TRUE
  }, SENS_OUTCOMES)

  if (!length(todo)) {
    say("=== Coverage bands: ", REGION_LABELS[[reg]],
        " - all outcomes already done")
    next
  }

  say("=== Coverage bands: ", REGION_LABELS[[reg]], " ===")

  dt_reg <- prepare_region(panel, reg, k_max = LAG_MAX) %>%
    mutate(fhs_band = band_of(cob_esf))

  # How much of the sample sits in each band, and over which years - needed to
  # interpret the comparison.
  share <- dt_reg %>%
    group_by(fhs_band) %>%
    summarise(muni_days = n(), n_muni = n_distinct(code_muni),
              mean_year = mean(anio), .groups = "drop") %>%
    mutate(region = REGION_LABELS[[reg]], share = muni_days / sum(muni_days))
  shares[[reg]] <- share

  for (out in todo) {

    say("  outcome: ", out)
    t0 <- Sys.time()
    per_band <- list()

    for (bnd in BAND_LEVELS) {
      dsub <- dt_reg %>% filter(fhs_band == bnd)
      if (nrow(dsub) < 1000) {
        say("    skipping band ", bnd, " (only ", nrow(dsub), " rows)")
        next
      }
      say("    band: ", bnd, " (", format(nrow(dsub), big.mark = " "), " rows)")

      fit <- cache_vcov(fit_dlm(dsub, out, lag_vars, ref_bin, interact = FALSE))

      per_band[[bnd]] <-
        purrr::map_dfr(BIN_LEVELS, ~ cum_effect(fit, lag_vars, .x, ref_bin)) %>%
        mutate(region = REGION_LABELS[[reg]], outcome = out, fhs_band = bnd,
               ref_bin = ref_bin, n_obs = model_nobs(fit),
               n_muni = n_distinct(dsub$code_muni))

      rm(dsub, fit); gc()
    }

    res <- bind_rows(per_band)
    if (!nrow(res)) {
      warning("No band had enough rows for ", out, " in ",
              REGION_LABELS[[reg]], "; nothing saved for this cell.",
              call. = FALSE)
      next
    }
    saveRDS(list(res = res, share = share),
            file.path(dir_out("models"), sprintf("sens_band_%s_%s.rds", reg, out)))
    results[[paste(reg, out)]] <- res
    say(sprintf("    done in %.1f minutes",
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  rm(dt_reg); gc()
}

if (!length(results) || !nrow(bind_rows(results)))
  stop("No results were produced. Check the messages above: every model for 10_sens_coverage_bands.R was skipped or failed.", call. = FALSE)

res <- bind_rows(results) %>%
  mutate(temp_bin  = factor(as.character(temp_bin), levels = BIN_LEVELS),
         fhs_band  = factor(fhs_band, levels = BAND_LEVELS),
         region    = factor(region, levels = unname(REGION_LABELS)),
         age_group = factor(age_label(outcome), levels = age_label(SENS_OUTCOMES)))

write_out(res, "csv", "appendixF_coverage_bands.csv")
write_out(bind_rows(shares), "csv", "appendixF_coverage_band_composition.csv")

# The comparison the Results paragraph quotes: the >30 C bin in the lowest and
# highest bands, by age group and region.
say("--- >30 C bin by coverage band ---")
print(as.data.frame(
  res %>% filter(temp_bin == ">30") %>%
    mutate(cell = sprintf("%.3f (%.3f to %.3f)", estimate, lwr, upr)) %>%
    select(region, age_group, fhs_band, cell) %>%
    tidyr::pivot_wider(names_from = fhs_band, values_from = cell)))

for (r in levels(droplevels(res$region))) {
  df <- res %>% filter(region == r)
  p <- ggplot(df, aes(temp_bin, estimate, colour = fhs_band)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(position = position_dodge(0.5), size = 2) +
    geom_errorbar(aes(ymin = lwr, ymax = upr),
                  position = position_dodge(0.5), width = 0.2) +
    scale_colour_manual(values = setNames(c("red", "black", "blue"), BAND_LEVELS)) +
    facet_wrap(~ age_group, scales = "free_y", ncol = 1) +
    theme_minimal() +
    labs(x = "Temperature bin (ºC)",
         y = paste0(LAG_MAX, "-day cumulative impact (deaths per 100 000)"),
         colour = "Coverage band",
         caption = paste0("Fixed cut-points at 25% and 75%, not empirical ",
                          "terciles. Bands are split at municipality-year ",
                          "level and so differ in calendar composition; see ",
                          "appendixF_coverage_band_composition.csv."))
  save_fig(p, sprintf("appendixF_coverage_bands_%s.png",
                      gsub("[^A-Za-z]", "", r)), width = 9, height = 9)
}

say("10_sens_coverage_bands.R done")
