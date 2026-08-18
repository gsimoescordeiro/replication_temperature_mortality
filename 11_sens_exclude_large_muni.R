# ==============================================================================
# 11_sens_exclude_large_muni.R
# ------------------------------------------------------------------------------
# Appendix F(iv), Figures F7 (Southeast) and F8 (South): excluding large
# municipalities.  Population-weighted estimates are dominated by a handful of
# very large cities, where a single centroid also represents the temperature of
# a large and heterogeneous area, and where the urban heat island makes the ERA5
# grid value a poorer measure of what residents experienced.  Dropping them
# tests whether the result is a big-city artefact.
#
# DEFINITION: a municipality is excluded if its total population EVER exceeded
# the threshold during 2000-2019, so a municipality is either in or out for the
# whole period.  Excluding it only in the years it exceeded the threshold would
# make the panel composition change with city growth.
#
# The size filter requires a FINITE maximum population.  max(pop_total,
# na.rm = TRUE) returns -Inf for a municipality whose population is missing in
# every year, and -Inf <= 500000 is TRUE, so a plain threshold would silently
# keep exactly the municipalities it cannot evaluate.
#
# STRUCTURE: one figure per region, panels = age group, matching the appendix.  Scope comes from SENS_REGIONS and SENS_OUTCOMES in 00_config.R.
#
# Inputs : panel (via get_panel())
# Outputs: csv/appendixF_exclude_large_muni.csv
#          figures/appendixF_exclude_large_muni_<region>.png
# ==============================================================================

POP_THRESHOLD <- 500000

panel <- get_panel()

pop_max <- panel %>%
  group_by(code_muni) %>%
  summarise(max_pop = suppressWarnings(max(pop_total, na.rm = TRUE)),
            .groups = "drop")

n_no_pop <- sum(!is.finite(pop_max$max_pop))
if (n_no_pop > 0)
  say(sprintf("%d municipalities have no population data at all and are ",
              n_no_pop), "excluded from this sensitivity analysis.")

keep <- pop_max %>% filter(is.finite(max_pop), max_pop <= POP_THRESHOLD)

say(sprintf("Keeping %d of %d municipalities (population never above %s)",
            nrow(keep), nrow(pop_max), format(POP_THRESHOLD, big.mark = " ")))

panel_small <- panel %>% semi_join(keep, by = "code_muni")

results <- list()

for (reg in SENS_REGIONS) {

  ref_bin  <- REF_BIN[[reg]]
  lag_vars <- lag_vars_for(LAG_MAX)

  todo <- Filter(function(out) {
    f <- file.path(dir_out("models"),
                   sprintf("sens_nolarge_%s_%s.rds", reg, out))
    if (SKIP_EXISTING && file.exists(f) && file.size(f) > 0) {
      results[[paste(reg, out)]] <<- readRDS(f); FALSE
    } else TRUE
  }, SENS_OUTCOMES)

  if (!length(todo)) {
    say("=== Excluding large municipalities: ", REGION_LABELS[[reg]],
        " - all outcomes already done")
    next
  }

  say("=== Excluding large municipalities: ", REGION_LABELS[[reg]], " ===")
  dt_reg <- prepare_region(panel_small, reg, k_max = LAG_MAX)
  if (!nrow(dt_reg)) { rm(dt_reg); next }

  for (out in todo) {
    say("  outcome: ", out)
    t0  <- Sys.time()
    fit <- cache_vcov(fit_dlm(dt_reg, out, lag_vars, ref_bin, interact = TRUE))

    res <- purrr::map_dfr(FHS_LEVELS, function(cv)
      purrr::map_dfr(BIN_LEVELS, ~ pred_at_coverage(fit, lag_vars, .x,
                                                    ref_bin, cv))) %>%
      mutate(region = REGION_LABELS[[reg]], outcome = out, ref_bin = ref_bin,
             n_obs = model_nobs(fit), n_muni = n_distinct(dt_reg$code_muni),
             pop_threshold = POP_THRESHOLD)

    saveRDS(res, file.path(dir_out("models"),
                           sprintf("sens_nolarge_%s_%s.rds", reg, out)))
    results[[paste(reg, out)]] <- res
    say(sprintf("    done in %.1f minutes",
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    rm(fit); gc()
  }
  rm(dt_reg); gc()
}

if (!length(results) || !nrow(bind_rows(results)))
  stop("No results were produced. Check the messages above: every model for 11_sens_exclude_large_muni.R was skipped or failed.", call. = FALSE)

res <- bind_rows(results) %>%
  mutate(temp_bin = factor(as.character(temp_bin), levels = BIN_LEVELS),
         coverage_label = factor(paste0("FHS = ", coverage),
                                 levels = paste0("FHS = ", FHS_LEVELS)),
         region = factor(region, levels = unname(REGION_LABELS)),
         age_group = factor(age_label(outcome),
                            levels = age_label(SENS_OUTCOMES)))

write_out(res, "csv", "appendixF_exclude_large_muni.csv")

for (r in levels(droplevels(res$region))) {
  df <- res %>% filter(region == r)
  p <- ggplot(df, aes(temp_bin, estimate, colour = coverage_label)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(position = position_dodge(0.5), size = 2) +
    geom_errorbar(aes(ymin = lwr, ymax = upr),
                  position = position_dodge(0.5), width = 0.2) +
    scale_colour_manual(values = FHS_COLOURS) +
    facet_wrap(~ age_group, scales = "free_y", ncol = 1) +
    theme_minimal() +
    labs(x = "Temperature bin (ºC)",
         y = paste0(LAG_MAX, "-day cumulative impact (deaths per 100 000)"),
         colour = "FHS coverage",
         caption = paste0("Municipalities whose population ever exceeded ",
                          format(POP_THRESHOLD, big.mark = " "),
                          " are excluded throughout the period."))
  save_fig(p, sprintf("appendixF_exclude_large_muni_%s.png",
                      gsub("[^A-Za-z]", "", r)), width = 9, height = 9)
}

say("11_sens_exclude_large_muni.R done")
