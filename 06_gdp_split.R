# ==============================================================================
# 06_gdp_split.R
# ------------------------------------------------------------------------------
# Appendix D, Figures D1 (Southeast) and D2 (South): socioeconomic
# heterogeneity.  The interacted model is re-estimated on the municipalities
# below and above the median GDP per capita within each macro-region.
#
# Municipalities are assigned to a group ONCE, using their average GDP per
# capita over the study period, rather than year by year.  A year-by-year split
# would let a municipality move between groups as the economy grows, so the two
# subsamples would differ in calendar composition as well as in income - the
# comparison would then confound income with period.
#
# STRUCTURE: like the submitted appendix, one figure per region laid out as a
# grid of GDP group (columns) by age group (rows).  Scope comes from
# SENS_REGIONS and SENS_OUTCOMES in 00_config.R.
#
# The submitted manuscript describes these results as underpowered.  The script
# writes every estimate and interval so that statement can be checked rather
# than assumed - on the corrected all-cause run it does not hold for the
# Southeast, where the above-median group shows a clear gradient.
#
# Inputs : panel (via get_panel()), GDP file from 00_config.R DATA_FILES$gdp
# Outputs: csv/appendixD_gdp_split_<region>.csv
#          csv/appendixD_gdp_split_ALL.csv
#          figures/appendixD_gdp_split_<region>.png
# ==============================================================================

gdp_file <- data_path("gdp")
if (!file.exists(gdp_file))
  stop("GDP file not found: ", gdp_file,
       "\nSet DATA_FILES$gdp in 00_config.R, or skip this script.")

gdp_raw <- haven::read_dta(gdp_file)

# The GDP table must be keyed on municipality-YEAR.  Joining an annual table on
# municipality alone silently multiplies every panel row by the number of years
# present, which corrupts every statistic computed afterwards.
year_col <- intersect(c("anio", "ano", "year"), names(gdp_raw))
gdp_col  <- intersect(c("gdp", "pib", "gdp_total"), names(gdp_raw))
if (!length(gdp_col)) stop("No GDP column found in ", gdp_file)

if (length(year_col)) {
  gdp <- gdp_raw %>%
    transmute(code_muni = as.numeric(code_muni),
              anio = as.numeric(.data[[year_col[1]]]),
              gdp  = as.numeric(.data[[gdp_col[1]]]))
  assert_unique_key(gdp, c("code_muni", "anio"), "GDP file")
} else {
  warning("GDP file has no year column; treating it as time-invariant.",
          call. = FALSE)
  gdp <- gdp_raw %>%
    transmute(code_muni = as.numeric(code_muni),
              gdp = as.numeric(.data[[gdp_col[1]]]))
  assert_unique_key(gdp, "code_muni", "GDP file")
}

panel <- get_panel()

GDP_GROUPS <- c("above median", "at or below median")

results <- list()

for (reg in SENS_REGIONS) {

  ref_bin  <- REF_BIN[[reg]]
  lag_vars <- lag_vars_for(LAG_MAX)

  todo <- Filter(function(out) {
    f <- file.path(dir_out("models"), sprintf("sens_gdp_%s_%s.rds", reg, out))
    if (SKIP_EXISTING && file.exists(f) && file.size(f) > 0) {
      results[[paste(reg, out)]] <<- readRDS(f); FALSE
    } else TRUE
  }, SENS_OUTCOMES)

  if (!length(todo)) {
    say("=== GDP split: ", REGION_LABELS[[reg]],
        " - all outcomes already done")
    next
  }

  say("=== GDP split: ", REGION_LABELS[[reg]], " ===")
  dt_reg <- prepare_region(panel, reg, k_max = LAG_MAX)

  # One GDP per capita per municipality: average over the study period.
  gdp_muni <- dt_reg %>%
    distinct(code_muni, anio, pop_20_plus) %>%
    { if (length(year_col)) left_join(., gdp, by = c("code_muni", "anio"))
      else left_join(., gdp, by = "code_muni") } %>%
    mutate(gdp_pc = gdp / pmax(pop_20_plus, 1)) %>%
    group_by(code_muni) %>%
    summarise(gdp_pc = mean(gdp_pc, na.rm = TRUE), .groups = "drop") %>%
    filter(is.finite(gdp_pc))

  cut_point <- median(gdp_muni$gdp_pc)
  say(sprintf("  median GDP per capita: %.1f (%d municipalities classified)",
              cut_point, nrow(gdp_muni)))

  gdp_muni <- gdp_muni %>%
    mutate(gdp_group = if_else(gdp_pc <= cut_point, "at or below median",
                               "above median"))

  dt_reg <- dt_reg %>% inner_join(gdp_muni, by = "code_muni")

  for (out in todo) {

    say("  outcome: ", out)
    t0 <- Sys.time()
    per_grp <- list()

    for (grp in intersect(GDP_GROUPS, unique(gdp_muni$gdp_group))) {
      dsub <- dt_reg %>% filter(gdp_group == grp)
      if (!nrow(dsub)) next
      say("    group: ", grp, " (", n_distinct(dsub$code_muni),
          " municipalities)")

      fit <- cache_vcov(fit_dlm(dsub, out, lag_vars, ref_bin, interact = TRUE))

      per_grp[[grp]] <- purrr::map_dfr(FHS_LEVELS, function(cv)
        purrr::map_dfr(BIN_LEVELS, ~ pred_at_coverage(fit, lag_vars, .x,
                                                      ref_bin, cv))) %>%
        mutate(region = REGION_LABELS[[reg]], outcome = out, gdp_group = grp,
               ref_bin = ref_bin, n_muni = n_distinct(dsub$code_muni),
               median_gdp_pc = cut_point)
      rm(dsub, fit); gc()
    }

    res <- bind_rows(per_grp)
    if (!nrow(res)) {
      warning("No GDP group had any rows for ", out, " in ",
              REGION_LABELS[[reg]], "; nothing saved for this cell.",
              call. = FALSE)
      next
    }
    saveRDS(res, file.path(dir_out("models"),
                           sprintf("sens_gdp_%s_%s.rds", reg, out)))
    results[[paste(reg, out)]] <- res
    say(sprintf("    done in %.1f minutes",
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }

  write_out(bind_rows(results[grepl(paste0("^", reg, " "), names(results))]),
            "csv", sprintf("appendixD_gdp_split_%s.csv", reg))
  rm(dt_reg); gc()
}

if (!length(results) || !nrow(bind_rows(results)))
  stop("No results were produced. Check the messages above: every model for 06_gdp_split.R was skipped or failed.", call. = FALSE)

res <- bind_rows(results) %>%
  mutate(temp_bin = factor(as.character(temp_bin), levels = BIN_LEVELS),
         coverage_label = factor(paste0("FHS = ", coverage),
                                 levels = paste0("FHS = ", FHS_LEVELS)),
         region = factor(region, levels = unname(REGION_LABELS)),
         gdp_panel = factor(ifelse(gdp_group == "above median",
                                   "GDP per capita Above Median",
                                   "GDP per capita Below Median"),
                            levels = c("GDP per capita Above Median",
                                       "GDP per capita Below Median")),
         age_group = factor(age_label(outcome),
                            levels = age_label(SENS_OUTCOMES)))

write_out(res, "csv", "appendixD_gdp_split_ALL.csv")

# The comparison the Results paragraph refers to: the >30 C bin at 0% and 90%
# coverage, by GDP group and age.
say("--- >30 C bin, 0% vs 90% coverage, by GDP group ---")
print(as.data.frame(
  res %>% filter(temp_bin == ">30", coverage %in% c(0, FHS_HIGH)) %>%
    mutate(cell = sprintf("%.3f (%.3f to %.3f)", estimate, lwr, upr)) %>%
    select(region, age_group, gdp_panel, coverage, cell) %>%
    tidyr::pivot_wider(names_from = coverage, values_from = cell,
                       names_prefix = "FHS ")))

for (r in levels(droplevels(res$region))) {
  df <- res %>% filter(region == r)
  p <- ggplot(df, aes(temp_bin, estimate, colour = coverage_label)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(position = position_dodge(0.5), size = 2) +
    geom_errorbar(aes(ymin = lwr, ymax = upr),
                  position = position_dodge(0.5), width = 0.2) +
    scale_colour_manual(values = FHS_COLOURS) +
    facet_grid(age_group ~ gdp_panel, scales = "free_y") +
    theme_minimal() +
    labs(x = "Temperature bin (ºC)",
         y = paste0(LAG_MAX, "-day cumulative impact (deaths per 100 000)"),
         colour = "FHS coverage",
         caption = paste0("Municipalities are assigned to a GDP group once, on ",
                          "their average GDP per capita over 2000-2019, so the ",
                          "two subsamples do not differ in calendar ",
                          "composition."))
  save_fig(p, sprintf("appendixD_gdp_split_%s.png", gsub("[^A-Za-z]", "", r)),
           width = 10, height = 9)
}

say("06_gdp_split.R done")
