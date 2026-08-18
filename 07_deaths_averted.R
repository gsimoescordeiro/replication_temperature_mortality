# ==============================================================================
# 07_deaths_averted.R
# ------------------------------------------------------------------------------
# Appendix E: deaths already averted by the observed expansion of FHS coverage,
# and the further deaths that near-universal coverage would have avoided,
# 2000-2019, with delta-method 95% uncertainty intervals.
#
# ------------------------------------------------------------------------------
# THE ARITHMETIC
# ------------------------------------------------------------------------------
# The interacted model gives, for temperature bin k, a cumulative effect that is
# linear in coverage c:
#
#     excess mortality rate(k, c) = B_k + c * D_k        (per 100 000, and
#                                                         already cumulated over
#                                                         the 0..30 day lags)
#
# B_k and D_k are each the SUM of the 31 lag coefficients for that bin.  Because
# they are already cumulative, one exposure day contributes its entire 31-day
# effect once; multiplying by the number of exposure days and summing over days
# therefore does NOT double count the lagged effects.
#
# Excess deaths in municipality i, year t:
#
#     E(i,t,c) = sum_k  n(i,t,k) * (B_k + c * D_k) * pop(i,t) / 100000
#
# where n(i,t,k) is the number of days that municipality-year spent in bin k.
# The three reported quantities are differences of this expression, so the main
# effects B_k cancel and each is a linear combination of the D_k alone:
#
#     averted to date      = E(0) - E(observed) = sum_k [ -c * n * pop/1e5 ] D_k
#     additional avoidable = E(observed) - E(90) = sum_k [ (c-90) * n * pop/1e5 ] D_k
#     total avertable      = E(0) - E(90)        = sum_k [ -90 * n * pop/1e5 ] D_k
#
# and by construction averted + additional = total.
#
# Uncertainty: each quantity is w'D for a known weight vector w, so
#     SE = sqrt(w' V_D w)
# with V_D the CUMULATIVE cluster-robust covariance matrix of the D_k, that is,
# V_D[k1,k2] = sum over the lag block of the full covariance matrix.  Cross-bin
# covariances matter here because the totals sum across bins, so the diagonal
# alone would understate (or overstate) the interval.
#
# ------------------------------------------------------------------------------
# THE 50+ AGE GROUP
# ------------------------------------------------------------------------------
# Adults aged 50 and over are covered by two separately fitted models, 50-69 and
# 70+, combined as
#     Q_50+        = Q_50-69 + Q_70+
#     Var(Q_50+)   = Var(Q_50-69) + Var(Q_70+)
#
# The point estimate is exact: each model is scaled by its own population
# (pop_50_69, pop_70_plus), so the two are counts of deaths among different
# people and simply add.
#
# The variance addition assumes the two fits are uncorrelated.  They are not,
# strictly - both use the same municipalities, the same realised weather and the
# same clusters, so Cov(Q_50-69, Q_70+) is not zero and is most likely positive,
# which would make these intervals slightly too narrow.  Two extra columns make
# the size of that concern visible without changing the reported numbers:
#
#     se_worst_case  = se_50-69 + se_70+   - the widest interval any correlation
#                                            can produce (Cauchy-Schwarz, rho=1)
#     rho_critical   - the correlation at which the interval would just touch
#                      zero.  Above 1 means no correlation whatsoever can put
#                      zero inside the interval, so the independence assumption
#                      cannot have affected the conclusion.
#
# ------------------------------------------------------------------------------
# CHANGE FROM THE ORIGINAL SCRIPT
# ------------------------------------------------------------------------------
# The model now includes the `date` fixed effect, matching the main model and
# Table 1.  The original omitted it, so its estimates came from a different
# specification than the one the paper describes.
#
# NOTE ON SIGN: municipality-years already above 90% coverage contribute
# negatively to "additional avoidable", because reaching exactly 90% everywhere
# would mean reducing their coverage.  The counterfactual is coverage set to 90
# everywhere, not "at least 90".  The share of such municipality-years is
# reported below so the size of this is visible.
#
# Inputs : models/fit_<region>_<outcome>.rds (from 03_main_dlm.R; refitted here
#          if absent), panel (via get_panel()) for the exposure days
# Outputs: csv/averted_deaths_grand_total.csv
#          csv/averted_deaths_by_year.csv
#          csv/averted_deaths_by_bin.csv
#          csv/averted_deaths_by_age.csv
#          csv/averted_deaths_decomposition_check.csv
#          csv/cumulative_interaction_coefficients_<region>_<outcome>.csv
#          figures/appendixE_averted_deaths_by_year_<region>.png
# ==============================================================================

panel <- get_panel()

per_age <- list()

for (reg in REGIONS_AVERTED) {

  ref_bin  <- REF_BIN[[reg]]
  lag_vars <- lag_vars_for(LAG_MAX)
  say("=== Averted deaths: ", REGION_LABELS[[reg]], " ===")

  dt_reg <- prepare_region(panel, reg, k_max = LAG_MAX)

  for (out in OUTCOMES_AVERTED) {

    say("  outcome: ", out)

    # -- The fitted model ------------------------------------------------------
    # 03_main_dlm.R already fits exactly this specification and saves it slim
    # (coefficients + cluster-robust covariance).  Reuse it; refitting is a
    # couple of hours per region-outcome for an identical answer.
    f <- file.path(dir_out("models"), sprintf("fit_%s_%s.rds", reg, out))
    if (file.exists(f)) {
      say("    reusing fit from 03_main_dlm.R")
      fit <- readRDS(f)$fit
    } else {
      say("    no saved fit found - estimating now")
      fit <- cache_vcov(fit_dlm(dt_reg, out, lag_vars, ref_bin, interact = TRUE))
    }

    # -- Cumulative interaction coefficients and their joint covariance --------
    ci <- cum_interaction_matrix(fit, lag_vars, BIN_LEVELS, ref_bin)

    write_out(tibble(temp_bin = ci$bins, D_cum = ci$D,
                     se_D = sqrt(diag(ci$V)), n_terms = ci$n_terms,
                     region = REGION_LABELS[[reg]], outcome = out),
              "csv", sprintf("cumulative_interaction_coefficients_%s_%s.csv",
                             reg, out))

    # -- Exposure days, restricted to the estimation sample --------------------
    # Rows with an incomplete lag history (the first LAG_MAX days of each
    # municipality) or missing precipitation contribute nothing to the fit, so
    # they must not contribute exposure days to the attribution either.
    pop_var <- weight_var_for(out)
    est_sample <- dt_reg %>%
      filter(!is.na(.data[[paste0("temp_L", LAG_MAX)]]),
             !is.na(PREC), !is.na(temp_bin), .data[[pop_var]] > 0)

    say(sprintf("    estimation sample: %s of %s municipality-days",
                format(nrow(est_sample), big.mark = " "),
                format(nrow(dt_reg), big.mark = " ")))

    pct_above <- 100 * mean(est_sample$cob_esf > FHS_HIGH)
    say(sprintf("    municipality-days already above %d%% coverage: %.1f%%",
                FHS_HIGH, pct_above))

    # Aggregate to municipality-year-bin: days in the bin and the population and
    # coverage that apply (both are annual, so they are constant within cell).
    cells <- est_sample %>%
      group_by(code_muni, anio, temp_bin) %>%
      summarise(days = n(),
                pop  = dplyr::first(.data[[pop_var]]),
                cov  = dplyr::first(cob_esf),
                .groups = "drop") %>%
      filter(temp_bin != ref_bin) %>%
      mutate(temp_bin = as.character(temp_bin),
             exposure = days * pop / 1e5)

    # -- Weight vectors over bins ---------------------------------------------
    # Each quantity is sum_k w_k * D_k; build w by bin, for the grand total, by
    # year and by bin.
    make_weights <- function(df, scenario) {
      mult <- switch(scenario,
                     averted    = -df$cov,
                     additional =  df$cov - FHS_HIGH,
                     total      = rep(-FHS_HIGH, nrow(df)))
      df %>% mutate(w = exposure * mult)
    }

    # Collapse the per-cell weights to one weight per bin, then apply the delta
    # method to w'D.  Bins with no exposure get weight zero, so the vector is
    # always in the same order as ci$D and ci$V.
    lincom_over_bins <- function(bin_w) {
      w <- setNames(rep(0, length(ci$bins)), ci$bins)
      w[bin_w$temp_bin] <- bin_w$w
      d <- delta_lincom(ci$D, ci$V, unname(w[ci$bins]))
      tibble(deaths = d$estimate, se = d$se)
    }

    quantity <- function(df, scenario, group_vars = NULL) {
      agg <- make_weights(df, scenario) %>%
        group_by(across(all_of(c(group_vars, "temp_bin")))) %>%
        summarise(w = sum(w), .groups = "drop")

      if (is.null(group_vars))
        return(bind_cols(tibble(scenario = scenario), lincom_over_bins(agg)))

      keys <- dplyr::distinct(agg[group_vars])
      purrr::map_dfr(seq_len(nrow(keys)), function(i) {
        g <- dplyr::semi_join(agg, keys[i, , drop = FALSE], by = group_vars)
        bind_cols(keys[i, , drop = FALSE], tibble(scenario = scenario),
                  lincom_over_bins(g))
      })
    }

    scen <- c("averted", "additional", "total")

    res <- bind_rows(
      purrr::map_dfr(scen, ~ quantity(cells, .x)) %>%
        mutate(grouping = "grand", anio = NA_integer_, temp_bin = NA_character_),
      purrr::map_dfr(scen, ~ quantity(cells, .x, group_vars = "anio")) %>%
        mutate(grouping = "year", temp_bin = NA_character_),
      purrr::map_dfr(scen, ~ quantity(cells, .x, group_vars = "temp_bin")) %>%
        mutate(grouping = "bin", anio = NA_integer_)
    ) %>%
      mutate(region = REGION_LABELS[[reg]], outcome = out,
             pct_muni_days_above_target = pct_above)

    per_age[[paste(reg, out)]] <- res

    rm(est_sample, cells); gc()
  }

  rm(dt_reg); gc()
}

Z <- qnorm(0.975)

age_df <- bind_rows(per_age) %>%
  mutate(lwr = deaths - Z * se, upr = deaths + Z * se) %>%
  select(region, outcome, grouping, scenario, anio, temp_bin,
         deaths, se, lwr, upr, pct_muni_days_above_target)

write_out(age_df, "csv", "averted_deaths_by_age.csv")

# ==============================================================================
# COMBINE THE TWO AGE BANDS INTO THE REPORTED 50+ FIGURES
# ==============================================================================
combined <- age_df %>%
  mutate(age = ifelse(outcome == "rate_age_50_69", "a69", "a70")) %>%
  select(region, grouping, scenario, anio, temp_bin, age, deaths, se) %>%
  tidyr::pivot_wider(names_from = age, values_from = c(deaths, se)) %>%
  filter(!is.na(deaths_a69), !is.na(deaths_a70)) %>%
  mutate(
    outcome = "age_50_plus",
    deaths  = deaths_a69 + deaths_a70,
    se      = sqrt(se_a69^2 + se_a70^2),
    lwr     = deaths - Z * se,
    upr     = deaths + Z * se,

    # Sensitivity to the independence assumption; see the header.
    se_worst_case  = se_a69 + se_a70,
    lwr_worst_case = deaths - Z * se_worst_case,
    upr_worst_case = deaths + Z * se_worst_case,
    rho_critical   = ((deaths / Z)^2 - se_a69^2 - se_a70^2) /
                     (2 * se_a69 * se_a70),
    excludes_zero_for_any_rho = rho_critical > 1
  ) %>%
  select(region, outcome, grouping, scenario, anio, temp_bin,
         deaths, se, lwr, upr,
         deaths_50_69 = deaths_a69, deaths_70_plus = deaths_a70,
         se_50_69 = se_a69, se_70_plus = se_a70,
         se_worst_case, lwr_worst_case, upr_worst_case,
         rho_critical, excludes_zero_for_any_rho)

grand_df <- combined %>% filter(grouping == "grand") %>% select(-anio, -temp_bin)
year_df  <- combined %>% filter(grouping == "year")  %>% select(-temp_bin)
bin_df   <- combined %>% filter(grouping == "bin")   %>% select(-anio)

write_out(grand_df, "csv", "averted_deaths_grand_total.csv")
write_out(year_df,  "csv", "averted_deaths_by_year.csv")
write_out(bin_df,   "csv", "averted_deaths_by_bin.csv")

# Internal check: averted + additional should equal total, to rounding.
chk <- grand_df %>%
  select(region, scenario, deaths) %>%
  tidyr::pivot_wider(names_from = scenario, values_from = deaths) %>%
  mutate(discrepancy = averted + additional - total)
if (any(abs(chk$discrepancy) > 1e-6 * pmax(abs(chk$total), 1)))
  warning("Decomposition does not add up - check the weight construction.",
          call. = FALSE)
write_out(chk, "csv", "averted_deaths_decomposition_check.csv")

for (r in unique(year_df$region)) {
  df  <- year_df %>% filter(region == r, scenario != "total")
  tot <- year_df %>% filter(region == r, scenario == "total")
  p <- ggplot(df, aes(anio, deaths, fill = scenario)) +
    geom_col() +
    geom_errorbar(data = tot, aes(x = anio, ymin = lwr, ymax = upr),
                  inherit.aes = FALSE, width = 0.3) +
    theme_minimal() +
    labs(x = NULL, y = "Deaths", fill = NULL,
         caption = paste0("Bars: deaths averted by observed coverage and ",
                          "further deaths avoidable at ", FHS_HIGH,
                          "% coverage. Error bars are 95% uncertainty ",
                          "intervals on the total, which are perfectly ",
                          "correlated across years because they derive from ",
                          "the same coefficients."))
  save_fig(p, sprintf("appendixE_averted_deaths_by_year_%s.png",
                      gsub("[^A-Za-z]", "", r)), width = 9, height = 5)
}

# ------------------------------------------------------------------------------
# Console summary, including whether the independence assumption could matter.
# ------------------------------------------------------------------------------
for (i in seq_len(nrow(grand_df))) {
  r <- grand_df[i, ]
  say(sprintf("%s / %s: %.0f deaths (%.0f to %.0f); %s",
              r$region, r$scenario, r$deaths, r$lwr, r$upr,
              if (isTRUE(r$excludes_zero_for_any_rho))
                "excludes zero for ANY age-model correlation"
              else if (r$lwr > 0 || r$upr < 0)
                sprintf("significant only if the age-model correlation is below %.2f",
                        r$rho_critical)
              else "includes zero"))
}

say("07_deaths_averted.R done")
