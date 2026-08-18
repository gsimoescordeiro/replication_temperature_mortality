# ==============================================================================
# 08_sens_placebo_future_change.R
# ------------------------------------------------------------------------------
# Appendix F(i), falsification test.
#
# ------------------------------------------------------------------------------
# WHY NOT A LEAD OF COVERAGE
# ------------------------------------------------------------------------------
# The obvious placebo for a rollout study is to substitute a lead of the
# treatment for the treatment itself.  It cannot discriminate here.  FHS
# coverage is a slow-moving stock, so corr(coverage_t, coverage_t+2) is very
# high and the lead is little more than a noisy copy of the treatment: it
# reproduces the main result whether or not the effect is causal, and a "pass"
# carries no information.
#
# A lead is also undefined for the last years of the panel.  Filling those rows
# with coverage = 0 - the highest-coverage years in the sample, about a tenth of
# it, recoded as "no FHS" - is measurement error that attenuates the placebo
# interaction toward zero.  That alone can manufacture a null placebo, so the
# test would appear to pass because the placebo variable was corrupted rather
# than because future coverage fails to predict today's heat sensitivity.
#
# ------------------------------------------------------------------------------
# THE TEST RUN HERE
# ------------------------------------------------------------------------------
# Put the placebo on the CHANGE rather than the level, and estimate it in the
# SAME model as actual coverage:
#
#     future growth  D_it = coverage_i,t+3 - coverage_i,t+1
#
# Everything in D is strictly future relative to t, and future growth is only
# weakly correlated with the current level, so the two interactions are
# separately identified.  The model is
#
#   rate ~ sum_l i(temp_L l)                          (temperature bins)
#        + sum_l i(temp_L l) x coverage_t             (the effect of interest)
#        + sum_l i(temp_L l) x D_it                   (the placebo)
#        + coverage_t + D_it + PREC | the usual three fixed effects
#
# Reading the result:
#   * contemporaneous interaction negative, future-change interaction null
#       -> what a causal effect of coverage implies;
#   * both negative
#       -> municipalities on a coverage-expansion path were already becoming
#          less heat-sensitive, and the main estimate is picking up that trend;
#   * contemporaneous null once the placebo is included
#       -> the main estimate was the trend.
#
# Because the two enter one model, they are estimated on one sample by
# construction - no separate comparison fit is needed.  The last three years
# have no future growth defined and are dropped from both.
#
# The correlation between coverage_t and D_it on the estimation sample is
# reported: it is the evidence that this test has power where the level-lead
# test does not.
#
# Inputs : panel (via get_panel())
# Outputs: csv/appendixF_placebo_future_change.csv
#          csv/appendixF_placebo_future_change_diagnostics.csv
#          figures/appendixF_placebo_future_change_<outcome>.png
# ==============================================================================

# Future growth is measured between t+FUT_FROM and t+FUT_TO.  Both endpoints are
# in the future, so the current level never enters the placebo.
if (!exists("FUT_FROM")) FUT_FROM <- 1
if (!exists("FUT_TO"))   FUT_TO   <- 3

# Figures F1 and F2 in the submitted appendix are the Southeast and the South,
# each broken out by the three age groups, so that is the default scope; all-
# cause is kept as well because the Results paragraph quotes it.  Models already
# on disk from a wider run are picked up below, so nothing computed is lost.
if (!exists("PLACEBO_OUTCOMES"))
  PLACEBO_OUTCOMES <- unique(c("rate_total", SENS_OUTCOMES))
if (!exists("PLACEBO_REGIONS")) PLACEBO_REGIONS <- SENS_REGIONS

# Coefficients are reported for a 90 point rise in current coverage (the paper's
# counterfactual) and for a 10 point future rise, which is nearer the scale
# actually observed.  Both are also reported on a common 90 point scale so the
# two terms can be compared directly.
SCALE_NOW <- 90
SCALE_FUT <- 10

panel <- get_panel()

# ------------------------------------------------------------------------------
# Future coverage growth, built by joining on the year rather than by row
# position.  dplyr::lead() would shift by position and silently misalign any
# municipality whose annual series has a gap.
# ------------------------------------------------------------------------------
cov_annual <- panel %>% distinct(code_muni, anio, cob_esf)
assert_unique_key(cov_annual, c("code_muni", "anio"), "annual coverage series")

shifted <- function(k, nm) {
  cov_annual %>%
    mutate(anio = anio - k) %>%
    select(code_muni, anio, !!nm := cob_esf)
}

cov_fut <- cov_annual %>%
  select(code_muni, anio) %>%
  left_join(shifted(FUT_FROM, "cob_from"), by = c("code_muni", "anio")) %>%
  left_join(shifted(FUT_TO,   "cob_to"),   by = c("code_muni", "anio")) %>%
  mutate(cob_fut_chg = cob_to - cob_from) %>%
  select(code_muni, anio, cob_fut_chg)

assert_unique_key(cov_fut, c("code_muni", "anio"), "future-growth lookup")
n_before <- nrow(panel)
panel_f  <- panel %>% left_join(cov_fut, by = c("code_muni", "anio"))
assert_no_fanout(n_before, nrow(panel_f), "future-growth join")

n_drop <- sum(is.na(panel_f$cob_fut_chg))
say(sprintf("Dropping %s of %s municipality-days with no future growth defined (%.1f%%)",
            format(n_drop, big.mark = " "), format(nrow(panel_f), big.mark = " "),
            100 * n_drop / nrow(panel_f)))
panel_f <- panel_f %>% filter(!is.na(cob_fut_chg))

results <- list(); diags <- list()

for (reg in PLACEBO_REGIONS) {

  ref_bin  <- REF_BIN[[reg]]
  lag_vars <- lag_vars_for(LAG_MAX)
  say("=== Future-change placebo [", sprintf("t+%d minus t+%d", FUT_TO, FUT_FROM),
    "]: ", REGION_LABELS[[reg]], " ===")

  # Lags are rebuilt on the restricted panel; add_lags() re-segments each
  # municipality, so dropping the final years cannot misalign them.
  dt_reg <- prepare_region(panel_f, reg, k_max = LAG_MAX)

  for (out in PLACEBO_OUTCOMES) {

    # The horizon is part of the cache key.  Without it, changing FUT_TO and
    # re-running would silently reuse the previous horizon's models and label
    # them with the new one.
    f_out <- file.path(dir_out("models"),
                       sprintf("placebo_futchg_h%d%d_%s_%s.rds",
                               FUT_FROM, FUT_TO, reg, out))
    # Models saved before the horizon was added to the name are the 1-to-3
    # horizon; adopt them rather than refitting.
    f_legacy <- file.path(dir_out("models"),
                          sprintf("placebo_futchg_%s_%s.rds", reg, out))
    if (identical(c(FUT_FROM, FUT_TO), c(1, 3)) &&
        !file.exists(f_out) && file.exists(f_legacy)) f_out <- f_legacy

    if (SKIP_EXISTING && file.exists(f_out) && file.size(f_out) > 0) {
      say("  ", out, " - already done, reusing")
      st <- readRDS(f_out)
      results[[paste(reg, out)]] <- st$res
      diags[[paste(reg, out)]]   <- st$diag
      next
    }

    say("  outcome: ", out)
    t0 <- Sys.time()

    pop_var <- weight_var_for(out)
    est <- dt_reg %>%
      filter(!is.na(.data[[paste0("temp_L", LAG_MAX)]]),
             !is.na(PREC), !is.na(temp_bin), .data[[pop_var]] > 0)

    # Does the placebo carry independent variation?  If this correlation were
    # near 1 the test would be as uninformative as the level-lead version.
    rho <- suppressWarnings(stats::cor(est$cob_esf, est$cob_fut_chg,
                                       use = "complete.obs"))
    say(sprintf("    corr(coverage, future growth) = %+.3f ; sd(future growth) = %.1f pp",
                rho, stats::sd(est$cob_fut_chg, na.rm = TRUE)))

    bin_terms <- paste0("i(", lag_vars, ", ref='", ref_bin, "')", collapse = " + ")
    inter_now <- paste0("i(", lag_vars, ", cob_esf, ref='", ref_bin, "')",
                        collapse = " + ")
    inter_fut <- paste0("i(", lag_vars, ", cob_fut_chg, ref='", ref_bin, "')",
                        collapse = " + ")

    f <- as.formula(paste0(
      out, " ~ ", bin_terms, " + ", inter_now, " + ", inter_fut,
      " + cob_esf + cob_fut_chg + PREC",
      " | code_muni^day_of_year + code_muni^anio + date"))

    fit <- cache_vcov(fixest::feols(f, data = est, weights = est[[pop_var]],
                                    cluster = ~code_muni))

    # Cumulative interaction coefficients for each term.  bin_coef_index()
    # matches coefficient names exactly, so "…::<15:cob_esf" and
    # "…::<15:cob_fut_chg" cannot be confused with one another.
    grab <- function(v, scale, label) {
      ci <- cum_interaction_matrix(fit, lag_vars, BIN_LEVELS, ref_bin,
                                   interact_var = v)
      se <- sqrt(diag(ci$V))
      tibble(temp_bin = ci$bins, term = label,
             D_cum = ci$D, se_D = se, n_terms = ci$n_terms,
             scale = scale,
             estimate = scale * ci$D,
             lwr = scale * (ci$D - qnorm(0.975) * se),
             upr = scale * (ci$D + qnorm(0.975) * se),
             # common 90 point scale, for a like-for-like comparison
             est_90 = 90 * ci$D,
             lwr_90 = 90 * (ci$D - qnorm(0.975) * se),
             upr_90 = 90 * (ci$D + qnorm(0.975) * se))
    }

    res <- bind_rows(grab("cob_esf", SCALE_NOW, "Current coverage"),
                     grab("cob_fut_chg", SCALE_FUT, "Future growth (placebo)")) %>%
      mutate(region = REGION_LABELS[[reg]], outcome = out, ref_bin = ref_bin,
             fut_from = FUT_FROM, fut_to = FUT_TO,
             horizon = sprintf("t+%d minus t+%d", FUT_TO, FUT_FROM),
             n_obs = model_nobs(fit))

    dg <- tibble(region = REGION_LABELS[[reg]], outcome = out,
                 fut_from = FUT_FROM, fut_to = FUT_TO,
                 horizon = sprintf("t+%d minus t+%d", FUT_TO, FUT_FROM),
                 corr_coverage_futgrowth = rho,
                 sd_future_growth = stats::sd(est$cob_fut_chg, na.rm = TRUE),
                 mean_future_growth = mean(est$cob_fut_chg, na.rm = TRUE),
                 n_obs = model_nobs(fit),
                 minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))

    saveRDS(list(res = res, diag = dg), f_out)
    results[[paste(reg, out)]] <- res
    diags[[paste(reg, out)]]   <- dg

    say(sprintf("    done in %.1f minutes", dg$minutes))
    rm(fit, est); gc()
  }
  rm(dt_reg); gc()
}

# Any placebo model left on disk by an earlier, wider run is folded in, so
# narrowing PLACEBO_REGIONS or PLACEBO_OUTCOMES never discards work already done.
# Every placebo model on disk is folded in, across ALL horizons, so the CSV
# accumulates the horizon panel as further horizons are run.
for (f in list.files(dir_out("models"), pattern = "^placebo_futchg_.*\\.rds$",
                     full.names = TRUE)) {
  k <- paste("ondisk", basename(f))
  if (is.null(results[[k]])) {
    st <- readRDS(f)
    if (is.null(st$res$fut_to)) {          # pre-horizon files are the 1-to-3 run
      st$res$fut_from  <- 1; st$res$fut_to  <- 3
      st$res$horizon   <- "t+3 minus t+1"
      st$diag$fut_from <- 1; st$diag$fut_to <- 3
      st$diag$horizon  <- "t+3 minus t+1"
    }
    results[[k]] <- st$res; diags[[k]] <- st$diag
  }
}

res <- bind_rows(results) %>%
  distinct(region, outcome, temp_bin, term, horizon, .keep_all = TRUE) %>%
  mutate(temp_bin  = factor(as.character(temp_bin), levels = BIN_LEVELS),
         region    = factor(region, levels = unname(REGION_LABELS)),
         age_group = factor(age_label(outcome),
                            levels = age_label(unique(outcome))))
dgs <- bind_rows(diags) %>% distinct(region, outcome, horizon, .keep_all = TRUE)

write_out(res, "csv", "appendixF_placebo_future_change.csv")
write_out(dgs, "csv", "appendixF_placebo_future_change_diagnostics.csv")

# ------------------------------------------------------------------------------
# The headline read: the >30 C bin, both terms on the common 90 point scale.
# ------------------------------------------------------------------------------
say("--- >30 C bin, both terms on a common 90 point scale ---")
hot <- res %>% filter(temp_bin == ">30") %>%
  mutate(cell = sprintf("%+.3f (%+.3f to %+.3f)", est_90, lwr_90, upr_90)) %>%
  select(region, outcome, horizon, term, cell) %>%
  tidyr::pivot_wider(names_from = term, values_from = cell) %>%
  arrange(region, outcome, horizon)
print(as.data.frame(hot))

# Figures F1 and F2: one per region, panels = age group, matching the appendix.
# The figure shows the horizon just run; the CSV keeps all of them.
res_fig <- res %>% filter(fut_from == FUT_FROM, fut_to == FUT_TO)
for (r in levels(droplevels(res_fig$region))) {
  df <- res_fig %>% filter(region == r)
  p <- ggplot(df, aes(temp_bin, est_90, colour = term)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(position = position_dodge(0.5), size = 2) +
    geom_errorbar(aes(ymin = lwr_90, ymax = upr_90),
                  position = position_dodge(0.5), width = 0.2) +
    scale_colour_manual(values = c("Current coverage" = "blue",
                                   "Future growth (placebo)" = "grey45")) +
    facet_wrap(~ age_group, scales = "free_y", ncol = 1) +
    theme_minimal() +
    labs(x = "Temperature bin (ºC)", colour = NULL,
         y = "Change in 30-day cumulative impact per 90 points (deaths per 100 000)",
         caption = paste0("Both terms come from the same model and the same ",
                          "sample. Future growth is coverage at t+", FUT_TO,
                          " minus coverage at t+", FUT_FROM,
                          ", so it contains no information about coverage at t. ",
                          "A causal effect of coverage implies the blue series ",
                          "is negative on hot days and the grey series is not."))
  save_fig(p, sprintf("appendixF_placebo_future_change_%s.png",
                      gsub("[^A-Za-z]", "", r)), width = 9, height = 10)
}

say("08_sens_placebo_future_change.R done")
