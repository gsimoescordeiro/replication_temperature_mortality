# ==============================================================================
# 05_table1.R
# ------------------------------------------------------------------------------
# Table 1: the difference in cumulative excess mortality on very hot days
# (>30 ºC, relative to the region's minimum-mortality bin) between FHS coverage
# of 90% and 0%, by cause of death and by age group, for each macro-region.
#
# The model is linear in coverage, so for bin k
#     effect(k, c) = B_k + c * D_k
# and the 90-vs-0 difference is exactly 90 * D_k: the main-effect terms cancel,
# and the difference depends only on the cumulative interaction coefficients.
# Its standard error follows from the same linear combination,
#     SE = 90 * sqrt(1' V_DD 1)
# where V_DD is the cluster-robust covariance block of those coefficients.
#
# A negative value means higher coverage is associated with LESS excess
# mortality on very hot days.  Units are deaths per 100 000 population.
#
# Reuses the models fitted by 03_main_dlm.R rather than refitting them, so the
# table is guaranteed to describe the same estimates as Figure 2.
#
# Inputs : models/fit_<region>_<outcome>.rds  (from 03_main_dlm.R)
# Outputs: tables/table1_fhs90_vs_fhs0_very_hot_days.csv
#          tables/table1_formatted.csv
# ==============================================================================

HOT_BIN <- ">30"

TABLE1_OUTCOMES <- c(
  "rate_circulatory", "rate_respiratory", "rate_external",
  "rate_cancer", "rate_nutrition", "rate_infectious",
  "rate_age_20_49", "rate_age_50_69", "rate_age_70_plus"
)

# Partial results are written after every model.  Reading 45 fitted models back
# from disk is slow when they were saved as full fixest objects, so an
# interruption part-way through should not throw away the work already done.
partial_path <- file.path(dir_out("models"), "table1_partial.rds")
rows <- if (SKIP_EXISTING && file.exists(partial_path)) readRDS(partial_path) else list()
if (length(rows)) say("Resuming with ", length(rows), " cells already computed")

n_total <- length(REGIONS) * length(TABLE1_OUTCOMES)
i_cell  <- 0

for (reg in REGIONS) {
  for (out in TABLE1_OUTCOMES) {

    i_cell <- i_cell + 1
    key <- paste(reg, out)
    if (!is.null(rows[[key]])) next

    f <- file.path(dir_out("models"), sprintf("fit_%s_%s.rds", reg, out))
    if (!file.exists(f)) {
      warning("Missing model: ", f, " - run 03_main_dlm.R first.", call. = FALSE)
      next
    }

    t0 <- Sys.time()
    m <- readRDS(f)

    d <- cov_difference(m$fit, m$lag_vars, HOT_BIN, m$ref_bin,
                        c_hi = FHS_HIGH, c_lo = 0)

    rows[[key]] <- d %>%
      mutate(region = REGION_LABELS[[reg]], outcome = out,
             ref_bin = m$ref_bin, n_obs = m$n_obs)

    say(sprintf("  [%d/%d] %s / %s : %+.4f  (%.1f s)", i_cell, n_total,
                REGION_LABELS[[reg]], out, d$estimate,
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))

    saveRDS(rows, partial_path)
    rm(m); gc()
  }
}

if (!length(rows))
  stop("No fitted models found in ", dir_out("models"),
       ". Run 03_main_dlm.R first - it saves the fits this table reads.",
       call. = FALSE)

tab <- bind_rows(rows) %>%
  select(region, outcome, temp_bin, estimate, se, lwr, upr, n_terms, n_obs, ref_bin)

write_out(tab, "tables", "table1_fhs90_vs_fhs0_very_hot_days.csv")

# ------------------------------------------------------------------------------
# Publication layout: rows = outcome, columns = region, cells = estimate (CI).
# Confidence intervals are written with "to" rather than a hyphen, which is the
# BMJ convention and avoids ambiguity with negative estimates.
# ------------------------------------------------------------------------------
fmt <- function(e, l, u) sprintf("%.3f (%.3f to %.3f)", e, l, u)

formatted <- tab %>%
  mutate(cell = fmt(estimate, lwr, upr)) %>%
  select(outcome, region, cell) %>%
  tidyr::pivot_wider(names_from = region, values_from = cell) %>%
  mutate(outcome = recode(outcome,
    rate_circulatory = "Circulatory", rate_respiratory = "Respiratory",
    rate_external = "External", rate_cancer = "Cancer",
    rate_nutrition = "Nutritional", rate_infectious = "Infectious",
    rate_age_20_49 = "20-49 years", rate_age_50_69 = "50-69 years",
    rate_age_70_plus = "70+ years"))

write_out(formatted, "tables", "table1_formatted.csv")

unlink(partial_path)   # the table is complete; the partial cache is no longer needed

say("05_table1.R done")
