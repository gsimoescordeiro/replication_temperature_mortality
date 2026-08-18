get_V <- function(fit) {

  if (!is.null(fit$.V_cluster)) return(fit$.V_cluster)

  stored <- fit$cov.scaled
  if (is.null(stored)) return(vcov(fit, cluster = ~code_muni))

  state <- get0(".vcov_state", envir = globalenv(), ifnotfound = NULL)

  if (is.null(state)) {
    # First call this session: work out how the stored matrix relates to the
    # clustered one.  Three possibilities, and only the first two are usable.
    fresh <- vcov(fit, cluster = ~code_muni)
    S <- unname(as.matrix(stored)); F <- unname(as.matrix(fresh))

    if (isTRUE(all.equal(S, F, tolerance = 1e-8))) {
      state <- list(use = "stored", k = 1)
      say("Stored covariance IS the clustered matrix - reusing it (fast path).")

    } else {
      # A small-sample correction in fixest is a SCALAR, so if the two differ
      # only by an ssc factor the ratio is constant across every element and the
      # stored matrix can be rescaled exactly rather than recomputed.
      big <- abs(S) > .Machine$double.eps
      ratio <- F[big] / S[big]
      k <- stats::median(ratio)
      constant <- length(ratio) > 0 &&
        max(abs(ratio - k)) < 1e-6 * max(1, abs(k))

      if (constant) {
        state <- list(use = "scaled", k = k)
        say(sprintf(paste0("Stored covariance differs from the clustered one by ",
                           "a constant factor of %.6f (a small-sample ",
                           "correction) - rescaling it instead of recomputing ",
                           "(fast path)."), k))
      } else {
        state <- list(use = "recompute", k = NA_real_)
        say("Stored covariance is not a rescaling of the clustered one ",
            "(probably the unclustered matrix) - recomputing for every model. ",
            "This is correct but slow; models saved by the current 03_main_dlm.R ",
            "carry the clustered matrix and avoid it entirely.")
      }
    }
    assign(".vcov_state", state, envir = globalenv())
    return(fresh)
  }

  switch(state$use,
         stored    = stored,
         scaled    = stored * state$k,
         recompute = vcov(fit, cluster = ~code_muni))
}

# ==============================================================================
# 01_helpers.R
# ------------------------------------------------------------------------------
# Every function the pipeline uses: panel construction, lag construction, model
# fitting, coefficient extraction and the delta-method calculations.  Sourced by
# each numbered script after 00_config.R.
#
# Nothing in this file runs on its own.
#
# Why one shared file: the eight model scripts are near-identical in structure,
# and coefficient selection is the step most easily got wrong in a way that
# still produces plausible numbers.  Keeping the extraction code here means the
# definition of "cumulative effect" exists in exactly one place and cannot drift
# between scripts.
# ==============================================================================

# ==============================================================================
# SECTION 1 - PANEL CONSTRUCTION
# ==============================================================================

#' Read the three source files and assemble the municipality-day panel.
#'
#' Integrity checks are deliberate and loud: a silent type mismatch on the join
#' key would set FHS coverage to zero everywhere and produce a well-behaved but
#' meaningless set of interaction estimates.
load_panel <- function() {

  say("Reading mortality/temperature panel ...")
  dt <- haven::read_dta(data_path("panel")) %>%
    mutate(anio      = as.numeric(anio),
           date      = as.Date(date),
           code_muni = as.numeric(code_muni))   # force numeric: the join keys
                                                # below are numeric, and Stata
                                                # files often arrive labelled

  say("Reading population ...")
  pop <- haven::read_dta(data_path("pop")) %>%
    mutate(anio = as.numeric(anio), code_muni = as.numeric(code_muni))

  # A duplicated municipality-year in a lookup table silently multiplies the
  # panel rows in a left_join, so check before joining rather than after.
  assert_unique_key(pop, c("code_muni", "anio"), "population file")

  n_before <- nrow(dt)
  dt <- dt %>% left_join(pop, by = c("code_muni", "anio"))
  assert_no_fanout(n_before, nrow(dt), "population join")

  say("Reading FHS coverage ...")
  esf <- haven::read_dta(data_path("esf")) %>%
    mutate(anio      = as.numeric(ano),
           code_muni = as.numeric(substr(id_municipio, 1, 6))) %>%
    select(code_muni, anio, cob_esf)
  assert_unique_key(esf, c("code_muni", "anio"), "FHS coverage file")

  n_before <- nrow(dt)
  dt <- dt %>% left_join(esf, by = c("code_muni", "anio"))
  assert_no_fanout(n_before, nrow(dt), "FHS coverage join")

  # --- Missing coverage --------------------------------------------------------
  # Report before imputing, so the size of the assumption is on the record.
  n_missing <- sum(is.na(dt$cob_esf))
  say(sprintf("FHS coverage missing for %s of %s municipality-days (%.1f%%)",
              format(n_missing, big.mark = " "),
              format(nrow(dt), big.mark = " "),
              100 * n_missing / nrow(dt)))
  if (TREAT_MISSING_COVERAGE_AS_ZERO) {
    dt <- dt %>% mutate(cob_esf = tidyr::replace_na(cob_esf, 0))
  } else {
    dt <- dt %>% filter(!is.na(cob_esf))
  }

  # --- Cap ---------------------------------------------------------------------
  # Applied after the missing-value fill so both decisions are visible in one
  # place.  The count is reported every run: it is a substantive choice about the
  # paper's key moderator, not a data-cleaning detail.
  if (CAP_COVERAGE_AT_100) {
    n_capped <- sum(dt$cob_esf > 100, na.rm = TRUE)
    say(sprintf(paste0("Capping coverage at 100: %s of %s municipality-days ",
                       "(%.1f%%) exceeded it, maximum %.1f"),
                format(n_capped, big.mark = " "),
                format(nrow(dt), big.mark = " "),
                100 * n_capped / nrow(dt), max(dt$cob_esf)))
    dt <- dt %>% mutate(cob_esf = pmin(cob_esf, 100))
  } else {
    say("Coverage NOT capped (CAP_COVERAGE_AT_100 = FALSE) - values above 100 ",
        "are retained; see tests/diagnose_coverage.R.")
  }

  if (mean(dt$cob_esf > 0) < 0.01)
    fail("FHS coverage is zero for essentially every row - the coverage join ",
         "almost certainly failed (check the type and width of code_muni).")

  # Both means are printed because they answer different questions and differ a
  # lot here: the unweighted mean counts every municipality-day equally, while
  # the population-weighted mean is the one the paper reports.  Small
  # municipalities have much higher coverage, so the unweighted figure is far
  # higher.  Coverage above 100 is possible because e-Gestor's indicator is a
  # capacity calculation (teams x 3450 / population x 100); the Ministry caps
  # its published indicator at 100, so a maximum above 100 means this file holds
  # the uncapped version.  See tests/diagnose_coverage.R.
  say(sprintf(paste0("Coverage: mean %.1f (unweighted), %.1f (population ",
                     "weighted); max %.1f; %.1f%% of rows above %d, %.1f%% ",
                     "above 100"),
              mean(dt$cob_esf),
              wmean(dt$cob_esf, dt[["pop_total"]]),
              max(dt$cob_esf),
              100 * mean(dt$cob_esf > FHS_HIGH), FHS_HIGH,
              100 * mean(dt$cob_esf > 100)))

  dt
}

#' Population-weighted mean that tolerates missing values in EITHER argument.
#'
#' stats::weighted.mean(x, w, na.rm = TRUE) drops missing x but NOT missing w,
#' so a single NA weight returns NA for the whole mean.  Municipal population is
#' missing for some municipality-years here, which is why the load message
#' printed NA.
wmean <- function(x, w) {
  if (is.null(w)) return(NA_real_)
  ok <- !is.na(x) & !is.na(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

assert_unique_key <- function(df, keys, what) {
  n_dup <- df %>% count(across(all_of(keys))) %>% filter(n > 1) %>% nrow()
  if (n_dup > 0)
    fail(sprintf("%s has %d duplicated %s combinations; a left_join would ",
                 what, n_dup, paste(keys, collapse = "-")),
         "multiply panel rows and corrupt every estimate.")
  invisible(TRUE)
}

assert_no_fanout <- function(n_before, n_after, what) {
  if (n_after != n_before)
    fail(sprintf("%s changed the row count (%s -> %s): the lookup table is not ",
                 what, format(n_before, big.mark = " "),
                 format(n_after, big.mark = " ")),
         "unique on the join keys.")
  invisible(TRUE)
}

#' Restrict to one macro-region (first digit of the IBGE municipality code).
subset_region <- function(dt, reg) {
  dt %>% filter(substr(as.character(code_muni), 1, 1) == reg)
}

#' Build the age aggregates, the population denominators and the mortality rates.
#'
#' Deaths in the source panel are already restricted to ages 20+.  Rates are
#' deaths per 100 000 of the matching population: age-specific rates use the
#' matching age band, everything else uses the adult (20+) population.
add_outcomes <- function(dt) {
  dt %>%
    mutate(across(starts_with("deaths"), ~ tidyr::replace_na(., 0))) %>%
    mutate(
      deaths_50_69   = deaths_50_59 + deaths_60_69,
      pop_50_69      = pop_50_59    + pop_60_69,
      deaths_50_plus = deaths_50_69 + deaths_70_plus,
      pop_50_plus    = pop_50_69    + pop_70_plus,
      deaths_20_plus = deaths_20_49 + deaths_50_69 + deaths_70_plus,
      pop_20_plus    = pop_20_49    + pop_50_69    + pop_70_plus
    ) %>%
    # Municipality-days with no adult population carry no information and would
    # divide by zero.  Filtering here (before lag construction) is safe only
    # because add_lags() re-validates date continuity afterwards.
    filter(!is.na(pop_20_plus), pop_20_plus > 0) %>%
    mutate(
      rate_total        = deaths_20_plus  / pop_20_plus * 1e5,
      rate_respiratory  = deaths_respiratory / pop_20_plus * 1e5,
      rate_circulatory  = deaths_circulatory / pop_20_plus * 1e5,
      rate_external     = deaths_external    / pop_20_plus * 1e5,
      rate_cancer       = deaths_cancer      / pop_20_plus * 1e5,
      rate_nutrition    = deaths_nutrition   / pop_20_plus * 1e5,
      rate_infectious   = deaths_infectious  / pop_20_plus * 1e5,
      rate_age_20_49    = deaths_20_49   / pmax(pop_20_49,   1) * 1e5,
      rate_age_50_69    = deaths_50_69   / pmax(pop_50_69,   1) * 1e5,
      rate_age_70_plus  = deaths_70_plus / pmax(pop_70_plus, 1) * 1e5,
      rate_age_50_plus  = deaths_50_plus / pmax(pop_50_plus, 1) * 1e5,
      day_of_year       = as.integer(format(date, "%j"))
    )
}

#' Population weight matching an outcome (age-specific outcomes use their own
#' age band, so that a municipality-day contributes in proportion to the
#' population actually at risk).
weight_var_for <- function(outcome) {
  switch(outcome,
         rate_age_20_49   = "pop_20_49",
         rate_age_50_69   = "pop_50_69",
         rate_age_70_plus = "pop_70_plus",
         rate_age_50_plus = "pop_50_plus",
         "pop_20_plus")
}

# ==============================================================================
# SECTION 2 - TEMPERATURE BINS
# ==============================================================================

#' Absolute temperature bins.
#'
#' The final arm is an explicit DAT >= 30 test rather than TRUE: with a
#' catch-all, a missing temperature falls through every comparison (NA < 15 is
#' NA, not FALSE) and would be silently coded as the hottest bin - the paper's
#' headline exposure.
add_temp_bins <- function(dt) {
  dt %>%
    mutate(temp_bin = factor(case_when(
      is.na(DAT)  ~ NA_character_,
      DAT <  15   ~ "<15",
      DAT <  20   ~ "15-20",
      DAT <  25   ~ "20-25",
      DAT <  30   ~ "25-30",
      DAT >= 30   ~ ">30"
    ), levels = BIN_LEVELS))
}

#' Region-specific percentile bins (Appendix F(v)).
#'
#' Cut-points are computed on the region's own distribution of daily mean
#' temperature, so "extreme heat" means locally extreme rather than absolutely
#' hot.  Labels are chosen so that no label is a substring of another; the
#' coefficient matcher below is exact regardless, but readable labels help.
PCT_BIN_LEVELS <- c("p00_p10", "p10_p25", "p25_p75", "p75_p90", "p90_p100")
PCT_REF_BIN    <- "p25_p75"

add_percentile_bins <- function(dt) {
  qs <- quantile(dt$DAT, probs = c(.10, .25, .75, .90), na.rm = TRUE)
  say(sprintf("Percentile cut-points: p10=%.1f p25=%.1f p75=%.1f p90=%.1f",
              qs[1], qs[2], qs[3], qs[4]))
  dt %>%
    mutate(temp_bin = factor(case_when(
      is.na(DAT)   ~ NA_character_,
      DAT <  qs[1] ~ "p00_p10",
      DAT <  qs[2] ~ "p10_p25",
      DAT <  qs[3] ~ "p25_p75",
      DAT <  qs[4] ~ "p75_p90",
      DAT >= qs[4] ~ "p90_p100"
    ), levels = PCT_BIN_LEVELS))
}

# ==============================================================================
# SECTION 3 - DISTRIBUTED LAGS
# ==============================================================================

#' Create temp_L0 ... temp_Lk within municipality.
#'
#' dplyr::lag() shifts by ROW POSITION, not by date, so two things have to hold
#' and both are enforced here rather than assumed:
#'
#'   (a) rows are sorted by municipality and date - handled by the arrange();
#'   (b) lags never cross a break in the daily series.
#'
#' Some municipalities have gaps: a municipality created part-way through the
#' period, or a year with no population estimate, leaves a hole after
#' add_outcomes() drops those rows.  Lagging straight through a hole would make
#' temp_L1 mean "the day before the gap", silently misaligning every cumulative
#' effect for that municipality.
#'
#' Rather than stopping, the series is split into CONTIGUOUS SEGMENTS and lags
#' are built within each segment.  The first k_max days of a segment then have
#' missing lags and are dropped by the estimator - exactly what completing the
#' daily grid and dropping incomplete rows would do, without materialising the
#' missing rows.  Duplicated dates are still a hard error: they mean the panel
#' is malformed, and no lag structure can be correct.
add_lags <- function(dt, k_max = LAG_MAX) {

  dt <- dt %>% arrange(code_muni, date)

  dup <- dt %>% count(code_muni, date) %>% filter(n > 1)
  if (nrow(dup) > 0)
    fail(sprintf(paste0("%d municipality-date combinations appear more than ",
                        "once. The panel must have one row per municipality ",
                        "per day."), nrow(dup)))

  # Segment id increments whenever the gap to the previous day is not one day.
  dt <- dt %>%
    group_by(code_muni) %>%
    mutate(.seg = cumsum(c(TRUE, as.numeric(diff(date)) != 1))) %>%
    ungroup()

  n_gapped <- dt %>% group_by(code_muni) %>%
    summarise(segs = max(.seg), .groups = "drop") %>%
    filter(segs > 1) %>% nrow()

  if (n_gapped > 0) {
    lost <- dt %>% group_by(code_muni, .seg) %>%
      summarise(n = n(), .groups = "drop") %>%
      summarise(lost = sum(pmin(n, k_max))) %>% pull(lost)
    say(sprintf(paste0("%d municipalities have breaks in their daily series; ",
                       "lags are built within contiguous segments, so about %s ",
                       "municipality-days at segment starts have incomplete ",
                       "lag histories and will be dropped by the estimator."),
                n_gapped, format(lost, big.mark = " ")))
  }

  say(sprintf("Building lags 0..%d ...", k_max))
  dt <- dt %>% group_by(code_muni, .seg)
  for (l in 0:k_max) {
    dt <- dt %>% mutate(!!paste0("temp_L", l) := dplyr::lag(temp_bin, l))
  }
  dt %>% ungroup() %>% select(-.seg)
}

#' Number of observations actually used by a fixest fit.
#'
#' fixest registers an S3 method for stats::nobs but does not export a function
#' called nobs, so model_nobs() fails.  The stored value is used first and the
#' generic is the fallback.
model_nobs <- function(fit) {
  if (!is.null(fit$nobs)) return(as.integer(fit$nobs))
  tryCatch(as.integer(stats::nobs(fit)), error = function(e) NA_integer_)
}

lag_vars_for <- function(k_max = LAG_MAX) paste0("temp_L", 0:k_max)

# ==============================================================================
# SECTION 4 - ESTIMATION
# ==============================================================================

#' Fit the distributed-lag model for one outcome.
#'
#' Fixed effects:
#'   code_muni^day_of_year - municipality-by-calendar-day: absorbs each place's
#'                           own seasonal pattern (influenza season, school
#'                           calendar, fixed holidays)
#'   code_muni^anio        - municipality-by-year: absorbs slow-moving local
#'                           conditions, including the level of FHS coverage
#'                           itself, so identification comes from the
#'                           interaction rather than from the main effect
#'   date                  - day-by-year: absorbs national daily shocks common
#'                           to all municipalities
#' Standard errors are clustered by municipality; observations are weighted by
#' the population at risk.
fit_dlm <- function(dt, outcome, lag_vars, ref_bin, interact = TRUE,
                    interact_var = "cob_esf") {

  bin_terms <- paste0("i(", lag_vars, ", ref='", ref_bin, "')", collapse = " + ")
  rhs <- bin_terms

  if (interact) {
    inter_terms <- paste0("i(", lag_vars, ", ", interact_var,
                          ", ref='", ref_bin, "')", collapse = " + ")
    rhs <- paste0(rhs, " + ", inter_terms, " + ", interact_var)
  }

  f <- as.formula(paste0(
    outcome, " ~ ", rhs,
    " + PREC | code_muni^day_of_year + code_muni^anio + date"))

  w <- dt[[weight_var_for(outcome)]]

  fixest::feols(f, data = dt, weights = w, cluster = ~code_muni)
}

# ==============================================================================
# SECTION 5 - COEFFICIENT EXTRACTION
# ==============================================================================

#' Compute the cluster-robust covariance matrix ONCE and attach it to the fit.
#'
#' Every extractor below needs this matrix.  Computing it per call meant roughly
#' twenty full cluster-robust covariance computations per fitted model - on a
#' 300-parameter model over ten million rows, that dominated the runtime by an
#' order of magnitude over the fit itself.  Call this immediately after fitting;
#' the extractors then reuse the stored matrix.
cache_vcov <- function(fit) {
  fit$.V_cluster <- vcov(fit, cluster = ~code_muni)
  fit
}

#' Strip a fitted model down to what the extractors actually need.
#'
#' A fixest object carries residuals and fitted values - one value per
#' observation, so hundreds of megabytes on a ten-million-row regression.
#' Writing that to a network drive can take longer than estimating the model.
#' Every function in this file uses only the coefficients and the cluster-robust
#' covariance matrix, so that is all that is saved.  stats::coef() reads
#' $coefficients from a plain list, and get_V() reads $.V_cluster, so the slim
#' object is a drop-in replacement.
slim_fit <- function(fit) {
  list(coefficients = coef(fit),
       .V_cluster   = get_V(fit),
       nobs         = model_nobs(fit))
}

#' Cluster-robust covariance for a fit, using the cheapest valid source.
#'
#' Three sources, in order:
#'   1. the copy attached by cache_vcov() or slim_fit()  - free;
#'   2. the matrix fixest already stored at fit time in $cov.scaled - free, but
#'      only valid if it is the CLUSTERED one, so the first call in a session
#'      recomputes and compares before trusting it;
#'   3. a fresh vcov() call - correct but expensive, minutes on a large model.
#'
#' Source 2 matters when reading fits back from disk: models saved before
#' slim_fit() existed carry no cached matrix, and recomputing for each of the
#' 45 region-outcome fits turns a two-minute table into a two-hour one.
get_V <- function(fit) {

  if (!is.null(fit$.V_cluster)) return(fit$.V_cluster)

  stored <- fit$cov.scaled
  if (is.null(stored)) return(vcov(fit, cluster = ~code_muni))

  tested <- get0(".vcov_stored_ok", envir = globalenv(), ifnotfound = NULL)

  if (is.null(tested)) {
    # First call this session: verify the stored matrix IS the clustered one.
    fresh <- vcov(fit, cluster = ~code_muni)
    same  <- isTRUE(all.equal(unname(as.matrix(stored)), unname(as.matrix(fresh)),
                              tolerance = 1e-8))
    assign(".vcov_stored_ok", same, envir = globalenv())
    say("Stored covariance matrix ",
        if (same) "matches the clustered one and will be reused (fast path)."
        else "is NOT the clustered one; recomputing for every model (slow).")
    return(fresh)
  }

  if (isTRUE(tested)) return(stored)
  vcov(fit, cluster = ~code_muni)
}

#' Exact positions of the coefficients for one temperature bin.
#'
#' This replaces grep(label, names(coef(fit))).  Substring matching is unsafe
#' whenever one bin label contains another ("cold" also matches "ext_cold"), in
#' which case the "cumulative" effect silently sums two bins and its standard
#' error is the SE of the wrong linear combination.
#'
#' fixest names interaction coefficients "temp_L0::<15" and, with a continuous
#' interaction, "temp_L0::<15:cob_esf".  Both orderings of the interaction name
#' are checked so the function is robust across fixest versions; whichever
#' exists in the fit is used.
bin_coef_index <- function(fit, lag_vars, label, interact_var = NULL) {

  nms <- names(coef(fit))
  main <- paste0(lag_vars, "::", label)

  targets <- if (is.null(interact_var)) main else
    c(paste0(main, ":", interact_var), paste0(interact_var, ":", main))

  idx <- match(intersect(targets, nms), nms)
  idx <- sort(idx[!is.na(idx)])

  if (length(idx) == 0)
    fail(sprintf(paste0("No coefficients found for bin '%s'%s.\nFirst few ",
                        "coefficient names are: %s"),
                 label,
                 if (is.null(interact_var)) "" else paste0(" x ", interact_var),
                 paste(utils::head(nms, 6), collapse = ", ")))

  if (length(idx) < length(lag_vars))
    warning(sprintf(paste0("Bin '%s': %d of %d lag terms present - fixest ",
                           "dropped %d for collinearity, so the cumulative ",
                           "effect covers fewer than %d days."),
                    label, length(idx), length(lag_vars),
                    length(lag_vars) - length(idx), length(lag_vars)),
            call. = FALSE)
  idx
}

#' Delta method for a linear combination L'b of estimated coefficients.
#' Returns the estimate, its standard error and a normal-approximation 95% CI.
delta_lincom <- function(b, V, L, level = 0.95) {
  est <- sum(L * b)
  se  <- sqrt(as.numeric(t(L) %*% V %*% L))
  z   <- qnorm(1 - (1 - level) / 2)
  list(estimate = est, se = se, lwr = est - z * se, upr = est + z * se)
}

#' Cumulative effect of one temperature bin, summed over the lag window, from a
#' model WITHOUT the coverage interaction.  This is the excess mortality on a
#' day in that bin relative to the region's reference bin, accumulated over the
#' following k_max days.
cum_effect <- function(fit, lag_vars, label, ref_bin) {

  if (identical(label, ref_bin))
    return(tibble(temp_bin = label, estimate = 0, se = 0, lwr = 0, upr = 0,
                  n_terms = NA_integer_))

  idx <- bin_coef_index(fit, lag_vars, label)
  V   <- get_V(fit)[idx, idx, drop = FALSE]
  d   <- delta_lincom(coef(fit)[idx], V, rep(1, length(idx)))

  tibble(temp_bin = label, estimate = d$estimate, se = d$se,
         lwr = d$lwr, upr = d$upr, n_terms = length(idx))
}

#' Cumulative effect of one temperature bin evaluated at a given level of FHS
#' coverage, from the INTERACTED model.
#'
#' The model is linear in coverage, so for bin k
#'     effect(k, c) = sum_l beta_kl + c * sum_l delta_kl = B_k + c * D_k
#' and the loading vector is (1,...,1, c,...,c).  The standard error must be
#' computed from the joint covariance of B and D - it cannot be interpolated
#' between the SEs at two other coverage levels, because
#'     SE(c) = sqrt(V_BB + 2c V_BD + c^2 V_DD)
#' is not linear in c.
pred_at_coverage <- function(fit, lag_vars, label, ref_bin, coverage,
                             interact_var = "cob_esf") {

  if (identical(label, ref_bin))
    return(tibble(temp_bin = label, coverage = coverage,
                  estimate = 0, se = 0, lwr = 0, upr = 0))

  idx_b <- bin_coef_index(fit, lag_vars, label)
  idx_d <- bin_coef_index(fit, lag_vars, label, interact_var = interact_var)
  idx   <- c(idx_b, idx_d)

  V <- get_V(fit)[idx, idx, drop = FALSE]
  L <- c(rep(1, length(idx_b)), rep(coverage, length(idx_d)))
  d <- delta_lincom(coef(fit)[idx], V, L)

  tibble(temp_bin = label, coverage = coverage, estimate = d$estimate,
         se = d$se, lwr = d$lwr, upr = d$upr)
}

#' Difference in cumulative excess mortality between two coverage levels.
#'
#' The main-effect terms cancel, so the difference is (c_hi - c_lo) * D_k and
#' depends only on the interaction coefficients:
#'     Var = (c_hi - c_lo)^2 * L' V_DD L,  L = vector of ones
cov_difference <- function(fit, lag_vars, label, ref_bin,
                           c_hi = FHS_HIGH, c_lo = 0, interact_var = "cob_esf") {

  if (identical(label, ref_bin))
    return(tibble(temp_bin = label, estimate = 0, se = 0, lwr = 0, upr = 0,
                  n_terms = NA_integer_))

  idx <- bin_coef_index(fit, lag_vars, label, interact_var = interact_var)
  V   <- get_V(fit)[idx, idx, drop = FALSE]
  d   <- delta_lincom(coef(fit)[idx], V, rep(c_hi - c_lo, length(idx)))

  tibble(temp_bin = label, estimate = d$estimate, se = d$se,
         lwr = d$lwr, upr = d$upr, n_terms = length(idx))
}

#' Cumulative interaction coefficients D_k for every non-reference bin, together
#' with their joint cumulative covariance matrix.  Used by the deaths-averted
#' calculation, where the quantities of interest are linear combinations across
#' bins and the cross-bin covariances matter.
cum_interaction_matrix <- function(fit, lag_vars, bins, ref_bin,
                                   interact_var = "cob_esf") {

  bins_nr <- setdiff(bins, ref_bin)
  idx_list <- lapply(bins_nr, function(b)
    bin_coef_index(fit, lag_vars, b, interact_var = interact_var))
  names(idx_list) <- bins_nr

  b_full <- coef(fit)
  V_full <- get_V(fit)

  D <- vapply(idx_list, function(i) sum(b_full[i]), numeric(1))

  # Cumulative covariance between bins k1 and k2 is the sum of the full
  # covariance block: Cov(sum_l d_k1l, sum_m d_k2m) = 1' V[k1, k2] 1
  V_D <- outer(seq_along(bins_nr), seq_along(bins_nr),
               Vectorize(function(i, j) sum(V_full[idx_list[[i]], idx_list[[j]]])))
  dimnames(V_D) <- list(bins_nr, bins_nr)

  list(D = D, V = V_D, bins = bins_nr,
       n_terms = vapply(idx_list, length, integer(1)))
}

# ==============================================================================
# SECTION 6 - OUTPUT
# ==============================================================================

#' Write a CSV, retrying on transient failures.
#'
#' Outputs often live on a network or cloud-synced drive, where a momentary
#' hiccup raises "error writing to connection" and would otherwise discard hours
#' of completed computation.  Three attempts, then a clear error.
write_out <- function(df, subdir, filename, attempts = 3) {
  path <- file.path(dir_out(subdir), filename)
  for (i in seq_len(attempts)) {
    okw <- tryCatch({ readr::write_csv(df, path); TRUE },
                    error = function(e) { say("  write failed (attempt ", i, "): ",
                                              conditionMessage(e)); FALSE })
    if (okw) { say("wrote ", path); return(invisible(path)) }
    Sys.sleep(5 * i)
  }
  stop("Could not write ", path, " after ", attempts, " attempts.", call. = FALSE)
}

save_fig <- function(plot, filename, width = 8, height = 5, attempts = 3) {
  path <- file.path(dir_out("figures"), filename)
  for (i in seq_len(attempts)) {
    okf <- tryCatch({ ggplot2::ggsave(path, plot, width = width, height = height,
                                      dpi = 300); TRUE },
                    error = function(e) { say("  figure write failed (attempt ", i,
                                              "): ", conditionMessage(e)); FALSE })
    if (okf) { say("wrote ", path); return(invisible(path)) }
    Sys.sleep(5 * i)
  }
  stop("Could not write ", path, " after ", attempts, " attempts.", call. = FALSE)
}

FHS_COLOURS <- c("FHS = 0" = "red", "FHS = 50" = "black", "FHS = 90" = "blue")

#' Standard coefficient plot: cumulative effect by temperature bin, one series
#' per coverage level.
plot_by_coverage <- function(df, ref_bin, ylab = NULL) {
  ggplot(df, aes(temp_bin, estimate, colour = coverage_label)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(position = position_dodge(0.5), size = 2.5) +
    geom_errorbar(aes(ymin = lwr, ymax = upr),
                  position = position_dodge(0.5), width = 0.2) +
    scale_colour_manual(values = FHS_COLOURS) +
    theme_minimal() +
    labs(x = "Temperature bin (ºC)",
         y = ylab %||% paste0(LAG_MAX,
              "-day cumulative impact on mortality rate (deaths per 100 000)"),
         colour = "FHS coverage",
         caption = paste0("Reference bin: ", ref_bin,
                          " (region-specific minimum-mortality bin)"))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

message("Helpers loaded.")

# ==============================================================================
# SECTION 7 - PANEL CACHE
# ==============================================================================
#' Build the panel once per session and reuse it.
#'
#' Reading and joining the full 2000-2019 municipality-day panel is the slowest
#' step in the pipeline, and every numbered script needs it.  The cache is
#' invalidated automatically if any source file is newer than the cache.
get_panel <- function(refresh = FALSE) {

  if (exists(".panel_cache", envir = globalenv()) && !refresh)
    return(get(".panel_cache", envir = globalenv()))

  # The cache is keyed on the capping choice, so flipping CAP_COVERAGE_AT_100
  # cannot silently reuse a panel built under the other setting.
  cache_dir <- if (!is.null(CACHE_DIR)) {
    if (!dir.exists(CACHE_DIR)) dir.create(CACHE_DIR, recursive = TRUE)
    CACHE_DIR
  } else dir_out("models")

  cache_file <- file.path(cache_dir,
                          sprintf("panel_cache_%s.rds",
                                  if (CAP_COVERAGE_AT_100) "capped100" else "uncapped"))
  src_time <- max(file.mtime(c(data_path("panel"), data_path("pop"),
                               data_path("esf"))), na.rm = TRUE)

  if (file.exists(cache_file) && file.mtime(cache_file) > src_time && !refresh) {
    say("Loading cached panel ...")
    dt <- readRDS(cache_file)
  } else {
    dt <- load_panel()
    say("Caching panel for reuse ...")
    saveRDS(dt, cache_file, compress = FALSE)
  }

  assign(".panel_cache", dt, envir = globalenv())
  dt
}

#' Region-level analysis frame: subset, outcomes, bins, lags.  One call gives a
#' frame ready for fit_dlm().
prepare_region <- function(dt, reg, k_max = LAG_MAX, bins = c("absolute", "percentile")) {
  bins <- match.arg(bins)
  out <- dt %>% subset_region(reg) %>% add_outcomes()
  out <- if (bins == "absolute") add_temp_bins(out) else add_percentile_bins(out)
  add_lags(out, k_max)
}

# ==============================================================================
# SECTION 8 - RESUMING AN INTERRUPTED RUN
# ==============================================================================

#' TRUE if every file exists and is non-empty.
outputs_exist <- function(paths) {
  length(paths) > 0 && all(file.exists(paths)) && all(file.size(paths) > 0, na.rm = TRUE)
}

#' The three files 03_main_dlm.R produces for one region and outcome.
main_dlm_outputs <- function(reg, out) {
  c(file.path(dir_out("csv"),    sprintf("cumulative_nointeraction_%s_%s.csv", reg, out)),
    file.path(dir_out("csv"),    sprintf("predicted_by_coverage_%s_%s.csv",    reg, out)),
    file.path(dir_out("models"), sprintf("fit_%s_%s.rds",                      reg, out)))
}
