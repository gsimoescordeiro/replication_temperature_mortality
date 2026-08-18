# ==============================================================================
# 02_descriptives.R
# ------------------------------------------------------------------------------
# Descriptive statistics reported in the Results ("Population-weighted mean
# daily mortality ... was X per 100 000 person-days") and in Appendix A,
# Table A1, plus the coverage series behind Figure 1.
#
# Unit of analysis is the MUNICIPALITY-DAY.  Means are population weighted, so
# they answer "what did the average person experience"; standard deviations are
# the unweighted spread across municipality-days, so they describe variation
# across places and days rather than variation of an already-averaged series.
# Both are labelled as such in the output.
#
# Table A1 uses a different unit again - the region-day series - which is why
# its standard deviations are an order of magnitude smaller than the
# municipality-day ones.  Section 5 builds it, and the header there explains the
# difference.
#
# Inputs : panel (via get_panel()), GDP file from 00_config.R DATA_FILES$gdp
# Outputs: csv/descriptives_by_region.csv
#          csv/descriptives_national.csv
#          csv/sample_size.csv
#          csv/exposure_days_per_year_by_region.csv
#          csv/fhs_coverage_by_region_year.csv
#          csv/tableA1_summary_stats_long.csv
#          csv/fig1_coverage_by_municipality.csv
#          csv/fig1_coverage_weighted_by_year.csv
#          tables/tableA1_summary_stats.csv
#          figures/fig1_fhs_coverage_trend.png
#          figures/fig1_map_fhs_coverage_{continuous,binned}_{2000,2010,2019}.png
#          figures/fig1_map_fhs_coverage_{continuous,binned}_panel.png
# ==============================================================================

dt <- get_panel() %>% add_outcomes() %>% add_temp_bins()

RATE_VARS <- c("rate_total", "rate_respiratory", "rate_circulatory",
               "rate_external", "rate_cancer", "rate_nutrition",
               "rate_infectious", "rate_age_20_49", "rate_age_50_69",
               "rate_age_70_plus", "rate_age_50_plus")

# ------------------------------------------------------------------------------
# 1. Rates and coverage by region
# ------------------------------------------------------------------------------
describe <- function(df, group_label) {
  rates <- purrr::map_dfr(RATE_VARS, function(v) {
    w <- df[[weight_var_for(v)]]
    x <- df[[v]]
    ok <- !is.na(x) & !is.na(w) & w > 0
    tibble(
      group          = group_label,
      variable       = v,
      mean_weighted  = weighted.mean(x[ok], w[ok]),
      sd_unweighted  = sd(x[ok]),
      p50            = median(x[ok]),
      n_muni_days    = sum(ok)
    )
  })

  extra <- tibble(
    group = group_label,
    variable = c("cob_esf", "DAT", "PREC"),
    mean_weighted = c(weighted.mean(df$cob_esf, df$pop_20_plus, na.rm = TRUE),
                      weighted.mean(df$DAT,     df$pop_20_plus, na.rm = TRUE),
                      weighted.mean(df$PREC,    df$pop_20_plus, na.rm = TRUE)),
    sd_unweighted = c(sd(df$cob_esf, na.rm = TRUE), sd(df$DAT, na.rm = TRUE),
                      sd(df$PREC, na.rm = TRUE)),
    p50 = c(median(df$cob_esf, na.rm = TRUE), median(df$DAT, na.rm = TRUE),
            median(df$PREC, na.rm = TRUE)),
    n_muni_days = nrow(df)
  )

  bind_rows(rates, extra)
}

by_region <- purrr::map_dfr(REGIONS, function(r) {
  say("Describing region ", REGION_LABELS[[r]])
  describe(subset_region(dt, r), REGION_LABELS[[r]])
})

national <- describe(dt, "Brazil")

# Sample size statement for the Methods
sample_stats <- tibble(
  n_municipalities = n_distinct(dt$code_muni),
  n_municipality_days = nrow(dt),
  first_date = min(dt$date), last_date = max(dt$date),
  total_deaths_20_plus = sum(dt$deaths_20_plus)
)

write_out(by_region, "csv", "descriptives_by_region.csv")
write_out(national,  "csv", "descriptives_national.csv")
write_out(sample_stats, "csv", "sample_size.csv")

say(sprintf("Panel: %s municipalities, %s municipality-days, %s deaths (20+)",
            format(sample_stats$n_municipalities, big.mark = " "),
            format(sample_stats$n_municipality_days, big.mark = " "),
            format(sample_stats$total_deaths_20_plus, big.mark = " ")))

# ------------------------------------------------------------------------------
# 2. Days per year in each temperature bin (the exposure histogram under Fig 2)
# ------------------------------------------------------------------------------
exposure <- dt %>%
  mutate(region = REGION_LABELS[substr(as.character(code_muni), 1, 1)]) %>%
  filter(!is.na(temp_bin)) %>%
  group_by(region, temp_bin) %>%
  summarise(pop_days = sum(pop_20_plus), .groups = "drop_last") %>%
  mutate(days_per_year = pop_days / sum(pop_days) * 365) %>%
  ungroup()

write_out(exposure, "csv", "exposure_days_per_year_by_region.csv")

# ------------------------------------------------------------------------------
# 3. FHS coverage by region and year (the trend behind Figure 1)
# ------------------------------------------------------------------------------
coverage <- dt %>%
  mutate(region = REGION_LABELS[substr(as.character(code_muni), 1, 1)]) %>%
  group_by(region, anio) %>%
  summarise(coverage_pop_weighted = weighted.mean(cob_esf, pop_20_plus),
            coverage_unweighted   = mean(cob_esf),
            n_municipalities      = n_distinct(code_muni),
            .groups = "drop")

national_coverage <- dt %>%
  group_by(anio) %>%
  summarise(region = "Brazil",
            coverage_pop_weighted = weighted.mean(cob_esf, pop_20_plus),
            coverage_unweighted   = mean(cob_esf),
            n_municipalities      = n_distinct(code_muni),
            .groups = "drop")

coverage <- bind_rows(coverage, national_coverage)
write_out(coverage, "csv", "fhs_coverage_by_region_year.csv")

p <- ggplot(coverage, aes(anio, coverage_pop_weighted, colour = region)) +
  geom_line(linewidth = 0.8) +
  theme_minimal() +
  labs(x = NULL, y = "FHS coverage (% of population, population weighted)",
       colour = NULL)
save_fig(p, "fig1_fhs_coverage_trend.png", width = 8, height = 5)

# ==============================================================================
# 4. Figure 1 - municipal coverage maps for 2000, 2010 and 2019
# ------------------------------------------------------------------------------
# Figure 1 in the paper is a map, not the trend line above; the trend is the
# series behind it and is quoted in the text.  Both the continuous and the
# binned version are drawn; the submitted figure is the binned panel.
#
# Coverage here is read straight from the FHS file rather than from the
# estimation panel, because the panel is municipality-DAY and this figure is
# municipality-YEAR.  It receives the same cap as everywhere else, so the map
# and the estimates describe the same variable.
# ==============================================================================

MAP_YEARS <- c(2000, 2010, 2019)

# Boundary vintage.  2020 is the last year of the study period plus one, so
# every municipality that existed during the panel appears; codes are truncated
# to six digits to match the panel.
MAP_BOUNDARY_YEAR <- 2020

# ------------------------------------------------------------------------------
# 4a. Coverage and population, straight from the source files
# ------------------------------------------------------------------------------
# The estimation panel is not used here: it is municipality-DAY and this figure
# is municipality-YEAR, so reading the two annual files directly is both faster
# and clearer about what is being plotted.
pop <- haven::read_dta(data_path("pop")) %>%
  mutate(anio = as.numeric(anio), code_muni = as.numeric(code_muni))
assert_unique_key(pop, c("code_muni", "anio"), "population file")

pop_col <- if ("pop_total" %in% names(pop)) "pop_total" else
  stop("No pop_total column in ", data_path("pop"))

pop <- pop %>%
  transmute(code_muni, anio, pop = as.numeric(.data[[pop_col]]))

esf <- haven::read_dta(data_path("esf")) %>%
  mutate(anio      = as.numeric(ano),
         code_muni = as.numeric(substr(id_municipio, 1, 6)),
         cob_esf   = as.numeric(cob_esf)) %>%
  select(code_muni, anio, cob_esf)
assert_unique_key(esf, c("code_muni", "anio"), "FHS coverage file")

if (CAP_COVERAGE_AT_100) {
  n_capped <- sum(esf$cob_esf > 100, na.rm = TRUE)
  say(sprintf("Map: capping coverage at 100 (%s of %s municipality-years)",
              format(n_capped, big.mark = " "),
              format(nrow(esf), big.mark = " ")))
  esf <- esf %>% mutate(cob_esf = pmin(cob_esf, 100))
}

cov_year <- esf %>%
  filter(anio %in% MAP_YEARS) %>%
  left_join(pop, by = c("code_muni", "anio"))

write_out(cov_year, "csv", "fig1_coverage_by_municipality.csv")

# The population-weighted national figures quoted in the text.
cov_weighted <- cov_year %>%
  group_by(anio) %>%
  summarise(coverage_pop_weighted = wmean(cob_esf, pop),
            coverage_unweighted   = mean(cob_esf, na.rm = TRUE),
            n_municipalities      = n_distinct(code_muni),
            .groups = "drop")

write_out(cov_weighted, "csv", "fig1_coverage_weighted_by_year.csv")

for (i in seq_len(nrow(cov_weighted)))
  say(sprintf("  %d: %.2f%% population weighted, %.2f%% unweighted",
              cov_weighted$anio[i], cov_weighted$coverage_pop_weighted[i],
              cov_weighted$coverage_unweighted[i]))

# ------------------------------------------------------------------------------
# 4b. Bins and theme
# ------------------------------------------------------------------------------
# Zero is its own category rather than the bottom of a continuous ramp: a
# municipality with no FHS team is a different object from one with a small
# programme, and the rollout story is largely the disappearance of that
# category.
BIN_LABELS <- c("0", ">0-25", "25-50", "50-75", "75-100")

coverage_bin <- function(x) {
  out <- rep(NA_character_, length(x))
  out[is.na(x) | x <= 0]     <- "0"
  out[!is.na(x) & x > 0  & x < 25] <- ">0-25"
  out[!is.na(x) & x >= 25 & x < 50] <- "25-50"
  out[!is.na(x) & x >= 50 & x < 75] <- "50-75"
  out[!is.na(x) & x >= 75]          <- "75-100"
  factor(out, levels = BIN_LABELS)
}

theme_map <- function() {
  ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title    = ggplot2::element_text(size = 10),
      legend.text     = ggplot2::element_text(size = 9),
      plot.title      = ggplot2::element_text(face = "bold", hjust = 0.5,
                                              size = 13)
    )
}

# ------------------------------------------------------------------------------
# 4c. Boundaries
# ------------------------------------------------------------------------------
# geobr and sf are not in REQUIRED_PKGS because no estimate depends on them, and
# sf in particular pulls in system libraries (GDAL, GEOS, PROJ) that a machine
# running only the models does not need.  So they are installed here, on demand,
# rather than being made a condition of running the pipeline at all.  If the
# install fails - no network, no compiler, no GDAL - the coverage data behind
# Figure 1 has already been written and the run continues without the drawing.
# Set INSTALL_MAP_PKGS <- FALSE before sourcing on a machine where installing
# from CRAN is not wanted; the maps are then skipped rather than attempted.
if (!exists("INSTALL_MAP_PKGS")) INSTALL_MAP_PKGS <- TRUE

MAP_PKGS <- c("geobr", "sf")

map_pkgs_ready <- function()
  all(vapply(MAP_PKGS, requireNamespace, logical(1), quietly = TRUE))

if (INSTALL_MAP_PKGS && !map_pkgs_ready()) {
  missing_map <- MAP_PKGS[!vapply(MAP_PKGS, requireNamespace, logical(1),
                                  quietly = TRUE)]
  say("Installing the mapping packages: ", paste(missing_map, collapse = ", "))
  repo <- getOption("repos")
  if (is.null(repo[["CRAN"]]) || repo[["CRAN"]] == "@CRAN@")
    repo <- c(CRAN = "https://cloud.r-project.org")
  try(utils::install.packages(missing_map, repos = repo), silent = TRUE)
}

if (!map_pkgs_ready()) {
  say("Skipping the maps: ",
      paste(MAP_PKGS[!vapply(MAP_PKGS, requireNamespace, logical(1),
                             quietly = TRUE)], collapse = ", "),
      " could not be installed. The coverage data behind Figure 1 has been ",
      "written; install the packages by hand and re-run to draw it.")
} else {

br_muni <- geobr::read_municipality(year = MAP_BOUNDARY_YEAR,
                                    showProgress = FALSE) %>%
  mutate(code_muni = as.numeric(substr(as.character(code_muni), 1, 6)))

# ------------------------------------------------------------------------------
# 4d. Draw
# ------------------------------------------------------------------------------
maps_cont <- list()
maps_bin  <- list()

for (yy in MAP_YEARS) {

  year_data <- cov_year %>%
    filter(anio == yy) %>%
    group_by(code_muni) %>%
    summarise(cob_esf = mean(cob_esf, na.rm = TRUE), .groups = "drop")

  map_sf <- br_muni %>%
    left_join(year_data, by = "code_muni") %>%
    mutate(cob_esf     = ifelse(is.na(cob_esf), 0, cob_esf),
           cob_esf_bin = coverage_bin(cob_esf))

  n_unmatched <- sum(!br_muni$code_muni %in% year_data$code_muni)
  if (n_unmatched)
    say(sprintf(paste0("  %d: %d of %d municipalities have no coverage record ",
                       "and are drawn as zero"),
                yy, n_unmatched, nrow(br_muni)))

  maps_cont[[as.character(yy)]] <- ggplot(map_sf) +
    geom_sf(aes(fill = cob_esf), colour = NA, linewidth = 0) +
    scale_fill_viridis_c(option = "magma", direction = -1,
                         limits = c(0, 100), oob = scales::squish,
                         name = "FHS coverage (%)") +
    labs(title = yy) +
    theme_map()

  maps_bin[[as.character(yy)]] <- ggplot(map_sf) +
    geom_sf(aes(fill = cob_esf_bin), colour = NA, linewidth = 0) +
    scale_fill_brewer(palette = "YlOrRd", drop = FALSE,
                      name = "FHS coverage (%)") +
    labs(title = yy) +
    theme_map()

  save_fig(maps_cont[[as.character(yy)]],
           sprintf("fig1_map_fhs_coverage_continuous_%d.png", yy),
           width = 7, height = 7)
  save_fig(maps_bin[[as.character(yy)]],
           sprintf("fig1_map_fhs_coverage_binned_%d.png", yy),
           width = 7, height = 7)
}

# ------------------------------------------------------------------------------
# 4e. Panels
# ------------------------------------------------------------------------------
# One shared legend, so the three years are read against the same scale.
collect <- function(lst) {
  p <- lst[[1]] | lst[[2]] | lst[[3]]
  p + patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
}

save_fig(collect(maps_cont), "fig1_map_fhs_coverage_continuous_panel.png",
         width = 16, height = 7)
save_fig(collect(maps_bin),  "fig1_map_fhs_coverage_binned_panel.png",
         width = 16, height = 7)

}  # end of the geobr/sf branch


# ==============================================================================
# 5. Appendix A, Table A1 - the submitted layout
# ------------------------------------------------------------------------------
# Table A1 is built on a different unit of analysis from section 1 above, which
# is why its standard deviations are so much smaller.  Section 1 describes
# MUNICIPALITY-DAYS: how much rates vary across places and days.
# Table A1 first collapses to a REGION-DAY series - one observation per region
# per day, population weighted - and then reports the mean and standard
# deviation of that series over time.  Averaging across ~1600 Southeast
# municipalities removes almost all the cross-sectional spread, so the standard
# deviation there measures day-to-day movement of a regional average and is an
# order of magnitude smaller.  Neither is wrong; they answer different
# questions, and both are written out so the table can be labelled honestly.
#
# This section reproduces the submitted table.  Section 1 is retained because
# the Methods and Results quote municipality-day figures.
# ==============================================================================

# GDP per capita uses total population.  The Appendix D split in
# 06_gdp_split.R uses the 20+ denominator instead, so the two are not directly
# comparable and are not meant to be.
gdp_pc_lookup <- NULL
gdp_file <- data_path("gdp")
if (file.exists(gdp_file)) {
  gdp_raw  <- haven::read_dta(gdp_file)
  gdp_col  <- intersect(c("gdp", "pib", "gdp_total"), names(gdp_raw))
  year_col <- intersect(c("anio", "ano", "year"), names(gdp_raw))
  if (length(gdp_col)) {
    gdp_pc_lookup <- gdp_raw %>%
      mutate(code_muni = as.numeric(code_muni),
             gdp = as.numeric(.data[[gdp_col[1]]]))
    if (length(year_col)) {
      gdp_pc_lookup <- gdp_pc_lookup %>%
        mutate(anio = as.numeric(.data[[year_col[1]]])) %>%
        select(code_muni, anio, gdp)
      assert_unique_key(gdp_pc_lookup, c("code_muni", "anio"), "GDP file")
    } else {
      gdp_pc_lookup <- gdp_pc_lookup %>% select(code_muni, gdp)
      assert_unique_key(gdp_pc_lookup, "code_muni", "GDP file")
    }
  }
} else {
  say("GDP file not found - Table A1 will omit the GDP per capita row.")
}

dt_a1 <- dt %>%
  mutate(region = REGION_LABELS[substr(as.character(code_muni), 1, 1)])

if (!is.null(gdp_pc_lookup)) {
  n_before <- nrow(dt_a1)
  dt_a1 <- if ("anio" %in% names(gdp_pc_lookup))
    left_join(dt_a1, gdp_pc_lookup, by = c("code_muni", "anio"))
  else
    left_join(dt_a1, gdp_pc_lookup, by = "code_muni")
  assert_no_fanout(n_before, nrow(dt_a1), "Table A1 GDP join")
  dt_a1 <- dt_a1 %>%
    mutate(gdp_pc = gdp / pmax(.data[["pop_total"]], 1))
}

# Heat-related deaths are in some versions of the panel and not others.
A1_RATES <- c("rate_total", "rate_respiratory", "rate_circulatory",
              "rate_external", "rate_nutrition", "rate_cancer",
              "rate_infectious",
              "rate_age_20_49", "rate_age_50_69", "rate_age_70_plus")
if ("deaths_heat_related" %in% names(dt_a1)) {
  dt_a1 <- dt_a1 %>%
    mutate(rate_heat_related = tidyr::replace_na(deaths_heat_related, 0) /
             pop_20_plus * 1e5)
  A1_RATES <- append(A1_RATES, "rate_heat_related", after = 3)
} else {
  say("deaths_heat_related not in the panel - Table A1 omits that row.")
}

A1_WEIGHTED <- intersect(c("DAT", "PREC", "cob_esf", "gdp_pc"), names(dt_a1))

# --- region-day series --------------------------------------------------------
# Rates are recomputed from summed deaths and summed population rather than
# averaged across municipalities, so each region-day rate is the rate an
# inhabitant of that region faced that day.
a1_region_day <- dt_a1 %>%
  group_by(region, date) %>%
  summarise(
    across(all_of(A1_WEIGHTED), ~ wmean(.x, pop_20_plus)),
    deaths_20_plus     = sum(deaths_20_plus,     na.rm = TRUE),
    deaths_respiratory = sum(deaths_respiratory, na.rm = TRUE),
    deaths_circulatory = sum(deaths_circulatory, na.rm = TRUE),
    deaths_external    = sum(deaths_external,    na.rm = TRUE),
    deaths_nutrition   = sum(deaths_nutrition,   na.rm = TRUE),
    deaths_cancer      = sum(deaths_cancer,      na.rm = TRUE),
    deaths_infectious  = sum(deaths_infectious,  na.rm = TRUE),
    deaths_20_49       = sum(deaths_20_49,       na.rm = TRUE),
    deaths_50_69       = sum(deaths_50_69,       na.rm = TRUE),
    deaths_70_plus     = sum(deaths_70_plus,     na.rm = TRUE),
    deaths_heat_related = if ("deaths_heat_related" %in% names(dt_a1))
      sum(deaths_heat_related, na.rm = TRUE) else NA_real_,
    pop_20_plus = sum(pop_20_plus, na.rm = TRUE),
    pop_20_49   = sum(pop_20_49,   na.rm = TRUE),
    pop_50_69   = sum(pop_50_69,   na.rm = TRUE),
    pop_70_plus = sum(pop_70_plus, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    rate_total         = deaths_20_plus     / pop_20_plus * 1e5,
    rate_respiratory   = deaths_respiratory / pop_20_plus * 1e5,
    rate_circulatory   = deaths_circulatory / pop_20_plus * 1e5,
    rate_external      = deaths_external    / pop_20_plus * 1e5,
    rate_nutrition     = deaths_nutrition   / pop_20_plus * 1e5,
    rate_cancer        = deaths_cancer      / pop_20_plus * 1e5,
    rate_infectious    = deaths_infectious  / pop_20_plus * 1e5,
    rate_heat_related  = deaths_heat_related / pop_20_plus * 1e5,
    rate_age_20_49     = deaths_20_49   / pmax(pop_20_49,   1) * 1e5,
    rate_age_50_69     = deaths_50_69   / pmax(pop_50_69,   1) * 1e5,
    rate_age_70_plus   = deaths_70_plus / pmax(pop_70_plus, 1) * 1e5
  )

# --- national series ----------------------------------------------------------
# Built from the region-day series so the national column is the population
# weighted combination of exactly the numbers in the other columns.
a1_brazil <- a1_region_day %>%
  group_by(date) %>%
  summarise(
    across(all_of(A1_WEIGHTED), ~ wmean(.x, pop_20_plus)),
    across(all_of(setdiff(A1_RATES,
                          c("rate_age_20_49", "rate_age_50_69",
                            "rate_age_70_plus"))),
           ~ wmean(.x, pop_20_plus)),
    rate_age_20_49   = wmean(rate_age_20_49,   pop_20_49),
    rate_age_50_69   = wmean(rate_age_50_69,   pop_50_69),
    rate_age_70_plus = wmean(rate_age_70_plus, pop_70_plus),
    .groups = "drop"
  ) %>%
  mutate(region = "Brazil")

a1_series <- bind_rows(a1_brazil, a1_region_day)

A1_LABELS <- c(
  rate_total        = "Mortality rate, total (20+)",
  rate_respiratory  = "Mortality rate, respiratory (20+)",
  rate_circulatory  = "Mortality rate, circulatory (20+)",
  rate_heat_related = "Mortality rate, heat-related (20+)",
  rate_external     = "Mortality rate, external (20+)",
  rate_nutrition    = "Mortality rate, nutrition (20+)",
  rate_cancer       = "Mortality rate, cancer (20+)",
  rate_infectious   = "Mortality rate, infectious (20+)",
  rate_age_20_49    = "Mortality rate, age 20-49",
  rate_age_50_69    = "Mortality rate, age 50-69",
  rate_age_70_plus  = "Mortality rate, age 70+",
  DAT               = "Temperature",
  PREC              = "Precipitation",
  cob_esf           = "FHS coverage",
  gdp_pc            = "GDP per capita"
)

a1_long <- purrr::map_dfr(c(A1_RATES, A1_WEIGHTED), function(v) {
  a1_series %>%
    group_by(region) %>%
    summarise(variable = v,
              mean = mean(.data[[v]], na.rm = TRUE),
              sd   = sd(.data[[v]],   na.rm = TRUE),
              n_days = sum(!is.na(.data[[v]])),
              .groups = "drop")
}) %>%
  mutate(label = ifelse(variable %in% names(A1_LABELS),
                        A1_LABELS[variable], variable))

write_out(a1_long, "csv", "tableA1_summary_stats_long.csv")

# Wide layout, Brazil first, as the table appears in the appendix.
a1_wide <- a1_long %>%
  select(label, variable, region, mean, sd) %>%
  tidyr::pivot_longer(c(mean, sd), names_to = "stat", values_to = "value") %>%
  mutate(region_stat = paste0(region, "_", stat)) %>%
  select(label, variable, region_stat, value) %>%
  tidyr::pivot_wider(names_from = region_stat, values_from = value)

a1_wide <- a1_wide %>%
  select(all_of(c("label", "variable",
                  intersect(c("Brazil_mean", "Brazil_sd"), names(a1_wide)),
                  sort(setdiff(names(a1_wide),
                               c("label", "variable",
                                 "Brazil_mean", "Brazil_sd")))))) %>%
  mutate(across(-c(label, variable), ~ round(., 3)))

# Row order follows the labels above rather than the alphabetical order the
# pivot produces.
a1_wide <- a1_wide[order(match(a1_wide$variable, names(A1_LABELS))), ]

write_out(a1_wide, "tables", "tableA1_summary_stats.csv")

say(sprintf("Table A1: %d rows over %d region-day series",
            nrow(a1_wide), n_distinct(a1_series$region)))

say("02_descriptives.R done")
