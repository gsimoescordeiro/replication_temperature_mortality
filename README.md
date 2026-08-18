# Reproduction package

Analysis code for **"Effects of primary health care coverage on temperature-related mortality in
Brazil: national quasi-experimental study"** (Cordeiro, Lagarde, Pereda, Millett, Gonçalves).

Running `run_all.R` reproduces every number, table and figure in the paper and its appendices,
in the order they appear in the text.

---

## Quick start

1. Open `00_config.R` and set three paths:

   ```r
   BASE_DIR <- "…/scripts_reproduction"          # this folder
   DATA_DIR <- "…"                                # where the .dta inputs live
   OUT_DIR  <- "…/reproduction_output"            # everything is written here
   ```

   Nothing else in the pipeline contains a path.

2. From that folder:

   ```r
   source("run_all.R")        # or: Rscript run_all.R
   ```

Outputs land in `OUT_DIR` under `csv/`, `tables/`, `figures/`, `models/` and `logs/`. Nothing is
ever written to `DATA_DIR`.

## Software

R 4.3.2, fixest 0.11.2 (`remotes::install_version("fixest", version = "0.11.2")`), plus haven,
dplyr, tidyr, stringr, readr, purrr, ggplot2, patchwork and scales. `13_session_info.R` writes the
exact environment to `logs/session_info.txt`. `00_config.R` stops with a clear message if a package
is missing and warns if the fixest version differs.

## Data

Not redistributed. All inputs are public: mortality from SIM (DATASUS), population from IBGE, FHS
coverage from e-Gestor Atenção Básica, temperature and precipitation from ERA5. The pipeline reads
three prepared Stata files named in `DATA_FILES`:

| Key | Contents |
|---|---|
| `panel` | Municipality-day panel: deaths (total, by cause, by age band, ages 20+), daily mean temperature `DAT`, precipitation `PREC` |
| `pop` | Annual municipal population by age band |
| `esf` | Annual municipal FHS coverage `cob_esf`, percentage of population registered with a team |
| `gdp` | Municipal GDP, used only by `06_gdp_split.R` |

## Files, in run order

| Script | Produces | Where it appears in the paper |
|---|---|---|
| `00_config.R` | Paths, packages, design constants | — |
| `01_helpers.R` | All functions: panel build, lags, estimation, coefficient extraction, delta method | Methods |
| `02_descriptives.R` | Mortality rates and coverage by region; exposure days per bin; municipal coverage maps; Table A1 | Results ¶1; **Figure 1**; **Appendix A, Table A1** |
| `03_main_dlm.R` | Main models with and without the coverage interaction; predictions at 0/50/90% coverage | **Figure 2**; Appendix B and C |
| `04_non_interacted_figure.R` | Temperature-mortality relationship pooled across coverage levels | Appendix B |
| `05_table1.R` | 90%-vs-0% difference on days above 30 °C, by cause and age | **Table 1**; Appendix C |
| `06_gdp_split.R` | Split at the regional median GDP per capita, by age group | Appendix D, Figures D1-D2 |
| `07_deaths_averted.R` | Deaths averted, further avoidable and total avertable, with delta-method intervals | Results, "Deaths averted"; Appendix E |
| `08_sens_placebo_future_change.R` | Falsification test: future coverage growth against current coverage, in one model | **Appendix F(i)**, Figures F1-F2 |
| `09_sens_lag_windows.R` | 7-day and 14-day lag windows | Appendix F(ii), Figures F3-F4 |
| `10_sens_coverage_bands.R` | Split-sample by coverage band | Appendix F(iii), Figures F5-F6 |
| `11_sens_exclude_large_muni.R` | Excluding municipalities over 500 000 | Appendix F(iv), Figures F7-F8 |
| `12_sens_percentile_bins.R` | Region-specific percentile bins | Appendix F(v), Figures F9-F10 |
| `13_session_info.R` | Software record | Methods |
| `run_all.R` | Runs 02-13, times each step, logs failures | — |

**Figure 1 needs `geobr` and `sf`.** They are not in `REQUIRED_PKGS` because no estimate
depends on them, and `sf` pulls in GDAL, GEOS and PROJ, which a machine running only the
models does not need. Section 4 of `02_descriptives.R` installs them on demand; if that
fails it writes the coverage data behind Figure 1 and skips the drawing rather than
failing the run. Set `INSTALL_MAP_PKGS <- FALSE` to stop it trying.

**Two units of analysis in Appendix A.** Section 1 of `02_descriptives.R` describes
municipality-days: the standard deviation there is variation across places and days. Table A1
collapses to a region-day series first and reports the standard deviation of that series over
time, which is an order of magnitude smaller, because averaging across a region removes the
cross-sectional spread. Both are written out and labelled: Table A1 is the appendix table, and
the municipality-day figures are the ones quoted in the Methods and Results.

## Scope of Appendices D, E and F

The appendix is in two halves. Appendices A, B and C cover all five macro-regions; Appendices D, E
and F cover the **Southeast and South only**, and every figure in Appendix F is broken out by the
three **age groups** (20-49, 50-69, 70+) rather than shown for all-cause mortality. Scripts `06`
and `08`-`12` follow that structure, driven by two constants in `00_config.R`:

```r
SENS_REGIONS  <- c("3", "4")                       # Southeast, South
SENS_OUTCOMES <- c("rate_age_20_49", "rate_age_50_69", "rate_age_70_plus")
```

Each script writes one figure per region with the age groups as panels, and one CSV carrying an
`outcome` column. Widen either constant to cover more ground; each addition is one more model per
script. Every fit is cached as an `.rds` under `models/`, so narrowing the scope later never
discards work and re-running only computes what is missing.

## The falsification test

The natural placebo for a rollout study is a lead of the treatment, but a lead cannot discriminate
here. FHS coverage is a slow-moving stock, so its own two-year lead is a near-perfect proxy for
current coverage and would reproduce the main result whether or not the effect is causal. A lead is
also undefined for the final years of the panel, and filling those rows with zero coverage — the
highest-coverage years in the sample, recoded as "no FHS" — attenuates the placebo interaction
toward zero, so a null placebo would be an artefact of the missing-value rule rather than evidence.

`08_sens_placebo_future_change.R` therefore puts the placebo on the **change** rather than the
level, and estimates it in the **same model** as actual coverage:

```
D_it = coverage(t+3) − coverage(t+1)

rate ~ Σ i(temp_L) + Σ i(temp_L)×coverage_t + Σ i(temp_L)×D_it + … | FE
```

Both endpoints of `D` are in the future, so current coverage never enters the placebo, and future
growth is only weakly correlated with the current level (−0.31 to −0.52 across regions), so the two
interactions are separately identified and the test has power. Because they enter one model they
are estimated on one sample by construction; no separate comparison fit is needed. The final three
years have no future growth defined and are dropped from both terms.

Reading the result:

| Contemporaneous term | Placebo term | Interpretation |
|---|---|---|
| negative on hot days | null | what a causal effect of coverage implies |
| negative | negative | municipalities on an expansion path were already becoming less heat-sensitive |
| null | either | the contemporaneous estimate was the trend |

`appendixF_placebo_future_change_diagnostics.csv` reports the correlation between coverage and
future growth for each region. Report it alongside the result: it is the evidence that the test can
discriminate at all.

The horizon is a parameter. `FUT_FROM` and `FUT_TO` are part of the model cache key, so running a
second horizon adds to the panel of results rather than overwriting the first.

## Resuming an interrupted run

A full run is many hours, so it is built to be restarted.

- `SKIP_EXISTING <- TRUE` (in `00_config.R`) makes `run_all.R` skip any step whose output is
  already on disk, and makes `03_main_dlm.R` skip individual region-outcome combinations it has
  already written — so an interruption costs only the piece it stopped on. Set it to `FALSE` to
  force a clean recomputation.
- Writes retry three times before failing, so a momentary network or cloud-sync error cannot
  discard hours of completed work.
- `14_check_outputs.R` compares what is on disk against the full expected manifest and reports
  exactly what is missing. `run_all.R` runs it automatically at the end; run it on its own at any
  time with `source("14_check_outputs.R")`. It writes `logs/output_check.csv`.

**Performance note.** The cluster-robust covariance matrix is computed once per fitted model
(`cache_vcov()`) and reused by every extractor. Each region-outcome calls the extractors about
twenty times, and the covariance computation costs far more than the fit itself, so recomputing it
per call would turn a couple of minutes into close to an hour.

## Definitions that are easy to misread

- **Coverage bands are not terciles.** The split-sample analysis uses fixed cut-points at 25% and
  75% (`FHS_BANDS`), which do not produce equal groups: in the Southeast the three bands hold 21%,
  24% and 55% of municipality-days. `10_sens_coverage_bands.R` writes the composition alongside the
  estimates.
- **The high-coverage counterfactual is coverage set to exactly 90% everywhere**, not "at least
  90%". Municipality-years already above 90% are moved *down* to it, so they contribute negatively
  to "additional avoidable deaths". `07_deaths_averted.R` reports what share of municipality-days
  that affects.
- **Percentile bins are a different exposure from the absolute bins.** The region-specific 90th
  percentile is cooler than 30 °C (26.7 °C in the Southeast, 25.6 °C in the South), so a weaker
  gradient in Appendix F(v) is expected rather than a failure to replicate.
- **Totals across all non-reference bins are dominated by the 25-30 °C bin**, whose exposure days
  far outnumber those above 30 °C. The `>30 °C` figures are the ones that correspond to the
  heat-mitigation estimates in the main text; `07` reports both.

## Assumptions worth knowing about

- **Coverage is capped at 100%** (`CAP_COVERAGE_AT_100`). e-Gestor's indicator is a capacity
  calculation — `(n_eSF × 3450 + n_eAB × 3000) / population × 100` — so a single team in a small
  municipality can exceed 100%, and the Ministry caps its published indicator for that reason. The
  source file holds the uncapped version: 7285 of about 115 000 municipality-years exceed 100
  (6.3%), up to 672. Capping matches the published indicator and keeps the 90% counterfactual
  interpretable. Set the switch to `FALSE` for the uncapped specification; the two panels are
  cached separately, so flipping it cannot silently reuse the other one.
  `tests/diagnose_coverage.R` reports the distribution and what share of municipality-days sits
  above the 90% counterfactual.
- **Municipality-years with no FHS record are treated as zero coverage**, not missing
  (`TREAT_MISSING_COVERAGE_AS_ZERO`). Defensible for years before a team was registered; the share
  affected is printed when the panel loads.
- **The 50+ figures come from two separately fitted models (50-69 and 70+) added together.** The
  point estimates add exactly, because each model is scaled by its own population (`pop_50_69`,
  `pop_70_plus`), so the two are counts of deaths among different people. The variances are also
  added, which assumes the two fits are uncorrelated; they share municipalities, weather and
  clusters, so the true covariance is not zero and is most likely positive, making these intervals
  slightly too narrow. `07` quantifies that without changing the reported numbers: `se_worst_case`
  is the widest interval any correlation can produce (`SE_50-69 + SE_70+`, i.e. rho = 1) and
  `rho_critical` is the correlation at which the interval would just touch zero. **`rho_critical >
  1` means no correlation whatsoever can put zero inside the interval**, so the assumption cannot
  have affected that conclusion.
- **Lags are built by date, not by row position.** `add_lags()` sorts by municipality and date,
  splits each municipality into contiguous segments and lags within them, so a break in the series
  — a municipality created mid-period, a year with no population estimate — can never make
  `temp_L1` mean "the day before the gap". The first `LAG_MAX` days of each segment have
  incomplete lag histories and are dropped by the estimator; the number affected is reported.
  Duplicated municipality-dates are a hard error.
- **A missing temperature is not a hot day.** The final arm of the binning tests `DAT >= 30`
  explicitly rather than falling through a catch-all, because `NA < 15` is `NA` rather than `FALSE`
  and a catch-all would silently code missing temperatures as the paper's headline exposure.
- **Coefficients are selected by exact name.** Percentile bin labels (`cold`, `ext_cold`, `hot`,
  `ext_hot`) are substrings of one another, so pattern matching would sum two bins and return the
  standard error of the wrong linear combination. A selection that matches nothing is an error, not
  an empty result — as is a missing bin, which returned as `estimate = 0, se = 0` would plot as a
  precisely estimated null.
- **Incomplete lag windows are flagged.** If fixest drops a lag term for collinearity, a "30-day
  cumulative effect" silently covers fewer days; the count is reported and warned on.
- **Joins are checked for uniqueness and fan-out** before and after they happen, so a duplicated
  municipality-year cannot silently multiply the panel.
- `STRICT = TRUE` makes data integrity checks stop the run. Keep it on for anything reported.
- The panel is cached as an `.rds` in `models/` after the first build and reused; delete it or call
  `get_panel(refresh = TRUE)` after changing the inputs.

## Tests

`tests/` holds harnesses that exercise the pipeline's logic against small synthetic inputs, with
the data loaders and writers stubbed, so they run in seconds and need no data:

```r
Rscript tests/test_helpers.R                 # lags, bins, coefficient selection, delta method
Rscript tests/test_descriptives.R            # Figure 1 maps and Appendix A Table A1
Rscript tests/test_placebo_future_change.R   # falsification test construction
Rscript tests/test_placebo_horizon.R         # the horizon is part of the cache key
Rscript tests/test_sens_scope.R              # scope constants and the resume logic
```

`tests/diagnose_coverage.R` is different: it reads the real coverage file and reports its
distribution. It is the right place to look before changing `CAP_COVERAGE_AT_100`.

See `GITHUB.md` for publishing this package and minting a citable DOI.
