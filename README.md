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
| `02_descriptives.R` | Mortality rates and coverage by region; exposure days per bin; municipal coverage maps; Table A1 in the submitted layout | Results ¶1; **Figure 1**; **Appendix A, Table A1** |
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
municipality-days: the standard deviation there is variation across places and days.
Table A1 as submitted collapses to a region-day series first and reports the standard
deviation of that series over time, which is an order of magnitude smaller because
averaging across a region removes the cross-sectional spread. Section 4 reproduces the
submitted table; section 1 is kept because the Methods and Results quote
municipality-day figures. Both are written out, labelled.

## What changed relative to the code that produced the submitted manuscript

These are corrections, not cosmetic edits. Numbers will move; how much has to be checked by
re-running. Each is annotated at the point in the code where it applies.

**Changes that alter reported estimates**

1. **`date` fixed effect restored in the deaths-averted models** (`07`). The original
   averted-deaths script omitted the day-by-year fixed effect that all other models use, so it
   estimated a different specification from the one the Methods describe. All Appendix E numbers
   need re-running.
2. **Exact coefficient selection** (`01_helpers.R::bin_coef_index`). Coefficients were previously
   selected with `grep(label, …)`. With percentile bin labels (`cold`, `ext_cold`, `hot`,
   `ext_hot`) that matched two bins at once, so the reported percentile-bin effects summed two
   bins and their standard errors were those of the wrong linear combination. Selection is now by
   exact name and errors out if nothing matches.
3. **The falsification test was replaced** (`08`). The two-year-lead test cannot discriminate for a
   slowly-moving treatment, and its apparent pass came from filling 2018-19 with zero coverage.
   See "The falsification test was replaced" below.
4. **Figure 2 is estimated at 0/50/90% coverage** (`03`). The models previously predicted at
   0/10/90 and the 50% series shown in the paper came from a post-processing step that is not in
   the original code folder. It is now estimated directly, so its interval is `sqrt(L'VL)` at
   loading 50 — which is not what interpolating between the 0% and 90% intervals would give,
   because the standard error is a square-root quadratic in coverage. **No post-processing step is
   needed any more**; `post_process_FHS_50_all_runs.R` is superseded.

5. **FHS coverage is capped at 100%** (`CAP_COVERAGE_AT_100` in `00_config.R`). e-Gestor's
   indicator is a capacity calculation — `(n_eSF x 3450 + n_eAB x 3000) / population x 100` — so a
   single team in a small municipality can exceed 100%, and the Ministry caps its published
   indicator for that reason. The source file holds the uncapped version: 7285 of about 115 000
   municipality-years exceed 100 (6.3%), up to 672. Capping matches the published indicator and the
   paper's own definition of the variable, and it keeps the 90% counterfactual interpretable. Set
   the switch to `FALSE` to reproduce the uncapped specification; the two are cached separately.

**Changes that prevent silent errors**

6. **Missing temperatures no longer become the hottest bin.** `case_when(…, TRUE ~ ">30")` sent
   any `NA` temperature into the `>30` bin, because `NA < 15` is `NA`, not `FALSE`. The final arm
   now tests `DAT >= 30` explicitly.
7. **Lag construction is gap-aware.** `dplyr::lag()` shifts by row position, not by date.
   `add_lags()` sorts by municipality and date, then splits each municipality into contiguous
   segments and lags within them, so a break in the series (a municipality created mid-period, or
   a year with no population estimate) can never make `temp_L1` mean "the day before the gap". The
   first `LAG_MAX` days of each segment have incomplete lag histories and are dropped by the
   estimator, exactly as completing the daily grid would do; the number affected is reported.
   Duplicated municipality-dates remain a hard error.
8. **Joins are checked for uniqueness and fan-out** before and after they happen, so a duplicated
   municipality-year cannot silently multiply the panel.
9. **Municipalities with no population data are excluded** from the large-municipality sensitivity
   analysis. `max(pop_total, na.rm = TRUE)` returns `-Inf` when every value is missing, and
   `-Inf <= 500000` is `TRUE`, so they were previously retained.
10. **Missing bins are an error, not a zero.** Extractors previously returned `estimate = 0,
    se = 0` when no coefficient matched, which plots as a precisely estimated null.
11. **Incomplete lag windows are flagged.** If fixest drops a lag term for collinearity, the
    "30-day cumulative effect" silently covers fewer days; the count is now reported and warned on.
12. **Exposure days are counted on the estimation sample** in `07`, excluding the lag burn-in and
    rows with missing precipitation, which previously contributed exposure but no identification.

## Scope of Appendices D, E and F

The submitted appendix is in two halves. Appendices A, B and C cover all five macro-regions;
Appendices D, E and F cover the **Southeast and South only**, and every figure in Appendix F is
broken out by the three **age groups** (20-49, 50-69, 70+) rather than shown for all-cause
mortality. Scripts `06` and `08`-`12` follow that structure, driven by two constants in `00_config.R`:

```r
SENS_REGIONS  <- c("3", "4")                       # Southeast, South
SENS_OUTCOMES <- c("rate_age_20_49", "rate_age_50_69", "rate_age_70_plus")
```

Each script writes one figure per region with the age groups as panels, and one CSV carrying an
`outcome` column. Widen either constant to cover more ground; each addition is one more model per
script. Every fit is cached as an `.rds` under `models/`, so narrowing the scope later never
discards work and re-running only computes what is missing.

## The falsification test was replaced

The previous `08_sens_placebo_lead.R` substituted the two-year lead of coverage for coverage itself. That test
cannot work here. FHS coverage is a slow-moving stock, so the lead is a near-perfect proxy for the
treatment and reproduces the main result whether or not the effect is causal.

In the code behind the submitted manuscript the lead was undefined for 2018-19 and those rows were
filled with coverage `= 0` — the two highest-coverage years in the panel, about a tenth of the
sample, recoded as "no FHS". That measurement error attenuated the placebo interaction toward zero,
which is why the test appeared to pass. It passed because the placebo variable was corrupted.

`08_sens_placebo_future_change.R` puts the placebo on the **change** instead, and estimates it in
the **same model** as actual coverage:

```
D_it = coverage(t+3) − coverage(t+1)
rate ~ Σ i(temp_L) + Σ i(temp_L)×coverage_t + Σ i(temp_L)×D_it + … | FE
```

Both endpoints of `D` are in the future, and future growth is only weakly correlated with the
current level (−0.31 to −0.52 across regions), so the two interactions are separately identified
and the test has power. A causal effect implies the contemporaneous term is negative on hot days
and the placebo term is not. `appendixF_placebo_future_change_diagnostics.csv` reports the
correlation per region — report it alongside the result, because it is what distinguishes this test
from the one it replaces.

The old `08_sens_placebo_lead.R` has been removed from the package; this file replaces it.

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
(`cache_vcov()`) and reused by every extractor. Computing it per call meant about twenty full
covariance computations per region-outcome, which dominated runtime — the fits themselves take a
couple of minutes, the extraction was taking close to an hour.

## Two things to state in the paper

- **Coverage bands are not terciles.** The split uses fixed cut-points at 25% and 75%
  (`FHS_BANDS`), which do not produce equal groups. The Methods currently say "tercile".
- **The high-coverage counterfactual is coverage set to exactly 90% everywhere**, so
  municipality-years already above 90% contribute negatively to "additional avoidable deaths".
  `07` reports what share of municipality-days that affects.

## Assumptions worth knowing about

- Coverage is **capped at 100%** by default; `tests/diagnose_coverage.R` reports the distribution,
  how much of the variance sits above 100, and what share of municipality-days sits above the 90%
  counterfactual.
- **The 50+ figures come from two separately fitted models (50-69 and 70+) added together**, as in
  the submitted manuscript. The point estimates add exactly, because each model is scaled by its own
  population (`pop_50_69`, `pop_70_plus`), so the two are counts of deaths among different people.
  The variances are also added, which assumes the two fits are uncorrelated; they share
  municipalities, weather and clusters, so the true covariance is not zero and is most likely
  positive, making these intervals slightly too narrow. `07` quantifies that without changing the
  reported numbers: `se_worst_case` is the widest interval any correlation can produce
  (`SE_50-69 + SE_70+`, i.e. rho = 1) and `rho_critical` is the correlation at which the interval
  would just touch zero. **`rho_critical > 1` means no correlation whatsoever can put zero inside
  the interval**, so the assumption cannot have affected that conclusion.
- Municipality-years with no FHS record are treated as **zero coverage**, not missing
  (`TREAT_MISSING_COVERAGE_AS_ZERO` in `00_config.R`). Defensible for years before a team was
  registered; the share affected is printed when the panel loads.
- `STRICT = TRUE` makes data integrity checks stop the run. Keep it on for anything reported.
- The panel is cached as an `.rds` in `models/` after the first build and reused; delete it or call
  `get_panel(refresh = TRUE)` after changing the inputs.
