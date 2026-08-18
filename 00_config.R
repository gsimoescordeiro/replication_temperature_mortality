# ==============================================================================
# 00_config.R
# ------------------------------------------------------------------------------
# Single place where every path, package and global setting lives.  Every other
# script starts by sourcing this file; nothing else in the pipeline contains a
# hard-coded path.
#
# TO RUN THE PIPELINE ON A NEW MACHINE: edit BASE_DIR, DATA_DIR and OUT_DIR
# below (and the three file names in DATA_FILES if yours differ).  Nothing else
# needs to change.
#
# Paper: Cordeiro G, Lagarde M, Pereda P, Millett C, Goncalves J.  Effects of
# primary health care coverage on temperature-related mortality in Brazil:
# national quasi-experimental study.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. PATHS
# ------------------------------------------------------------------------------
# BASE_DIR  : root of this reproduction package (the folder holding these scripts)
# DATA_DIR  : where the three prepared .dta inputs live (never written to)
# OUT_DIR   : where everything produced by the pipeline is written
BASE_DIR <- "Z:/Projects/temp_mortality/scripts_reproduction"
DATA_DIR <- "Z:/Projects/temp_mortality"
OUT_DIR  <- "Z:/Projects/temp_mortality/reproduction_output"

DATA_FILES <- list(
  panel = "final_merged_data_with_EHI_A_20plus.dta",  # municipality-day panel, deaths 20+
  pop   = "pop_race.dta",                             # annual municipal population by age band
  esf   = "#final/data/ESF_coverage.dta",             # annual municipal FHS coverage
  gdp   = "gdp.dta"                                   # municipal GDP (Appendix D split only)
)

# Where the assembled panel is cached between sessions.  Building it takes a few
# minutes and writing it takes longer again on a network drive, so point this at
# a local disk if OUT_DIR is on a mapped share (Z:, OneDrive, etc.).
# Example: CACHE_DIR <- "C:/temp/temp_mortality_cache"
CACHE_DIR <- NULL          # NULL = keep the cache under OUT_DIR/models

# Output subfolders, created on demand by dir_out()
OUT_SUBDIRS <- c("tables", "figures", "csv", "models", "logs")

# ------------------------------------------------------------------------------
# 2. PACKAGES
# ------------------------------------------------------------------------------
# Versions used for the published results: R 4.3.2, fixest 0.11.2.
# Install the pinned estimation package with:
#   remotes::install_version("fixest", version = "0.11.2")
REQUIRED_PKGS <- c("haven", "dplyr", "tidyr", "stringr", "readr",
                   "purrr", "ggplot2", "patchwork", "scales", "fixest")

missing_pkgs <- REQUIRED_PKGS[!vapply(REQUIRED_PKGS, requireNamespace,
                                      logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "),
       "\nInstall them, then re-run.", call. = FALSE)
}
invisible(lapply(REQUIRED_PKGS, function(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))))

if (getRversion() < "4.0.0")
  warning("Results were produced under R 4.3.2; you are on ", getRversion())
if (packageVersion("fixest") != "0.11.2")
  message("NOTE: results were produced with fixest 0.11.2; you have ",
          packageVersion("fixest"),
          ". Coefficients should be identical, but check before publishing.")

# ------------------------------------------------------------------------------
# 3. STUDY DESIGN CONSTANTS
# ------------------------------------------------------------------------------
STUDY_YEARS <- c(2000, 2019)

# Macro-regions are identified by the first digit of the 6-digit IBGE municipality
# code.  Models are estimated separately by region because the temperature
# distribution, the FHS implementation history and the minimum-mortality
# temperature all differ markedly across them.
REGIONS <- c("1", "2", "3", "4", "5")
REGION_LABELS <- c("1" = "North", "2" = "Northeast", "3" = "Southeast",
                   "4" = "South", "5" = "Centre-West")

# Absolute daily mean temperature bins (degrees C).
BIN_LEVELS <- c("<15", "15-20", "20-25", "25-30", ">30")

# Reference bin = the region's minimum-mortality temperature bin.  All reported
# effects are excess mortality relative to this bin, so it differs by region.
REF_BIN <- c("1" = "25-30", "2" = "25-30", "3" = "20-25",
             "4" = "20-25", "5" = "25-30")

# Distributed-lag window.  The main analysis cumulates over lags 0..30, i.e. the
# 31 days beginning with the day of exposure.
LAG_MAX      <- 30
LAG_MAX_SENS <- c(7, 14)          # Appendix F(ii)

# FHS coverage levels (percentage points) at which the interacted model is
# evaluated for Figure 2, and the high-coverage counterfactual used throughout.
FHS_LEVELS <- c(0, 50, 90)
FHS_HIGH   <- 90

# Coverage bands for the split-sample analysis, Appendix F(iii).  NOTE these are
# fixed cut-points, not empirical terciles - describe them as coverage bands.
FHS_BANDS <- list(low = c(-Inf, 25), mid = c(25, 75), high = c(75, Inf))

# Outcomes.  All rates are deaths per 100 000 population.  Numerators are
# already restricted to deaths at age 20+ in the source panel.
OUTCOMES <- c(
  total       = "rate_total",
  respiratory = "rate_respiratory",
  circulatory = "rate_circulatory",
  external    = "rate_external",
  cancer      = "rate_cancer",
  nutrition   = "rate_nutrition",
  infectious  = "rate_infectious",
  age_20_49   = "rate_age_20_49",
  age_50_69   = "rate_age_50_69",
  age_70_plus = "rate_age_70_plus"
)

# Scope of the robustness checks, Appendices D, E and F.  In the submitted
# manuscript these appendices cover the Southeast and South only, and every
# figure in Appendix F is broken out by the three age groups rather than shown
# for all-cause mortality.  Scripts 06 and 08-12 follow that structure.
# Add "rate_total" to SENS_OUTCOMES, or more regions to SENS_REGIONS, to widen
# it; each addition is one more model per script.
SENS_REGIONS  <- c("3", "4")
SENS_OUTCOMES <- c("rate_age_20_49", "rate_age_50_69", "rate_age_70_plus")

# Panel labels used by the Appendix F figures.
AGE_LABELS <- c(rate_total       = "All adults 20+",
                rate_age_20_49   = "20-49 years old",
                rate_age_50_69   = "50-69 years old",
                rate_age_70_plus = "70+ years old")

age_label <- function(x) unname(ifelse(x %in% names(AGE_LABELS), AGE_LABELS[x], x))

# Outcomes used for the deaths-averted scenarios (Appendix E).  The paper reports
# adults aged 50 and over, where the protective association is concentrated; that
# group is built from these two separately fitted age bands, whose death counts
# add because each is scaled by its own population.
OUTCOMES_AVERTED <- c("rate_age_50_69", "rate_age_70_plus")
REGIONS_AVERTED  <- c("3", "4")     # Southeast and South

# ------------------------------------------------------------------------------
# 4. BEHAVIOUR SWITCHES
# ------------------------------------------------------------------------------
# SKIP_EXISTING: when TRUE, work whose outputs already exist and are non-empty is
# not repeated.  A full run takes many hours, so an interruption - a network
# write failure, a reboot - should cost only the step it stopped on.  Set FALSE
# to force everything to be recomputed from scratch.
SKIP_EXISTING <- TRUE

# STRICT: stop on a failed data integrity check rather than warn.  Keep TRUE for
# anything whose results are reported; set FALSE only when exploring.
STRICT <- TRUE

# e-Gestor's coverage indicator is a capacity calculation, not a count of people
# enrolled: (n_eSF * 3450 + n_eAB * 3000) / population estimate * 100.  In a
# small municipality one team can therefore exceed 100%, and the Ministry's own
# methodological note caps the published indicator ("O indicador de cobertura
# nao deve passar de 100%").  The source file here holds the UNCAPPED version:
# 7285 of about 115 000 municipality-years exceed 100 (6.3%), reaching 672.
# Capping matches the published indicator, matches how the paper defines the
# variable, and keeps the 90% counterfactual meaningful.
# Set FALSE to reproduce the uncapped specification as a sensitivity analysis;
# the two versions are cached separately, so switching is safe.
CAP_COVERAGE_AT_100 <- TRUE

# Coverage recorded as missing after the join is treated as zero coverage.  This
# is deliberate: e-Gestor has no record for municipality-years before a team was
# registered, and coverage in those years genuinely was zero.  01_helpers.R
# reports how many municipality-years this affects so the assumption is visible.
TREAT_MISSING_COVERAGE_AS_ZERO <- TRUE

set.seed(20260806)   # nothing here is stochastic; set for safety

# ------------------------------------------------------------------------------
# 5. DERIVED PATHS AND SMALL UTILITIES
# ------------------------------------------------------------------------------
data_path <- function(key) file.path(DATA_DIR, DATA_FILES[[key]])

dir_out <- function(...) {
  p <- file.path(OUT_DIR, ...)
  if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  p
}
invisible(lapply(OUT_SUBDIRS, dir_out))

say <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

fail <- function(...) if (STRICT) stop(..., call. = FALSE) else warning(..., call. = FALSE)

message("Config loaded. Output root: ", OUT_DIR)
