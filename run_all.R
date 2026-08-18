# ==============================================================================
# run_all.R
# ------------------------------------------------------------------------------
# Reproduces every result in the paper, in the order the results appear in the
# text.  Run this file and nothing else:
#
#     Rscript run_all.R          (or source("run_all.R") in an R session)
#
# Before the first run, edit the three paths at the top of 00_config.R.
#
# Each step is timed and logged.  A step that fails is reported and the pipeline
# continues, so one missing input file does not stop everything; the summary at
# the end lists what succeeded and what did not.
#
# Runtime is dominated by the model fits: expect minutes per region x outcome on
# a typical workstation, and several hours for the full pipeline.
# ==============================================================================

t_start <- Sys.time()

source("00_config.R")
source("01_helpers.R")

STEPS <- c(
  "02_descriptives.R",            # Results, first paragraph; Figure 1; Appendix A
  "03_main_dlm.R",                # Figure 2; Appendix B and C inputs
  "04_non_interacted_figure.R",   # Appendix B
  "05_table1.R",                  # Table 1
  "06_gdp_split.R",               # Appendix D
  "07_deaths_averted.R",          # Appendix E
  "08_sens_placebo_future_change.R",   # Appendix F(i), Figures F1-F2
  "09_sens_lag_windows.R",        # Appendix F(ii)
  "10_sens_coverage_bands.R",     # Appendix F(iii)
  "11_sens_exclude_large_muni.R", # Appendix F(iv)
  "12_sens_percentile_bins.R",    # Appendix F(v)
  "13_session_info.R"             # software record
)

# ------------------------------------------------------------------------------
# Is a step already complete?
# ------------------------------------------------------------------------------
# A single sentinel file is not enough.  Several steps were restructured after
# they had already run once, and the new outputs sit alongside the old files of
# the same name - so a sentinel could be present while most of the step's work
# was still outstanding, and run_all.R would silently skip it.
#
# The manifest in 14_check_outputs.R already lists every file each step is
# expected to produce.  A step is skipped only when ALL of them are present and
# non-empty, which cannot drift out of step with what the scripts actually write.
# ------------------------------------------------------------------------------
.sourced_by_run_all <- TRUE
source("14_check_outputs.R")
MANIFEST <- expected_outputs()

STEP_KEY <- c(
  "02_descriptives.R"                = "02_descriptives",
  "03_main_dlm.R"                    = "03_main_dlm",
  "04_non_interacted_figure.R"       = "04_non_interacted",
  "05_table1.R"                      = "05_table1",
  "06_gdp_split.R"                   = "06_gdp_split",
  "07_deaths_averted.R"              = "07_deaths_averted",
  "08_sens_placebo_future_change.R"  = "08_placebo_future_change",
  "09_sens_lag_windows.R"            = "09_lag_windows",
  "10_sens_coverage_bands.R"         = "10_coverage_bands",
  "11_sens_exclude_large_muni.R"     = "11_exclude_large",
  "12_sens_percentile_bins.R"        = "12_percentile_bins",
  "13_session_info.R"                = "13_session_info"
)

step_status <- function(step) {
  key <- STEP_KEY[[step]]
  want <- MANIFEST[MANIFEST$step == key, , drop = FALSE]
  if (!nrow(want)) return(list(complete = FALSE, have = 0, need = 0))
  p  <- file.path(OUT_DIR, want$subdir, want$file)
  ok <- file.exists(p) & !is.na(file.size(p)) & file.size(p) > 0
  list(complete = all(ok), have = sum(ok), need = length(ok))
}

log_path <- file.path(dir_out("logs"),
                      format(Sys.time(), "run_all_%Y%m%d_%H%M%S.log"))
say("Logging to ", log_path)

status <- data.frame(step = STEPS, ok = NA, minutes = NA_real_,
                     message = NA_character_, stringsAsFactors = FALSE)

for (i in seq_along(STEPS)) {

  step <- STEPS[i]

  st <- step_status(step)
  if (SKIP_EXISTING && st$complete) {
    say("################ ", step, " - SKIPPED (all ", st$need,
        " outputs present) ############")
    status$ok[i] <- TRUE; status$minutes[i] <- 0
    status$message[i] <- sprintf("skipped: %d of %d outputs already present",
                                 st$have, st$need)
    next
  }
  if (st$have > 0 && st$have < st$need)
    say("(", step, ": ", st$have, " of ", st$need,
        " outputs present - running to produce the rest)")

  say("################ ", step, " ################")
  t0 <- Sys.time()

  res <- tryCatch({
    source(step, local = new.env(parent = globalenv()))
    list(ok = TRUE, msg = "")
  }, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))

  status$ok[i]      <- res$ok
  status$minutes[i] <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  status$message[i] <- res$msg

  if (!res$ok) say("FAILED: ", step, " - ", res$msg)
  gc()
}

status$minutes <- round(status$minutes, 1)
write.csv(status, file.path(dir_out("logs"), "run_all_status.csv"), row.names = FALSE)

say("=========================================")
say(sprintf("Finished in %.1f minutes: %d of %d steps succeeded",
            as.numeric(difftime(Sys.time(), t_start, units = "mins")),
            sum(status$ok, na.rm = TRUE), nrow(status)))
print(status[, c("step", "ok", "minutes")])
if (any(!status$ok, na.rm = TRUE)) {
  say("Failed steps:")
  print(status[!status$ok, c("step", "message")])
}
say("Outputs are in ", OUT_DIR)

# Whether the steps reported success is not the same question as whether every
# expected file is on disk, so verify the outputs themselves.
res <- check_outputs()
if (!res$complete)
  say("Run is INCOMPLETE - see the list above and logs/output_check.csv. ",
      "Re-running run_all.R will resume where it stopped.")
