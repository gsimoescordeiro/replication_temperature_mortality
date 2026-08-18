# ==============================================================================
# 13_session_info.R
# ------------------------------------------------------------------------------
# Records the exact software environment that produced the results, and writes
# it next to them.  The BMJ asks for software and package versions to be stated;
# this file is the machine-generated record behind the sentence in the Methods.
#
# Outputs: logs/session_info.txt
# ==============================================================================

path <- file.path(dir_out("logs"), "session_info.txt")
con  <- file(path, open = "wt")

writeLines(c(
  R.version.string,
  paste("Platform:", R.version$platform),
  paste("Run at:  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "Package versions:"), con)

for (p in REQUIRED_PKGS)
  writeLines(paste0("  ", p, " ", as.character(packageVersion(p))), con)

writeLines(c("", "Results in the paper were produced with R 4.3.2 and fixest 0.11.2.",
             "", "Full sessionInfo():", ""), con)
capture.output(sessionInfo(), file = con)
close(con)

say("wrote ", path)
say("13_session_info.R done")
