# ==============================================================================
# 04_non_interacted_figure.R
# ------------------------------------------------------------------------------
# Appendix B: the temperature-mortality relationship from the model WITHOUT the
# FHS interaction, i.e. pooled across all observed coverage levels.  This is the
# quantity most comparable with the existing literature.
#
# Reads the CSVs written by 03_main_dlm.R; fits nothing itself.
#
# The reference bin is plotted as an open symbol, not as an estimate: it is a
# normalisation (zero by construction), and a filled point with no interval
# reads as a precisely estimated null.
#
# Inputs : csv/cumulative_nointeraction_ALL.csv
# Outputs: figures/appendixB_temperature_mortality_no_interaction.png
# ==============================================================================

f <- file.path(dir_out("csv"), "cumulative_nointeraction_ALL.csv")
if (!file.exists(f)) stop("Run 03_main_dlm.R first: ", f, " not found.")

cum <- readr::read_csv(f, show_col_types = FALSE) %>%
  filter(outcome == "rate_total") %>%
  mutate(region   = factor(region, levels = unname(REGION_LABELS)),
         temp_bin = factor(temp_bin, levels = BIN_LEVELS),
         is_ref   = temp_bin == ref_bin)

p <- ggplot(cum, aes(temp_bin, estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(data = ~ dplyr::filter(.x, !is_ref),
                aes(ymin = lwr, ymax = upr), width = 0.2) +
  geom_point(aes(shape = is_ref, fill = is_ref), size = 2.4) +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 1), guide = "none") +
  scale_fill_manual(values = c(`FALSE` = "black", `TRUE` = NA), guide = "none") +
  facet_wrap(~ region, scales = "free_y") +
  theme_minimal() +
  labs(x = "Temperature bin (ºC)",
       y = paste0(LAG_MAX, "-day cumulative impact on mortality rate ",
                  "(deaths per 100 000)"),
       caption = paste0("Open symbol marks the region-specific reference bin, ",
                        "which is zero by construction. Vertical lines are 95% ",
                        "confidence intervals from municipality-clustered ",
                        "standard errors."))

save_fig(p, "appendixB_temperature_mortality_no_interaction.png",
         width = 11, height = 7)

say("04_non_interacted_figure.R done")
