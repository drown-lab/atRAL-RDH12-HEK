# atRAL dose-response: curve fitting, IC50 estimation, and two-factor ANOVA
#
# Input:  Cell_Viability/data/atRAL_cell_viability.csv
# Output: Cell_Viability/figures/
#   - atRAL_dose_response_curves.pdf/.png  (fitted curves, mean +/- SEM)
#   - atRAL_IC50_by_group.pdf/.png         (per-biorep IC50s with group means)
#   - atRAL_IC50_estimates.csv             (IC50 per biorep and per group w/ 95% CI)
#   - atRAL_IC50_anova.txt                 (two-way ANOVA on log10(IC50) + Tukey HSD)
#
# Run from the project root (working directory = repo root).

library(drc)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

data_file <- "Cell_Viability/data/atRAL_cell_viability.csv"
fig_dir   <- "Cell_Viability/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load and tidy -----------------------------------------------------------

viab <- read_csv(data_file, show_col_types = FALSE)
names(viab)[names(viab) == "ferrostatin-1"] <- "fer1"

viab <- viab |>
  dplyr::mutate(
    genotype = factor(genotype, levels = c("WT", "RDH12")),
    fer1     = factor(fer1, levels = c(0, 30),
                      labels = c("Vehicle", "Fer-1 (30 µM)")),
    biorep   = factor(biorep),
    group    = interaction(genotype, fer1, sep = " / ")
  )

# ---- Dose-response fits ------------------------------------------------------
# Log-logistic with variable slope, Bottom fixed at 0 and Top fixed at 100 --
# equivalent to Prism's "log(inhibitor) vs. normalized response -- variable
# slope", for consistency with prior fits of these data.
prism_model <- LL.4(fixed = c(NA, 0, 100, NA),
                    names = c("HillSlope", "Bottom", "Top", "IC50"))

# One fit per genotype x fer-1 group (all bioreps pooled) for plotting and
# group-level IC50 with delta-method 95% CI.
group_fits <- viab |>
  group_by(genotype, fer1) |>
  group_map(\(d, key) {
    fit <- drm(viability ~ atRAL, data = d, fct = prism_model)
    ed  <- ED(fit, 50, interval = "delta", display = FALSE)
    list(key = key, fit = fit,
         ic50 = tibble(genotype = key$genotype, fer1 = key$fer1,
                       IC50 = ed[1, "Estimate"], SE = ed[1, "Std. Error"],
                       CI_lower = ed[1, "Lower"], CI_upper = ed[1, "Upper"]))
  })

group_ic50 <- bind_rows(lapply(group_fits, `[[`, "ic50"))

# One fit per biorep (n = 3 per group) so each biological replicate yields an
# independent IC50 for the ANOVA.
biorep_ic50 <- viab |>
  group_by(genotype, fer1, biorep) |>
  group_modify(\(d, key) {
    fit <- drm(viability ~ atRAL, data = d, fct = prism_model)
    tibble(IC50 = ED(fit, 50, display = FALSE)[1, "Estimate"])
  }) |>
  ungroup()

write_csv(
  bind_rows(
    biorep_ic50 |> dplyr::mutate(level = "biorep"),
    group_ic50  |> dplyr::mutate(level = "group (pooled, 95% CI)")
  ),
  file.path(fig_dir, "atRAL_IC50_estimates.csv")
)

# ---- Two-factor ANOVA on log10(IC50) ----------------------------------------
# IC50s are log-normally distributed, so test on the log scale. Design is
# balanced (2 x 2 x 3), so Type I and Type III sums of squares coincide.

anova_fit <- aov(log10(IC50) ~ genotype * fer1, data = biorep_ic50)

anova_out <- file.path(fig_dir, "atRAL_IC50_anova.txt")
sink(anova_out)
cat("Two-factor ANOVA on log10(IC50)\n")
cat("Model: log10(IC50) ~ genotype * ferrostatin-1\n")
cat("n = 3 biological replicates per group; IC50 from per-biorep log-logistic\n")
cat("fits (variable slope, Top = 100, Bottom = 0; Prism 'log(inhibitor) vs.\n")
cat("normalized response -- variable slope' equivalent)\n\n")
print(summary(anova_fit))
cat("\nGroup means (IC50, uM; geometric mean across bioreps):\n")
biorep_ic50 |>
  group_by(genotype, fer1) |>
  dplyr::summarise(geomean_IC50 = 10^mean(log10(IC50)), .groups = "drop") |>
  as.data.frame() |>
  print()
cat("\nTukey HSD on genotype:fer1 cell means (log10 scale):\n")
print(TukeyHSD(anova_fit, "genotype:fer1"))
sink()

# ---- Dose-response figure ----------------------------------------------------

pal <- c(WT = "#0072B2", RDH12 = "#D55E00")  # Okabe-Ito, CVD-safe

# Mean +/- SEM per dose within each group
summ <- viab |>
  group_by(genotype, fer1, atRAL) |>
  dplyr::summarise(mean = mean(viability),
                   sem  = sd(viability) / sqrt(n()), .groups = "drop") |>
  dplyr::mutate(group = interaction(genotype, fer1, sep = " / "))

# Smooth prediction curves from the pooled group fits
pred <- bind_rows(lapply(group_fits, \(g) {
  doses <- exp(seq(log(25), log(400), length.out = 200))
  tibble(genotype = g$key$genotype, fer1 = g$key$fer1, atRAL = doses,
         viability = predict(g$fit, newdata = data.frame(atRAL = doses)))
}))

# IC50 (95% CI) annotation, bottom left under the curves
ic50_lab <- group_ic50 |>
  dplyr::mutate(
    label = sprintf("%s + %s: %.0f (%.0f–%.0f)",
                    genotype,
                    ifelse(fer1 == "Vehicle", "Vehicle", "Fer-1"),
                    IC50, CI_lower, CI_upper),
    x = 25,
    y = 26 - 6.5 * dplyr::row_number()
  )

p_curves <- ggplot(mapping = aes(atRAL, color = genotype, linetype = fer1)) +
  geom_line(data = pred, aes(y = viability), linewidth = 0.7) +
  geom_errorbar(data = summ,
                aes(ymin = mean - sem, ymax = mean + sem, group = group),
                width = 0.03, linewidth = 0.4, linetype = "solid",
                show.legend = FALSE) +
  geom_point(data = summ, aes(y = mean, shape = fer1), size = 2) +
  annotate("text", x = 25, y = 26, hjust = 0, size = 3.2, fontface = "bold",
           label = "IC[50]*', µM (95% CI)'", parse = TRUE) +
  geom_text(data = ic50_lab,
            aes(x = x, y = y, label = label, colour = genotype),
            inherit.aes = FALSE, hjust = 0, size = 3.2, show.legend = FALSE) +
  scale_x_log10(breaks = c(25, 50, 100, 200, 400)) +
  scale_color_manual(values = pal) +
  scale_shape_manual(values = c(16, 17)) +
  labs(x = "atRAL (µM)", y = "Percent Viability",
       color = "Genotype", shape = "Pretreatment", linetype = "Pretreatment") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "atRAL_dose_response_curves.pdf"), p_curves,
       width = 6.5, height = 4.5)
ggsave(file.path(fig_dir, "atRAL_dose_response_curves.png"), p_curves,
       width = 6.5, height = 4.5, dpi = 300)

# ---- IC50 figure -------------------------------------------------------------

p_ic50 <- ggplot(biorep_ic50, aes(fer1, IC50, colour = genotype)) +
  geom_point(position = position_dodge(width = 0.4), size = 2, alpha = 0.8) +
  stat_summary(fun = \(x) 10^mean(log10(x)), geom = "crossbar",
               position = position_dodge(width = 0.4),
               width = 0.3, linewidth = 0.4) +
  scale_y_log10(breaks = c(100, 150, 200, 300, 400)) +
  scale_color_manual(values = pal) +
  labs(x = NULL, y = expression(IC[50]~"(µM, log scale)"), colour = "Genotype") +
  theme_classic(base_size = 12)

ggsave(file.path(fig_dir, "atRAL_IC50_by_group.pdf"), p_ic50,
       width = 5, height = 4)
ggsave(file.path(fig_dir, "atRAL_IC50_by_group.png"), p_ic50,
       width = 5, height = 4, dpi = 300)
