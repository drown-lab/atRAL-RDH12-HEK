n_reps_per_condition <- meta %>%
  dplyr::count(condition, name = "n_reps")


run_fisher_contrast3 <- function(det_long, meta, contrast_name) {
  
  parts  <- strsplit(contrast_name, "_vs_")[[1]]
  group1 <- parts[1]
  group2 <- parts[2]
  
  n1 <- meta %>% filter(condition == group1) %>% nrow()
  n2 <- meta %>% filter(condition == group2) %>% nrow()
  
  counts <- det_long %>%
    filter(condition %in% c(group1, group2)) %>%
    group_by(protein, condition) %>%
    summarise(detected = sum(detected), .groups = "drop")
  
  g1 <- counts %>%
    filter(condition == group1) %>%
    transmute(protein,
              detected_1 = detected,
              not_detected_1 = n1 - detected)
  
  g2 <- counts %>%
    filter(condition == group2) %>%
    transmute(protein,
              detected_2 = detected,
              not_detected_2 = n2 - detected)
  
  df <- full_join(g1, g2, by = "protein") %>%
    replace_na(list(
      detected_1 = 0, not_detected_1 = n1,
      detected_2 = 0, not_detected_2 = n2
    ))
  
  res <- df %>%
    rowwise() %>%
    mutate(
      fisher_p = fisher.test(
        matrix(c(detected_1, not_detected_1,
                 detected_2, not_detected_2),
               nrow = 2, byrow = TRUE)
      )$p.value,
      delta_detect = (detected_1 / n1) - (detected_2 / n2)
    ) %>%
    ungroup() %>%
    mutate(
      fisher_p_adj = p.adjust(fisher_p, method = "BH"),
      contrast = contrast_name
    )
  
  res
}
fisher_results <- purrr::map_dfr(
  contrasts,
  ~run_fisher_contrast3(det_long, meta, .x)
)


fisher_results %>%
  filter(protein == "RAB15") %>%
  select(contrast,
         detected_1, not_detected_1,
         detected_2, not_detected_2,
         fisher_p, fisher_p_adj, delta_detect)


fisher_results %>%
  filter(protein == "RDH12") %>%
  select(contrast,
         detected_1, not_detected_1,
         detected_2, not_detected_2,
         fisher_p, fisher_p_adj, delta_detect)
#==============================================

library(dplyr)
library(ggplot2)
library(forcats)

library(dplyr)
library(tidyr)

delta_hist <- fisher_results %>%
  mutate(
    frac_1 = detected_1 / (detected_1 + not_detected_1),
    frac_2 = detected_2 / (detected_2 + not_detected_2),
    delta_detect = frac_1 - frac_2,
    # bin to exact thirds to avoid floating rounding issues
    delta_bin = round(delta_detect * 3) / 3,
    delta_bin = factor(
      delta_bin,
      levels = c(-1, -2/3, -1/3, 0, 1/3, 2/3, 1),
      labels = c("-1", "-2/3", "-1/3", "0", "1/3", "2/3", "1")
    )
  ) %>%
  dplyr::count(contrast, delta_bin, name = "n_proteins") %>%
  group_by(contrast) %>%
  tidyr::complete(delta_bin, fill = list(n_proteins = 0)) %>%
  ungroup()

ggplot(delta_hist, aes(x = delta_bin, y = n_proteins, fill = contrast)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8) +
  theme_bw(base_size = 12) +
  labs(
    title = "Detection frequency shift distribution (all contrasts)",
    x = expression(Delta~"detection fraction (group1 - group2)"),
    y = "# proteins",
    fill = "Contrast"
  ) +
  facet_wrap(~contrast, nrow=2)+
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )
