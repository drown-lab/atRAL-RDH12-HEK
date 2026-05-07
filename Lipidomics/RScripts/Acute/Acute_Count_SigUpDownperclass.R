# =========================================================
# Acute lipid class summary:
# number of significant up/down lipids per class
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(tidyr)
  library(forcats)
})

# ---------------------------------------------------------
# 1. INPUT
# ---------------------------------------------------------
acute_file <- "Lipidomics/output_txts/Acute_DEA_results_v1_corrected.csv"

out_dir <- "Lipidomics/Figures/Acute/Lipid_class_direction_counts/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# significance thresholds
padj_cutoff <- 0.1
fc_cutoff <- 1.3
lfc_cutoff <- log2(fc_cutoff)

# output files
out_counts_csv <- file.path(out_dir, "Acute_lipid_class_up_down_counts.csv")
out_long_csv   <- file.path(out_dir, "Acute_lipid_class_up_down_counts_long.csv")

out_stacked_pdf <- file.path(out_dir, "Acute_lipid_class_up_down_stacked_bar.pdf")
out_stacked_png <- file.path(out_dir, "Acute_lipid_class_up_down_stacked_bar.png")

out_grouped_pdf <- file.path(out_dir, "Acute_lipid_class_up_down_grouped_bar.pdf")
out_grouped_png <- file.path(out_dir, "Acute_lipid_class_up_down_grouped_bar.png")

out_lollipop_pdf <- file.path(out_dir, "Acute_lipid_class_up_down_lollipop.pdf")
out_lollipop_png <- file.path(out_dir, "Acute_lipid_class_up_down_lollipop.png")

# ---------------------------------------------------------
# 2. READ DATA
# ---------------------------------------------------------
lipid_df <- read.csv(
  acute_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_cols <- c(
  "name",
  "RDH12_100_vs_RDH12_Veh_p.val",
  "RDH12_200_vs_RDH12_Veh_p.val",
  "RDH12_100_vs_RDH12_Veh_ratio",
  "RDH12_200_vs_RDH12_Veh_ratio"
)

missing_cols <- setdiff(required_cols, colnames(lipid_df))
if (length(missing_cols) > 0) {
  stop(
    "These required columns are missing:\n",
    paste(missing_cols, collapse = "\n")
  )
}

# ---------------------------------------------------------
# 3. ASSIGN LIPID CLASS
# ---------------------------------------------------------
get_lipid_class <- function(x) {
  case_when(
    str_detect(x, "^LPC\\s")   ~ "LPC",
    str_detect(x, "^PC O-")    ~ "PC ether",
    str_detect(x, "^PC P-")    ~ "PC plasmalogen",
    str_detect(x, "^PC\\s")    ~ "PC",
    str_detect(x, "^LPE\\s")   ~ "LPE",
    str_detect(x, "^PE O-")    ~ "PE ether",
    str_detect(x, "^PE P-")    ~ "PE plasmalogen",
    str_detect(x, "^PE\\s")    ~ "PE",
    str_detect(x, "^PI O-")    ~ "PI ether",
    str_detect(x, "^PI\\s")    ~ "PI",
    str_detect(x, "^PS O-")    ~ "PS ether",
    str_detect(x, "^PS\\s")    ~ "PS",
    str_detect(x, "^PG\\s")    ~ "PG",
    str_detect(x, "^DG O-")    ~ "DG ether",
    str_detect(x, "^DG\\s")    ~ "DG",
    str_detect(x, "^TG\\s")    ~ "TG",
    str_detect(x, "^Cer\\s")   ~ "Cer",
    str_detect(x, "^dhCer\\s") ~ "dhCer",
    str_detect(x, "^SM\\s")    ~ "SM",
    str_detect(x, "^CE\\s")    ~ "CE",
    str_detect(x, "^CAR\\s")   ~ "CAR",
    TRUE ~ "Other"
  )
}

# ---------------------------------------------------------
# 4. FLAG UP / DOWN BY CONTRAST
# ---------------------------------------------------------
lipid_df2 <- lipid_df %>%
  mutate(
    name = str_trim(name),
    Class = get_lipid_class(name),
    
    direction_100 = case_when(
      !is.na(`RDH12_100_vs_RDH12_Veh_p.val`) &
        !is.na(`RDH12_100_vs_RDH12_Veh_ratio`) &
        `RDH12_100_vs_RDH12_Veh_p.val` <= padj_cutoff &
        `RDH12_100_vs_RDH12_Veh_ratio` >= lfc_cutoff ~ "Up",
      !is.na(`RDH12_100_vs_RDH12_Veh_p.val`) &
        !is.na(`RDH12_100_vs_RDH12_Veh_ratio`) &
        `RDH12_100_vs_RDH12_Veh_p.val` <= padj_cutoff &
        `RDH12_100_vs_RDH12_Veh_ratio` <= -lfc_cutoff ~ "Down",
      TRUE ~ "NS"
    ),
    
    direction_200 = case_when(
      !is.na(`RDH12_200_vs_RDH12_Veh_p.val`) &
        !is.na(`RDH12_200_vs_RDH12_Veh_ratio`) &
        `RDH12_200_vs_RDH12_Veh_p.val` <= padj_cutoff &
        `RDH12_200_vs_RDH12_Veh_ratio` >= lfc_cutoff ~ "Up",
      !is.na(`RDH12_200_vs_RDH12_Veh_p.val`) &
        !is.na(`RDH12_200_vs_RDH12_Veh_ratio`) &
        `RDH12_200_vs_RDH12_Veh_p.val` <= padj_cutoff &
        `RDH12_200_vs_RDH12_Veh_ratio` <= -lfc_cutoff ~ "Down",
      TRUE ~ "NS"
    )
  )

# ---------------------------------------------------------
# 5. COUNT UP / DOWN LIPIDS PER CLASS
# ---------------------------------------------------------
count_100 <- lipid_df2 %>%
  filter(direction_100 %in% c("Up", "Down")) %>%
  count(Class, Direction = direction_100, name = "Count") %>%
  mutate(Contrast = "100 vs Veh")

count_200 <- lipid_df2 %>%
  filter(direction_200 %in% c("Up", "Down")) %>%
  count(Class, Direction = direction_200, name = "Count") %>%
  mutate(Contrast = "200 vs Veh")

count_long <- bind_rows(count_100, count_200)

# fill in missing combinations with zero
all_classes <- sort(unique(lipid_df2$Class))
count_long <- count_long %>%
  complete(
    Contrast = c("100 vs Veh", "200 vs Veh"),
    Class = all_classes,
    Direction = c("Up", "Down"),
    fill = list(Count = 0)
  )

# signed counts for lollipop plot
count_long <- count_long %>%
  mutate(
    SignedCount = if_else(Direction == "Down", -Count, Count)
  )

# order classes by total significant lipids across both contrasts
class_order <- count_long %>%
  group_by(Class) %>%
  summarise(Total = sum(Count), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(Class)

count_long <- count_long %>%
  mutate(Class = factor(Class, levels = class_order))

count_wide <- count_long %>%
  select(Contrast, Class, Direction, Count) %>%
  tidyr::pivot_wider(
    names_from = c(Contrast, Direction),
    values_from = Count
  ) %>%
  arrange(Class)

write.csv(count_wide, out_counts_csv, row.names = FALSE)
write.csv(count_long, out_long_csv, row.names = FALSE)

# ---------------------------------------------------------
# 6. COLORS
# ---------------------------------------------------------
direction_colors <- c(
  "Up" = "#B40426",
  "Down" = "#3B4CC0"
)

# ---------------------------------------------------------
# 7. STACKED BAR PLOT
# ---------------------------------------------------------
p_stacked <- ggplot(count_long, aes(x = Class, y = Count, fill = Direction)) +
  geom_col() +
  facet_wrap(~ Contrast, ncol = 1) +
  scale_fill_manual(values = direction_colors) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Acute significant lipids per class",
    x = NULL,
    y = "Number of significant lipids",
    fill = "Direction"
  )

print(p_stacked)

ggsave(out_stacked_pdf, p_stacked, width = 9, height = 7)
ggsave(out_stacked_png, p_stacked, width = 9, height = 7, dpi = 300)

# ---------------------------------------------------------
# 8. GROUPED BAR PLOT
# ---------------------------------------------------------
p_grouped <- ggplot(count_long, aes(x = Class, y = Count, fill = Direction)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  facet_wrap(~ Contrast, ncol = 1) +
  scale_fill_manual(values = direction_colors) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Acute significant lipids per class",
    x = NULL,
    y = "Number of significant lipids",
    fill = "Direction"
  )

print(p_grouped)

ggsave(out_grouped_pdf, p_grouped, width = 9, height = 7)
ggsave(out_grouped_png, p_grouped, width = 9, height = 7, dpi = 300)

# ---------------------------------------------------------
# 9. LOLLIPOP PLOT
# ---------------------------------------------------------
p_lollipop <- ggplot(count_long, aes(x = Class, y = SignedCount, color = Direction)) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4) +
  geom_segment(aes(xend = Class, y = 0, yend = SignedCount), linewidth = 0.8) +
  geom_point(size = 3) +
  facet_wrap(~ Contrast, ncol = 1) +
  scale_color_manual(values = direction_colors) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Acute significant lipids per class",
    x = NULL,
    y = "number of significant lipids",
    color = "Direction"
  )

print(p_lollipop)

ggsave(out_lollipop_pdf, p_lollipop, width = 9, height = 7)
ggsave(out_lollipop_png, p_lollipop, width = 9, height = 7, dpi = 300)

# ---------------------------------------------------------
# 10. CONSOLE SUMMARY
# ---------------------------------------------------------
cat("\nSaved files:\n")
cat(out_counts_csv, "\n")
cat(out_long_csv, "\n")
cat(out_stacked_pdf, "\n")
cat(out_stacked_png, "\n")
cat(out_grouped_pdf, "\n")
cat(out_grouped_png, "\n")
cat(out_lollipop_pdf, "\n")
cat(out_lollipop_png, "\n")