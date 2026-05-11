# =========================================================
# Recovery lipid class summary:
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
recovery_file <- "Lipidomics/output_txts/Recovery_DEA_results_v1_corrected.csv"

out_dir <- "Lipidomics/Figures/Recovery/Lipid_class_direction_counts/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# significance thresholds
padj_cutoff <- 0.1
fc_cutoff <- 1.3
lfc_cutoff <- log2(fc_cutoff)

# output files
out_counts_csv <- file.path(out_dir, "Recovery_lipid_class_up_down_counts.csv")
out_long_csv   <- file.path(out_dir, "Recovery_lipid_class_up_down_counts_long.csv")

out_stacked_pdf <- file.path(out_dir, "Recovery_lipid_class_up_down_stacked_bar.pdf")
out_stacked_png <- file.path(out_dir, "Recovery_lipid_class_up_down_stacked_bar.png")

out_grouped_pdf <- file.path(out_dir, "Recovery_lipid_class_up_down_grouped_bar.pdf")
out_grouped_png <- file.path(out_dir, "Recovery_lipid_class_up_down_grouped_bar.png")

out_lollipop_pdf <- file.path(out_dir, "Recovery_lipid_class_up_down_lollipop.pdf")
out_lollipop_png <- file.path(out_dir, "Recovery_lipid_class_up_down_lollipop.png")

# ---------------------------------------------------------
# 2. READ DATA
# ---------------------------------------------------------
lipid_df <- read.csv(
  recovery_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------
# 3. DEFINE RECOVERY CONTRASTS
# ---------------------------------------------------------
# Positive direction is always the first group in the contrast name.
# Example:
# RDH12_200_vs_Control_200_ratio > 0 means enriched in RDH12 200
# RDH12_200_vs_Control_200_ratio < 0 means enriched in WT 200

contrast_tbl <- tibble::tribble(
  ~ContrastID,                       ~ContrastLabel,              ~UpLabel,             ~DownLabel,
  "Control_100_vs_Control_Veh",       "WT 100 vs WT Veh",          "WT 100",             "WT Veh",
  "Control_200_vs_Control_Veh",       "WT 200 vs WT Veh",          "WT 200",             "WT Veh",
  "RDH12_100_vs_RDH12_Veh",           "RDH12 100 vs RDH12 Veh",    "RDH12 100",          "RDH12 Veh",
  "RDH12_200_vs_RDH12_Veh",           "RDH12 200 vs RDH12 Veh",    "RDH12 200",          "RDH12 Veh",
  "RDH12_Veh_vs_Control_Veh",         "RDH12 Veh vs WT Veh",       "RDH12 Veh",          "WT Veh",
  "RDH12_100_vs_Control_100",         "RDH12 100 vs WT 100",       "RDH12 100",          "WT 100",
  "RDH12_200_vs_Control_200",         "RDH12 200 vs WT 200",       "RDH12 200",          "WT 200"
)

# Required columns based on contrast table
required_cols <- c(
  "name",
  paste0(contrast_tbl$ContrastID, "_p.val"),
  paste0(contrast_tbl$ContrastID, "_ratio")
)

missing_cols <- setdiff(required_cols, colnames(lipid_df))

if (length(missing_cols) > 0) {
  stop(
    "These required columns are missing:\n",
    paste(missing_cols, collapse = "\n")
  )
}

# ---------------------------------------------------------
# 4. ASSIGN LIPID CLASS
# ---------------------------------------------------------
get_lipid_class <- function(x) {
  x <- str_trim(x)
  
  case_when(
    # Oxidized PCs / truncated oxidized phospholipids
    str_detect(x, regex("^PC\\(|Aze|COOH|Azelaoyl", ignore_case = TRUE)) ~ "Oxidized PC",
    
    # Cholesteryl esters, including [CE 18:1 +NH4]
    str_detect(x, "^\\[?CE\\s") ~ "CE",
    
    # TG, including [TG44:5] NL 20:0 and TG52:2] NL 18:1
    str_detect(x, "^\\[?TG") ~ "TG",
    
    # DG, including [DG 38:1
    str_detect(x, "^\\[?DG\\s") ~ "DG",
    str_detect(x, "^DG O-") ~ "DG ether",
    
    # LPC / PC
    str_detect(x, "^LPC\\s") ~ "LPC",
    str_detect(x, "^PC O-") ~ "PC ether",
    str_detect(x, "^PC P-") ~ "PC plasmalogen",
    str_detect(x, "^PC\\s") ~ "PC",
    
    # LPE / PE
    str_detect(x, "^LPE\\s") ~ "LPE",
    str_detect(x, "^PE O-") ~ "PE ether",
    str_detect(x, "^PE P-") ~ "PE plasmalogen",
    str_detect(x, "^PE\\s") ~ "PE",
    
    # PI / PS / PG
    str_detect(x, "^PI O-") ~ "PI ether",
    str_detect(x, "^PI\\s") ~ "PI",
    str_detect(x, "^PS O-") ~ "PS ether",
    str_detect(x, "^PS\\s") ~ "PS",
    str_detect(x, "^PG\\s") ~ "PG",
    
    # Ceramides, including Cer(d14:1/16:0)
    str_detect(x, "^Cer\\(") ~ "Cer",
    str_detect(x, "^Cer\\s") ~ "Cer",
    str_detect(x, "^dhCer\\(") ~ "dhCer",
    str_detect(x, "^dhCer\\s") ~ "dhCer",
    
    # Sphingomyelins, including SM(34:1
    str_detect(x, "^SM\\(") ~ "SM",
    str_detect(x, "^SM\\s") ~ "SM",
    
    # Acyl-carnitines
    str_detect(x, "^CAR\\s") ~ "CAR",
    
    TRUE ~ "Other"
  )
}

# ---------------------------------------------------------
# 5. FLAG UP / DOWN FOR EACH CONTRAST
# ---------------------------------------------------------
lipid_df2 <- lipid_df %>%
  mutate(
    name = str_trim(name),
    Class = get_lipid_class(name)
  )

count_list <- list()
lipid_direction_list <- list()

for (i in seq_len(nrow(contrast_tbl))) {
  
  contrast_id <- contrast_tbl$ContrastID[i]
  contrast_label <- contrast_tbl$ContrastLabel[i]
  up_label <- contrast_tbl$UpLabel[i]
  down_label <- contrast_tbl$DownLabel[i]
  
  p_col <- paste0(contrast_id, "_p.val")
  ratio_col <- paste0(contrast_id, "_ratio")
  
  tmp <- lipid_df2 %>%
    mutate(
      Contrast = contrast_label,
      UpGroup = up_label,
      DownGroup = down_label,
      Direction = case_when(
        !is.na(.data[[p_col]]) &
          !is.na(.data[[ratio_col]]) &
          .data[[p_col]] <= padj_cutoff &
          .data[[ratio_col]] >= lfc_cutoff ~ "Up",
        
        !is.na(.data[[p_col]]) &
          !is.na(.data[[ratio_col]]) &
          .data[[p_col]] <= padj_cutoff &
          .data[[ratio_col]] <= -lfc_cutoff ~ "Down",
        
        TRUE ~ "NS"
      ),
      EnrichedIn = case_when(
        Direction == "Up" ~ up_label,
        Direction == "Down" ~ down_label,
        TRUE ~ "NS"
      ),
      p_value = .data[[p_col]],
      ratio = .data[[ratio_col]]
    )
  
  lipid_direction_list[[contrast_label]] <- tmp %>%
    select(
      name,
      Class,
      Contrast,
      UpGroup,
      DownGroup,
      Direction,
      EnrichedIn,
      p_value,
      ratio
    )
  
  count_list[[contrast_label]] <- tmp %>%
    filter(Direction %in% c("Up", "Down")) %>%
    count(Class, Direction, EnrichedIn, name = "Count") %>%
    mutate(
      Contrast = contrast_label,
      UpGroup = up_label,
      DownGroup = down_label
    )
}

count_long <- bind_rows(count_list)
lipid_direction_long <- bind_rows(lipid_direction_list)

# ---------------------------------------------------------
# 6. COMPLETE MISSING COMBINATIONS WITH ZERO
# ---------------------------------------------------------
all_classes <- sort(unique(lipid_df2$Class))

count_long <- count_long %>%
  complete(
    Contrast = contrast_tbl$ContrastLabel,
    Class = all_classes,
    Direction = c("Up", "Down"),
    fill = list(Count = 0)
  ) %>%
  left_join(
    contrast_tbl %>%
      select(ContrastLabel, UpLabel, DownLabel),
    by = c("Contrast" = "ContrastLabel")
  ) %>%
  mutate(
    EnrichedIn = case_when(
      Direction == "Up" ~ UpLabel,
      Direction == "Down" ~ DownLabel,
      TRUE ~ NA_character_
    ),
    SignedCount = if_else(Direction == "Down", -Count, Count)
  )

# ---------------------------------------------------------
# 7. ORDER CONTRASTS AND CLASSES
# ---------------------------------------------------------
contrast_order <- contrast_tbl$ContrastLabel

class_order <- count_long %>%
  group_by(Class) %>%
  summarise(Total = sum(Count), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(Class)

count_long <- count_long %>%
  mutate(
    Contrast = factor(Contrast, levels = contrast_order),
    Class = factor(Class, levels = class_order),
    Direction = factor(Direction, levels = c("Up", "Down"))
  )

count_wide <- count_long %>%
  select(Contrast, Class, Direction, Count) %>%
  tidyr::pivot_wider(
    names_from = c(Contrast, Direction),
    values_from = Count
  ) %>%
  arrange(Class)

write.csv(count_wide, out_counts_csv, row.names = FALSE)
write.csv(count_long, out_long_csv, row.names = FALSE)

write.csv(
  lipid_direction_long,
  file.path(out_dir, "Recovery_lipid_direction_calls_all_contrasts.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------
# 8. COLORS
# ---------------------------------------------------------
direction_colors <- c(
  "Up" = "#B40426",
  "Down" = "#3B4CC0"
)

# ---------------------------------------------------------
# 9. STACKED BAR PLOT
# ---------------------------------------------------------
p_stacked <- ggplot(count_long, aes(x = Class, y = Count, fill = Direction)) +
  geom_col() +
  facet_wrap(~ Contrast, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = direction_colors) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Recovery significant lipids per class",
    x = NULL,
    y = "Number of significant lipids",
    fill = "Direction"
  )

print(p_stacked)

ggsave(out_stacked_pdf, p_stacked, width = 10, height = 13)
ggsave(out_stacked_png, p_stacked, width = 10, height = 13, dpi = 300)

# ---------------------------------------------------------
# 10. GROUPED BAR PLOT
# ---------------------------------------------------------
p_grouped <- ggplot(count_long, aes(x = Class, y = Count, fill = Direction)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  facet_wrap(~ Contrast, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = direction_colors) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Recovery significant lipids per class",
    x = NULL,
    y = "Number of significant lipids",
    fill = "Direction"
  )

print(p_grouped)

ggsave(out_grouped_pdf, p_grouped, width = 10, height = 13)
ggsave(out_grouped_png, p_grouped, width = 10, height = 13, dpi = 300)

# ---------------------------------------------------------
# 11. LOLLIPOP PLOT
# ---------------------------------------------------------
p_lollipop <- ggplot(count_long, aes(x = Class, y = SignedCount, color = Direction)) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4) +
  geom_segment(aes(xend = Class, y = 0, yend = SignedCount), linewidth = 0.8) +
  geom_point(size = 3) +
  facet_wrap(~ Contrast, ncol = 1, scales = "free_y") +
  scale_color_manual(values = direction_colors) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Recovery significant lipids per class",
    x = NULL,
    y = "Number of significant lipids",
    color = "Direction"
  )

print(p_lollipop)

ggsave(out_lollipop_pdf, p_lollipop, width = 10, height = 13)
ggsave(out_lollipop_png, p_lollipop, width = 10, height = 13, dpi = 300)

# ---------------------------------------------------------
# 12. OPTIONAL: MORE COMPACT LOLLIPOP FOR RDH12 vs WT ONLY
# ---------------------------------------------------------
rdh12_vs_wt_contrasts <- c(
  "RDH12 Veh vs WT Veh",
  "RDH12 100 vs WT 100",
  "RDH12 200 vs WT 200"
)

count_rdh12_vs_wt <- count_long %>%
  filter(Contrast %in% rdh12_vs_wt_contrasts)

p_lollipop_rdh12_vs_wt <- ggplot(
  count_rdh12_vs_wt,
  aes(x = Class, y = SignedCount, color = Direction)
) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4) +
  geom_segment(aes(xend = Class, y = 0, yend = SignedCount), linewidth = 0.8) +
  geom_point(size = 3) +
  facet_wrap(~ Contrast, ncol = 1, scales = "free_y") +
  scale_color_manual(values = direction_colors) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Recovery RDH12 vs WT significant lipids per class",
    x = NULL,
    y = "Number of significant lipids",
    color = "Direction"
  )

print(p_lollipop_rdh12_vs_wt)

ggsave(
  file.path(out_dir, "Recovery_RDH12_vs_WT_lipid_class_lollipop.pdf"),
  p_lollipop_rdh12_vs_wt,
  width = 10,
  height = 7
)

ggsave(
  file.path(out_dir, "Recovery_RDH12_vs_WT_lipid_class_lollipop.png"),
  p_lollipop_rdh12_vs_wt,
  width = 10,
  height = 7,
  dpi = 300
)

# ---------------------------------------------------------
# 13. CONSOLE SUMMARY
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
cat(file.path(out_dir, "Recovery_RDH12_vs_WT_lipid_class_lollipop.pdf"), "\n")
cat(file.path(out_dir, "Recovery_RDH12_vs_WT_lipid_class_lollipop.png"), "\n")