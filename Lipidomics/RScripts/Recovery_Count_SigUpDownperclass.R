# =========================================================
# Recovery lipid class summary:
# number of significant up/down lipids per class
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(stringr)
  library(tidyr)
})

# ---------------------------------------------------------
# 1. INPUT
# ---------------------------------------------------------
recovery_file <- "Lipidomics/output_txts/Recovery_DEA_results_v1_corrected.csv"
out_dir <- "Lipidomics/Figures/Recovery/Lipid_class_direction_counts/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

padj_cutoff <- 0.1
lfc_cutoff <- log2(1.3)

out_counts_csv <- file.path(out_dir, "Recovery_lipid_class_up_down_counts.csv")
out_long_csv <- file.path(out_dir, "Recovery_lipid_class_up_down_counts_long.csv")
out_calls_csv <- file.path(out_dir, "Recovery_lipid_direction_calls_all_contrasts.csv")
out_main_pdf <- file.path(out_dir, "Recovery_lipid_class_direction_counts_within_genotype_main_text.pdf")
out_main_png <- file.path(out_dir, "Recovery_lipid_class_direction_counts_within_genotype_main_text.png")
out_supp_pdf <- file.path(out_dir, "Recovery_lipid_class_direction_counts_all_contrasts_supplement.pdf")
out_supp_png <- file.path(out_dir, "Recovery_lipid_class_direction_counts_all_contrasts_supplement.png")

# ---------------------------------------------------------
# 2. CONTRASTS
# ---------------------------------------------------------
contrast_tbl <- tibble::tribble(
  ~ContrastID, ~ContrastLabel, ~ComparisonType, ~UpLabel, ~DownLabel,
  "Control_100_vs_Control_Veh", "WT 100 vs WT EtOH", "Within genotype", "WT 100", "WT EtOH",
  "Control_200_vs_Control_Veh", "WT 200 vs WT EtOH", "Within genotype", "WT 200", "WT EtOH",
  "RDH12_100_vs_RDH12_Veh", "RDH12 100 vs RDH12 EtOH", "Within genotype", "RDH12 100", "RDH12 EtOH",
  "RDH12_200_vs_RDH12_Veh", "RDH12 200 vs RDH12 EtOH", "Within genotype", "RDH12 200", "RDH12 EtOH",
  "RDH12_Veh_vs_Control_Veh", "RDH12 EtOH vs WT EtOH", "RDH12 vs WT matched dose", "RDH12 EtOH", "WT EtOH",
  "RDH12_100_vs_Control_100", "RDH12 100 vs WT 100", "RDH12 vs WT matched dose", "RDH12 100", "WT 100",
  "RDH12_200_vs_Control_200", "RDH12 200 vs WT 200", "RDH12 vs WT matched dose", "RDH12 200", "WT 200",
  "Control_200_vs_Control_100", "WT 200 vs WT 100", "Dose-to-dose context", "WT 200", "WT 100",
  "RDH12_200_vs_RDH12_100", "RDH12 200 vs RDH12 100", "Dose-to-dose context", "RDH12 200", "RDH12 100"
)

main_contrasts <- contrast_tbl |>
  filter(ComparisonType == "Within genotype") |>
  pull(ContrastLabel)

# ---------------------------------------------------------
# 3. HELPERS
# ---------------------------------------------------------
get_lipid_class <- function(x) {
  x <- str_trim(x)
  
  case_when(
    str_detect(x, regex("^PC\\(|Aze|COOH|Azelaoyl", ignore_case = TRUE)) ~ "Oxidized PC",
    str_detect(x, "^\\[?CE\\s") ~ "CE",
    str_detect(x, "^\\[?TG") ~ "TG",
    str_detect(x, "^\\[?DG\\s") ~ "DG",
    str_detect(x, "^DG O-") ~ "DG ether",
    str_detect(x, "^LPC\\s") ~ "LPC",
    str_detect(x, "^PC O-") ~ "PC ether",
    str_detect(x, "^PC P-") ~ "PC plasmalogen",
    str_detect(x, "^PC\\s") ~ "PC",
    str_detect(x, "^LPE\\s") ~ "LPE",
    str_detect(x, "^PE O-") ~ "PE ether",
    str_detect(x, "^PE P-") ~ "PE plasmalogen",
    str_detect(x, "^PE\\s") ~ "PE",
    str_detect(x, "^PI O-") ~ "PI ether",
    str_detect(x, "^PI\\s") ~ "PI",
    str_detect(x, "^PS O-") ~ "PS ether",
    str_detect(x, "^PS\\s") ~ "PS",
    str_detect(x, "^PG\\s") ~ "PG",
    str_detect(x, "^Cer\\(") ~ "Cer",
    str_detect(x, "^Cer\\s") ~ "Cer",
    str_detect(x, "^dhCer\\(") ~ "dhCer",
    str_detect(x, "^dhCer\\s") ~ "dhCer",
    str_detect(x, "^SM\\(") ~ "SM",
    str_detect(x, "^SM\\s") ~ "SM",
    str_detect(x, "^CAR\\s") ~ "CAR",
    TRUE ~ "Other"
  )
}

make_direction_calls <- function(lipid_df, contrast_tbl, padj_cutoff, lfc_cutoff) {
  lipid_df2 <- lipid_df |>
    mutate(name = str_trim(name), Class = get_lipid_class(name))
  
  call_list <- vector("list", nrow(contrast_tbl))
  count_list <- vector("list", nrow(contrast_tbl))
  
  for (i in seq_len(nrow(contrast_tbl))) {
    contrast_id <- contrast_tbl$ContrastID[i]
    p_col <- paste0(contrast_id, "_p.val")
    ratio_col <- paste0(contrast_id, "_ratio")
    
    tmp <- lipid_df2 |>
      mutate(
        Contrast = contrast_tbl$ContrastLabel[i],
        ContrastID = contrast_id,
        ComparisonType = contrast_tbl$ComparisonType[i],
        UpGroup = contrast_tbl$UpLabel[i],
        DownGroup = contrast_tbl$DownLabel[i],
        Direction = case_when(
          !is.na(.data[[p_col]]) & !is.na(.data[[ratio_col]]) &
            .data[[p_col]] <= padj_cutoff & .data[[ratio_col]] >= lfc_cutoff ~ "Up",
          !is.na(.data[[p_col]]) & !is.na(.data[[ratio_col]]) &
            .data[[p_col]] <= padj_cutoff & .data[[ratio_col]] <= -lfc_cutoff ~ "Down",
          TRUE ~ "NS"
        ),
        EnrichedIn = case_when(
          Direction == "Up" ~ UpGroup,
          Direction == "Down" ~ DownGroup,
          TRUE ~ "NS"
        ),
        p_value = .data[[p_col]],
        ratio = .data[[ratio_col]]
      )
    
    call_list[[i]] <- tmp |>
      select(name, Class, ContrastID, Contrast, ComparisonType, UpGroup, DownGroup, Direction, EnrichedIn, p_value, ratio)
    
    count_list[[i]] <- tmp |>
      filter(Direction %in% c("Up", "Down")) |>
      count(ContrastID, Contrast, ComparisonType, Class, Direction, EnrichedIn, name = "Count")
  }
  
  list(
    calls = bind_rows(call_list),
    counts = bind_rows(count_list),
    classes = sort(unique(lipid_df2$Class))
  )
}

make_diverging_plot <- function(count_df, title, subtitle = NULL, ncol = 2) {
  plot_df <- count_df |>
    filter(Count > 0) |>
    group_by(Class) |>
    mutate(ClassTotal = sum(Count)) |>
    ungroup()
  
  if (nrow(plot_df) == 0) {
    stop("No significant up/down lipids were found for this plot.")
  }
  
  x_limit <- max(abs(plot_df$SignedCount), na.rm = TRUE) + 0.8
  
  ggplot(plot_df, aes(x = SignedCount, y = fct_reorder(Class, ClassTotal), fill = Direction)) +
    geom_vline(xintercept = 0, color = "grey35", linewidth = 0.45) +
    geom_col(width = 0.62, color = "white", linewidth = 0.25) +
    geom_text(
      aes(
        x = SignedCount + if_else(SignedCount > 0, 0.18, -0.18),
        label = Count,
        hjust = if_else(SignedCount > 0, 0, 1)
      ),
      size = 4,
      color = "grey15"
    ) +
    facet_wrap(~ Contrast, ncol = ncol) +
    scale_fill_manual(values = c("Up" = "#B40426", "Down" = "#3B4CC0")) +
    scale_x_continuous(
      limits = c(-x_limit, x_limit),
      labels = function(x) abs(x),
      expand = expansion(mult = c(0, 0))
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 10),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(color = "black"),
      axis.text.x = element_text(color = "grey20"),
      panel.spacing = unit(0.7, "lines"),
      plot.title = element_text(face = "bold", size = 12, hjust = 0),
      plot.subtitle = element_text(size = 8.5, color = "grey25"),
      plot.margin = margin(6, 16, 6, 6)
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Number of significant lipids",
      y = NULL
    )
}

# ---------------------------------------------------------
# 4. BUILD TABLES AND FIGURES
# ---------------------------------------------------------
lipid_df <- read.csv(recovery_file, check.names = FALSE, stringsAsFactors = FALSE)

required_cols <- c(
  "name",
  paste0(contrast_tbl$ContrastID, "_p.val"),
  paste0(contrast_tbl$ContrastID, "_ratio")
)
missing_cols <- setdiff(required_cols, colnames(lipid_df))
if (length(missing_cols) > 0) {
  stop("Missing required columns:\n", paste(missing_cols, collapse = "\n"))
}

direction_results <- make_direction_calls(lipid_df, contrast_tbl, padj_cutoff, lfc_cutoff)

count_long <- direction_results$counts |>
  select(ContrastID, Class, Direction, Count) |>
  complete(
    ContrastID = contrast_tbl$ContrastID,
    Class = direction_results$classes,
    Direction = c("Up", "Down"),
    fill = list(Count = 0)
  ) |>
  left_join(contrast_tbl, by = "ContrastID") |>
  mutate(
    EnrichedIn = case_when(
      Direction == "Up" ~ UpLabel,
      Direction == "Down" ~ DownLabel,
      TRUE ~ NA_character_
    ),
    SignedCount = if_else(Direction == "Down", -Count, Count),
    Contrast = factor(ContrastLabel, levels = contrast_tbl$ContrastLabel),
    ComparisonType = factor(ComparisonType, levels = c("Within genotype", "RDH12 vs WT matched dose", "Dose-to-dose context")),
    Direction = factor(Direction, levels = c("Up", "Down"))
  ) |>
  select(ContrastID, Contrast, ComparisonType, Class, Direction, EnrichedIn, Count, SignedCount)

count_wide <- count_long |>
  select(Contrast, Class, Direction, Count) |>
  pivot_wider(names_from = c(Contrast, Direction), values_from = Count) |>
  arrange(Class)

write.csv(count_wide, out_counts_csv, row.names = FALSE)
write.csv(count_long, out_long_csv, row.names = FALSE)
write.csv(direction_results$calls, out_calls_csv, row.names = FALSE)

main_count <- count_long |>
  filter(Contrast %in% main_contrasts) |>
  mutate(Contrast = fct_drop(Contrast))

supp_count <- count_long |>
  mutate(Contrast = fct_drop(Contrast))

p_main <- make_diverging_plot(
  main_count,
  title = "Recovery treatment response within genotype",
  subtitle = "WT and RDH12 cells are each compared with their own EtOH vehicle.",
  ncol = 2
)

p_supp <- make_diverging_plot(
  supp_count,
  title = "Recovery significant lipid classes across all exported contrasts",
  subtitle = "Includes within-genotype, RDH12-vs-WT, and dose-to-dose context contrasts.",
  ncol = 3
)

ggsave(out_main_pdf, p_main, width = 8, height = 6, useDingbats = FALSE)
ggsave(out_main_png, p_main, width = 8, height = 6, dpi = 900)
ggsave(out_supp_pdf, p_supp, width = 9.2, height = 7.0, useDingbats = FALSE)
ggsave(out_supp_png, p_supp, width = 9.2, height = 7.0, dpi = 600)

cat("\nSaved files:\n")
cat(out_counts_csv, "\n")
cat(out_long_csv, "\n")
cat(out_calls_csv, "\n")
cat(out_main_pdf, "\n")
cat(out_main_png, "\n")
cat(out_supp_pdf, "\n")
cat(out_supp_png, "\n")
