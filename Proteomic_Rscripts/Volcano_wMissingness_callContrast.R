# =========================================================
# Volcano plot with missingness-aware point annotation
# Works on wide DEP-style result table like DEPresults_v2
# =========================================================

# -------------------------
# Packages
# -------------------------
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(forcats)

# =========================================================
# 1. USER INPUTS
# =========================================================

# Main table
# Assumes your data frame is already in memory as DEPresults_v2

contrast <- "RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr"

favorite_IDs <- c("RDH12", "RDH11", "AGR2", "POLD4", "SBSN")

padj_cutoff  <- 0.05
lfc_cutoff   <- log2(1.5)
n_top_labels <- 3

# Optional behavior
drop_mnar_vs_mnar <- FALSE     # TRUE if you want to remove these points
deemphasize_mnar_vs_mnar <- TRUE  # FALSE if you want same alpha for all points

# =========================================================
# 2. HELPER FUNCTIONS
# =========================================================

# Parse raw missclass values like:
# Present_3of3, MAR_2of3, MAR_1of3, MNAR_0of3
parse_missclass_type <- function(x) {
  case_when(
    is.na(x) ~ "Unknown",
    str_detect(x, "^Present") ~ "Present",
    str_detect(x, "^MNAR") ~ "MNAR",
    str_detect(x, "^MAR") ~ "MAR",
    TRUE ~ "Other"
  )
}

# Extract the "3of3", "2of3", "1of3", "0of3" part if present
parse_missclass_count <- function(x) {
  out <- str_extract(x, "[0-9]+of[0-9]+")
  if_else(is.na(out), "Unknown", out)
}

# Build a pairwise class from two parsed miss classes
make_missing_pair <- function(miss_1, miss_2) {
  case_when(
    miss_1 == "Present" & miss_2 == "Present" ~ "Present vs Present",
    
    miss_1 == "MNAR" & miss_2 == "Present" ~ "MNAR vs Present",
    miss_1 == "Present" & miss_2 == "MNAR" ~ "Present vs MNAR",
    
    miss_1 == "MAR" & miss_2 == "Present" ~ "MAR vs Present",
    miss_1 == "Present" & miss_2 == "MAR" ~ "Present vs MAR",
    
    miss_1 == "MAR" & miss_2 == "MAR" ~ "MAR vs MAR",
    miss_1 == "MNAR" & miss_2 == "MNAR" ~ "MNAR vs MNAR",
    
    miss_1 == "MNAR" & miss_2 == "MAR" ~ "MNAR vs MAR",
    miss_1 == "MAR" & miss_2 == "MNAR" ~ "MAR vs MNAR",
    
    TRUE ~ "Mixed / Other"
  )
}

# =========================================================
# 3. MAIN PREP FUNCTION
# =========================================================

make_volcano_df <- function(df,
                            contrast,
                            padj_cutoff = 0.05,
                            lfc_cutoff = log2(2),
                            favorite_IDs = character(),
                            n_top_labels = 3) {
  
  # Derive column names from contrast string
  ratio_col <- paste0(contrast, "_ratio")
  padj_col  <- paste0(contrast, "_p.adj")
  
  # Split contrast into the two compared groups
  groups <- str_split(contrast, "_vs_", simplify = TRUE)
  
  if (ncol(groups) != 2) {
    stop("Contrast must look like 'group1_vs_group2'")
  }
  
  grp1 <- groups[1]
  grp2 <- groups[2]
  
  miss1_col <- paste0("missclass_", grp1)
  miss2_col <- paste0("missclass_", grp2)
  
  needed_cols <- c("Gene", ratio_col, padj_col, miss1_col, miss2_col)
  missing_cols <- setdiff(needed_cols, colnames(df))
  
  if (length(missing_cols) > 0) {
    stop(
      "These required columns are missing:\n",
      paste(missing_cols, collapse = "\n")
    )
  }
  
  # Build plotting table
  out <- df %>%
    transmute(
      Gene = Gene,
      ratio = .data[[ratio_col]],
      padj  = .data[[padj_col]],
      miss_1_raw = .data[[miss1_col]],
      miss_2_raw = .data[[miss2_col]]
    ) %>%
    mutate(
      # protect against padj == 0 causing Inf on -log10 scale
      padj_plot = pmax(padj, 1e-300),
      
      sig = !is.na(padj) & padj < padj_cutoff,
      
      direction = case_when(
        !is.na(ratio) & ratio >  lfc_cutoff & sig ~ "Up",
        !is.na(ratio) & ratio < -lfc_cutoff & sig ~ "Down",
        TRUE ~ "NS"
      ),
      
      miss_1 = parse_missclass_type(miss_1_raw),
      miss_2 = parse_missclass_type(miss_2_raw),
      
      miss_1_count = parse_missclass_count(miss_1_raw),
      miss_2_count = parse_missclass_count(miss_2_raw),
      
      miss_pair = make_missing_pair(miss_1, miss_2),
      
      miss_pair = factor(
        miss_pair,
        levels = c(
          "Present vs Present",
          "MNAR vs Present",
          "Present vs MNAR",
          "MAR vs Present",
          "Present vs MAR",
          "MAR vs MAR",
          "MNAR vs MAR",
          "MAR vs MNAR",
          "MNAR vs MNAR",
          "Mixed / Other"
        )
      )
    )
  
  # Label top significant hits + favorites
  top_hits <- out %>%
    filter(!is.na(padj)) %>%
    arrange(padj) %>%
    slice_head(n = n_top_labels) %>%
    pull(Gene)
  
  out <- out %>%
    mutate(
      label = if_else(Gene %in% c(top_hits, favorite_IDs), Gene, NA_character_)
    )
  
  out
}

# =========================================================
# 4. PLOTTING FUNCTION
# =========================================================

plot_volcano_missingness <- function(volcano_df,
                                     contrast,
                                     padj_cutoff = 0.05,
                                     lfc_cutoff = log2(2),
                                     drop_mnar_vs_mnar = FALSE,
                                     deemphasize_mnar_vs_mnar = TRUE) {
  
  plot_df <- volcano_df
  
  if (drop_mnar_vs_mnar) {
    plot_df <- plot_df %>%
      filter(miss_pair != "MNAR vs MNAR")
  }
  
  plot_df <- plot_df %>%
    mutate(
      point_alpha = case_when(
        deemphasize_mnar_vs_mnar & miss_pair == "MNAR vs MNAR" ~ 0.20,
        TRUE ~ 0.75
      )
    )
  
  ggplot(plot_df, aes(x = ratio, y = -log10(padj_plot))) +
    
    # Main points
    geom_point(
      aes(
        color = direction,
        shape = miss_pair,
        alpha = point_alpha
      ),
      size = 2.6
    ) +
    
    # Outline labeled points
    geom_point(
      data = subset(plot_df, !is.na(label)),
      shape = 21,
      size = 3.3,
      stroke = 0.8,
      color = "black",
      fill = NA,
      inherit.aes = FALSE,
      aes(x = ratio, y = -log10(padj_plot))
    ) +
    
    # Fold-change cutoffs
    geom_vline(
      xintercept = c(-lfc_cutoff, lfc_cutoff),
      linetype = "dashed",
      color = "grey40"
    ) +
    
    # Adjusted p-value cutoff
    geom_hline(
      yintercept = -log10(padj_cutoff),
      linetype = "dashed",
      color = "grey40"
    ) +
    
    # Labels
    geom_text_repel(
      aes(label = label),
      max.overlaps = 25,
      size = 3.5,
      color = "black",
      segment.color = "grey20",
      segment.size = 0.3,
      box.padding = 0.4,
      point.padding = 0.3,
      force = 2,
      min.segment.length = 0,
      na.rm = TRUE
    ) +
    
    theme_bw() +
    labs(
      title = paste("Volcano:", contrast),
      x = "log2 Fold Change",
      y = "-log10(adjusted p-value)",
      color = "Regulation",
      shape = "Missingness class"
    ) +
    
    scale_color_manual(
      values = c(
        "Up" = "purple",
        "Down" = "#00B8B8",
        "NS" = "grey70"
      ),
      drop = FALSE
    ) +
    
    scale_shape_manual(
      values = c(
        "Present vs Present" = 16,  # filled circle
        "MNAR vs Present"    = 17,  # triangle
        "Present vs MNAR"    = 17,  # triangle
        "MAR vs Present"     = 16,   # filled
        "Present vs MAR"     = 16,   # filled
        "MAR vs MAR"         = 16,  # filled
        "MNAR vs MAR"        = 1,   # x 
        "MAR vs MNAR"        = 1,   # x
        "MNAR vs MNAR"       = 4,   # open circle
        "Mixed / Other"      = 13
      ),
      drop = FALSE
    ) +
    
    scale_alpha_identity() +
    
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 12),
      axis.title = element_text(size = 12)
    )
}

# =========================================================
# 5. RUN IT
# =========================================================

res_volcano <- make_volcano_df(
  df = DEPresults_v2,
  contrast = contrast,
  padj_cutoff = padj_cutoff,
  lfc_cutoff = lfc_cutoff,
  favorite_IDs = favorite_IDs,
  n_top_labels = n_top_labels
)

p <- plot_volcano_missingness(
  volcano_df = res_volcano,
  contrast = contrast,
  padj_cutoff = padj_cutoff,
  lfc_cutoff = lfc_cutoff,
  drop_mnar_vs_mnar = drop_mnar_vs_mnar,
  deemphasize_mnar_vs_mnar = deemphasize_mnar_vs_mnar
)

print(p)