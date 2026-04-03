

#File Needed
DEPresults_v2<- read.csv("Proteomic_output_txts/data_results_w_missingclass_BH_readjusted.csv")

# =========================================================
# Batch volcano plots with missingness-aware point annotation
# - loops over all contrasts automatically
# - uses ratio vs adjusted p-value
# - adds missingness class as point shape
# - saves each plot as SVG
# - also writes a per-contrast table with rounded display values
# =========================================================

# -------------------------
# Packages
# -------------------------
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)

# =========================================================
# 1. USER INPUTS
# =========================================================

# Assumes your main table is already in memory as:
# DEPresults_v2

favorite_IDs <- c("RDH12", "RDH11")

padj_cutoff  <- 0.055
lfc_cutoff   <- log2(1.45)
n_top_labels <- 3

drop_mnar_vs_mnar <- FALSE
deemphasize_mnar_vs_mnar <- TRUE

# output folder
out_dir <- "Proteomic_Figs/VolcanoPlots/wMissingessInfo/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# optional png output too
save_png <- FALSE

# plot size
plot_width  <- 12
plot_height <- 10

# whether labels should include rounded values
include_values_in_labels <- FALSE

# =========================================================
# 2. HELPER FUNCTIONS
# =========================================================

parse_missclass_type <- function(x) {
  case_when(
    is.na(x) ~ "Unknown",
    str_detect(x, "^Present") ~ "Present",
    str_detect(x, "^MNAR") ~ "MNAR",
    str_detect(x, "^MAR") ~ "MAR",
    TRUE ~ "Other"
  )
}

parse_missclass_count <- function(x) {
  out <- str_extract(x, "[0-9]+of[0-9]+")
  if_else(is.na(out), "Unknown", out)
}

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

sanitize_filename <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9_\\-]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "")
}

get_all_contrasts <- function(df) {
  colnames(df) %>%
    .[str_detect(., "_ratio$")] %>%
    str_remove("_ratio$")
}

format_padj_display <- function(x) {
  case_when(
    is.na(x) ~ NA_character_,
    x < 0.01 ~ "<0.01",
    TRUE ~ sprintf("%.2f", round(x, 2))
  )
}

format_ratio_display <- function(x) {
  ifelse(is.na(x), NA_character_, sprintf("%.1f", round(x, 1)))
}

# =========================================================
# 3. MAIN PREP FUNCTION
# =========================================================

make_volcano_df <- function(df,
                            contrast,
                            padj_cutoff = 0.05,
                            lfc_cutoff = log2(1.5),
                            favorite_IDs = character(),
                            n_top_labels = 3,
                            include_values_in_labels = FALSE) {
  
  ratio_col <- paste0(contrast, "_ratio")
  padj_col  <- paste0(contrast, "_p.adj")
  
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
      "These required columns are missing for contrast '", contrast, "':\n",
      paste(missing_cols, collapse = "\n")
    )
  }
  
  out <- df %>%
    transmute(
      Gene = Gene,
      ratio = .data[[ratio_col]],
      padj  = .data[[padj_col]],
      miss_1_raw = .data[[miss1_col]],
      miss_2_raw = .data[[miss2_col]]
    ) %>%
    mutate(
      # keep full precision for actual plotting
      padj_plot = pmax(padj, 1e-300),
      
      # rounded/formatted values for display only
      ratio_round = round(ratio, 1),
      padj_round  = round(padj, 2),
      
      ratio_display = format_ratio_display(ratio),
      padj_display  = format_padj_display(padj),
      
      sig = !is.na(padj) & padj <= padj_cutoff,
      
      direction = case_when(
        !is.na(ratio) & ratio >=  lfc_cutoff & sig ~ "Up",
        !is.na(ratio) & ratio <= -lfc_cutoff & sig ~ "Down",
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
  
  top_hits <- out %>%
    filter(!is.na(padj)) %>%
    arrange(padj) %>%
    slice_head(n = n_top_labels) %>%
    pull(Gene)
  
  out <- out %>%
    mutate(
      label_base = if_else(Gene %in% c(top_hits, favorite_IDs), Gene, NA_character_),
      label = case_when(
        is.na(label_base) ~ NA_character_,
        include_values_in_labels ~ paste0(
          label_base, "\n",
          "ratio=", ratio_display, ", padj=", padj_display
        ),
        TRUE ~ label_base
      )
    )
  
  out
}

# =========================================================
# 4. PLOTTING FUNCTION
# =========================================================

plot_volcano_missingness <- function(volcano_df,
                                     contrast,
                                     padj_cutoff = 0.05,
                                     lfc_cutoff = log2(1.5),
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
    geom_point(
      aes(
        color = direction,
        shape = miss_pair,
        alpha = point_alpha
      ),
      size = 2.6
    ) +
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
    geom_vline(
      xintercept = c(-lfc_cutoff, lfc_cutoff),
      linetype = "dashed",
      color = "grey40"
    ) +
    geom_hline(
      yintercept = -log10(padj_cutoff),
      linetype = "dashed",
      color = "grey40"
    ) +
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
        "Present vs Present" = 16,
        "MNAR vs Present"    = 17,
        "Present vs MNAR"    = 17,
        "MAR vs Present"     = 16,
        "Present vs MAR"     = 16,
        "MAR vs MAR"         = 16,
        "MNAR vs MAR"        = 1,
        "MAR vs MNAR"        = 1,
        "MNAR vs MNAR"       = 4,
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
# 5. BATCH SAVE FUNCTION
# =========================================================

save_all_volcanoes <- function(df,
                               out_dir = "volcano_plots_svg",
                               favorite_IDs = character(),
                               padj_cutoff = 0.05,
                               lfc_cutoff = log2(1.5),
                               n_top_labels = 3,
                               drop_mnar_vs_mnar = FALSE,
                               deemphasize_mnar_vs_mnar = TRUE,
                               width = 12,
                               height = 7,
                               save_png = FALSE,
                               include_values_in_labels = FALSE) {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  contrasts <- get_all_contrasts(df)
  
  if (length(contrasts) == 0) {
    stop("No contrast columns ending in '_ratio' were found.")
  }
  
  message("Found ", length(contrasts), " contrasts.")
  
  summary_tbl <- vector("list", length(contrasts))
  
  for (i in seq_along(contrasts)) {
    contrast <- contrasts[i]
    message("[", i, "/", length(contrasts), "] Processing: ", contrast)
    
    result <- tryCatch({
      volcano_df <- make_volcano_df(
        df = df,
        contrast = contrast,
        padj_cutoff = padj_cutoff,
        lfc_cutoff = lfc_cutoff,
        favorite_IDs = favorite_IDs,
        n_top_labels = n_top_labels,
        include_values_in_labels = include_values_in_labels
      )
      
      p <- plot_volcano_missingness(
        volcano_df = volcano_df,
        contrast = contrast,
        padj_cutoff = padj_cutoff,
        lfc_cutoff = lfc_cutoff,
        drop_mnar_vs_mnar = drop_mnar_vs_mnar,
        deemphasize_mnar_vs_mnar = deemphasize_mnar_vs_mnar
      )
      
      base_name <- sanitize_filename(paste0("volcano_", contrast))
      
      svg_file <- file.path(out_dir, paste0(base_name, ".svg"))
      
      ggsave(
        filename = svg_file,
        plot = p,
        width = width,
        height = height,
        units = "in"
      )
      
      if (save_png) {
        png_file <- file.path(out_dir, paste0(base_name, ".png"))
        ggsave(
          filename = png_file,
          plot = p,
          width = width,
          height = height,
          units = "in",
          dpi = 300
        )
      }
      
      # save rounded display table
      table_file <- file.path(out_dir, paste0(base_name, "_table.csv"))
      
      write.csv(
        volcano_df %>%
          select(
            Gene,
            ratio,
            padj,
            ratio_round,
            padj_round,
            ratio_display,
            padj_display,
            direction,
            miss_1_raw,
            miss_2_raw,
            miss_pair,
            label
          ),
        file = table_file,
        row.names = FALSE
      )
      
      data.frame(
        contrast = contrast,
        status = "saved",
        output_svg = svg_file,
        output_table = table_file,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      message("  Skipped: ", contrast)
      message("  Reason: ", e$message)
      
      data.frame(
        contrast = contrast,
        status = paste("failed:", e$message),
        output_svg = NA_character_,
        output_table = NA_character_,
        stringsAsFactors = FALSE
      )
    })
    
    summary_tbl[[i]] <- result
  }
  
  bind_rows(summary_tbl)
}

# =========================================================
# 6. RUN BATCH EXPORT
# =========================================================

volcano_export_summary <- save_all_volcanoes(
  df = DEPresults_v2,
  out_dir = out_dir,
  favorite_IDs = favorite_IDs,
  padj_cutoff = padj_cutoff,
  lfc_cutoff = lfc_cutoff,
  n_top_labels = n_top_labels,
  drop_mnar_vs_mnar = drop_mnar_vs_mnar,
  deemphasize_mnar_vs_mnar = deemphasize_mnar_vs_mnar,
  width = plot_width,
  height = plot_height,
  save_png = save_png,
  include_values_in_labels = include_values_in_labels
)

print(volcano_export_summary)

write.csv(
  volcano_export_summary,
  file = file.path(out_dir, "volcano_export_summary.csv"),
  row.names = FALSE
)