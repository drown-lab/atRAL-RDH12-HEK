# =========================================================
# Recovery lipidomics standard volcano plots
# - Exports every recovery contrast
# - Labels only colored significant points
# - Uses shared x-axes within biologically matched contrast families
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(stringr)
})

# ---------------------------------------------------------
# 1. USER INPUTS
# ---------------------------------------------------------
recovery_file <- "Lipidomics/output_txts/Recovery_DEA_results_v1_corrected.csv"
out_dir <- "Lipidomics/Figures/Recovery/VolcanoPlots/StandardVolcano/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

favorite_IDs <- c()
focal_label_contrasts <- c(
  "RDH12_Veh_vs_Control_Veh",
  "RDH12_100_vs_Control_100",
  "RDH12_200_vs_Control_200"
)

padj_cutoff <- 0.1
lfc_cutoff <- log2(1.3)
n_top_labels <- 6
save_png <- FALSE
plot_width <- 7
plot_height <- 7
include_values_in_labels <- FALSE

shared_axis_families <- list(
  WT_dose_vs_EtOH = c(
    "Control_100_vs_Control_Veh",
    "Control_200_vs_Control_Veh"
  ),
  RDH12_dose_vs_EtOH = c(
    "RDH12_100_vs_RDH12_Veh",
    "RDH12_200_vs_RDH12_Veh"
  ),
  RDH12_vs_WT_matched_dose = c(
    "RDH12_Veh_vs_Control_Veh",
    "RDH12_100_vs_Control_100",
    "RDH12_200_vs_Control_200"
  )
)

# ---------------------------------------------------------
# 2. HELPER FUNCTIONS
# ---------------------------------------------------------
sanitize_filename <- function(x) {
  x |>
    str_replace_all("[^A-Za-z0-9_\\-]+", "_") |>
    str_replace_all("_+", "_") |>
    str_replace_all("^_|_$", "")
}

get_all_contrasts <- function(df) {
  ratio_cols <- colnames(df)[str_detect(colnames(df), "_ratio$")]
  str_remove(ratio_cols, "_ratio$")
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

is_story_relevant_lipid <- function(x) {
  str_detect(x, "^PC\\s") |
    str_detect(x, "^PC O-") |
    str_detect(x, "^PC P-") |
    str_detect(x, regex("^PC\\(|Aze|COOH|Azelaoyl", ignore_case = TRUE)) |
    str_detect(x, "^\\[?DG\\s")
}

get_symmetric_x_limits <- function(df, contrasts) {
  ratio_cols <- paste0(contrasts, "_ratio")
  missing_cols <- setdiff(ratio_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Missing ratio columns for shared x-axis:\n", paste(missing_cols, collapse = "\n"))
  }
  
  ratio_values <- unlist(df[ratio_cols], use.names = FALSE)
  ratio_values <- ratio_values[is.finite(ratio_values)]
  if (length(ratio_values) == 0) {
    return(NULL)
  }
  
  limit <- max(abs(ratio_values), na.rm = TRUE)
  c(-limit, limit)
}

build_shared_axis_lookup <- function(df, contrasts, shared_axis_families) {
  lookup <- list()
  
  for (family_name in names(shared_axis_families)) {
    family_contrasts <- intersect(shared_axis_families[[family_name]], contrasts)
    if (length(family_contrasts) > 1) {
      family_limits <- get_symmetric_x_limits(df, family_contrasts)
      for (contrast in family_contrasts) {
        lookup[[contrast]] <- family_limits
      }
    }
  }
  
  lookup
}

make_volcano_df_standard <- function(df,
                                     contrast,
                                     padj_cutoff = 0.1,
                                     lfc_cutoff = log2(1.3),
                                     favorite_IDs = character(),
                                     n_top_labels = 10,
                                     include_values_in_labels = FALSE,
                                     focal_label_contrasts = character()) {
  ratio_col <- paste0(contrast, "_ratio")
  padj_col <- paste0(contrast, "_p.val")
  
  needed_cols <- c("name", ratio_col, padj_col)
  missing_cols <- setdiff(needed_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Missing columns for contrast ", contrast, ":\n", paste(missing_cols, collapse = "\n"))
  }
  
  out <- df |>
    transmute(
      name = name,
      ratio = .data[[ratio_col]],
      padj = .data[[padj_col]]
    ) |>
    mutate(
      padj_plot = pmax(padj, 1e-300),
      ratio_round = round(ratio, 1),
      padj_round = round(padj, 2),
      ratio_display = format_ratio_display(ratio),
      padj_display = format_padj_display(padj),
      sig = !is.na(padj) & padj <= padj_cutoff,
      direction = case_when(
        !is.na(ratio) & ratio >= lfc_cutoff & sig ~ "Up",
        !is.na(ratio) & ratio <= -lfc_cutoff & sig ~ "Down",
        TRUE ~ "NS"
      ),
      story_relevant_lipid = is_story_relevant_lipid(name)
    )
  
  top_hits <- out |>
    filter(
      direction %in% c("Up", "Down"),
      !is.na(padj),
      !is.na(ratio)
    ) |>
    arrange(padj) |>
    slice_head(n = n_top_labels) |>
    pull(name)
  
  out |>
    mutate(
      can_label = direction %in% c("Up", "Down") & !is.na(padj) & !is.na(ratio),
      label_base = if_else(
        can_label &
          (
            name %in% c(top_hits, favorite_IDs) |
              (contrast %in% focal_label_contrasts & story_relevant_lipid)
          ),
        name,
        NA_character_
      ),
      label = case_when(
        is.na(label_base) ~ NA_character_,
        include_values_in_labels ~ paste0(
          label_base, "\n",
          "ratio=", ratio_display, ", p=", padj_display
        ),
        TRUE ~ label_base
      )
    )
}

plot_volcano_standard <- function(volcano_df,
                                  contrast,
                                  padj_cutoff = 0.1,
                                  lfc_cutoff = log2(1.3),
                                  x_limits = NULL,
                                  focal_label_contrasts = character()) {
  is_focal_label_plot <- contrast %in% focal_label_contrasts
  label_size <- if (is_focal_label_plot) 3 else 3.6
  label_max_overlaps <- if (is_focal_label_plot) Inf else 25
  
  p <- ggplot(volcano_df, aes(x = ratio, y = -log10(padj_plot))) +
    geom_point(aes(color = direction), alpha = 0.75, size = 2.6) +
    geom_point(
      data = subset(volcano_df, !is.na(label)),
      aes(x = ratio, y = -log10(padj_plot)),
      inherit.aes = FALSE,
      shape = 21,
      size = 3.3,
      stroke = 0.8,
      color = "black",
      fill = NA
    ) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "grey40") +
    geom_text_repel(
      aes(label = label),
      max.overlaps = label_max_overlaps,
      size = label_size,
      color = "black",
      segment.color = "grey20",
      segment.size = 0.3,
      box.padding = 0.4,
      point.padding = 0.3,
      force = 2,
      min.segment.length = 0,
      na.rm = TRUE
    ) +
    scale_color_manual(
      values = c("Up" = "#D73027", "Down" = "#4575B4", "NS" = "grey70"),
      drop = FALSE
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 12),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = paste("Volcano:", contrast),
      x = "log2 fold change",
      y = "-log10(p-value)",
      color = "Regulation"
    )
  
  if (!is.null(x_limits)) {
    p <- p + coord_cartesian(xlim = x_limits)
  }
  
  p
}

save_all_standard_volcanoes <- function(df,
                                        out_dir,
                                        favorite_IDs = character(),
                                        padj_cutoff = 0.1,
                                        lfc_cutoff = log2(1.3),
                                        n_top_labels = 10,
                                        width = 7,
                                        height = 7,
                                        save_png = FALSE,
                                        include_values_in_labels = FALSE,
                                        shared_axis_families = list(),
                                        focal_label_contrasts = character()) {
  contrasts <- get_all_contrasts(df)
  if (length(contrasts) == 0) {
    stop("No contrast columns ending in '_ratio' were found.")
  }
  
  x_limit_lookup <- build_shared_axis_lookup(df, contrasts, shared_axis_families)
  message("Found ", length(contrasts), " contrasts.")
  
  summary_tbl <- vector("list", length(contrasts))
  
  for (i in seq_along(contrasts)) {
    contrast <- contrasts[i]
    message("[", i, "/", length(contrasts), "] Processing: ", contrast)
    
    result <- tryCatch({
      volcano_df <- make_volcano_df_standard(
        df = df,
        contrast = contrast,
        padj_cutoff = padj_cutoff,
        lfc_cutoff = lfc_cutoff,
        favorite_IDs = favorite_IDs,
        n_top_labels = n_top_labels,
        include_values_in_labels = include_values_in_labels,
        focal_label_contrasts = focal_label_contrasts
      )
      
      p <- plot_volcano_standard(
        volcano_df = volcano_df,
        contrast = contrast,
        padj_cutoff = padj_cutoff,
        lfc_cutoff = lfc_cutoff,
        x_limits = x_limit_lookup[[contrast]],
        focal_label_contrasts = focal_label_contrasts
      )
      
      base_name <- sanitize_filename(paste0("volcano_", contrast))
      svg_file <- file.path(out_dir, paste0(base_name, ".svg"))
      table_file <- file.path(out_dir, paste0(base_name, "_table.csv"))
      
      ggsave(svg_file, p, width = width, height = height, units = "in")
      
      if (save_png) {
        ggsave(file.path(out_dir, paste0(base_name, ".png")), p, width = width, height = height, units = "in", dpi = 300)
      }
      
      write.csv(
        volcano_df |>
          select(name, ratio, padj, ratio_round, padj_round, ratio_display, padj_display, direction, story_relevant_lipid, label),
        table_file,
        row.names = FALSE
      )
      
      data.frame(
        contrast = contrast,
        status = "saved",
        shared_x_min = if (!is.null(x_limit_lookup[[contrast]])) x_limit_lookup[[contrast]][1] else NA_real_,
        shared_x_max = if (!is.null(x_limit_lookup[[contrast]])) x_limit_lookup[[contrast]][2] else NA_real_,
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
        shared_x_min = NA_real_,
        shared_x_max = NA_real_,
        output_svg = NA_character_,
        output_table = NA_character_,
        stringsAsFactors = FALSE
      )
    })
    
    summary_tbl[[i]] <- result
  }
  
  bind_rows(summary_tbl)
}

# ---------------------------------------------------------
# 3. RUN BATCH EXPORT
# ---------------------------------------------------------
DEPresults_v2 <- read.csv(recovery_file, check.names = FALSE, stringsAsFactors = FALSE)

volcano_export_summary <- save_all_standard_volcanoes(
  df = DEPresults_v2,
  out_dir = out_dir,
  favorite_IDs = favorite_IDs,
  padj_cutoff = padj_cutoff,
  lfc_cutoff = lfc_cutoff,
  n_top_labels = n_top_labels,
  width = plot_width,
  height = plot_height,
  save_png = save_png,
  include_values_in_labels = include_values_in_labels,
  shared_axis_families = shared_axis_families,
  focal_label_contrasts = focal_label_contrasts
)

print(volcano_export_summary)

write.csv(
  volcano_export_summary,
  file.path(out_dir, "volcano_export_summary.csv"),
  row.names = FALSE
)
