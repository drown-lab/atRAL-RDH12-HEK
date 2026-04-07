#File Needed
DEPresults_v2 <- read.csv("Lipidomics/output_txts/Recovery_DEA_results_v1_corrected.csv")

# =========================================================
# Batch standard volcano plots
# - loops over all contrasts automatically
# - uses ratio vs adjusted p-value
# - no missingness-class annotation
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

favorite_IDs <- c( )

padj_cutoff  <- 0.1
lfc_cutoff   <- log2(1.45)
n_top_labels <- 10

# output folder
out_dir <- "Lipidomics/Figures/Recovery/VolcanoPlots/StandardVolcano/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# optional png output too
save_png <- FALSE

# plot size
plot_width  <- 7
plot_height <- 7

# whether labels should include rounded values
include_values_in_labels <- FALSE

# =========================================================
# 2. HELPER FUNCTIONS
# =========================================================

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

make_volcano_df_standard <- function(df,
                                     contrast,
                                     padj_cutoff = 0.1,
                                     lfc_cutoff = log2(1.5),
                                     favorite_IDs = character(),
                                     n_top_labels = 10,
                                     include_values_in_labels = FALSE) {
  
  ratio_col <- paste0(contrast, "_ratio")
  padj_col  <- paste0(contrast, "_p.adj")
  
  needed_cols <- c("name", ratio_col, padj_col)
  missing_cols <- setdiff(needed_cols, colnames(df))
  
  if (length(missing_cols) > 0) {
    stop(
      "These required columns are missing for contrast '", contrast, "':\n",
      paste(missing_cols, collapse = "\n")
    )
  }
  
  out <- df %>%
    transmute(
      name = name,
      ratio = .data[[ratio_col]],
      padj  = .data[[padj_col]]
    ) %>%
    mutate(
      # keep full precision for plotting
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
      )
    )
  
  top_hits <- out %>%
    filter(!is.na(padj)) %>%
    arrange(padj) %>%
    slice_head(n = n_top_labels) %>%
    pull(name)
  
  out <- out %>%
    mutate(
      label_base = if_else(name %in% c(top_hits, favorite_IDs), name, NA_character_),
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

plot_volcano_standard <- function(volcano_df,
                                  contrast,
                                  padj_cutoff = 0.1,
                                  lfc_cutoff = log2(1.5)) {
  
  ggplot(volcano_df, aes(x = ratio, y = -log10(padj_plot))) +
    geom_point(
      aes(color = direction),
      alpha = 0.75,
      size = 2.6
    ) +
    geom_point(
      data = subset(volcano_df, !is.na(label)),
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
      color = "Regulation"
    ) +
    scale_color_manual(
      values = c(
        "Up" = "purple",
        "Down" = "#00B8B8",
        "NS" = "grey70"
      ),
      drop = FALSE
    ) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 12),
      axis.title = element_text(size = 12)
    )
}

# =========================================================
# 5. BATCH SAVE FUNCTION
# =========================================================

save_all_standard_volcanoes <- function(df,
                                        out_dir = "Lipidomics/Figures/Recovery/VolcanoPlots/StandardVolcano/",
                                        favorite_IDs = character(),
                                        padj_cutoff = 0.1,
                                        lfc_cutoff = log2(1.3),
                                        n_top_labels = 10,
                                        width = 12,
                                        height = 10,
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
      volcano_df <- make_volcano_df_standard(
        df = df,
        contrast = contrast,
        padj_cutoff = padj_cutoff,
        lfc_cutoff = lfc_cutoff,
        favorite_IDs = favorite_IDs,
        n_top_labels = n_top_labels,
        include_values_in_labels = include_values_in_labels
      )
      
      p <- plot_volcano_standard(
        volcano_df = volcano_df,
        contrast = contrast,
        padj_cutoff = padj_cutoff,
        lfc_cutoff = lfc_cutoff
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
            name,
            ratio,
            padj,
            ratio_round,
            padj_round,
            ratio_display,
            padj_display,
            direction,
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
  include_values_in_labels = include_values_in_labels
)

print(volcano_export_summary)

write.csv(
  volcano_export_summary,
  file = file.path(out_dir, "volcano_export_summary.csv"),
  row.names = FALSE
)