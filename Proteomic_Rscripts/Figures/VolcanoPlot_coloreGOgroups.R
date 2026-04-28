# =========================================================
# Batch volcano plots with GO-term/pathway highlighting
# - loops over all contrasts automatically
# - no grid lines
# - significant proteins (padj <= 0.05 and abs(FC) >= 2) are black
# - if a significant protein belongs to one of the selected GO/pathway groups,
#   it is color-coded by that group instead of black
# - only selected genes are labeled
# - saves SVG + optional PNG + per-contrast table
# =========================================================

# -------------------------
# Packages
# -------------------------
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(tibble)

# =========================================================
# 1. INPUT
# =========================================================

DEPresults_v2 <- read.csv(
  "Proteomic_output_txts/data_results_w_missingclass_BH_readjusted.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
colnames(DEPresults_v2)

# =========================================================
# 2. USER SETTINGS
# =========================================================

padj_cutoff <- 0.05
fc_cutoff   <- 2
lfc_cutoff  <- log2(fc_cutoff)

# Only these genes get text labels
label_genes <- c("ALG1", "STT3A", "ZDHHC13", "ZZDHHC13", "HSPA6", "GPAT4", "EIF3C", "JUNB", "TUBA4A")

# Output folder
out_dir <- "Proteomic_Figs/VolcanoPlots/GOterm_highlighted_blackSig/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Optional PNG output
save_png <- FALSE

# Plot size
plot_width  <- 8
plot_height <- 7

# Leave NULL to auto-detect all contrasts
selected_contrasts <- NULL

# =========================================================
# 3. PATHWAY / GO-TERM GENE SETS
# =========================================================

# ---- CCT/TRiC complex / actin-tubulin folding related ----
CCT_TRiC <- unique(c(
  "CCT5","CCT4","CCT3","CCT8","CCT2","CCT6A","CCT7","TCP1",
  "TUBB4B","TUBA1C","TUBB2A","TUBAL3","TUBB4A","TUBB6","TUBA4A",
  "GNB2","GNB4","GNA11"
))

# ---- N-glycan biosynthesis / N-linked glycosylation ----
N_Glycans <- unique(c(
  "ALG3","DPAGT1","ALG1","ALG2","ALG11","ALG5","ALG8","ALG9",
  "RPN2","KRTCAP2","STT3B","MAN1A1","RPN1","DDOST","MGAT4B",
  "DPM1","STT3A","PGM3","RFT1","GMPPA","FPGT","ST3GAL6","HK1"
))

# ---- SLC transporter / ion-metabolite transport ----
SLC_transport <- unique(c(
  "SLC39A7","SLC30A5","SLC9A6","SLC31A1","SLC39A14","SLC12A2",
  "SLC33A1","SLC25A10","SLC25A1","SLC7A8","SLC6A9","SLC25A11",
  "SLC38A2","SLC15A4","SLC16A1","SLC4A7","SLC29A1","SLC7A1",
  "SLC7A6","SLC6A15","SLC20A1","SLC2A1","SLC7A5","SLC7A11"
))

# ---- Lipid metabolic processing / lipid synthesis ----
Lipid_synthesis <- unique(c(
  "GGCX","PRKAG1","KBTBD2","HACD3","TYSND1","OCRL","DDX20","DPM1",
  "BLVRA","ECHDC1","DHCR24","TM9SF2","CHKA","INPP1","ERLIN1","SACM1L",
  "ACSS3","MAPK1","AGK","MAPK14","FASN","LTA4H","SGPL1","HPS6","GPAT4",
  "WDTC1","DDHD1","EPHX1","CLN6","AGPAT4","FDPS","DDHD2","MGST2",
  "EFR3B","LDLRAP1","ABHD12","RDH11","ABHD6","RDH14","PCCA","B4GALT4",
  "ACOT9","PTGR3","HSD17B12","SLC30A5","SIRT3","HSD17B11","ZMPSTE24",
  "ATP1A1","NSDHL","ETNK1","SOAT1","INPP4A","PCK2","KPNB1","SPTLC2",
  "BCKDK","CYP51A1","ECHS1","PRKG1","SCD","ACSS2","PECR","MMUT",
  "PLCG1","PLAA","OSBPL3","SESN2","CYB5R3","CSNK1G2","OGT","PIK3C3",
  "ERLIN2","VAC14","NUDT8","GNPAT","ABCA3","POR","HADHA","PC","SLC27A3",
  "ABHD11","ACSF3","TRIB3","KAT5","ABHD16A","LPCAT1","MVK","PTDSS1",
  "NR3C1","EFR3A","GBA2","GPCPD1","DGKE","PDK1","MTMR1","MFSD8",
  "PRKAA2","PDSS1","ACAD9","PHB2","ST3GAL6","PIGU","MGST3","LDAH",
  "CERS2","PLPP3","HSD17B7","IMPA2","CHST10","PTEN","SNX17","PRKACA",
  "SPTLC1","SERINC1","DHCR7","ACSL3","ATP5F1A","ACOX3","TMEM43","SCD5",
  "HDHD5","PIP5K1A","EPHX2","UGCG","AGPAT3","FDFT1","PRKD2","LPGAT1",
  "PI4KB","SCAP","NUDT19","INPP5K","PLPP6","DHRS7","FADS1","SLC16A1",
  "MLYCD","GPR180","PEDS1","PIGG","CREM","CTDNEP1","PRKCE","PTPMT1",
  "MARCHF6","TMEM38B","HDLBP","PI4K2B","CDIPT","GNAI1","DAGLB","PDK2",
  "ACLY"
))

# ---- Translation / ribosome-related ----
Translation_Ribosome <- unique(c(
  "RPL35","RPL35A","RPLP2","RPS8","RPS15A","RPS24","APEH","RPL7A",
  "RPL19","RPL27A","RPS27L","RPL37A","RPS3A","RPS11","RPS18","RPS27",
  "RPL3","RPL12","RPL23A","RPL14","RPL32","RPL10L","RPLP0","RPS6",
  "RPS14","RPS21","RPL6","RPL17","RPL27","RPL13A","RPL36AL","RPS2",
  "RPS9","RPS16","RPS25","FAU","RPL8","RPL21","RPL26L1","RPL38","RPS4X",
  "RPS12","RPS19","RPS28","RPL29","RPL4","RPL13","RPL24","RPL23","RPL34",
  "RPL22L1","RPLP1","RPS7","RPS15","RPS23","RPL7","RPL18","RPL30",
  "RPL36","RPL37","RPS3","RPS10","RPS17","RPS26","RPSA","RPL10","RPL22",
  "TRMT112","RPL36A","RPS5","RPS13","RPS20","RPS29","RPL31","RPL5",
  "RPL15","RPL26", "EIF3C"
))

# Build pathway table
pathway_tbl <- bind_rows(
  tibble(Gene = CCT_TRiC,             Pathway = "CCT/TRiC"),
  tibble(Gene = N_Glycans,            Pathway = "N-glycans"),
  tibble(Gene = SLC_transport,        Pathway = "SLC transport"),
  tibble(Gene = Lipid_synthesis,      Pathway = "Lipid synthesis"),
  tibble(Gene = Translation_Ribosome, Pathway = "Translation/Ribosome")
) %>%
  distinct(Gene, .keep_all = TRUE)

# Optional duplicate check
dup_check <- bind_rows(
  tibble(Gene = CCT_TRiC,             Pathway = "CCT/TRiC"),
  tibble(Gene = N_Glycans,            Pathway = "N-glycans"),
  tibble(Gene = SLC_transport,        Pathway = "SLC transport"),
  tibble(Gene = Lipid_synthesis,      Pathway = "Lipid synthesis"),
  tibble(Gene = Translation_Ribosome, Pathway = "Translation/Ribosome")
) %>%
  count(Gene, sort = TRUE) %>%
  filter(n > 1)

if (nrow(dup_check) > 0) {
  message("Genes found in more than one pathway list; first assignment kept:")
  print(dup_check)
}

# =========================================================
# 4. HELPER FUNCTIONS
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
    x < 0.001 ~ "<0.001",
    x < 0.01  ~ "<0.01",
    TRUE ~ sprintf("%.3f", x)
  )
}

format_ratio_display <- function(x) {
  ifelse(is.na(x), NA_character_, sprintf("%.2f", x))
}

# =========================================================
# 5. PREP DATA FOR ONE VOLCANO
# =========================================================

make_volcano_df_pathway <- function(df,
                                    contrast,
                                    pathway_tbl,
                                    label_genes,
                                    padj_cutoff = 0.05,
                                    lfc_cutoff = log2(2)) {
  
  ratio_col <- paste0(contrast, "_ratio")
  padj_col  <- paste0(contrast, "_p.adj")
  
  needed_cols <- c("Gene", "name", ratio_col, padj_col)
  missing_cols <- setdiff(needed_cols, colnames(df))
  
  if (length(missing_cols) > 0) {
    stop(
      "These required columns are missing for contrast '", contrast, "':\n",
      paste(missing_cols, collapse = "\n")
    )
  }
  
  out <- df %>%
    mutate(
      Gene_use = dplyr::coalesce(na_if(Gene, ""), na_if(name, ""))
    ) %>%
    filter(!is.na(Gene_use), Gene_use != "") %>%
    distinct(Gene_use, .keep_all = TRUE) %>%
    transmute(
      Gene  = Gene_use,
      ratio = .data[[ratio_col]],
      padj  = .data[[padj_col]]
    ) %>%
    mutate(
      padj_plot = pmax(padj, 1e-300),
      sig = !is.na(padj) & padj <= padj_cutoff,
      pass_fc = !is.na(ratio) & abs(ratio) >= lfc_cutoff
    ) %>%
    left_join(pathway_tbl, by = "Gene") %>%
    mutate(
      # significant + FC + in pathway group
      pathway_highlight = sig & pass_fc & !is.na(Pathway),
      
      # significant + FC but not in selected pathway groups
      sig_other = sig & pass_fc & is.na(Pathway),
      
      # plot category drives point color
      plot_group = case_when(
        pathway_highlight & Pathway == "CCT/TRiC"             ~ "CCT/TRiC",
        pathway_highlight & Pathway == "N-glycans"            ~ "N-glycans",
        pathway_highlight & Pathway == "SLC transport"        ~ "SLC transport",
        pathway_highlight & Pathway == "Lipid synthesis"      ~ "Lipid synthesis",
        pathway_highlight & Pathway == "Translation/Ribosome" ~ "Translation/Ribosome",
        sig_other                                              ~ "Significant other",
        TRUE                                                   ~ "NS"
      ),
      
      label = case_when(
        Gene %in% label_genes & sig & pass_fc ~ Gene,
        TRUE ~ NA_character_
      ),
      
      ratio_display = format_ratio_display(ratio),
      padj_display  = format_padj_display(padj)
    )
  
  out
}

# =========================================================
# 6. PLOT FUNCTION
# =========================================================

plot_volcano_pathway <- function(volcano_df,
                                 contrast,
                                 padj_cutoff = 0.05,
                                 lfc_cutoff = log2(2)) {
  
  plot_colors <- c(
    "NS" = "grey78",
    "Significant other" = "grey25",
    "CCT/TRiC" = "#1b9e77",
    "N-glycans" = "#d95f02",
    "SLC transport" = "#7570b3",
    "Lipid synthesis" = "#e7298a",
    "Translation/Ribosome" = "#66a61e"
  )
  
  # Plot in layers so highlighted GO/pathway proteins are always on top
  df_ns <- volcano_df %>%
    filter(plot_group == "NS")
  
  df_sig_other <- volcano_df %>%
    filter(plot_group == "Significant other")
  
  df_pathway <- volcano_df %>%
    filter(plot_group %in% c(
      "CCT/TRiC",
      "N-glycans",
      "SLC transport",
      "Lipid synthesis",
      "Translation/Ribosome"
    ))
  
  df_labels <- volcano_df %>%
    filter(!is.na(label))
  
  ggplot(mapping = aes(x = ratio, y = -log10(padj_plot))) +
    
    # 1. background NS points
    geom_point(
      data = df_ns,
      aes(color = plot_group),
      alpha = 0.45,
      size = 2.0
    ) +
    
    # 2. significant points not in GO groups
    geom_point(
      data = df_sig_other,
      aes(color = plot_group),
      alpha = 0.65,
      size = 2.4
    ) +
    
    # 3. GO/pathway-highlighted points ON TOP
    geom_point(
      data = df_pathway,
      aes(color = plot_group),
      alpha = 1,
      size = 3.2
    ) +
    
    # thresholds
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
    
    # labels only for selected genes
    geom_text_repel(
      data = df_labels,
      aes(label = label),
      size = 3.4,
      color = "black",
      box.padding = 0.35,
      point.padding = 0.25,
      segment.color = "grey30",
      segment.size = 0.3,
      max.overlaps = Inf,
      min.segment.length = 0,
      na.rm = TRUE
    ) +
    
    scale_color_manual(
      values = plot_colors,
      breaks = c(
        "CCT/TRiC",
        "N-glycans",
        "SLC transport",
        "Lipid synthesis",
        "Translation/Ribosome",
        "Significant other",
        "NS"
      ),
      drop = FALSE,
      name = "Protein group"
    ) +
    
    labs(
      title = paste("Volcano:", contrast),
      x = "log2 Fold Change",
      y = "-log10(adjusted p-value)"
    ) +
    
    theme_classic() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
}
# =========================================================
# 7. BATCH SAVE FUNCTION
# =========================================================

save_all_pathway_volcanoes <- function(df,
                                       pathway_tbl,
                                       label_genes,
                                       out_dir,
                                       padj_cutoff = 0.05,
                                       lfc_cutoff = log2(2),
                                       width = 8,
                                       height = 7,
                                       save_png = FALSE,
                                       selected_contrasts = NULL) {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  contrasts <- get_all_contrasts(df)
  
  if (!is.null(selected_contrasts)) {
    contrasts <- intersect(contrasts, selected_contrasts)
  }
  
  if (length(contrasts) == 0) {
    stop("No contrasts found to plot.")
  }
  
  message("Found ", length(contrasts), " contrasts.")
  
  summary_tbl <- vector("list", length(contrasts))
  
  for (i in seq_along(contrasts)) {
    contrast <- contrasts[i]
    message("[", i, "/", length(contrasts), "] Processing: ", contrast)
    
    result <- tryCatch({
      volcano_df <- make_volcano_df_pathway(
        df = df,
        contrast = contrast,
        pathway_tbl = pathway_tbl,
        label_genes = label_genes,
        padj_cutoff = padj_cutoff,
        lfc_cutoff = lfc_cutoff
      )
      
      p <- plot_volcano_pathway(
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
          dpi = 900
        )
      }
      
      table_file <- file.path(out_dir, paste0(base_name, "_table.csv"))
      write.csv(
        volcano_df %>%
          select(
            Gene,
            Pathway,
            ratio,
            padj,
            ratio_display,
            padj_display,
            sig,
            pass_fc,
            pathway_highlight,
            sig_other,
            plot_group,
            label
          ),
        file = table_file,
        row.names = FALSE
      )
      
      data.frame(
        contrast = contrast,
        n_pathway_highlighted = sum(volcano_df$pathway_highlight, na.rm = TRUE),
        n_sig_other = sum(volcano_df$sig_other, na.rm = TRUE),
        n_labeled = sum(!is.na(volcano_df$label)),
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
        n_pathway_highlighted = NA_integer_,
        n_sig_other = NA_integer_,
        n_labeled = NA_integer_,
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
# 8. RUN
# =========================================================

volcano_export_summary <- save_all_pathway_volcanoes(
  df = DEPresults_v2,
  pathway_tbl = pathway_tbl,
  label_genes = label_genes,
  out_dir = out_dir,
  padj_cutoff = padj_cutoff,
  lfc_cutoff = lfc_cutoff,
  width = plot_width,
  height = plot_height,
  save_png = save_png,
  selected_contrasts = selected_contrasts
)

print(volcano_export_summary)

write.csv(
  volcano_export_summary,
  file = file.path(out_dir, "volcano_export_summary.csv"),
  row.names = FALSE
)