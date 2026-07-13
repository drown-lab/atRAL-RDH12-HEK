library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(forcats)

# ---------------------------------------------------------
# 1. INPUT
# ---------------------------------------------------------

enrich_df <- read.csv("Proteomic_output_txts/PANGEA_results/Acute_100vsVeh_100uMenrichment_2026-04-06 135232.csv",
                      check.names = FALSE)

outdir <- "Proteomic_Figs/PANGEA_top_terms_plots/Acute_100vsVeh_100uMenrichment/dotplots/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

bh_cutoff <- 0.5
colnames(enrich_df)
# ---------------------------------------------------------
# 2. CLEAN COLUMN TYPES
# ---------------------------------------------------------

enrich_df <- enrich_df %>%
  mutate(
    `Fold Enrichment` = as.numeric(`log2 fold`),
    `Benjamini & Hochberg` = as.numeric(`Benjamini & Hochberg`)
  )

# ---------------------------------------------------------
# 3. HELPER FUNCTION TO MAKE PLOTS
# ---------------------------------------------------------

make_enrichment_plot <- function(df, n_terms, title_text, out_file) {
  
  if (nrow(df) == 0) {
    message("No rows found for: ", title_text)
    return(NULL)
  }
  
  plot_df <- df %>%
    filter(!is.na(`Fold Enrichment`),
           !is.na(`Benjamini & Hochberg`),
           `Benjamini & Hochberg` <= bh_cutoff) %>%
    arrange(desc(`Fold Enrichment`), `Benjamini & Hochberg`) %>%
    distinct(`Gene Set Name`, .keep_all = TRUE) %>%
    slice_head(n = n_terms) %>%
    mutate(
      `Gene Set Name` = fct_reorder(`Gene Set Name`, `Fold Enrichment`)
    )
  
  if (nrow(plot_df) == 0) {
    message("No significant rows after BH filtering for: ", title_text)
    return(NULL)
  }
  
  library(viridis)
  
  p <- ggplot(plot_df,
         aes(x = `Fold Enrichment`,
             y = fct_reorder(`Gene Set Name`, `Fold Enrichment`),
             size = `Count Overlap Gene`,
             color = -log10(`Benjamini & Hochberg`))) +
    geom_point() +
    scale_color_viridis_c(option = "cividis") +
    labs(
      x = "Fold Enrichment",
      y = NULL,
      color = "-log10(BH)",
      size = "Gene Count"
    ) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(size = 14)
    )
  print(p)
  
  ggsave(
    filename = file.path(outdir, out_file),
    plot = p,
    width = 10,
    height = 7
  )
  
  return(plot_df)
}

# ---------------------------------------------------------
# 4. DEFINE CATEGORIES
# ---------------------------------------------------------

pathway_categories <- c(
  "REACTOME pathway",
  "PANTHER Pathway H.sap",
  "KEGG Pathway H.sap"
)

bp_categories <- c(
  "GO Biological Processes",
  "Direct GO Biological Processes",
  "SLIM2 GO BP"
)

cc_categories <- c(
  "GO Cellular Component",
  "Direct GO Cellular Component",
  "SLIM2 GO CC"
  
)

# ---------------------------------------------------------
# 5. SUBSET TABLES
# ---------------------------------------------------------

pathway_df <- enrich_df %>%
  filter(`Gene Set Category` %in% pathway_categories)

bp_df <- enrich_df %>%
  filter(`Gene Set Category` %in% bp_categories)

cc_df <- enrich_df %>%
  filter(`Gene Set Category` %in% cc_categories)

# ---------------------------------------------------------
# 6. MAKE PLOTS
# ---------------------------------------------------------

top_pathways <- make_enrichment_plot(
  df = pathway_df,
  n_terms = 20,
  title_text = "Top 20 Enriched Pathways (BH <= 0.1)",
  out_file = "top20_pathways_BH0.1.pdf"

)

top_bp <- make_enrichment_plot(
  df = bp_df,
  n_terms = 20,
  title_text = "Top 20 Enriched Biological Processes (BH <= 0.1)",
  out_file = "top20_biological_processes_BH0.1.pdf"

)

top_cc <- make_enrichment_plot(
  df = cc_df,
  n_terms = 10,
  title_text = "Top 10 Enriched Cellular Components (BH <= 0.1)",
  out_file = "top10_cellular_components_BH0.1.pdf"
 
)

# ---------------------------------------------------------
# 7. WRITE TABLES USED FOR EACH PLOT
# ---------------------------------------------------------

if (!is.null(top_pathways)) {
  write.csv(top_pathways,
            file.path(outdir, "top20_pathways_BH0.1.csv"),
            row.names = FALSE)
}

if (!is.null(top_bp)) {
  write.csv(top_bp,
            file.path(outdir, "top20_biological_processes_BH0.1.csv"),
            row.names = FALSE)
}

if (!is.null(top_cc)) {
  write.csv(top_cc,
            file.path(outdir, "top10_cellular_components_BH0.1.csv"),
            row.names = FALSE)
}
