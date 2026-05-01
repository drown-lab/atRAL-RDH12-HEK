library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(stringr)
library(enrichplot)
library(ggplot2)
library(forcats)
library(readr)
library(scales)
library(paletteer)
library(ReactomePA)
# =========================================================
# 1. USER INPUTS
# =========================================================

up_file <- "Proteomic_output_txts/PANGEA_exports/RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr_up.csv"
all_file <- "Proteomic_output_txts/PANGEA_exports/RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr_all.csv"

contrast_name <- tools::file_path_sans_ext(basename(up_file))

out_dir <- file.path(
  "PathwayAnalysis/output/Reactome",
  contrast_name
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pvalue_cutoff <- 0.2
qvalue_cutoff <- 0.2
min_gs_size <- 10
max_gs_size <- 700

use_manual_filter <- TRUE
min_count_keep <- 3
max_terms_to_plot <- 15
max_padj_keep <- 0.10
terms_to_plot <- 25

terms_to_remove_exact <- character(0)
patterns_to_remove <- character(0)

# =========================================================
# 2. READ FILES
# =========================================================

df_up <- read.csv(up_file, stringsAsFactors = FALSE, check.names = FALSE)
df_all <- read.csv(all_file, stringsAsFactors = FALSE, check.names = FALSE)

required_cols <- c("ID", "Gene")
missing_up <- setdiff(required_cols, colnames(df_up))
missing_all <- setdiff(required_cols, colnames(df_all))

if (length(missing_up) > 0) {
  stop("df_up is missing required columns:\n", paste(missing_up, collapse = "\n"))
}
if (length(missing_all) > 0) {
  stop("df_all is missing required columns:\n", paste(missing_all, collapse = "\n"))
}

# =========================================================
# 3. CLEAN UNIPROT IDS
# =========================================================

clean_uniprot_vector <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_split(";") %>%
    unlist() %>%
    str_trim() %>%
    str_remove("-\\d+$") %>%
    .[. != ""] %>%
    unique()
}

sig_up_uniprot <- clean_uniprot_vector(df_up$ID)
universe_uniprot <- clean_uniprot_vector(df_all$ID)

# =========================================================
# 4. MAP UNIPROT -> ENTREZID / SYMBOL
# =========================================================

id_map <- bitr(
  geneID   = universe_uniprot,
  fromType = "UNIPROT",
  toType   = c("ENTREZID", "SYMBOL"),
  OrgDb    = org.Hs.eg.db,
  drop     = TRUE
)

write.csv(
  id_map,
  file.path(out_dir, "uniprot_to_entrez_symbol_map.csv"),
  row.names = FALSE
)

sig_up_entrez <- id_map %>%
  filter(UNIPROT %in% sig_up_uniprot) %>%
  pull(ENTREZID) %>%
  unique()

universe_entrez <- id_map %>%
  pull(ENTREZID) %>%
  unique()

if (length(sig_up_entrez) == 0) {
  stop("No significant up UniProt IDs mapped to ENTREZID.")
}

if (length(universe_entrez) == 0) {
  stop("No universe UniProt IDs mapped to ENTREZID.")
}

# =========================================================
# 5. RUN REACTOME PATHWAY ENRICHMENT
# =========================================================

er_reactome <- enrichPathway(
  gene          = sig_up_entrez,
  universe      = universe_entrez,
  organism      = "human",
  pvalueCutoff  = pvalue_cutoff,
  pAdjustMethod = "BH",
  qvalueCutoff  = qvalue_cutoff,
  minGSSize     = min_gs_size,
  maxGSSize     = max_gs_size,
  readable      = TRUE
)

if (is.null(er_reactome) || nrow(as.data.frame(er_reactome)) == 0) {
  cat("No Reactome pathways were enriched.\n")
} else {
  
  reactome_df <- as.data.frame(er_reactome)
  
  write.csv(
    reactome_df,
    file.path(out_dir, "Reactome_raw.csv"),
    row.names = FALSE
  )
  
  # =======================================================
  # 6. OPTIONAL MANUAL FILTER / CURATION
  # =======================================================
  
  if (use_manual_filter) {
    reactome_manual_df <- reactome_df %>%
      mutate(
        Description_lower = str_to_lower(Description),
        GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
      ) %>%
      filter(Count >= min_count_keep) %>%
      filter(p.adjust <= max_padj_keep) %>%
      filter(!Description %in% terms_to_remove_exact)
    
    if (length(patterns_to_remove) > 0) {
      combined_pattern <- paste(patterns_to_remove, collapse = "|")
      reactome_manual_df <- reactome_manual_df %>%
        filter(!str_detect(Description_lower, regex(combined_pattern, ignore_case = TRUE)))
    }
    
    reactome_manual_df <- reactome_manual_df %>%
      arrange(p.adjust, desc(Count), desc(GeneRatio_num)) %>%
      slice_head(n = max_terms_to_plot) %>%
      mutate(
        Description = fct_reorder(Description, GeneRatio_num)
      )
    
    removed_terms_df <- reactome_df %>%
      mutate(
        Description_lower = str_to_lower(Description),
        removed_exact = Description %in% terms_to_remove_exact,
        removed_pattern = if (length(patterns_to_remove) > 0) {
          str_detect(Description_lower, regex(paste(patterns_to_remove, collapse = "|"), ignore_case = TRUE))
        } else {
          FALSE
        },
        removed_low_count = Count < min_count_keep,
        removed_high_padj = p.adjust > max_padj_keep
      ) %>%
      filter(removed_exact | removed_pattern | removed_low_count | removed_high_padj)
    
    write.csv(
      reactome_manual_df,
      file.path(out_dir, "Reactome_manual_filtered.csv"),
      row.names = FALSE
    )
    
    write.csv(
      removed_terms_df,
      file.path(out_dir, "Reactome_removed_terms.csv"),
      row.names = FALSE
    )
  } else {
    reactome_manual_df <- reactome_df %>%
      mutate(
        GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
        Description = fct_reorder(Description, GeneRatio_num)
      ) %>%
      arrange(p.adjust, desc(Count), desc(GeneRatio_num)) %>%
      slice_head(n = max_terms_to_plot)
    
    write.csv(
      reactome_manual_df,
      file.path(out_dir, "Reactome_manual_filtered.csv"),
      row.names = FALSE
    )
  }
  
  # =======================================================
  # 7. PLOTS
  # =======================================================
  
  p_dot_reactome <- dotplot(er_reactome, showCategory = terms_to_plot) +
    ggtitle(paste0(contrast_name, ": Reactome pathways")) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  print(p_dot_reactome)
  
  ggsave(
    filename = file.path(out_dir, "Reactome_dotplot.pdf"),
    plot = p_dot_reactome,
    width = 10,
    height = 10
  )
  
  ggsave(
    filename = file.path(out_dir, "Reactome_dotplot.png"),
    plot = p_dot_reactome,
    width = 10,
    height = 10,
    dpi = 300
  )
  
  if (nrow(reactome_manual_df) > 0) {
    
    p_dot_reactome_manual <- ggplot(
      reactome_manual_df,
      aes(x = GeneRatio_num, y = Description, size = Count, color = p.adjust)
    ) +
      geom_point() +
      theme_bw() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 12)
      ) +
      labs(
        title = paste0(contrast_name, ": Reactome pathways (manual filtered)"),
        x = "GeneRatio",
        y = NULL,
        color = "p.adjust"
      ) +
      scale_x_continuous(labels = percent_format(accuracy = 1)) +
      scale_color_paletteer_c("grDevices::Plasma", direction = -1)
    
    print(p_dot_reactome_manual)
    
    ggsave(
      filename = file.path(out_dir, "Reactome_manual_filtered_dotplot.pdf"),
      plot = p_dot_reactome_manual,
      width = 10,
      height = 10
    )
    
    ggsave(
      filename = file.path(out_dir, "Reactome_manual_filtered_dotplot.png"),
      plot = p_dot_reactome_manual,
      width = 10,
      height = 10,
      dpi = 300
    )
  } else {
    cat("No Reactome pathways remained after manual filtering.\n")
  }
  
  cat("\nTop Reactome pathways:\n")
  print(
    reactome_manual_df %>%
      select(ID, Description, GeneRatio, Count, p.adjust)
  )
}