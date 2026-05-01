library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(stringr)
library(enrichplot)
library(ggplot2)
library(forcats)
library(GOSemSim)

# =========================================================
# 1. READ FILES
# =========================================================

df_up <- read.csv(
  "Cell_Perturbed/Paper1/output_txts/PANGEA_exports/redo_BHadj_sepchnl/X6h_MG132_vs_X6h_UT_Heavy_up.csv",
  stringsAsFactors = FALSE
)

df_all <- read.csv(
  "Cell_Perturbed/Paper1/output_txts/PANGEA_exports/redo_BHadj_sepchnl/X6h_MG132_vs_X6h_UT_Heavy_all.csv",
  stringsAsFactors = FALSE
)

out_dir <- "Cell_Perturbed/Paper1/output_txts/clusterProfiler/X6h_MG132_vs_X6h_UT_Heavy_up/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# 2. CLEAN UNIPROT IDS
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

cat("Number of unique sig_up UniProt IDs:", length(sig_up_uniprot), "\n")
cat("Number of unique universe UniProt IDs:", length(universe_uniprot), "\n")

# =========================================================
# 3. MAP UNIPROT -> ENTREZID / SYMBOL
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

cat("Mapped sig_up Entrez IDs:", length(sig_up_entrez), "\n")
cat("Mapped universe Entrez IDs:", length(universe_entrez), "\n")

if (length(sig_up_entrez) == 0) {
  stop("No significant up UniProt IDs mapped to ENTREZID.")
}

if (length(universe_entrez) == 0) {
  stop("No universe UniProt IDs mapped to ENTREZID.")
}

# =========================================================
# 4. RUN enrichGO
# =========================================================

ego_bp <- enrichGO(
  gene          = sig_up_entrez,
  universe      = universe_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.1,
  qvalueCutoff  = 0.2,
  minGSSize     = 10,
  maxGSSize     = 500,
  readable      = TRUE
)

if (is.null(ego_bp) || nrow(as.data.frame(ego_bp)) == 0) {
  cat("No GO BP terms were enriched.\n")
} else {
  
  # =======================================================
  # 5. SIMPLIFY
  # =======================================================
  
  hsGO <- godata(annoDb = "org.Hs.eg.db", ont = "BP")
  
  ego_bp_s <- simplify(
    ego_bp,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min,
    semData = hsGO
  )
  
  ego_bp_df <- as.data.frame(ego_bp)
  ego_bp_s_df <- as.data.frame(ego_bp_s)
  
  write.csv(
    ego_bp_df,
    file.path(out_dir, "enrichGO_BP_raw.csv"),
    row.names = FALSE
  )
  
  write.csv(
    ego_bp_s_df,
    file.path(out_dir, "enrichGO_BP_simplified.csv"),
    row.names = FALSE
  )
  
  # =======================================================
  # 6. OPTIONAL: GO LEVEL FILTER
  # =======================================================
  
  use_go_level_filter <- FALSE
  go_level_to_keep <- 3
  
  ego_bp_for_manual <- ego_bp_s
  
  if (use_go_level_filter) {
    ego_bp_for_manual <- gofilter(ego_bp_for_manual, level = go_level_to_keep)
  }
  
  ego_bp_for_manual_df <- as.data.frame(ego_bp_for_manual)
  
  # =======================================================
  # 7. MANUAL GO FILTER / CURATION
  # =======================================================
  
  min_count_keep <- 3
  max_terms_to_plot <- 20
  max_padj_keep <- 0.10
  
  terms_to_remove_exact <- c(
    "heart looping",
    "embryonic heart tube development",
    "determination of heart left/right asymmetry",
    "left/right pattern formation",
    "determination of left/right symmetry",
    "tube closure",
    "neural tube closure",
    "camera-type eye morphogenesis",
    "muscle cell apoptotic process",
    "regulation of muscle cell apoptotic process"
  )
  
  patterns_to_remove <- c(
    "heart",
    "embryonic",
    "left/right",
    "pattern formation",
    "smoothened",
    "tube closure",
    "eye morphogenesis",
    "muscle"
  )
  
  ego_bp_manual_df <- ego_bp_for_manual_df %>%
    mutate(
      Description_lower = str_to_lower(Description),
      GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
    ) %>%
    filter(Count >= min_count_keep) %>%
    filter(p.adjust <= max_padj_keep) %>%
    filter(!Description %in% terms_to_remove_exact)
  
  if (length(patterns_to_remove) > 0) {
    combined_pattern <- paste(patterns_to_remove, collapse = "|")
    ego_bp_manual_df <- ego_bp_manual_df %>%
      filter(!str_detect(Description_lower, regex(combined_pattern, ignore_case = TRUE)))
  }
  
  ego_bp_manual_df <- ego_bp_manual_df %>%
    arrange(p.adjust, desc(Count), desc(GeneRatio_num)) %>%
    slice_head(n = max_terms_to_plot) %>%
    mutate(
      Description = fct_reorder(Description, GeneRatio_num)
    )
  
  write.csv(
    ego_bp_manual_df,
    file.path(out_dir, "enrichGO_BP_manual_filtered.csv"),
    row.names = FALSE
  )
  
  removed_terms_df <- ego_bp_for_manual_df %>%
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
    removed_terms_df,
    file.path(out_dir, "enrichGO_BP_removed_terms.csv"),
    row.names = FALSE
  )
  
  # =======================================================
  # 8. PLOTS
  # =======================================================
  
  p_dot_simplified <- dotplot(ego_bp_s, showCategory = 20) +
    ggtitle("X6h MG132 vs X6h UT Heavy up: GO BP (simplified)")
  
  print(p_dot_simplified)
  
  ggsave(
    filename = file.path(out_dir, "enrichGO_BP_simplified_dotplot.pdf"),
    plot = p_dot_simplified,
    width = 10,
    height = 7
  )
  
  ggsave(
    filename = file.path(out_dir, "enrichGO_BP_simplified_dotplot.png"),
    plot = p_dot_simplified,
    width = 10,
    height = 7,
    dpi = 300
  )
  
  if (nrow(ego_bp_manual_df) > 0) {
    
    p_dot_manual <- ggplot(
      ego_bp_manual_df,
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
        title = "X6h MG132 vs X6h UT Heavy up: GO BP (manual filtered)",
        x = "GeneRatio",
        y = NULL,
        color = "p.adjust"
      ) +
      scale_x_continuous(labels = scales::percent_format(accuracy = 1))
    
    print(p_dot_manual)
    
    ggsave(
      filename = file.path(out_dir, "enrichGO_BP_manual_filtered_dotplot.pdf"),
      plot = p_dot_manual,
      width = 10,
      height = 7
    )
    
    ggsave(
      filename = file.path(out_dir, "enrichGO_BP_manual_filtered_dotplot.png"),
      plot = p_dot_manual,
      width = 10,
      height = 7,
      dpi = 300
    )
  } else {
    cat("No GO terms remained after manual filtering.\n")
  }
  
  # =======================================================
  # 9. CONSOLE OUTPUT
  # =======================================================
  
  cat("\nTop simplified GO BP terms:\n")
  print(ego_bp_s_df %>% select(ID, Description, GeneRatio, Count, p.adjust) %>% head(20))
  
  cat("\nTop manual-filtered GO BP terms:\n")
  print(ego_bp_manual_df %>% select(ID, Description, GeneRatio, Count, p.adjust))
}