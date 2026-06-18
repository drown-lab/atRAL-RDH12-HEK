library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(stringr)
library(enrichplot)
library(ggplot2)
library(forcats)
library(GOSemSim)
library(readr)
library(scales)
library(paletteer)

# =========================================================
# 1. USER INPUTS
# =========================================================

up_file  <- "Proteomic_Figs/Overlap_100uM_vs_200uM/shared_down_100uM_200uM.csv"
all_file <- "Proteomic_output_txts/PANGEA_exports/RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr_all.csv"

contrast_name <- "shared_down_100uM_200uM"

out_dir <- file.path(
  "GOanalysis/output/clusterProfiler/OverlapSharedproteins/",
  contrast_name
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Enrichment settings
pvalue_cutoff <- 0.3
qvalue_cutoff <- 0.3
min_gs_size <- 3
max_gs_size <- 700

# Simplify settings
simplify_cutoff <- 0.9

# Optional manual filtering
use_manual_filter <- TRUE
min_count_keep <- 3
max_terms_to_plot <- 10
max_padj_keep <- 0.20
terms_to_plot <- 10

terms_to_remove_exact <- character(0)
patterns_to_remove <- character(0)

use_go_level_filter <- FALSE
go_level_to_keep <- 4

# =========================================================
# 2. READ FILES
# =========================================================

df_up <- read.csv(up_file, stringsAsFactors = FALSE, check.names = FALSE)
df_all <- read.csv(all_file, stringsAsFactors = FALSE, check.names = FALSE)

cat("Read up file:", up_file, "\n")
cat("Read all file:", all_file, "\n")
cat("df_up rows:", nrow(df_up), "\n")
cat("df_all rows:", nrow(df_all), "\n")

# Check expected columns
if (!"Protein" %in% colnames(df_up)) {
  stop("df_up must contain a column named 'Protein'")
}

if (!"Gene" %in% colnames(df_all)) {
  stop("df_all must contain a column named 'Gene'")
}

# =========================================================
# 3. CLEAN GENE SYMBOLS
# =========================================================

clean_symbol_vector <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_split(";") %>%
    unlist() %>%
    str_trim() %>%
    .[. != ""] %>%
    unique()
}

sig_up_symbols <- clean_symbol_vector(df_up$Protein)
universe_symbols <- clean_symbol_vector(df_all$Gene)

cat("Number of unique sig_up symbols:", length(sig_up_symbols), "\n")
cat("Number of unique universe symbols:", length(universe_symbols), "\n")

# Keep only sig genes that are also in universe
sig_up_symbols <- intersect(sig_up_symbols, universe_symbols)

cat("Number of sig_up symbols present in universe:", length(sig_up_symbols), "\n")

if (length(sig_up_symbols) == 0) {
  stop("No significant symbols from up_file were found in the universe file.")
}

# =========================================================
# 4. MAP SYMBOL -> ENTREZID
# =========================================================

id_map <- bitr(
  geneID   = universe_symbols,
  fromType = "SYMBOL",
  toType   = c("ENTREZID", "SYMBOL"),
  OrgDb    = org.Hs.eg.db,
  drop     = TRUE
)

write.csv(
  id_map,
  file.path(out_dir, "symbol_to_entrez_map.csv"),
  row.names = FALSE
)

sig_up_entrez <- id_map %>%
  filter(SYMBOL %in% sig_up_symbols) %>%
  pull(ENTREZID) %>%
  unique()

universe_entrez <- id_map %>%
  pull(ENTREZID) %>%
  unique()

cat("Mapped sig_up Entrez IDs:", length(sig_up_entrez), "\n")
cat("Mapped universe Entrez IDs:", length(universe_entrez), "\n")

if (length(sig_up_entrez) == 0) {
  stop("No significant symbols mapped to ENTREZID.")
}

if (length(universe_entrez) == 0) {
  stop("No universe symbols mapped to ENTREZID.")
}

mapped_symbols <- unique(id_map$SYMBOL)

write.csv(
  data.frame(SYMBOL = setdiff(sig_up_symbols, mapped_symbols)),
  file.path(out_dir, "unmapped_sig_up_symbols.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(SYMBOL = setdiff(universe_symbols, mapped_symbols)),
  file.path(out_dir, "unmapped_universe_symbols.csv"),
  row.names = FALSE
)

# =========================================================
# 5. RUN enrichGO
# =========================================================

ego_bp <- enrichGO(
  gene          = sig_up_entrez,
  universe      = universe_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = pvalue_cutoff,
  qvalueCutoff  = qvalue_cutoff,
  minGSSize     = min_gs_size,
  maxGSSize     = max_gs_size,
  readable      = TRUE
)

if (is.null(ego_bp) || nrow(as.data.frame(ego_bp)) == 0) {
  cat("No GO BP terms were enriched.\n")
} else {
  
  # =======================================================
  # 6. SIMPLIFY
  # =======================================================
  
  hsGO <- godata(annoDb = "org.Hs.eg.db", ont = "BP")
  
  ego_bp_s <- simplify(
    ego_bp,
    cutoff = simplify_cutoff,
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
  # 7. OPTIONAL GO LEVEL FILTER
  # =======================================================
  
  ego_bp_for_manual <- ego_bp_s
  
  if (use_go_level_filter) {
    ego_bp_for_manual <- gofilter(ego_bp_for_manual, level = go_level_to_keep)
  }
  
  ego_bp_for_manual_df <- as.data.frame(ego_bp_for_manual)
  
  # =======================================================
  # 8. OPTIONAL MANUAL FILTER / CURATION
  # =======================================================
  
  if (use_manual_filter) {
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
      ego_bp_manual_df,
      file.path(out_dir, "enrichGO_BP_manual_filtered.csv"),
      row.names = FALSE
    )
    
    write.csv(
      removed_terms_df,
      file.path(out_dir, "enrichGO_BP_removed_terms.csv"),
      row.names = FALSE
    )
  } else {
    ego_bp_manual_df <- ego_bp_for_manual_df %>%
      mutate(
        GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
        Description = fct_reorder(Description, GeneRatio_num)
      ) %>%
      arrange(p.adjust, desc(Count), desc(GeneRatio_num)) %>%
      slice_head(n = max_terms_to_plot)
    
    write.csv(
      ego_bp_manual_df,
      file.path(out_dir, "enrichGO_BP_manual_filtered.csv"),
      row.names = FALSE
    )
  }
  
  # =======================================================
  # 9. PLOTS
  # =======================================================
  
  p_dot_simplified <- dotplot(ego_bp_s, showCategory = terms_to_plot) +
    ggtitle(paste0(contrast_name, ": GO BP (simplified)")) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  print(p_dot_simplified)
  
  ggsave(
    filename = file.path(out_dir, "enrichGO_BP_simplified_dotplot.pdf"),
    plot = p_dot_simplified,
    width = 9,
    height = 5
  )
  
  ggsave(
    filename = file.path(out_dir, "enrichGO_BP_simplified_dotplot.png"),
    plot = p_dot_simplified,
    width = 9,
    height = 5,
    dpi = 900
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
        axis.text.y = element_text(size = 14)
      ) +
      labs(
        title = paste0(contrast_name, ": GO BP (manual filtered)"),
        x = "GeneRatio",
        y = NULL,
        color = "p.adjust"
      ) +
      scale_x_continuous(labels = percent_format(accuracy = 1)) +
      scale_color_paletteer_c("grDevices::Plasma", direction = -1)
    
    print(p_dot_manual)
    
    ggsave(
      filename = file.path(out_dir, "enrichGO_BP_manual_filtered_dotplot.pdf"),
      plot = p_dot_manual,
      width = 9,
      height = 5
    )
    
    ggsave(
      filename = file.path(out_dir, "enrichGO_BP_manual_filtered_dotplot.png"),
      plot = p_dot_manual,
      width = 9,
      height = 5,
      dpi = 900
    )
  } else {
    cat("No GO terms remained after manual filtering.\n")
  }
  
  # =======================================================
  # 10. CONSOLE OUTPUT
  # =======================================================
  
  cat("\nTop simplified GO BP terms:\n")
  print(
    ego_bp_s_df %>%
      select(ID, Description, GeneRatio, Count, p.adjust) %>%
      head(30)
  )
  
  cat("\nTop manual-filtered GO BP terms:\n")
  print(
    ego_bp_manual_df %>%
      select(ID, Description, GeneRatio, Count, p.adjust)
  )
}