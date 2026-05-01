# =========================================================
# Shared proteins between 100 uM and 200 uM atRAL
# - makes Venn diagrams for UP and DOWN proteins
# - writes shared/unique/discordant protein tables
# - uses your exported CSV files directly
# =========================================================

# -------------------------
# Packages
# -------------------------
packages_needed <- c(
  "dplyr",
  "stringr",
  "ggplot2",
  "ggVennDiagram",
  "patchwork",
  "readr"
)

for (pkg in packages_needed) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(dplyr)
library(stringr)
library(ggplot2)
library(ggVennDiagram)
library(patchwork)
library(readr)

# =========================================================
# 1. INPUT FILES
# =========================================================

file_100_up   <- "Proteomic_output_txts/PANGEA_exports/RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr_up.csv"
file_100_down <- "Proteomic_output_txts/PANGEA_exports/RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr_down.csv"
file_200_up   <- "Proteomic_output_txts/PANGEA_exports/RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr_up.csv"
file_200_down <- "Proteomic_output_txts/PANGEA_exports/RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr_down.csv"

# optional "all" files
file_100_all  <- "Proteomic_output_txts/PANGEA_exports/RDH12_100_atRAL5hr_vs_RDH12_control_atRAL5hr_all.csv"
file_200_all  <- "Proteomic_output_txts/PANGEA_exports/RDH12_200_atRAL5hr_vs_RDH12_control_atRAL5hr_all.csv"

out_dir <- "Proteomic_Figs/Overlap_100uM_vs_200uM/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# 2. READ FILES
# =========================================================

df_100_up   <- read.csv(file_100_up,   check.names = FALSE, stringsAsFactors = FALSE)
df_100_down <- read.csv(file_100_down, check.names = FALSE, stringsAsFactors = FALSE)
df_200_up   <- read.csv(file_200_up,   check.names = FALSE, stringsAsFactors = FALSE)
df_200_down <- read.csv(file_200_down, check.names = FALSE, stringsAsFactors = FALSE)

df_100_all  <- read.csv(file_100_all,  check.names = FALSE, stringsAsFactors = FALSE)
df_200_all  <- read.csv(file_200_all,  check.names = FALSE, stringsAsFactors = FALSE)

# =========================================================
# 3. HELPER FUNCTIONS
# =========================================================

# Choose a paper-friendly identifier.
# Preference order:
#   1. Gene
#   2. name
#   3. ID
extract_protein_ids <- function(df) {
  candidate_cols <- intersect(c("Gene", "name", "ID"), colnames(df))
  
  if (length(candidate_cols) == 0) {
    stop("No usable identifier columns found. Expected one of: Gene, name, ID")
  }
  
  tmp <- df[, candidate_cols, drop = FALSE]
  
  tmp[] <- lapply(tmp, function(x) {
    x <- as.character(x)
    x <- str_trim(x)
    x[x == ""] <- NA_character_
    x
  })
  
  out <- tmp[[1]]
  
  if (ncol(tmp) > 1) {
    for (i in 2:ncol(tmp)) {
      out <- ifelse(is.na(out), tmp[[i]], out)
    }
  }
  
  out <- unique(out[!is.na(out)])
  out
}

write_vector_csv <- function(x, file, column_name = "Protein") {
  out_df <- data.frame(sort(unique(x)), stringsAsFactors = FALSE)
  colnames(out_df) <- column_name
  write.csv(out_df, file = file, row.names = FALSE)
}

# =========================================================
# 4. EXTRACT PROTEIN SETS
# =========================================================

up_100   <- extract_protein_ids(df_100_up)
down_100 <- extract_protein_ids(df_100_down)
up_200   <- extract_protein_ids(df_200_up)
down_200 <- extract_protein_ids(df_200_down)

all_100  <- extract_protein_ids(df_100_all)
all_200  <- extract_protein_ids(df_200_all)

# =========================================================
# 5. COMPUTE OVERLAPS
# =========================================================

# Shared in same direction
shared_up   <- intersect(up_100, up_200)
shared_down <- intersect(down_100, down_200)

# Unique in same direction
up_100_only   <- setdiff(up_100, up_200)
up_200_only   <- setdiff(up_200, up_100)
down_100_only <- setdiff(down_100, down_200)
down_200_only <- setdiff(down_200, down_100)

# Direction flips
up100_down200 <- intersect(up_100, down_200)
down100_up200 <- intersect(down_100, up_200)

# Any overlap regardless of direction
all_sig_100 <- union(up_100, down_100)
all_sig_200 <- union(up_200, down_200)
shared_any  <- intersect(all_sig_100, all_sig_200)

# =========================================================
# 6. SUMMARY TABLE
# =========================================================

summary_tbl <- data.frame(
  Metric = c(
    "100uM up",
    "200uM up",
    "Shared up",
    "100uM down",
    "200uM down",
    "Shared down",
    "100uM up & 200uM down",
    "100uM down & 200uM up",
    "Any significant at 100uM",
    "Any significant at 200uM",
    "Shared significant regardless of direction",
    "All proteins in 100uM all-file",
    "All proteins in 200uM all-file"
  ),
  Count = c(
    length(up_100),
    length(up_200),
    length(shared_up),
    length(down_100),
    length(down_200),
    length(shared_down),
    length(up100_down200),
    length(down100_up200),
    length(all_sig_100),
    length(all_sig_200),
    length(shared_any),
    length(all_100),
    length(all_200)
  ),
  stringsAsFactors = FALSE
)

print(summary_tbl)

write.csv(
  summary_tbl,
  file = file.path(out_dir, "summary_counts_100uM_vs_200uM.csv"),
  row.names = FALSE
)

# =========================================================
# 7. WRITE OVERLAP TABLES
# =========================================================

write_vector_csv(shared_up,       file.path(out_dir, "shared_up_100uM_200uM.csv"))
write_vector_csv(shared_down,     file.path(out_dir, "shared_down_100uM_200uM.csv"))
write_vector_csv(up_100_only,     file.path(out_dir, "up_100uM_only.csv"))
write_vector_csv(up_200_only,     file.path(out_dir, "up_200uM_only.csv"))
write_vector_csv(down_100_only,   file.path(out_dir, "down_100uM_only.csv"))
write_vector_csv(down_200_only,   file.path(out_dir, "down_200uM_only.csv"))
write_vector_csv(up100_down200,   file.path(out_dir, "up_100uM_down_200uM.csv"))
write_vector_csv(down100_up200,   file.path(out_dir, "down_100uM_up_200uM.csv"))
write_vector_csv(shared_any,      file.path(out_dir, "shared_any_direction_100uM_200uM.csv"))

# =========================================================
# 8. MAKE PAPER-FRIENDLY VENN DIAGRAMS
# =========================================================

# =========================================================
# 8. MAKE PAPER-FRIENDLY VENN DIAGRAMS
# =========================================================

venn_up_list <- list(
  "100 uM Up" = up_100,
  "200 uM Up" = up_200
)

venn_down_list <- list(
  "100 uM Down" = down_100,
  "200 uM Down" = down_200
)

p_up <- ggVennDiagram(
  venn_up_list,
  label_alpha = 0
) +
  scale_fill_gradient(low = "white", high = "#6A3D9A") +
  labs(title = "Shared up-regulated proteins") +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "none",
    plot.margin = margin(20, 35, 20, 35)
  )

p_down <- ggVennDiagram(
  venn_down_list,
  label_alpha = 0
) +
  scale_fill_gradient(low = "white", high = "#1F78B4") +
  labs(title = "Shared down-regulated proteins") +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "none",
    plot.margin = margin(20, 35, 20, 35)
  )

combined_plot <- (p_up + p_down) +
  plot_layout(ncol = 2)

print(combined_plot)

ggsave(
  filename = file.path(out_dir, "Venn_100uM_vs_200uM_up_down.pdf"),
  plot = combined_plot,
  width = 12,
  height = 5.5
)

ggsave(
  filename = file.path(out_dir, "Venn_100uM_vs_200uM_up_down.png"),
  plot = combined_plot,
  width = 12,
  height = 5.5,
  dpi = 900
)
# =========================================================
# 9. OPTIONAL: BAR PLOT SUMMARY FOR PAPER/SUPPLEMENT
# =========================================================

summary_plot_df <- summary_tbl %>%
  filter(Metric %in% c(
    "100uM up", "200uM up", "Shared up",
    "100uM down", "200uM down", "Shared down",
    "100uM up & 200uM down", "100uM down & 200uM up"
  )) %>%
  mutate(Metric = factor(Metric, levels = Metric))

p_bar <- ggplot(summary_plot_df, aes(x = Metric, y = Count)) +
  geom_col() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Overlap summary between 100 uM and 200 uM atRAL",
    x = NULL,
    y = "Protein count"
  )

print(p_bar)

ggsave(
  filename = file.path(out_dir, "Overlap_summary_barplot.pdf"),
  plot = p_bar,
  width = 12,
  height = 5
)

ggsave(
  filename = file.path(out_dir, "Overlap_summary_barplot.png"),
  plot = p_bar,
  width = 15,
  height = 5,
  dpi = 900
)

# =========================================================
# 10. CONSOLE OUTPUT
# =========================================================

cat("\nShared UP proteins:\n")
print(sort(shared_up))

cat("\nShared DOWN proteins:\n")
print(sort(shared_down))

cat("\n100uM UP and 200uM DOWN proteins:\n")
print(sort(up100_down200))

cat("\n100uM DOWN and 200uM UP proteins:\n")
print(sort(down100_up200))

cat("\nFiles written to:\n", out_dir, "\n")