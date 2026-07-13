# =========================================================
# RDH vs GFP two-circle Venn diagram for PPF protein IDs
# From PPFtable_v4 created by 01_Prep_data.R
#
# Input:
#   PPFtable_v4
#
# Outputs:
#   1. Two-circle Venn diagram of unique Protein.Group IDs
#   2. CSV tables for RDH-only, GFP-only, shared, and membership counts
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(readr)
  library(tibble)
})

# ---------------------------------------------------------
# 1. INPUT / OUTPUT
# ---------------------------------------------------------

input_dir <- if (dir.exists("PPF_datasets")) "PPF_datasets" else "."
prep_script_path <- file.path(input_dir, "01_Prep_data.R")
out_dir <- file.path(input_dir, "output_venn")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Detection rule
# A protein is counted as present in RDH or GFP if it is detected in any
# replicate for that group.
min_intensity_for_detected <- 0


# ---------------------------------------------------------
# 2. READ PREPARED PPF DATA
# ---------------------------------------------------------

if (!file.exists(prep_script_path)) {
  prep_script_path <- "01_Prep_data.R"
}

if (!file.exists(prep_script_path)) {
  stop("Could not find 01_Prep_data.R. Run this script from the project root or PPF_datasets folder.")
}

source(prep_script_path)

if (!exists("PPFtable_v4")) {
  stop("PPFtable_v4 was not created by 01_Prep_data.R.")
}

required_cols <- c("sample_name", "Protein.Group", "Genes", "PG.MaxLFQ")
missing_cols <- setdiff(required_cols, colnames(PPFtable_v4))

if (length(missing_cols) > 0) {
  stop("PPFtable_v4 is missing required columns: ", paste(missing_cols, collapse = ", "))
}

if (!is.numeric(PPFtable_v4$PG.MaxLFQ)) {
  stop("PG.MaxLFQ must be numeric before calculating RDH/GFP protein overlap.")
}

cat("\nPPFtable_v4 dimensions:\n")
print(dim(PPFtable_v4))

cat("\nSamples found:\n")
print(sort(unique(PPFtable_v4$sample_name)))


# ---------------------------------------------------------
# 3. BUILD RDH/GFP PROTEIN MEMBERSHIP
# ---------------------------------------------------------

protein_detection_long <- PPFtable_v4 %>%
  mutate(
    SampleName = as.character(sample_name),
    Group = str_extract(SampleName, "^(RDH|GFP)"),
    ProteinID = Protein.Group,
    Gene = if_else(is.na(Genes) | Genes == "", Protein.Group, Genes),
    Detected = !is.na(PG.MaxLFQ) & PG.MaxLFQ > min_intensity_for_detected
  ) %>%
  filter(Group %in% c("RDH", "GFP"))

if (!all(c("RDH", "GFP") %in% unique(protein_detection_long$Group))) {
  stop("Could not find both RDH and GFP sample groups in sample_name.")
}

protein_membership <- protein_detection_long %>%
  group_by(ProteinID, Gene, Group) %>%
  summarise(
    n_replicates = n_distinct(SampleName),
    n_detected_replicates = n_distinct(SampleName[Detected]),
    present = any(Detected),
    .groups = "drop"
  ) %>%
  select(ProteinID, Gene, Group, present, n_replicates, n_detected_replicates) %>%
  pivot_wider(
    names_from = Group,
    values_from = c(present, n_replicates, n_detected_replicates),
    values_fill = list(
      present = FALSE,
      n_replicates = 0,
      n_detected_replicates = 0
    )
  ) %>%
  mutate(
    present_RDH = coalesce(present_RDH, FALSE),
    present_GFP = coalesce(present_GFP, FALSE),
    membership = case_when(
      present_RDH & present_GFP ~ "Shared",
      present_RDH & !present_GFP ~ "RDH_only",
      !present_RDH & present_GFP ~ "GFP_only",
      TRUE ~ "Not_detected"
    )
  )

protein_membership_detected <- protein_membership %>%
  filter(membership != "Not_detected")

rdh_only_proteins <- protein_membership_detected %>%
  filter(membership == "RDH_only")

gfp_only_proteins <- protein_membership_detected %>%
  filter(membership == "GFP_only")

shared_proteins <- protein_membership_detected %>%
  filter(membership == "Shared")

venn_counts <- tibble(
  category = c("RDH_only", "GFP_only", "Shared", "RDH_total", "GFP_total"),
  n_proteins = c(
    nrow(rdh_only_proteins),
    nrow(gfp_only_proteins),
    nrow(shared_proteins),
    sum(protein_membership_detected$present_RDH),
    sum(protein_membership_detected$present_GFP)
  )
)

cat("\nRDH/GFP Venn counts:\n")
print(venn_counts)


# ---------------------------------------------------------
# 4. WRITE TABLES
# ---------------------------------------------------------

write.csv(
  venn_counts,
  file.path(out_dir, "RDH_GFP_venn_counts.csv"),
  row.names = FALSE
)

write.csv(
  protein_membership_detected,
  file.path(out_dir, "RDH_GFP_protein_membership.csv"),
  row.names = FALSE
)

write.csv(
  rdh_only_proteins,
  file.path(out_dir, "RDH_only_proteins.csv"),
  row.names = FALSE
)

write.csv(
  gfp_only_proteins,
  file.path(out_dir, "GFP_only_proteins.csv"),
  row.names = FALSE
)

write.csv(
  shared_proteins,
  file.path(out_dir, "RDH_GFP_shared_proteins.csv"),
  row.names = FALSE
)


# ---------------------------------------------------------
# 5. TWO-CIRCLE VENN DIAGRAM
# ---------------------------------------------------------

circle_points <- function(center_x, center_y, radius, label, n = 360) {
  theta <- seq(0, 2 * pi, length.out = n)
  tibble(
    x = center_x + radius * cos(theta),
    y = center_y + radius * sin(theta),
    circle = label
  )
}

venn_circles <- bind_rows(
  circle_points(center_x = -0.75, center_y = 0, radius = 1.4, label = "RDH"),
  circle_points(center_x = 0.75, center_y = 0, radius = 1.4, label = "GFP")
)

rdh_total <- venn_counts %>%
  filter(category == "RDH_total") %>%
  pull(n_proteins)

gfp_total <- venn_counts %>%
  filter(category == "GFP_total") %>%
  pull(n_proteins)

rdh_only <- venn_counts %>%
  filter(category == "RDH_only") %>%
  pull(n_proteins)

gfp_only <- venn_counts %>%
  filter(category == "GFP_only") %>%
  pull(n_proteins)

shared <- venn_counts %>%
  filter(category == "Shared") %>%
  pull(n_proteins)

p_venn <- ggplot() +
  geom_polygon(
    data = venn_circles,
    aes(x = x, y = y, fill = circle, group = circle),
    alpha = 0.35,
    color = "grey20",
    size = 0.8
  ) +
  annotate("text", x = -1.35, y = 0, label = rdh_only, size = 7, fontface = "bold") +
  annotate("text", x = 0, y = 0, label = shared, size = 7, fontface = "bold") +
  annotate("text", x = 1.35, y = 0, label = gfp_only, size = 7, fontface = "bold") +
  annotate("text", x = -1.2, y = 1.55, label = paste0("RDH\nn = ", rdh_total), size = 5, fontface = "bold") +
  annotate("text", x = 1.2, y = 1.55, label = paste0("GFP\nn = ", gfp_total), size = 5, fontface = "bold") +
  scale_fill_manual(values = c("RDH" = "#2B8CBE", "GFP" = "#F03B20")) +
  coord_equal(xlim = c(-2.6, 2.6), ylim = c(-1.7, 1.9), expand = FALSE) +
  labs(
    title = "Unique Protein.Group IDs Detected in RDH vs GFP",
    subtitle = paste0("Detected in any replicate; PG.MaxLFQ > ", min_intensity_for_detected),
    caption = "Left = RDH only, center = shared, right = GFP only"
  ) +
  theme_void(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 0.5, color = "grey35")
  )

print(p_venn)

ggsave(
  file.path(out_dir, "RDH_GFP_unique_protein_venn.pdf"),
  p_venn,
  width = 7,
  height = 5
)

ggsave(
  file.path(out_dir, "RDH_GFP_unique_protein_venn.png"),
  p_venn,
  width = 7,
  height = 5,
  dpi = 300
)


# ---------------------------------------------------------
# 6. FINAL MESSAGE
# ---------------------------------------------------------

cat("\nSaved output directory:\n")
cat(out_dir, "\n")

cat("\nMain output files:\n")
cat(file.path(out_dir, "RDH_GFP_unique_protein_venn.pdf"), "\n")
cat(file.path(out_dir, "RDH_GFP_unique_protein_venn.png"), "\n")
cat(file.path(out_dir, "RDH_GFP_venn_counts.csv"), "\n")
cat(file.path(out_dir, "RDH_GFP_protein_membership.csv"), "\n")
cat(file.path(out_dir, "RDH_only_proteins.csv"), "\n")
cat(file.path(out_dir, "GFP_only_proteins.csv"), "\n")
cat(file.path(out_dir, "RDH_GFP_shared_proteins.csv"), "\n")
