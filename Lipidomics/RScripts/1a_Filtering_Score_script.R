library(dplyr)
library(stringr)
library(tidyr)

df <- Recovery_lipidsall_v7

df<-df


# -----------------------------------------------------------------------------
# Function: score_candidate
# Purpose:
#   Assign a heuristic score to a lipid candidate parsed from a name string.
#   Higher scores indicate more biologically plausible annotations given
#   common lipidomics patterns and expected chemistry for PC/PE lipids.
#
# Inputs:
#   C_total  : total number of carbons parsed from candidate (e.g., 38 from 38:4)
#   DB_total : total number of double bonds parsed from candidate
#   text     : full candidate annotation string (used to detect lipid class
#              or special prefixes like O-, dO-, etc.)
#
# Output:
#   numeric score used to rank competing lipid candidates
# -----------------------------------------------------------------------------

score_candidate <- function(C_total, DB_total, text) {
  
  # Initialize score accumulator
  score <- 0
  
  # ---------------------------------------------------------------------------
  # 1) Carbon parity rule
  # Most mammalian lipids have even total carbons due to fatty acid synthesis
  # occurring via 2-carbon elongation steps (acetyl-CoA).
  #
  # +2 points if even carbon total
  # -1 point if odd carbon total or missing
  # ---------------------------------------------------------------------------
  score <- score + ifelse(!is.na(C_total) & C_total %% 2 == 0, 2, -1)
  
  
  # ---------------------------------------------------------------------------
  # 2) Plausible carbon total range
  # Ensures the lipid total carbon count falls within a realistic biological
  # range for most phospholipids measured in typical lipidomics panels.
  #
  # Range chosen: 14–60 carbons
  #
  # +1 if within range
  # -1 if outside range or missing
  # ---------------------------------------------------------------------------
  score <- score + ifelse(!is.na(C_total) & C_total >= 14 & C_total <= 60, 1, -1)
  
  
  # ---------------------------------------------------------------------------
  # 3) Plausible double bond range
  # Double bonds in phospholipids typically fall within 0–20 total DB.
  # This rule just ensures the parsed value isn't unrealistic.
  #
  # +1 if DB_total within range
  # -1 if outside range or missing
  # ---------------------------------------------------------------------------
  score <- score + ifelse(!is.na(DB_total) & DB_total >= 0 & DB_total <= 20, 1, -1)
  
  
  # ---------------------------------------------------------------------------
  # 4) Detect lipid class or special prefixes from annotation string
  #
  # These flags allow class-specific scoring rules later in the function.
  # ---------------------------------------------------------------------------
  is_pc <- stringr::str_detect(text, "\\bPC\\b")        # Phosphatidylcholine
  is_pe <- stringr::str_detect(text, "\\bPE\\b")        # Phosphatidylethanolamine
  is_dO <- stringr::str_detect(text, "\\bPC\\s+dO-")    # unusual "dO-" PC notation
  is_pi <- stringr::str_detect(text, "\\bPE\\b")        # Phosphatidylethanolamine
  
  
  # ---------------------------------------------------------------------------
  # 5) Penalize extreme unsaturation for PCs
  #
  # HG-only transitions (e.g., m/z 184 for PC) cannot distinguish many species,
  # so highly unsaturated PCs (≥6 DB) are less likely and often artifacts of
  # ambiguous naming. These are penalized.
  #
  # ≥6 DB → moderate penalty
  # ≥8 DB → stronger penalty
  # ---------------------------------------------------------------------------
  score <- score + ifelse(is_pc & !is.na(DB_total) & DB_total >= 7, -3, 0)
  score <- score + ifelse(is_pc & !is.na(DB_total) & DB_total >= 8, -2, 0)
  
  
  # ---------------------------------------------------------------------------
  # 6) Penalize extreme unsaturation for PE
  #
  # Similar logic as PC but slightly relaxed. Extremely high unsaturation in
  # PE species can occur but is less common; thus a mild penalty is applied
  # only at ≥8 DB.
  # ---------------------------------------------------------------------------
  score <- score + ifelse(is_pe & !is.na(DB_total) & DB_total >= 8, -3, 0)
  score <- score + ifelse(is_pi & !is.na(DB_total) & DB_total >= 7, -3, 0)
  
  
  
  # ---------------------------------------------------------------------------
  # 7) Penalize unusual PC dO- prefixes
  #
  # The "dO-" prefix is uncommon for PC annotations and often appears due to
  # naming inconsistencies or mis-parsing. These are discouraged.
  # ---------------------------------------------------------------------------
  score <- score + ifelse(is_dO, -3, 0)
  
  
  # ---------------------------------------------------------------------------
  # 8) Strong penalty for internal standards
  #
  # Candidates containing "STD" likely represent internal standards rather
  # than endogenous lipids, so they are strongly penalized.
  # ---------------------------------------------------------------------------
  score <- score + ifelse(stringr::str_detect(text, "\\bSTD\\b"), -10, 0)
  
  
  # ---------------------------------------------------------------------------
  # Return final score
  # Higher score = more plausible candidate
  # ---------------------------------------------------------------------------
  score
}

parsed <- df %>%
  mutate(
    lipid_name_raw = lipid_name,
    lipid_name = str_squish(lipid_name),
    
    # Split into candidates if multiple are separated by "_"
    # e.g., "PC 34:1_PC O-35:1" -> 2 candidates
    candidate_list = str_split(lipid_name, "_")
  ) %>%
  unnest(candidate_list) %>%
  mutate(candidate = str_squish(candidate_list)) %>%
  select(-candidate_list) %>%
  mutate(
    # Extract FIRST "C:DB" occurrence inside each candidate (treats d18:1 as 18:1, etc.)
    C_total  = as.integer(str_match(candidate, "(\\d{1,2})\\s*:\\s*(\\d{1,2})")[,2]),
    DB_total = as.integer(str_match(candidate, "(\\d{1,2})\\s*:\\s*(\\d{1,2})")[,3]),
    
    # Extract NL chain when present: "NL 18:1"
    NL_C  = as.integer(str_match(candidate, "NL\\s*(\\d{1,2})\\s*:\\s*(\\d{1,2})")[,2]),
    NL_DB = as.integer(str_match(candidate, "NL\\s*(\\d{1,2})\\s*:\\s*(\\d{1,2})")[,3]),
    
    # Also capture patterns like "]_22:6" or "_15:0" used in some TG naming
    NL_C_alt  = as.integer(str_match(candidate, "(?<=\\]_)(\\d{1,2})\\s*:\\s*(\\d{1,2})")[,2]),
    NL_DB_alt = as.integer(str_match(candidate, "(?<=\\]_)(\\d{1,2})\\s*:\\s*(\\d{1,2})")[,3]),
    
    NL_C_alt2  = as.integer(str_match(candidate, "(?<=\\]_)(\\d{1,2})\\s*:\\s*(\\d{1,2})")[,2]),
    NL_DB_alt2 = as.integer(str_match(candidate, "(?<=\\]_)(\\d{1,2})\\s*:\\s*(\\d{1,2})")[,3])
  ) %>%
  mutate(
    # If NL not found via "NL 18:1", try "_15:0" at end (common in your TG strings)
    NL_C2  = as.integer(str_match(candidate, "(?:\\]|\\b|_)\\s*(\\d{1,2})\\s*:\\s*(\\d{1,2})\\s*$")[,2]),
    NL_DB2 = as.integer(str_match(candidate, "(?:\\]|\\b|_)\\s*(\\d{1,2})\\s*:\\s*(\\d{1,2})\\s*$")[,3]),
    
    # Final NL picks (first non-NA among patterns)
    NL_C_final  = coalesce(NL_C, NL_C_alt, NL_C2),
    NL_DB_final = coalesce(NL_DB, NL_DB_alt, NL_DB2),
    
    # Score each candidate so we can pick one per original lipid_name row
    cand_score = score_candidate(C_total, DB_total, candidate)
  )

# ---- pick the best candidate per original row (lipid_name_raw + mrm typically define the feature) ----
picked <- parsed %>%
  group_by(lipid_name_raw, mrm) %>%
  arrange(desc(cand_score), .by_group = TRUE) %>%
 dplyr:: slice(1) %>%
  ungroup() %>%
  transmute(
    lipid_name_raw,
    mrm,
    picked_candidate = candidate,
    C_total,
    DB_total,
    NL_C = NL_C_final,
    NL_DB = NL_DB_final,
    cand_score
  )

# Join back to your df (so every row gets the picked annotation + unsaturation)
df_with_unsat <- df %>%
  left_join(picked, by = c("lipid_name" = "lipid_name_raw", "mrm"))

# Look at the result
df_with_unsat %>% select(lipid_name, mrm, picked_candidate, C_total, DB_total, NL_C, NL_DB, cand_score) %>% head(20)

#===================================================
library(dplyr)
library(stringr)

df <- df

cer_parsed <- df %>%
  mutate(lipid_name_clean = str_squish(lipid_name)) %>%
  # keep only ceramides that have the (d..../..:..) pattern
  filter(str_detect(lipid_name_clean, "^Cer\\(d\\d{1,2}:\\d{1,2}/\\d{1,2}:\\d{1,2}")) %>%
  mutate(
    # Extract sn1 (sphingoid base): d15:2  -> C=15, DB=2
    sn1_C  = as.integer(str_match(lipid_name_clean, "Cer\\(d(\\d{1,2}):(\\d{1,2})/")[,2]),
    sn1_DB = as.integer(str_match(lipid_name_clean, "Cer\\(d(\\d{1,2}):(\\d{1,2})/")[,3]),
    
    # Extract sn2 (acyl chain): 22:0(2OH) -> C=22, DB=0
    sn2_C  = as.integer(str_match(lipid_name_clean, "/(\\d{1,2}):(\\d{1,2})")[,2]),
    sn2_DB = as.integer(str_match(lipid_name_clean, "/(\\d{1,2}):(\\d{1,2})")[,3]),
    
    # Extract hydroxyl count on sn2 if present: (2OH) -> 2
    sn2_OH = as.integer(str_match(lipid_name_clean, "\\((\\d+)OH\\)")[,2]),
    
    # Sum composition
    Cer_C_total  = sn1_C + sn2_C,
    Cer_DB_total = sn1_DB + sn2_DB
  )

# Look at a few examples
cer_parsed %>%
  select(lipid_name, sn1_C, sn1_DB, sn2_C, sn2_DB, sn2_OH, Cer_C_total, Cer_DB_total) %>%
  head(20)

