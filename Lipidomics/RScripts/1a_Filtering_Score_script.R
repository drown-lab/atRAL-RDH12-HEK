library(dplyr)
library(stringr)
library(tidyr)


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
  
  # 1) Prefer even total carbons
  score <- score + ifelse(!is.na(C_total) & C_total %% 2 == 0, 2, -1)
  
  # 2) Plausible carbon total range
  score <- score + ifelse(!is.na(C_total) & C_total >= 14 & C_total <= 60, 1, -1)
  
  # 3) Plausible double bond range
  score <- score + ifelse(!is.na(DB_total) & DB_total >= 0 & DB_total <= 20, 1, -1)
  
  # 4) Detect lipid classes / prefixes
  is_pc <- str_detect(text, "\\bPC\\b")
  is_pe <- str_detect(text, "\\bPE\\b")
  is_pi <- str_detect(text, "\\bPI\\b")
  is_dO <- str_detect(text, "\\bPC\\s+dO-")
  
  # 5) PC-specific unsaturation logic
  # High DB PCs can be real, especially at higher total carbons,
  # so use milder penalties than before.
  score <- score + ifelse(is_pc & !is.na(DB_total) & DB_total >= 8, -1, 0)
  
  # Small bonus for large, highly unsaturated even-carbon PCs
  score <- score + ifelse(
    is_pc & !is.na(C_total) & !is.na(DB_total) &
      C_total >= 42 & DB_total >= 7 & C_total %% 2 == 0,
    2, 0
  )
  
  # 6) PE / PI extreme unsaturation penalties
  score <- score + ifelse(is_pe & !is.na(DB_total) & DB_total >= 8, -2, 0)
  score <- score + ifelse(is_pi & !is.na(DB_total) & DB_total >= 7, -2, 0)
  
  # 7) Penalize unusual PC dO- prefix
  score <- score + ifelse(is_dO, -3, 0)
  
  # 8) Penalize internal standards
  score <- score + ifelse(str_detect(text, "\\bSTD\\b"), -10, 0)
  
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


##merge tables
df_with_unsat_updated <- df_with_unsat %>%
  left_join(
    cer_parsed %>%
      select(lipid_name, mrm1, Cer_C_total, Cer_DB_total),
    by = c("lipid_name", "mrm1")
  ) %>%
  mutate(
    C_total  = ifelse(lipid_class == "Cer" & !is.na(Cer_C_total), Cer_C_total, C_total),
    DB_total = ifelse(lipid_class == "Cer" & !is.na(Cer_DB_total), Cer_DB_total, DB_total)
  ) %>%
  select(-Cer_C_total, -Cer_DB_total)