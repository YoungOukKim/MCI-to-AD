# ==============================================================================
# 03_ADNI_clinical_analysis.R
#
# Purpose : ADNI CSF proteomics analysis — clinical validation of the PTGDS
#           phase boundary (Fig. 5, Table S6, Table S7)
# Input   : ADNIMERGE2 .rda files  +  CruchagaLab SomaScan 7K CSV
#
# Outputs → Paper figures / tables
#   output/TableS6_Demographics.csv   →  Additional file 1: Table S6
#   output/TableS7_Robustness.csv     →  Additional file 1: Table S7
#   ROC analysis                      →  Fig. 5E  (AUC reported in Results)
#   Survival analysis (Cox + KM)      →  Fig. 5F  (HR reported in Abstract)
#
# Key results reproduced
#   Full-cohort AUC  ≈  0.741  (PTGDS/LCN2 Checkpoint Integrity Index, MCI vs AD)
#   Table S7 AUC range: 0.56–0.72  (MMSE threshold subgroups)
#   HR  ≈  3.2  (75th-percentile risk stratification, MCI → AD conversion)
#   Cohort N = 735  (CN=172, EMCI=183, LMCI=234, AD=146)
#
# Note on AUC values
#   The full-cohort AUC uses all MCI + AD participants with complete protein
#   and covariate data. Table S7 reports subgroup AUCs restricted to
#   participants with MMSE ≤ threshold, yielding lower N and AUC.
#
# ⚠  ADNI data are subject to the ADNI Data Use Agreement and cannot be
#    redistributed. Register and download at: https://adni.loni.usc.edu
#    See data/NOTE_restricted_data.md
#
# Paper reference
#   Kim Y†, Heo W†, Park S†, et al. (2025). bioRxiv.
#   https://github.com/YoungOukKim/MCI-to-AD
#
# Usage
#   Rscript analysis/Fig5_ADNI/03_ADNI_clinical_analysis.R
#
# Requirements
#   R >= 4.3.2  |  tidyverse, pROC, survival, survminer, lubridate
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(pROC)
  library(survival)
  library(survminer)
  library(lubridate)
})

set.seed(42)

# ── Paths (relative to repository root) ───────────────────────────────────────
# Set environment variables or edit the two lines below:
#   export ADNI_PATH=/path/to/ADNIMERGE2/data/
#   export SOMA_CSV=/path/to/CruchagaLab_CSF_SOMAscan7k.csv
ADNI_PATH <- Sys.getenv("ADNI_PATH", unset = "path/to/ADNIMERGE2/data/")
SOMA_CSV  <- Sys.getenv("SOMA_CSV",  unset = "path/to/CruchagaLab_CSF_SOMAscan7k.csv")
OUT_DIR   <- "output"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ══════════════════════════════════════════════════════════════════════════════
# 1. Load data
# ══════════════════════════════════════════════════════════════════════════════
message(">>> [1/9] Loading data...")
somascan <- read.csv(SOMA_CSV, stringsAsFactors = FALSE)
for (f in c("DXSUM", "PTDEMOG", "APOERES", "UPENNBIOMK_MASTER", "MMSE", "ARM")) {
  rda <- file.path(ADNI_PATH, paste0(f, ".rda"))
  if (!file.exists(rda)) stop("File not found: ", rda,
                               "\n  Download from https://adni.loni.usc.edu")
  load(rda)
}

# ══════════════════════════════════════════════════════════════════════════════
# 2. SomaScan protein extraction (baseline, 1 row per participant)
# ══════════════════════════════════════════════════════════════════════════════
message(">>> [2/9] Extracting SomaScan proteins (PTGDS, LCN2, NEFL, TAU)...")
protein_df <- somascan %>%
  mutate(RID = as.numeric(as.character(RID))) %>%
  group_by(RID) %>% slice(1) %>% ungroup() %>%
  transmute(
    RID,
    PTGDS = as.numeric(X10514.5),
    LCN2  = as.numeric(X2836.68),
    NEFL  = as.numeric(X10082.251),
    TAU   = as.numeric(X5854.60)
  )
soma_rids <- unique(protein_df$RID)
message("  SomaScan participants: n = ", length(soma_rids))

# ══════════════════════════════════════════════════════════════════════════════
# 3. Diagnosis (3-group: CN / MCI / AD; 4-group: CN / EMCI / LMCI / AD)
# ══════════════════════════════════════════════════════════════════════════════
message(">>> [3/9] Processing diagnosis codes...")

dx_col <- grep("DIAGNOSIS|DXCURREN|DXCHANGE", colnames(DXSUM), value = TRUE)[1]

clean_dx <- DXSUM %>%
  mutate(RID = as.numeric(as.character(RID)), DX_val = .[[dx_col]]) %>%
  filter(!is.na(DX_val)) %>%
  arrange(RID, EXAMDATE) %>%
  group_by(RID) %>% slice(1) %>% ungroup() %>%
  mutate(DX_3group = case_when(
    DX_val %in% c(3, "AD", "Dementia")      ~ "AD",
    DX_val %in% c(2, "MCI", "LMCI", "EMCI") ~ "MCI",
    DX_val %in% c(1, "CN", "Normal", "NL")  ~ "CN",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(DX_3group)) %>%
  select(RID, DX_3group)

# EMCI / LMCI subdivision from ARM table
arm_dx <- ARM %>%
  mutate(RID = as.numeric(as.character(RID)),
         DX_ARM = case_when(
           str_detect(ARM, "EMCI")           ~ "EMCI",
           str_detect(ARM, "LMCI")           ~ "LMCI",
           str_detect(ARM, "AD")             ~ "AD",
           str_detect(ARM, "NL|Normal|SMC")  ~ "CN",
           str_detect(ARM, "MCI")            ~ "LMCI",
           TRUE ~ NA_character_
         )) %>%
  filter(!is.na(DX_ARM)) %>%
  distinct(RID, .keep_all = TRUE) %>%
  select(RID, DX_ARM)

clean_dx <- clean_dx %>%
  left_join(arm_dx, by = "RID") %>%
  mutate(DX_4group = case_when(
    DX_3group == "AD"  ~ "AD",
    DX_3group == "CN"  ~ "CN",
    DX_3group == "MCI" & DX_ARM == "EMCI" ~ "EMCI",
    TRUE ~ "LMCI"
  )) %>%
  filter(RID %in% soma_rids)

# ══════════════════════════════════════════════════════════════════════════════
# 4. Demographics, APOE genotype, CSF biomarkers, MMSE
# ══════════════════════════════════════════════════════════════════════════════
message(">>> [4/9] Processing covariates...")

clean_demog <- PTDEMOG %>%
  mutate(RID = as.numeric(as.character(RID))) %>%
  filter(RID %in% soma_rids) %>%
  group_by(RID) %>% slice(1) %>% ungroup() %>%
  mutate(AGE = year(as.Date(VISDATE)) - PTDOBYY,
         SEX = ifelse(grepl("Female|2", PTGENDER, ignore.case = TRUE),
                      "Female", "Male"),
         EDU = as.numeric(PTEDUCAT)) %>%
  filter(!is.na(AGE), AGE > 40, AGE < 120) %>%
  select(RID, AGE, SEX, EDU)

clean_apoe <- APOERES %>%
  mutate(RID = as.numeric(as.character(RID))) %>%
  filter(RID %in% soma_rids) %>%
  group_by(RID) %>% slice(1) %>% ungroup() %>%
  mutate(APOE4 = as.integer(grepl("4", GENOTYPE))) %>%
  select(RID, APOE4)

clean_csf <- UPENNBIOMK_MASTER %>%
  mutate(RID = as.numeric(as.character(RID))) %>%
  filter(grepl("bl|sc", VISCODE, ignore.case = TRUE),
         RID %in% soma_rids) %>%
  mutate(ABETA_num = as.numeric(gsub("[<>]", "", ABETA)),
         PTAU_num  = as.numeric(gsub("[<>]", "", PTAU))) %>%
  group_by(RID) %>% slice(1) %>% ungroup() %>%
  select(RID, ABETA_num, PTAU_num)

clean_mmse <- MMSE %>%
  mutate(RID     = as.numeric(as.character(RID)),
         MMSCORE = as.numeric(as.character(MMSCORE))) %>%
  filter(grepl("^bl$|^sc$", VISCODE, ignore.case = TRUE),
         RID %in% soma_rids, !is.na(MMSCORE),
         MMSCORE > 0, MMSCORE <= 30) %>%
  group_by(RID) %>%
  summarise(MMSE_score = mean(MMSCORE), .groups = "drop")

# ══════════════════════════════════════════════════════════════════════════════
# 5. Merge final cohort
# ══════════════════════════════════════════════════════════════════════════════
message(">>> [5/9] Merging cohort...")

full_cohort <- protein_df %>%
  inner_join(clean_dx,    by = "RID") %>%
  left_join(clean_demog,  by = "RID") %>%
  left_join(clean_apoe,   by = "RID") %>%
  left_join(clean_csf,    by = "RID") %>%
  left_join(clean_mmse,   by = "RID") %>%
  filter(!is.na(AGE))

message("  Final cohort N = ", nrow(full_cohort))
print(table(full_cohort$DX_4group))

# ══════════════════════════════════════════════════════════════════════════════
# 6. TABLE S6 — Participant demographics
# ══════════════════════════════════════════════════════════════════════════════
message(">>> [6/9] Generating Table S6 (demographics)...")

table_s6 <- full_cohort %>%
  mutate(DX_4group = factor(DX_4group, levels = c("CN","EMCI","LMCI","AD"))) %>%
  group_by(DX_4group) %>%
  summarise(
    N             = n(),
    Age_mean      = round(mean(AGE),              1),
    Age_sd        = round(sd(AGE),                1),
    Female_pct    = round(mean(SEX == "Female") * 100, 1),
    Education_mean = round(mean(EDU,  na.rm = TRUE), 1),
    Education_sd   = round(sd(EDU,   na.rm = TRUE), 1),
    APOE4_pct     = round(mean(APOE4, na.rm = TRUE) * 100, 1),
    MMSE_mean     = round(mean(MMSE_score, na.rm = TRUE), 1),
    MMSE_sd       = round(sd(MMSE_score,  na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  arrange(DX_4group)

print(table_s6)
message("  Total N = ", sum(table_s6$N))
write.csv(table_s6, file.path(OUT_DIR, "TableS6_Demographics.csv"),
          row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# 7. ROC — Full-cohort AUC (PTGDS/LCN2 Checkpoint Integrity Index)
# ══════════════════════════════════════════════════════════════════════════════
message(">>> [7/9] ROC analysis (full cohort, MCI vs AD)...")

roc_df <- full_cohort %>%
  filter(DX_3group %in% c("MCI", "AD")) %>%
  mutate(
    outcome   = as.integer(DX_3group == "AD"),
    log_ratio = log2(PTGDS / (LCN2 + 1e-6)),
    log_nefl  = log2(NEFL  + 1e-6),
    log_tau   = log2(TAU   + 1e-6),
    inter_act = log_ratio * APOE4
  ) %>%
  drop_na(outcome, log_ratio, log_nefl, log_tau,
          AGE, APOE4, ABETA_num, PTAU_num)

final_model <- glm(
  outcome ~ log_ratio + log_nefl + log_tau + AGE + APOE4 +
    ABETA_num + PTAU_num + inter_act,
  data = roc_df, family = binomial)

roc_df$checkpoint_index <- predict(final_model, type = "response")
roc_obj <- roc(roc_df$outcome, roc_df$checkpoint_index, ci = TRUE, quiet = TRUE)

message(sprintf("  Full-cohort ROC: N = %d  |  AUC = %.4f (95%% CI: %.3f–%.3f)",
                nrow(roc_df), auc(roc_obj), ci(roc_obj)[1], ci(roc_obj)[3]))

# ══════════════════════════════════════════════════════════════════════════════
# 8. TABLE S7 — AUC robustness across MMSE inclusion thresholds
# ══════════════════════════════════════════════════════════════════════════════
message(">>> [8/9] Generating Table S7 (MMSE threshold robustness)...")

roc_mmse     <- roc_df %>% filter(!is.na(MMSE_score))
thresholds   <- c(24, 25, 26, 27, 28, 30)
table_s7_rows <- list()

for (th in thresholds) {
  sub   <- roc_mmse %>% filter(MMSE_score <= th)
  n_ad  <- sum(sub$outcome == 1L)
  n_mci <- sum(sub$outcome == 0L)
  if (n_ad >= 10L && n_mci >= 10L) {
    r  <- roc(sub$outcome, sub$checkpoint_index, quiet = TRUE)
    ci <- ci.auc(r)
    table_s7_rows[[as.character(th)]] <- data.frame(
      Inclusion_Criteria = paste0("MMSE <= ", th),
      N = nrow(sub), N_MCI = n_mci, N_AD = n_ad,
      AUC = round(auc(r), 3),
      CI_Lower = round(ci[1], 3), CI_Upper = round(ci[3], 3))
  }
}

table_s7 <- do.call(rbind, table_s7_rows)
print(table_s7)
write.csv(table_s7, file.path(OUT_DIR, "TableS7_Robustness.csv"),
          row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# 9. Survival analysis — Cox proportional hazards + Kaplan–Meier (Fig. 5F)
# ══════════════════════════════════════════════════════════════════════════════
message(">>> [9/9] Survival analysis (Cox + Kaplan-Meier, Fig. 5F)...")

surv_timeline <- DXSUM %>%
  mutate(RID      = as.numeric(as.character(RID)),
         DX_val   = .[[dx_col]],
         EXAMDATE = as.Date(EXAMDATE)) %>%
  filter(RID %in% roc_df$RID) %>%
  group_by(RID) %>%
  summarise(
    Start_Date = min(EXAMDATE),
    Last_Date  = max(EXAMDATE),
    Event      = as.integer(any(DX_val %in% c(3, "Dementia", "AD"))),
    AD_Date    = suppressWarnings(
      min(EXAMDATE[DX_val %in% c(3, "Dementia", "AD")], na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(Time_months = ifelse(
    Event == 1L & !is.infinite(as.numeric(AD_Date)),
    as.numeric(difftime(AD_Date,   Start_Date, units = "days")) / 30.44,
    as.numeric(difftime(Last_Date, Start_Date, units = "days")) / 30.44
  )) %>%
  filter(!is.na(Time_months), Time_months > 0)

cox_df <- roc_df %>% inner_join(surv_timeline, by = "RID")
cutoff_75      <- quantile(cox_df$checkpoint_index, 0.75)
cox_df$Risk_Group <- factor(
  ifelse(cox_df$checkpoint_index > cutoff_75, "High_Risk", "Low_Risk"),
  levels = c("Low_Risk", "High_Risk"))

cox_model <- coxph(Surv(Time_months, Event) ~ Risk_Group, data = cox_df)
hr    <- exp(coef(cox_model))
hr_ci <- exp(confint(cox_model))

message(sprintf("  Cox HR = %.3f  (95%% CI: %.3f–%.3f)",
                hr, hr_ci[1], hr_ci[2]))

# ══════════════════════════════════════════════════════════════════════════════
# 10. Summary
# ══════════════════════════════════════════════════════════════════════════════
message("\n", strrep("=", 65))
message("  ADNI Analysis Summary")
message(strrep("=", 65))
message(sprintf("  Cohort: N = %d  (CN=%d, EMCI=%d, LMCI=%d, AD=%d)",
                sum(table_s6$N),
                table_s6$N[table_s6$DX_4group == "CN"],
                table_s6$N[table_s6$DX_4group == "EMCI"],
                table_s6$N[table_s6$DX_4group == "LMCI"],
                table_s6$N[table_s6$DX_4group == "AD"]))
message(sprintf("  Full-cohort AUC = %.4f", auc(roc_obj)))
message(sprintf("  Cox HR = %.3f  (95%% CI: %.3f–%.3f)", hr, hr_ci[1], hr_ci[2]))
message(strrep("=", 65))
message("  Output files:")
message("    TableS6_Demographics.csv  →  Additional file 1: Table S6")
message("    TableS7_Robustness.csv    →  Additional file 1: Table S7")
message(">>> Script 03 complete.")
