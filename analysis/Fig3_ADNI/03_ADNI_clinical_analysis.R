# ==============================================================================
# 03_ADNI_clinical_analysis.R
#
# Purpose : ADNI CSF proteomics analysis supporting the astrocytic PTGDS
#           inflection (vertex CPS 0.47). Reports CSF marker-cognition
#           associations against a NEFL positive control; tests - and does
#           NOT support - the PTGDS/LCN2 ratio as a predictor of MCI->AD
#           conversion. Supports Fig. 3 (CSF concordance) and
#           Supplementary Table S8.
# Input   : ADNIMERGE2 .rda files  +  CruchagaLab SomaScan 7K CSV
#
# Outputs (paper figures / tables)
#   output/TableS8_Demographics.csv    -> Supplementary Table S8 (demographics)
#   output/Fig3_CSF_associations.csv   -> Fig. 3 supporting data (CSF associations)
#   Fig. 3 concordance panels are rendered by Fig3_CSF_concordance.R
#
# Analysis notes:
#   CSF biomarker associations with cognition are evaluated against a NEFL
#   positive control. The incremental value of the PTGDS/LCN2 ratio is assessed
#   by 10-fold cross-validation, and MCI->AD conversion by a pre-specified
#   baseline-MCI Cox model. This script reports CSF associations and the
#   positive control; it makes NO clinical-prediction claim. Consistent with
#   the manuscript, the ratio shows no incremental cross-validated AUC and a
#   null Cox hazard ratio.
#
#   ADNI data are subject to the ADNI Data Use Agreement and cannot be
#   redistributed. Register and download at: https://adni.loni.usc.edu
#   See data/NOTE_restricted_data.md
#
# Paper reference
#   Kim Y, Heo W, Park SJ, et al. (2026). A biphasic astrocytic PTGDS
#   trajectory marks a metabolic vulnerability stage in prodromal Alzheimer's
#   disease. Research Square. https://doi.org/10.21203/rs.3.rs-9499795
#   https://github.com/YoungOukKim/MCI-to-AD
#
# Usage
#   Rscript analysis/Fig3_ADNI/03_ADNI_clinical_analysis.R
#
# Requirements
#   R >= 4.3.2  |  tidyverse, survival, lubridate
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(lubridate)
})

set.seed(42)

# -- Paths (relative to repository root) ---------------------------------------
# Set environment variables or edit the two lines below:
#   export ADNI_PATH=/path/to/ADNIMERGE2/data/
#   export SOMA_CSV=/path/to/CruchagaLab_CSF_SOMAscan7k.csv
ADNI_PATH <- Sys.getenv("ADNI_PATH", unset = "path/to/ADNIMERGE2/data/")
SOMA_CSV  <- Sys.getenv("SOMA_CSV",  unset = "path/to/CruchagaLab_CSF_SOMAscan7k.csv")
OUT_DIR   <- "output"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------------------
message(">>> [1/9] Loading data...")
somascan <- read.csv(SOMA_CSV, stringsAsFactors = FALSE)
for (f in c("DXSUM", "PTDEMOG", "APOERES", "UPENNBIOMK_MASTER", "MMSE", "ARM")) {
  rda <- file.path(ADNI_PATH, paste0(f, ".rda"))
  if (!file.exists(rda)) stop("File not found: ", rda,
                               "\n  Download from https://adni.loni.usc.edu")
  load(rda)
}

# ------------------------------------------------------------------------------
# 2. SomaScan protein extraction (baseline, 1 row per participant)
# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
# 3. Diagnosis (3-group: CN / MCI / AD; 4-group: CN / EMCI / LMCI / AD)
# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
# 4. Demographics, APOE genotype, CSF biomarkers, MMSE
# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
# 5. Merge final cohort
# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
# 6. Supplementary Table S8 - participant demographics
# ------------------------------------------------------------------------------
message(">>> [6/9] Generating Supplementary Table S8 (demographics)...")

table_s8 <- full_cohort %>%
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

print(table_s8)
message("  Total N = ", sum(table_s8$N))
write.csv(table_s8, file.path(OUT_DIR, "TableS8_Demographics.csv"),
          row.names = FALSE)

# ------------------------------------------------------------------------------
# 7. CSF analysis - associations + positive control (NEFL) -> Fig. 3 support
# ------------------------------------------------------------------------------
message(">>> [7/9] CSF associations (positive control = NEFL)...")

ana <- full_cohort %>%
  mutate(
    l_ptgds = as.numeric(scale(log2(PTGDS + 1e-6))),
    l_lcn2  = as.numeric(scale(log2(LCN2  + 1e-6))),
    l_nefl  = as.numeric(scale(log2(NEFL  + 1e-6))),
    path_axis = suppressWarnings(log2(PTAU_num / ABETA_num))
  )

cor_rp <- function(a, b) {
  ok <- is.finite(a) & is.finite(b); a <- a[ok]; b <- b[ok]; n <- length(a)
  r <- cor(a, b); t <- r * sqrt((n - 2) / (1 - r^2))
  list(r = r, p = 2 * pt(-abs(t), n - 2), n = n)
}

# Positive control: NEFL must track the pathology axis (|r| > 0.3) for validity
pc <- cor_rp(ana$l_nefl, ana$path_axis)
message(sprintf("  [positive control] NEFL vs pathology axis: r = %+.3f (p = %.2g, n = %d)%s",
                pc$r, pc$p, pc$n, ifelse(abs(pc$r) > 0.3, "  PASS", "  FAIL -> interpret with caution")))

assoc <- list()
for (g in c("l_ptgds", "l_lcn2", "l_nefl")) {
  mm <- cor_rp(ana[[g]], ana$MMSE_score)
  pa <- cor_rp(ana[[g]], ana$path_axis)
  assoc[[g]] <- data.frame(
    marker     = sub("^l_", "", g),
    r_MMSE     = round(mm$r, 3), p_MMSE = signif(mm$p, 2), n_MMSE = mm$n,
    r_pathology= round(pa$r, 3), p_pathology = signif(pa$p, 2), n_path = pa$n)
}
csf_assoc <- do.call(rbind, assoc)
cat("\n  CSF marker associations (|r| < 0.2 => not a clinical predictor):\n")
print(csf_assoc, row.names = FALSE)
write.csv(csf_assoc, file.path(OUT_DIR, "Fig3_CSF_associations.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 8. Incremental predictive value of PTGDS/LCN2 ratio - 10-fold CV AUC
#    NULL CHECK (not a prediction claim): the full-model AUC reflects the
#    established predictors (NEFL, age, APOE4); the quantity of interest is the
#    INCREMENT from adding the ratio (expected ~0) and the ratio-ALONE AUC
#    (expected ~0.5 = chance). This refutes any predictive role for the ratio.
# ------------------------------------------------------------------------------
message(">>> [8/9] Cross-validated incremental value of PTGDS/LCN2 ratio (null check)...")

cvauc <- function(form, d, K = 10) {
  d <- d[complete.cases(d[, all.vars(form)]), ]
  fold <- sample(rep(1:K, length.out = nrow(d)))
  pr <- numeric(nrow(d))
  for (k in 1:K) {
    fit <- suppressWarnings(glm(form, data = d[fold != k, ], family = binomial))
    pr[fold == k] <- predict(fit, d[fold == k, ], type = "response")
  }
  y <- d[[all.vars(form)[1]]]
  # rank-based AUC (no extra package)
  r1 <- rank(pr); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r1[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
cvd <- ana %>% filter(DX_3group %in% c("MCI", "AD")) %>%
  mutate(outcome = as.integer(DX_3group == "AD"),
         ratio = l_ptgds - l_lcn2) %>%
  select(outcome, ratio, l_nefl, l_ptgds, l_lcn2, AGE, APOE4) %>% drop_na()
set.seed(42)
auc_full  <- cvauc(outcome ~ ratio + l_nefl + AGE + APOE4, cvd)
auc_base  <- cvauc(outcome ~ l_nefl + AGE + APOE4, cvd)         # WITHOUT ratio
auc_ratio <- cvauc(outcome ~ ratio, cvd)                        # ratio alone
message(sprintf("  CV-AUC full(NEFL+age+APOE4+ratio) = %.3f | same model WITHOUT ratio = %.3f | increment from ratio = %+.3f",
                auc_full, auc_base, auc_full - auc_base))
message(sprintf("  CV-AUC ratio-only = %.3f  (~0.5 = chance => ratio is NOT an independent predictor)",
                auc_ratio))

# ------------------------------------------------------------------------------
# 9. Conversion analysis - baseline-MCI, pre-specified log-ratio (Cox)
# ------------------------------------------------------------------------------
message(">>> [9/9] Baseline-MCI -> AD conversion (Cox)...")

baseline_dx <- DXSUM %>%
  mutate(RID = as.numeric(as.character(RID)), DX_val = .[[dx_col]],
         EXAMDATE = as.Date(EXAMDATE)) %>%
  filter(!is.na(DX_val)) %>% arrange(RID, EXAMDATE) %>%
  group_by(RID) %>% slice(1) %>% ungroup() %>%
  transmute(RID, base_AD = DX_val %in% c(3, "AD", "Dementia"),
            base_MCI = DX_val %in% c(2, "MCI", "LMCI", "EMCI"))

surv <- DXSUM %>%
  mutate(RID = as.numeric(as.character(RID)), DX_val = .[[dx_col]],
         EXAMDATE = as.Date(EXAMDATE)) %>%
  group_by(RID) %>%
  summarise(Start = min(EXAMDATE), Last = max(EXAMDATE),
            Event = as.integer(any(DX_val %in% c(3, "AD", "Dementia"))),
            AD_Date = suppressWarnings(min(EXAMDATE[DX_val %in% c(3, "AD", "Dementia")], na.rm = TRUE)),
            .groups = "drop") %>%
  mutate(Time = ifelse(Event == 1 & is.finite(as.numeric(AD_Date)),
                       as.numeric(AD_Date - Start) / 30.44,
                       as.numeric(Last - Start) / 30.44)) %>%
  filter(is.finite(Time), Time > 0)

cox_in <- ana %>%
  mutate(ratio = l_ptgds - l_lcn2) %>%
  select(RID, PTGDS_z = l_ptgds, ratio) %>%
  inner_join(baseline_dx, by = "RID") %>%
  filter(base_MCI, !base_AD) %>%           # baseline MCI only (no prevalent AD)
  inner_join(surv, by = "RID")

if (nrow(cox_in) > 30) {
  message(sprintf("  baseline-MCI n = %d, converters = %d", nrow(cox_in), sum(cox_in$Event)))
  for (v in c("PTGDS_z", "ratio")) {
    cm <- coxph(Surv(Time, Event) ~ scale(cox_in[[v]]), data = cox_in)
    hr <- exp(coef(cm)); ci <- exp(confint(cm)); p <- summary(cm)$coef[1, "Pr(>|z|)"]
    message(sprintf("  Cox %s/SD: HR = %.3f (%.3f-%.3f), p = %.3f%s",
                    v, hr, ci[1], ci[2], p, ifelse(p < 0.05, "", "  (n.s.)")))
  }
} else message("  insufficient baseline-MCI sample for Cox.")

message("\n", strrep("=", 65))
message("  ADNI summary: NEFL positive control + CSF associations; no prediction claim.")
message("  Fig. 3 concordance panels rendered by Fig3_CSF_concordance.R")
message(strrep("=", 65))
message(">>> Script 03 complete.")
