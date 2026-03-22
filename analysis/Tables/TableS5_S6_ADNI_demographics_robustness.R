# ==============================================================================
# TableS5_S6_ADNI_demographics_robustness.R
#
# Purpose : Generate Additional file 1: Table S5 and Table S6
#
# Outputs → Paper tables
#   output/tables/TableS5_ADNI_Demographics.csv   → Additional file 1: Table S5
#     ADNI participant demographics by diagnostic group
#     (CN=172 / EMCI=183 / LMCI=234 / AD=146; total N=735)
#     Columns: N | Age (mean±SD) | Female% | Education | APOE4% | MMSE
#
#   output/tables/TableS6_Robustness_AUC.csv      → Additional file 1: Table S6
#     Robustness of Phase Boundary Integrity Index across MMSE thresholds
#     Full-cohort AUC ~ 0.741 (MCI vs AD); Table S6 range 0.56-0.72
#
#   output/figures/Fig5E_ROC.png                  → Fig. 5E (ROC curve)
#   output/figures/Fig5F_KM_survival.png          → Fig. 5F (KM curve)
#
# NOTE on variable naming in original code:
#   The original script (Supp_Table_S5-S6-R-code.txt) used internal variable
#   names "table_s6" for demographics and "table_s7" for robustness — one
#   number offset from the actual paper table numbers. This corrected version
#   uses table_s5 / table_s6 to match the paper.
#
# Key results reproduced
#   Table S5: CN=172, EMCI=183, LMCI=234, AD=146 (total N=735)
#   Table S6: AUC range 0.56-0.72 across MMSE <= 24 to MMSE <= 30 thresholds
#   Fig. 5E: Full-cohort AUC ~ 0.741 (MCI vs AD)
#   Fig. 5F: HR ~ 2.89 (75th percentile risk stratification)
#
# Input
#   ADNI_PATH env var : folder with DXSUM.rda, PTDEMOG.rda, APOERES.rda,
#                       UPENNBIOMK_MASTER.rda, MMSE.rda, ARM.rda
#   SOMA_CSV  env var : CruchagaLab SomaScan 7K CSV
#
# ADNI data: register at https://adni.loni.usc.edu
# Paper: Kim Y et al. (2025). bioRxiv. https://github.com/YoungOukKim/MCI-to-AD
# Usage: Rscript analysis/Tables/TableS5_S6_ADNI_demographics_robustness.R
# Requirements: R >= 4.3.2 | tidyverse, pROC, survival, survminer, lubridate
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(pROC)
  library(survival)
  library(survminer)
  library(lubridate)
})
set.seed(42)

# Paths (relative to repository root)
ADNI_PATH <- Sys.getenv("ADNI_PATH", unset = "path/to/ADNIMERGE2/data/")
SOMA_CSV  <- Sys.getenv("SOMA_CSV",
  unset = "path/to/CruchagaLab_CSF_SOMAscan7k.csv")
OUT_TABLES  <- file.path("output", "tables")
OUT_FIGURES <- file.path("output", "figures")
dir.create(OUT_TABLES,  showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_FIGURES, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. Load data
# ==============================================================================
message(">>> [1/9] Loading data...")
somascan <- read.csv(SOMA_CSV, stringsAsFactors = FALSE)

for (f in c("DXSUM","PTDEMOG","APOERES","UPENNBIOMK_MASTER","MMSE","ARM")) {
  path <- file.path(ADNI_PATH, paste0(f, ".rda"))
  if (!file.exists(path)) stop("File not found: ", path,
                                "\n  Register and download at: https://adni.loni.usc.edu")
  load(path)
}

# ==============================================================================
# 2. SomaScan proteins (baseline, 1 row per participant)
#    PTGDS = X10514.5 | LCN2 = X2836.68 | NEFL = X10082.251 | TAU = X5854.60
# ==============================================================================
message(">>> [2/9] Extracting SomaScan proteins...")
protein_df <- somascan %>%
  mutate(RID = as.numeric(as.character(RID))) %>%
  group_by(RID) %>% slice(1) %>% ungroup() %>%
  transmute(RID,
    PTGDS = as.numeric(X10514.5),
    LCN2  = as.numeric(X2836.68),
    NEFL  = as.numeric(X10082.251),
    TAU   = as.numeric(X5854.60))
soma_rids <- unique(protein_df$RID)
message(sprintf("  SomaScan baseline participants: %d", length(soma_rids)))

# ==============================================================================
# 3. Diagnosis (3-group: CN/MCI/AD; 4-group: CN/EMCI/LMCI/AD)
#    Uses first non-missing DIAGNOSIS from DXSUM — matches AUC 0.741 logic
# ==============================================================================
message(">>> [3/9] Processing diagnosis (AUC-matched method)...")

dx_col <- grep("DIAGNOSIS|DXCURREN|DXCHANGE", colnames(DXSUM), value=TRUE)[1]

clean_dx <- DXSUM %>%
  mutate(RID=as.numeric(as.character(RID)), DX_val=.[[dx_col]]) %>%
  filter(!is.na(DX_val)) %>%
  arrange(RID, EXAMDATE) %>%
  group_by(RID) %>% slice(1) %>% ungroup() %>%
  mutate(DX_3group = case_when(
    DX_val %in% c(3,"AD","Dementia")      ~ "AD",
    DX_val %in% c(2,"MCI","LMCI","EMCI") ~ "MCI",
    DX_val %in% c(1,"CN","Normal","NL")  ~ "CN",
    TRUE ~ NA_character_)) %>%
  filter(!is.na(DX_3group)) %>%
  select(RID, DX_3group)

# EMCI/LMCI subdivision from ARM
arm_dx <- ARM %>%
  mutate(RID=as.numeric(as.character(RID)),
         DX_ARM=case_when(
           str_detect(ARM,"EMCI")          ~ "EMCI",
           str_detect(ARM,"LMCI")          ~ "LMCI",
           str_detect(ARM,"AD")            ~ "AD",
           str_detect(ARM,"NL|Normal|SMC") ~ "CN",
           str_detect(ARM,"MCI")           ~ "LMCI",
           TRUE ~ NA_character_)) %>%
  filter(!is.na(DX_ARM)) %>%
  distinct(RID, .keep_all=TRUE) %>%
  select(RID, DX_ARM)

clean_dx <- clean_dx %>%
  left_join(arm_dx, by="RID") %>%
  mutate(DX_4group = case_when(
    DX_3group=="AD"  ~ "AD",
    DX_3group=="CN"  ~ "CN",
    DX_3group=="MCI" & DX_ARM=="EMCI" ~ "EMCI",
    TRUE ~ "LMCI")) %>%
  filter(RID %in% soma_rids)

# ==============================================================================
# 4. Demographics
# ==============================================================================
message(">>> [4/9] Processing demographics...")

clean_demog <- PTDEMOG %>%
  mutate(RID=as.numeric(as.character(RID))) %>%
  filter(RID %in% soma_rids) %>%
  group_by(RID) %>% slice(1) %>% ungroup() %>%
  mutate(
    AGE = year(as.Date(VISDATE)) - PTDOBYY,
    SEX = ifelse(grepl("Female|2", PTGENDER, ignore.case=TRUE), "Female","Male"),
    EDU = as.numeric(PTEDUCAT)) %>%
  filter(!is.na(AGE), AGE>40, AGE<120) %>%
  select(RID, AGE, SEX, EDU)

# ==============================================================================
# 5. APOE, CSF biomarkers, MMSE
# ==============================================================================
message(">>> [5/9] Processing APOE, CSF, MMSE...")

clean_apoe <- APOERES %>%
  mutate(RID=as.numeric(as.character(RID))) %>%
  filter(RID %in% soma_rids) %>%
  group_by(RID) %>% slice(1) %>% ungroup() %>%
  mutate(APOE4=as.integer(grepl("4", GENOTYPE))) %>%
  select(RID, APOE4)

clean_csf <- UPENNBIOMK_MASTER %>%
  mutate(RID=as.numeric(as.character(RID))) %>%
  filter(grepl("bl|sc", VISCODE, ignore.case=TRUE), RID %in% soma_rids) %>%
  mutate(ABETA_num=as.numeric(gsub("[<>]","",ABETA)),
         PTAU_num =as.numeric(gsub("[<>]","",PTAU))) %>%
  group_by(RID) %>% slice(1) %>% ungroup() %>%
  select(RID, ABETA_num, PTAU_num)

clean_mmse <- MMSE %>%
  mutate(RID=as.numeric(as.character(RID)),
         MMSCORE=as.numeric(as.character(MMSCORE))) %>%
  filter(grepl("^bl$|^sc$", VISCODE, ignore.case=TRUE),
         RID %in% soma_rids, !is.na(MMSCORE), MMSCORE>0, MMSCORE<=30) %>%
  group_by(RID) %>% summarise(MMSE_score=mean(MMSCORE), .groups="drop")

# ==============================================================================
# 6. Merge full cohort
# ==============================================================================
message(">>> [6/9] Merging cohort...")

full_cohort <- protein_df %>%
  inner_join(clean_dx,    by="RID") %>%
  left_join(clean_demog,  by="RID") %>%
  left_join(clean_apoe,   by="RID") %>%
  left_join(clean_csf,    by="RID") %>%
  left_join(clean_mmse,   by="RID") %>%
  filter(!is.na(AGE))

message(sprintf("  Full cohort N = %d", nrow(full_cohort)))
print(table(full_cohort$DX_4group))

# ==============================================================================
# 7. TABLE S5 — ADNI participant demographics
#    Paper: Additional file 1: Table S5
#    Matches Supp header: "ADNI participant demographics"
#    Expected: CN=172, EMCI=183, LMCI=234, AD=146
# ==============================================================================
message(">>> [7/9] Generating Table S5 (demographics)...")

table_s5 <- full_cohort %>%
  mutate(DX_4group=factor(DX_4group, levels=c("CN","EMCI","LMCI","AD"))) %>%
  group_by(DX_4group) %>%
  summarise(
    N              = n(),
    Age_mean       = round(mean(AGE),1),
    Age_SD         = round(sd(AGE),1),
    Female_pct     = round(mean(SEX=="Female")*100,1),
    Education_mean = round(mean(EDU,na.rm=TRUE),1),
    Education_SD   = round(sd(EDU,na.rm=TRUE),1),
    APOE4_pct      = round(mean(APOE4,na.rm=TRUE)*100,1),
    MMSE_mean      = round(mean(MMSE_score,na.rm=TRUE),1),
    MMSE_SD        = round(sd(MMSE_score,na.rm=TRUE),1),
    .groups="drop") %>%
  arrange(DX_4group)

cat("\n  TABLE S5. ADNI Participant Demographics\n")
print(as.data.frame(table_s5))
cat(sprintf("  Total N = %d\n", sum(table_s5$N)))

# Verification
expected <- c(CN=172, EMCI=183, LMCI=234, AD=146)
actual   <- setNames(table_s5$N, table_s5$DX_4group)
for (g in names(expected)) {
  if (g %in% names(actual)) {
    diff <- abs(actual[g] - expected[g])
    if (diff > 5) message(sprintf(
      "  WARNING: %s N=%d, expected ~%d (diff=%d). Check AUC-matched DX logic.",
      g, actual[g], expected[g], diff))
  }
}

write.csv(table_s5,
          file.path(OUT_TABLES, "TableS5_ADNI_Demographics.csv"),
          row.names=FALSE)

# ==============================================================================
# 8. ROC analysis — Full-cohort AUC (Checkpoint Integrity Index)
#    Expected: AUC ~ 0.741
# ==============================================================================
message(">>> [8/9] ROC analysis (MCI vs AD)...")

roc_df <- full_cohort %>%
  filter(DX_3group %in% c("MCI","AD")) %>%
  mutate(
    outcome   = as.integer(DX_3group=="AD"),
    log_ratio = log2(PTGDS / (LCN2 + 1e-6)),
    log_nefl  = log2(NEFL  + 1e-6),
    log_tau   = log2(TAU   + 1e-6),
    inter_act = log_ratio * APOE4) %>%
  drop_na(outcome, log_ratio, log_nefl, log_tau, AGE, APOE4, ABETA_num, PTAU_num)

final_model <- glm(
  outcome ~ log_ratio + log_nefl + log_tau + AGE + APOE4 +
    ABETA_num + PTAU_num + inter_act,
  data=roc_df, family=binomial)

roc_df$checkpoint_index <- predict(final_model, type="response")
roc_obj <- roc(roc_df$outcome, roc_df$checkpoint_index, ci=TRUE, quiet=TRUE)
auc_val <- as.numeric(auc(roc_obj))

message(sprintf("  Full-cohort AUC = %.4f (95%% CI: %.3f-%.3f)",
                auc_val, ci(roc_obj)[1], ci(roc_obj)[3]))
if (abs(auc_val - 0.741) > 0.015)
  message("  WARNING: AUC deviates from expected 0.741. Check DX coding.")

# ROC plot → Fig. 5E
png(file.path(OUT_FIGURES,"Fig5E_ROC.png"), width=7, height=7, units="in", res=300)
plot(roc_obj, col="#0072B2", lwd=2, main="Phase Boundary Integrity Index (MCI vs AD)")
legend("bottomright",
       legend=sprintf("AUC = %.3f (%.3f-%.3f)",
                      auc_val, ci(roc_obj)[1], ci(roc_obj)[3]),
       col="#0072B2", lwd=2, bty="n")
dev.off()
message("  Saved: Fig5E_ROC.png")

# ==============================================================================
# TABLE S6 — Robustness analysis across MMSE thresholds
#    Paper: Additional file 1: Table S6
#    Matches Supp header: "Robustness analysis of the Phase Boundary
#    Integrity Index across varying MMSE inclusion thresholds"
#    Expected: AUC range 0.56-0.72
# ==============================================================================
message("  Generating Table S6 (robustness across MMSE thresholds)...")

roc_mmse     <- roc_df %>% filter(!is.na(MMSE_score))
thresholds   <- c(24, 25, 26, 27, 28, 30)
table_s6_rows <- list()

for (th in thresholds) {
  sub   <- roc_mmse %>% filter(MMSE_score <= th)
  n_mci <- sum(sub$outcome == 0)
  n_ad  <- sum(sub$outcome == 1)
  if (n_ad >= 10 && n_mci >= 10) {
    r_obj <- roc(sub$outcome, sub$checkpoint_index, quiet=TRUE)
    ci_r  <- ci.auc(r_obj)
    table_s6_rows[[as.character(th)]] <- data.frame(
      MMSE_Threshold  = paste0("MMSE <= ", th),
      N               = nrow(sub),
      N_MCI           = n_mci,
      N_AD            = n_ad,
      AUC             = round(as.numeric(auc(r_obj)), 3),
      CI_95_Lower     = round(ci_r[1], 3),
      CI_95_Upper     = round(ci_r[3], 3))
  } else {
    message(sprintf("  Skip MMSE <= %d: MCI=%d, AD=%d (< 10 each)", th, n_mci, n_ad))
  }
}

table_s6 <- do.call(rbind, table_s6_rows)
cat("\n  TABLE S6. Robustness Analysis (Phase Boundary Integrity Index)\n")
print(table_s6)

write.csv(table_s6,
          file.path(OUT_TABLES, "TableS6_Robustness_AUC.csv"),
          row.names=FALSE)

# ==============================================================================
# 9. Survival analysis — Cox + Kaplan-Meier (Fig. 5F)
#    Expected: HR ~ 2.89 (75th percentile threshold)
# ==============================================================================
message(">>> [9/9] Survival analysis (Cox + Kaplan-Meier)...")

surv_timeline <- DXSUM %>%
  mutate(RID=as.numeric(as.character(RID)),
         DX_val=.[[dx_col]], EXAMDATE=as.Date(EXAMDATE)) %>%
  filter(RID %in% roc_df$RID) %>%
  group_by(RID) %>%
  summarise(
    Start_Date = min(EXAMDATE, na.rm=TRUE),
    Last_Date  = max(EXAMDATE, na.rm=TRUE),
    Event      = as.integer(any(DX_val %in% c(3,"Dementia","AD"))),
    AD_Date    = suppressWarnings(
      min(EXAMDATE[DX_val %in% c(3,"Dementia","AD")], na.rm=TRUE)),
    .groups="drop") %>%
  mutate(Time_months = ifelse(
    Event==1 & !is.infinite(as.numeric(AD_Date)),
    as.numeric(difftime(AD_Date,   Start_Date, units="days"))/30.44,
    as.numeric(difftime(Last_Date, Start_Date, units="days"))/30.44)) %>%
  filter(!is.na(Time_months), Time_months > 0)

cox_df <- roc_df %>% inner_join(surv_timeline, by="RID")
threshold_75 <- quantile(cox_df$checkpoint_index, 0.75)
cox_df$Risk_Group <- factor(
  ifelse(cox_df$checkpoint_index > threshold_75, "High_Risk","Low_Risk"),
  levels=c("Low_Risk","High_Risk"))

cox_model <- coxph(Surv(Time_months, Event) ~ Risk_Group, data=cox_df)
hr    <- exp(coef(cox_model))
hr_ci <- exp(confint(cox_model))

message(sprintf("  Cox HR = %.3f (95%% CI: %.3f-%.3f)", hr, hr_ci[1], hr_ci[2]))
if (abs(hr - 2.89) > 0.3)
  message("  WARNING: HR deviates from expected 2.89. Check risk group cutoff.")

# KM plot → Fig. 5F
km_fit  <- survfit(Surv(Time_months, Event) ~ Risk_Group, data=cox_df)
km_plot <- ggsurvplot(
  km_fit, data=cox_df,
  pval=TRUE, conf.int=TRUE, risk.table=TRUE,
  palette=c("#377EB8","#E41A1C"),
  xlab="Time (months)", ylab="MCI-to-AD conversion-free probability",
  title="Phase Boundary Integrity Index — MCI-to-AD Conversion")

ggsave(file.path(OUT_FIGURES,"Fig5F_KM_survival.png"),
       km_plot$plot, width=10, height=7, dpi=300, bg="white")
message("  Saved: Fig5F_KM_survival.png")

# ==============================================================================
# 10. Summary
# ==============================================================================
cat("\n", strrep("=", 65), "\n")
cat("  ADNI Analysis Summary\n")
cat(strrep("=", 65), "\n")
cat(sprintf("  Cohort: N = %d (CN=%d, EMCI=%d, LMCI=%d, AD=%d)\n",
            sum(table_s5$N),
            table_s5$N[table_s5$DX_4group=="CN"],
            table_s5$N[table_s5$DX_4group=="EMCI"],
            table_s5$N[table_s5$DX_4group=="LMCI"],
            table_s5$N[table_s5$DX_4group=="AD"]))
cat(sprintf("  Full-cohort AUC = %.4f\n", auc_val))
cat(sprintf("  Cox HR = %.3f (95%% CI: %.3f-%.3f)\n", hr, hr_ci[1], hr_ci[2]))
cat(strrep("=", 65), "\n")
cat("  Output files:\n")
cat("    TableS5_ADNI_Demographics.csv  ->  Additional file 1: Table S5\n")
cat("    TableS6_Robustness_AUC.csv     ->  Additional file 1: Table S6\n")
cat("    Fig5E_ROC.png                  ->  Fig. 5E\n")
cat("    Fig5F_KM_survival.png          ->  Fig. 5F\n")
cat(strrep("=", 65), "\n")
message(">>> Script TableS5_S6 complete.")
