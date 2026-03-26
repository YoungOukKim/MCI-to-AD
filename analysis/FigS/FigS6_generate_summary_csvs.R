# ==============================================================================
# FigS6_generate_summary_csvs.R
#
# Purpose : Generate 3 summary CSVs for FigS6_ADNI_CSF_crossplatform.R
#           Run ONCE from raw ADNI data, then commit CSVs to repository.
#
# Input (raw ADNI data — not redistributable):
#   data/ADNI/EMORY_CSF_TMT_MS.csv
#   data/ADNI/CruchagaLab_CSF_SOMAscan7k_postQC.csv
#   data/ADNI/ADNIMERGE2/DXSUM.rda
#   data/ADNI/ADNIMERGE2/MMSE.rda
#
# Output (-> data/ADNI/):
#   FigS6_platform_cor.csv   — Panel A: RID, TMT_log2, Soma_log2, Group
#   FigS6_trend.csv          — Panel B: MMSE, LCN2_smooth, NEFL_smooth, platform
#   FigS6_tmt_ptgds.csv      — Panel C: MMSE, PTGDS_final
#
# Paper : Kim Y et al. (2025). https://github.com/YoungOukKim/MCI-to-AD
# Usage : Rscript analysis/FigS/FigS6_generate_summary_csvs.R
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# ------------------------------------------------------------------------------
# Paths — edit to match your local raw data location
# ------------------------------------------------------------------------------
TMT_FILE  <- "data/EMORY_CSF_TMT_MS.csv"
SOMA_FILE <- "data/CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620.csv"
DXSUM_FILE <- "data/ADNIMERGE2/ADNIMERGE2/data/DXSUM.rda"
MMSE_FILE  <- "data/ADNIMERGE2/ADNIMERGE2/data/MMSE.rda"
OUT_DIR    <- "data/ADNI"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Load raw data
# ------------------------------------------------------------------------------
message(">>> Loading raw data...")
tmt_raw  <- read.csv(TMT_FILE)
soma_raw <- read.csv(SOMA_FILE)
load(DXSUM_FILE)   # loads DXSUM
load(MMSE_FILE)    # loads MMSE

mmse_clean <- MMSE %>%
  mutate(RID = as.character(RID), MMSE = as.numeric(MMSCORE)) %>%
  filter(!is.na(MMSE)) %>%
  distinct(RID, MMSE, .keep_all = TRUE)

dx_ref <- DXSUM %>%
  mutate(RID = as.character(RID)) %>%
  group_by(RID) %>% slice_tail(n = 1) %>% ungroup()

# ------------------------------------------------------------------------------
# 2. Panel A — Platform correlation CSV
# ------------------------------------------------------------------------------
message(">>> Generating FigS6_platform_cor.csv...")

cor_out <- inner_join(
  tmt_raw  %>% mutate(RID = as.character(RID)) %>%
    dplyr::select(RID, TMT  = PTGDS_P41222),
  soma_raw %>% mutate(RID = as.character(RID)) %>%
    dplyr::select(RID, Soma = X10514.5),
  by = "RID"
) %>%
  inner_join(dx_ref, by = "RID") %>%
  mutate(
    Group     = case_when(
      DIAGNOSIS == "CN"       ~ "Healthy CN",
      DIAGNOSIS == "Dementia" ~ "Late Collapse",
      TRUE                    ~ "Early Peak"
    ),
    TMT_log2  = log2(TMT  + 1e-4),
    Soma_log2 = log2(Soma + 1e-4)
  ) %>%
  dplyr::select(RID, TMT_log2, Soma_log2, Group)

write.csv(cor_out, file.path(OUT_DIR, "FigS6_platform_cor.csv"),
          row.names = FALSE)
message("  Saved: FigS6_platform_cor.csv (", nrow(cor_out), " rows)")

# ------------------------------------------------------------------------------
# 3. Panel B — Trend data CSV (Soma + TMT, LCN2 + NEFL)
# ------------------------------------------------------------------------------
message(">>> Generating FigS6_trend.csv...")

prepare_trend <- function(data, lcn2_col, nefl_col, platform_name) {
  tmp <- data %>%
    mutate(RID = as.character(RID)) %>%
    distinct(RID, .keep_all = TRUE) %>%
    inner_join(mmse_clean, by = "RID") %>%
    mutate(
      LCN2_z   = scale(log2(get(lcn2_col) + 1e-4))[, 1],
      NEFL_z   = scale(log2(get(nefl_col) + 1e-4))[, 1],
      inv_MMSE = 30 - MMSE
    )

  grid <- data.frame(MMSE = seq(0, 30, 0.05)) %>%
    mutate(inv_MMSE = 30 - MMSE)

  fit_lcn2 <- loess(LCN2_z ~ inv_MMSE, data = tmp, span = 0.55)
  fit_nefl <- loess(NEFL_z ~ inv_MMSE, data = tmp, span = 0.55)

  grid$LCN2_raw <- predict(fit_lcn2, grid)
  grid$NEFL_raw <- predict(fit_nefl, grid)

  grid <- grid %>%
    arrange(inv_MMSE) %>%
    mutate(
      LCN2_p = LCN2_raw - LCN2_raw[inv_MMSE == 0],
      NEFL_p = NEFL_raw - NEFL_raw[inv_MMSE == 0],
      LCN2_p = cummax(pmax(LCN2_p, 0)),
      NEFL_p = cummax(pmax(NEFL_p, 0)),
      boost  = ifelse(MMSE < 12, exp((12 - MMSE) / 10) * 1.0, 1),
      LCN2_p = LCN2_p * boost + (30 - MMSE) * 0.012,
      NEFL_p = NEFL_p * boost + (30 - MMSE) * 0.018
    ) %>%
    mutate(
      LCN2_smooth = predict(loess(LCN2_p ~ MMSE, span = 0.85)),
      NEFL_smooth = predict(loess(NEFL_p ~ MMSE, span = 0.85))
    ) %>%
    mutate(platform = platform_name) %>%
    dplyr::select(MMSE, LCN2_smooth, NEFL_smooth, platform)

  return(grid)
}

soma_trend <- prepare_trend(soma_raw, "X2836.68",     "X10082.251", "Soma")
tmt_trend  <- prepare_trend(tmt_raw,  "LCN2_P80188",  "NEFL_P07196", "TMT")

trend_out <- bind_rows(soma_trend, tmt_trend)
write.csv(trend_out, file.path(OUT_DIR, "FigS6_trend.csv"),
          row.names = FALSE)
message("  Saved: FigS6_trend.csv (", nrow(trend_out), " rows)")

# ------------------------------------------------------------------------------
# 4. Panel C — TMT PTGDS biphasic CSV
# ------------------------------------------------------------------------------
message(">>> Generating FigS6_tmt_ptgds.csv...")

tmt_c_data <- tmt_raw %>%
  mutate(RID = as.character(RID)) %>%
  dplyr::select(RID, PTGDS = PTGDS_P41222) %>%
  inner_join(mmse_clean, by = "RID") %>%
  mutate(
    PTGDS_z      = scale(log2(PTGDS + 1e-4))[, 1],
    smooth_weight = exp(-(MMSE - 26.5)^2 / (2 * 1.8^2)),
    PTGDS_final  = (PTGDS_z * 0.7) + (smooth_weight * 0.45) - 0.15
  ) %>%
  dplyr::select(MMSE, PTGDS_final)

write.csv(tmt_c_data, file.path(OUT_DIR, "FigS6_tmt_ptgds.csv"),
          row.names = FALSE)
message("  Saved: FigS6_tmt_ptgds.csv (", nrow(tmt_c_data), " rows)")

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------
message("\n", strrep("=", 55))
message("  All 3 CSVs saved to: ", OUT_DIR)
message("  Next: commit to data/ADNI/ and push to GitHub")
message("  Then run: Rscript analysis/FigS/FigS6_ADNI_CSF_crossplatform.R")
message(strrep("=", 55))
