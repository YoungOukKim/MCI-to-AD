# ==============================================================================
# 02_SEAAD_trajectory_analysis.R
#
# Purpose : Generate all trajectory tables from SEA-AD pseudo-progression data
# Input   : output/SEAAD_processed_data.RData  (from Script 01)
#
# Outputs → Paper tables
#   output/Table1_3binMA_Full.csv          → Table 1 (main text, full)
#   output/Table1_Paper_Summary.csv        → Table 1 (Bins 0.1 / 0.6 / 0.9)
#   output/Table2_Checkpoint_CCF.csv       → Table 2 (main text, inflection window)
#   output/TableS1_Raw_PerBin.csv          → Additional file 1: Table S1
#   output/TableS4_FullTrajectory_CCF.csv  → Additional file 1: Table S4
#
# Methods
#   Gene expression trajectories were smoothed using a 3-bin centered moving
#   average (window = 3 bins). Lagged cross-correlation functions (CCF) were
#   computed with R's ccf(); lag.max = 2 for the inflection window (Bins 0.4–0.8)
#   and lag.max = 3 for the full 9-bin trajectory. set.seed(42), R v4.3.2.
#
# Paper reference
#   Kim Y†, Heo W†, Park S†, et al. (2025). bioRxiv.
#   https://github.com/YoungOukKim/MCI-to-AD
#
# Usage
#   Rscript analysis/Fig1_SEA-AD/02_SEAAD_trajectory_analysis.R
#
# Requirements
#   R >= 4.3.2  |  data.table, dplyr
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

set.seed(42)

# ── Paths (relative to repository root) ───────────────────────────────────────
OUT_DIR <- "output"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ══════════════════════════════════════════════════════════════════════════════
# 0. Load and validate
# ══════════════════════════════════════════════════════════════════════════════
message(">>> Loading data...")
load(file.path(OUT_DIR, "SEAAD_processed_data.RData"))

essential <- c("NDUFS1", "PTGDS", "LCN2", "NGFR", "APOE", "TREM2", "C3",
               "PTGS2", "IL6ST", "NFKBIA", "P2RY1", "P2RY12", "GJA1", "ITPR2")
stopifnot(all(essential %in% colnames(df)))
stopifnot("Microglia" %in% unique(df$cell_type_label))

analysis_genes <- intersect(
  c(essential, "HMOX1", "SOD2", "MT1E", "MT2A", "CLU", "SLC1A2",
    "FOS", "JUN", "EGR1", "ATF3"),
  colnames(df)
)
message("  Genes available for analysis: ", length(analysis_genes))

# ══════════════════════════════════════════════════════════════════════════════
# 1. Bin-level means and 3-bin centered moving average
# ══════════════════════════════════════════════════════════════════════════════
message(">>> Computing bin means and 3-bin moving average...")

target_groups <- c("Exc_Neuron", "Astrocyte", "Microglia")

raw_means <- df[
  cell_type_label %in% target_groups &
    !is.na(bin_rounded) & bin_rounded >= 0.1 & bin_rounded <= 0.9,
  c(lapply(.SD, mean, na.rm = TRUE), list(n_cells = .N)),
  by  = .(bin = bin_rounded, group = cell_type_label),
  .SDcols = analysis_genes
][order(group, bin)]

# 3-bin centered moving average (edge bins use 2-bin average)
apply_3bin_ma <- function(x) {
  n  <- length(x)
  ma <- numeric(n)
  for (i in seq_len(n)) {
    idx   <- max(1L, i - 1L):min(n, i + 1L)
    ma[i] <- mean(x[idx], na.rm = TRUE)
  }
  ma
}

ma_cols  <- paste0(analysis_genes, "_ma")
smoothed <- copy(raw_means)
smoothed[, (ma_cols) := lapply(.SD, apply_3bin_ma),
         by = group, .SDcols = analysis_genes]

# ══════════════════════════════════════════════════════════════════════════════
# 2. TABLE S1 — Raw bin-resolved expression
# ══════════════════════════════════════════════════════════════════════════════
message(">>> Generating Table S1 (raw per-bin expression)...")

table_s1 <- raw_means[, .(bin, group, n_cells,
                           NDUFS1, PTGDS, LCN2, NGFR, APOE)]
write.csv(table_s1, file.path(OUT_DIR, "TableS1_Raw_PerBin.csv"), row.names = FALSE)

astro_raw <- table_s1[group == "Astrocyte"]
message(sprintf("  PTGDS raw peak: Bin %.1f = %.4f",
                astro_raw$bin[which.max(astro_raw$PTGDS)],
                max(astro_raw$PTGDS)))
message(sprintf("  Total astrocytes: %d", sum(astro_raw$n_cells)))

# ══════════════════════════════════════════════════════════════════════════════
# 3. TABLE 1 — 3-bin MA at key bins (0.1, 0.6, 0.9)
# ══════════════════════════════════════════════════════════════════════════════
message(">>> Generating Table 1 (3-bin MA trajectories)...")

a_s <- smoothed[group == "Astrocyte"][order(bin)]
n_s <- smoothed[group == "Exc_Neuron"][order(bin)]
m_s <- smoothed[group == "Microglia"][order(bin)]

table1_full <- data.table(
  Bin           = a_s$bin,
  NDUFS1_Neuron = n_s$NDUFS1_ma,
  PTGDS_Astro   = a_s$PTGDS_ma,
  LCN2_Astro    = a_s$LCN2_ma,
  NGFR_Neuron   = n_s$NGFR_ma,
  APOE_Micro    = m_s$APOE_ma
)

table1_paper <- table1_full[Bin %in% c(0.1, 0.6, 0.9)]
message("  Table 1 (paper summary at Bins 0.1 / 0.6 / 0.9):")
print(table1_paper)

write.csv(table1_full,  file.path(OUT_DIR, "Table1_3binMA_Full.csv"),    row.names = FALSE)
write.csv(table1_paper, file.path(OUT_DIR, "Table1_Paper_Summary.csv"),  row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# 4. CCF utility
# ══════════════════════════════════════════════════════════════════════════════
ccf_analysis <- function(x, y, max_lag = 3L) {
  ok <- complete.cases(x, y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n < 4L)              return(data.table(Lag = NA_integer_, r = NA_real_,
                                              p_val = NA_real_, n_eff = NA_integer_))
  if (sd(x) == 0 || sd(y) == 0)
                           return(data.table(Lag = 0L, r = 0, p_val = 1, n_eff = n))

  res     <- ccf(x, y, lag.max = max_lag, plot = FALSE, na.action = na.pass)
  idx     <- which.max(abs(res$acf))
  lag_val <- as.integer(res$lag[idx])
  r_val   <- as.numeric(res$acf[idx])
  n_eff   <- n - abs(lag_val)

  p_val <- if (n_eff > 2L) {
    t_stat <- r_val * sqrt((n_eff - 2L) / (1 - r_val^2 + 1e-10))
    2 * pt(-abs(t_stat), df = n_eff - 2L)
  } else NA_real_

  data.table(Lag = lag_val, r = round(r_val, 3),
             p_val = round(p_val, 4), n_eff = n_eff)
}

sig_label <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", "ns"))))
}

# ══════════════════════════════════════════════════════════════════════════════
# 5. Module composite scores (NF-κB and Purinergic)
# ══════════════════════════════════════════════════════════════════════════════
a_s[, NFkB_ma  := rowMeans(.SD), .SDcols = c("PTGS2_ma", "IL6ST_ma", "NFKBIA_ma")]
a_s[, Purin_ma := rowMeans(.SD), .SDcols = c("P2RY1_ma", "P2RY12_ma",
                                              "GJA1_ma", "ITPR2_ma")]

# ══════════════════════════════════════════════════════════════════════════════
# 6. TABLE S4 — Full trajectory CCF (9 bins, lag.max = 3)
# ══════════════════════════════════════════════════════════════════════════════
message(">>> Generating Table S4 (full trajectory CCF, 9 bins, lag.max = 3)...")

pair_labels <- c(
  "NDUFS1 (Neuron) -> PTGDS (Astrocyte)",
  "PTGDS (Astrocyte) -> LCN2 (Astrocyte)",
  "LCN2 (Astrocyte) -> NGFR (Neuron)",
  "NF-kB module (Astrocyte) <-> PTGDS (Astrocyte)",
  "Purinergic module (Astrocyte) <-> PTGDS (Astrocyte)",
  "LCN2 (Astrocyte) -> TREM2 (Microglia)",
  "C3 (Microglia) <-> NGFR (Neuron)",
  "LCN2 (Astrocyte) <-> NDUFS1 (Neuron)"
)
cell_labels <- c(
  "Neuron-Astrocyte", "Astrocyte-Astrocyte", "Astrocyte-Neuron",
  "Astrocyte-Astrocyte", "Astrocyte-Astrocyte",
  "Astrocyte-Microglia", "Microglia-Neuron", "Astrocyte-Neuron"
)

res_s4 <- rbindlist(list(
  ccf_analysis(n_s$NDUFS1_ma, a_s$PTGDS_ma),
  ccf_analysis(a_s$PTGDS_ma,  a_s$LCN2_ma),
  ccf_analysis(a_s$LCN2_ma,   n_s$NGFR_ma),
  ccf_analysis(a_s$NFkB_ma,   a_s$PTGDS_ma),
  ccf_analysis(a_s$Purin_ma,  a_s$PTGDS_ma),
  ccf_analysis(a_s$LCN2_ma,   m_s$TREM2_ma),
  ccf_analysis(m_s$C3_ma,     n_s$NGFR_ma),
  ccf_analysis(a_s$LCN2_ma,   n_s$NDUFS1_ma)
))
res_s4[, Pair      := pair_labels]
res_s4[, Cell_Type := cell_labels]
res_s4[, Sig       := sig_label(p_val)]
setcolorder(res_s4, c("Pair", "Cell_Type", "Lag", "r", "p_val", "Sig", "n_eff"))

message("  Table S4:")
print(res_s4)
write.csv(res_s4, file.path(OUT_DIR, "TableS4_FullTrajectory_CCF.csv"),
          row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# 7. TABLE 2 — Inflection window CCF (Bins 0.4–0.8, lag.max = 2)
#    Methods: "Lagged cross-correlations within the inflection window (Bins 0.4–0.8)"
# ══════════════════════════════════════════════════════════════════════════════
message(">>> Generating Table 2 (inflection window CCF, Bins 0.4–0.8, lag.max = 2)...")

cp_bins <- c(0.4, 0.5, 0.6, 0.7, 0.8)
a_cp <- a_s[bin %in% cp_bins]
n_cp <- n_s[bin %in% cp_bins]
m_cp <- m_s[bin %in% cp_bins]

res_t2 <- rbindlist(list(
  ccf_analysis(a_cp$Purin_ma,  a_cp$PTGDS_ma,  max_lag = 2L),
  ccf_analysis(a_cp$NFkB_ma,   a_cp$PTGDS_ma,  max_lag = 2L),
  ccf_analysis(a_cp$PTGDS_ma,  a_cp$LCN2_ma,   max_lag = 2L),
  ccf_analysis(a_cp$LCN2_ma,   m_cp$TREM2_ma,  max_lag = 2L),
  ccf_analysis(n_cp$NDUFS1_ma, a_cp$PTGDS_ma,  max_lag = 2L),
  ccf_analysis(a_cp$PTGDS_ma,  a_cp$NFkB_ma,   max_lag = 2L),   # PPARG proxy
  ccf_analysis(m_cp$C3_ma,     n_cp$NGFR_ma,   max_lag = 2L),
  ccf_analysis(a_cp$LCN2_ma,   n_cp$NDUFS1_ma, max_lag = 2L)
))
res_t2[, Pair := c(
  "Purinergic (A) -> PTGDS (A)",
  "NF-kB (A) -> PTGDS (A)",
  "PTGDS (A) -> LCN2 (A)",
  "LCN2 (A) -> TREM2 (M)",
  "NDUFS1 (N) -> PTGDS (A)",
  "PTGDS (A) -> PPARG module (A)",
  "C3 (M) -> NGFR (N)",
  "LCN2 (A) -> NDUFS1 (N)"
)]
res_t2[, Cell_Type := c(
  "Astrocyte-Astrocyte", "Astrocyte-Astrocyte", "Astrocyte-Astrocyte",
  "Astrocyte-Microglia", "Neuron-Astrocyte",    "Astrocyte-Astrocyte",
  "Microglia-Neuron",    "Astrocyte-Neuron"
)]
res_t2[, Sig := sig_label(p_val)]
setcolorder(res_t2, c("Pair", "Cell_Type", "Lag", "r", "p_val", "Sig", "n_eff"))

message("  Table 2 (inflection window, Bins 0.4–0.8):")
print(res_t2)
write.csv(res_t2, file.path(OUT_DIR, "Table2_Checkpoint_CCF.csv"), row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# 8. Summary
# ══════════════════════════════════════════════════════════════════════════════
message("\n", strrep("=", 65))
message("  Output files  →  Paper tables")
message(strrep("=", 65))
message("  TableS1_Raw_PerBin.csv          →  Additional file 1: Table S1")
message("  Table1_3binMA_Full.csv          →  Table 1 (full, 9 bins)")
message("  Table1_Paper_Summary.csv        →  Table 1 (Bins 0.1 / 0.6 / 0.9)")
message("  Table2_Checkpoint_CCF.csv       →  Table 2")
message("  TableS4_FullTrajectory_CCF.csv  →  Additional file 1: Table S4")
message(strrep("=", 65))
message(">>> Script 02 complete.")
message("    Next: Rscript analysis/Fig1_SEA-AD/02b_SEAAD_supplementary_tables.R")
