# ==============================================================================
# TableS7_QC_Braak_mapping.R
#
# Purpose : Generate Additional file 1: Table S7
#           Quality control metrics and bin-to-Braak stage mapping
#
# Outputs (3 CSV files, all → Table S7 in paper)
#   output/tables/TableS7_A_Bin_CellType_Counts.csv
#     Part A: Bin-resolved cell type counts (Astrocyte / Exc.Neuron /
#             SST+Neuron / Microglia) across CPS Bins 0.1-0.9
#   output/tables/TableS7_B_QC_Criteria.csv
#     Part B: QC parameters applied during SEA-AD processing
#             (min genes, mito%, doublet removal, atlas version, etc.)
#   output/tables/TableS7_C_Braak_Summary.csv
#     Part C: Heuristic CPS-to-Braak stage mapping
#             Bin 0.1-0.3 ~ Braak I-II | 0.4-0.6 ~ III-IV | 0.7-0.9 ~ V-VI
#
# Matches Supp caption (Table S7):
#   "A. Bin-resolved cell type counts across pseudo-progression.
#    B. QC metrics. C. CPS-to-Braak stage mapping."
#
# Methods reference:
#   "Quality control metrics and bin-to-Braak stage mapping are detailed
#    in Additional file 1: Table S7."
#   "CPS bins were approximately aligned with Braak stages. This mapping
#    is heuristic and does not imply direct pathological equivalence."
#
# Input   : output/SEAAD_processed_data.RData  (from 01_SEAAD_data_extraction.R)
#           H5AD_PATH env var (optional, for reading live QC metrics)
# Output  : output/tables/TableS7_*.csv
#
# Paper: Kim Y et al. (2025). bioRxiv. https://github.com/YoungOukKim/MCI-to-AD
# Usage: Rscript analysis/Tables/TableS7_QC_Braak_mapping.R
# Requirements: R >= 4.3.2 | data.table, rhdf5 (optional for live QC)
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

set.seed(42)

# Paths (relative to repository root)
OUT_DIR  <- file.path("output", "tables")
H5AD_PATH <- Sys.getenv("H5AD_PATH",
  unset = "path/to/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. Load processed data
# ==============================================================================
message(">>> Loading SEAAD_processed_data.RData...")
load(file.path("output", "SEAAD_processed_data.RData"))
message(sprintf("  Total cells: %s", format(nrow(df), big.mark=",")))
message(sprintf("  Cell types: %s", paste(unique(df$cell_type_label), collapse=", ")))

# ==============================================================================
# 2. Part A: Bin-resolved cell type counts
#    Matches Table S7 columns: CPS_Bin | Astrocyte | Exc.Neuron |
#                               SST+Neuron | Microglia | Subtotal | Total | Braak
# ==============================================================================
message(">>> Part A: Bin-resolved cell type counts...")

target_types <- c("Astrocyte", "Exc_Neuron", "SST_Neuron", "Microglia")

counts_long <- df[
  cell_type_label %in% target_types &
    bin_rounded >= 0.1 & bin_rounded <= 0.9,
  .N, by = .(bin_rounded, cell_type_label)]

counts_wide <- dcast(counts_long,
                     bin_rounded ~ cell_type_label,
                     value.var = "N", fill = 0L)
counts_wide <- as.data.table(counts_wide)
setnames(counts_wide, "bin_rounded", "CPS_Bin")

# Ensure all target type columns exist
for (ct in target_types) {
  if (!ct %in% names(counts_wide)) counts_wide[, (ct) := 0L]
}
setcolorder(counts_wide, c("CPS_Bin", target_types))

# Subtotal (4 target types) and total (all cells)
counts_wide[, Subtotal_4types := rowSums(.SD), .SDcols = target_types]
all_per_bin <- df[bin_rounded >= 0.1 & bin_rounded <= 0.9,
                  .(Total_All = .N), by = bin_rounded]
counts_wide <- merge(counts_wide, all_per_bin,
                     by.x = "CPS_Bin", by.y = "bin_rounded", all.x = TRUE)

# Heuristic CPS-to-Braak mapping
# Note: tilde (~) = approximate; CPS is transcriptomic, not neuropathological
counts_wide[, Approx_Braak := fcase(
  CPS_Bin <= 0.3, "~I-II (Preclinical)",
  CPS_Bin <= 0.6, "~III-IV (Prodromal MCI)",
  CPS_Bin <= 0.9, "~V-VI (Dementia)"
)]

message("  Table S7-A:")
print(counts_wide)

write.csv(counts_wide,
          file.path(OUT_DIR, "TableS7_A_Bin_CellType_Counts.csv"),
          row.names = FALSE)

# ==============================================================================
# 3. Part B: QC criteria
#    Attempt to read live QC metrics from h5ad; fall back to documented values
# ==============================================================================
message("\n>>> Part B: QC metrics...")

qc_metrics_live <- list()

if (file.exists(H5AD_PATH)) {
  tryCatch({
    library(rhdf5)
    obs_keys   <- h5read(H5AD_PATH, "obs")
    obs_names  <- names(obs_keys)
    gene_field <- intersect(c("n_genes","n_genes_by_counts","nFeature_RNA"),
                             obs_names)[1]
    mito_field <- intersect(c("pct_counts_mt","percent.mt","percent_mito"),
                             obs_names)[1]
    if (!is.na(gene_field)) {
      ng <- as.numeric(h5read(H5AD_PATH, paste0("obs/", gene_field)))
      qc_metrics_live$median_genes <- median(ng)
      qc_metrics_live$pct_gte500   <- round(mean(ng >= 500)*100, 1)
    }
    if (!is.na(mito_field)) {
      mt <- as.numeric(h5read(H5AD_PATH, paste0("obs/", mito_field)))
      qc_metrics_live$median_mito  <- round(median(mt), 2)
      qc_metrics_live$pct_lt20mito <- round(mean(mt < 20)*100, 1)
    }
    h5closeAll()
    message("  Live QC metrics read from h5ad.")
  }, error = function(e) {
    message("  Could not read h5ad QC fields (", e$message, ").")
    message("  Using documented values from SEA-AD publication.")
  })
} else {
  message("  H5AD_PATH not set. Using documented QC values.")
}

qc_table <- data.table(
  Parameter = c(
    "Minimum genes per nucleus",
    "Maximum mitochondrial read fraction",
    "Doublet removal method",
    "Atlas version (file)",
    "Brain region",
    "Total donors",
    "Total nuclei (post-QC)",
    "Scrublet reference",
    "CPS range"
  ),
  Criterion = c(
    ">= 500 genes",
    "< 20%",
    "Scrublet-based [Wolock et al., 2019]",
    "SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad",
    "Middle temporal gyrus (MTG)",
    "84 donors",
    format(nrow(df), big.mark = ","),
    "Wolock SL et al. (2019) Cell Syst 8(4):329-337",
    "0.0 - 1.0 (continuous pseudo-progression score)"
  ),
  Note = c(
    if (!is.null(qc_metrics_live$pct_gte500))
      sprintf("%.1f%% of cells pass", qc_metrics_live$pct_gte500) else
      "Per SEA-AD QC protocol",
    if (!is.null(qc_metrics_live$pct_lt20mito))
      sprintf("%.1f%% of cells pass", qc_metrics_live$pct_lt20mito) else
      "Per SEA-AD QC protocol",
    "Scrublet score threshold = 0.25",
    "Download: portal.brain-map.org/atlases-and-data/rnaseq/sea-ad",
    "Gabitto et al. (2024) Nat Neurosci",
    "84 MTG donors, ages 44-90",
    "After QC, doublet removal, and CPS assignment",
    "https://doi.org/10.1016/j.cels.2019.02.002",
    "Rounded to 1 decimal for 9 discrete bins (0.1-0.9)"
  )
)

message("  Table S7-B (QC criteria):")
print(qc_table)
write.csv(qc_table,
          file.path(OUT_DIR, "TableS7_B_QC_Criteria.csv"),
          row.names = FALSE)

# ==============================================================================
# 4. Part C: Braak stage mapping summary
# ==============================================================================
message("\n>>> Part C: Braak stage summary...")

braak_summary <- counts_wide[, .(
  N_Bins       = .N,
  Total_Cells  = sum(Total_All, na.rm=TRUE),
  Astrocyte    = sum(Astrocyte),
  Exc_Neuron   = sum(Exc_Neuron),
  SST_Neuron   = sum(SST_Neuron),
  Microglia    = sum(Microglia)
), by = Approx_Braak][order(Approx_Braak)]

braak_summary[, Note := c(
  "Preclinical to early symptomatic; CPS 0.1-0.3",
  "Prodromal MCI spectrum; CPS 0.4-0.6",
  "Dementia; CPS 0.7-0.9"
)]

message("  Table S7-C (Braak summary):")
print(braak_summary)
write.csv(braak_summary,
          file.path(OUT_DIR, "TableS7_C_Braak_Summary.csv"),
          row.names = FALSE)

# ==============================================================================
# 5. Summary
# ==============================================================================
cat("\n", strrep("=", 65), "\n")
cat("  Table S7 — QC metrics and bin-to-Braak stage mapping\n")
cat(strrep("=", 65), "\n")
cat("\n  Part A: Bin-resolved cell type counts\n")
for (i in seq_len(nrow(counts_wide))) {
  row <- counts_wide[i]
  cat(sprintf("  Bin %.1f | Astro: %6s | Exc: %7s | SST: %5s | Micro: %5s | %s\n",
              row$CPS_Bin,
              format(row$Astrocyte,  big.mark=","),
              format(row$Exc_Neuron, big.mark=","),
              format(row$SST_Neuron, big.mark=","),
              format(row$Microglia,  big.mark=","),
              row$Approx_Braak))
}
cat(sprintf("\n  Total (4 types): %s\n",
            format(sum(counts_wide$Subtotal_4types), big.mark=",")))
cat(sprintf("  Total (all):     %s\n",
            format(sum(counts_wide$Total_All, na.rm=TRUE), big.mark=",")))
cat("\n  Part C: CPS-to-Braak heuristic mapping\n")
cat("  Bin 0.1-0.3  ->  ~Braak I-II   (Preclinical)\n")
cat("  Bin 0.4-0.6  ->  ~Braak III-IV (Prodromal MCI)\n")
cat("  Bin 0.7-0.9  ->  ~Braak V-VI   (Dementia)\n")
cat("  Note: ~ denotes approximate correspondence.\n")
cat("        CPS is transcriptomic, not neuropathological.\n")
cat(strrep("=", 65), "\n")
message(">>> Script TableS7 complete.")
