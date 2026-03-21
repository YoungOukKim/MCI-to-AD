# ==============================================================================
# 01_SEAAD_data_extraction.R
#
# Purpose : Extract gene expression data from the SEA-AD h5ad file (1.3M nuclei)
# Output  : output/SEAAD_processed_data.RData  (input for Script 02)
#
# Data source
#   SEA-AD atlas — Gabitto et al. (2024), Nature Neuroscience
#   Download : https://portal.brain-map.org/atlases-and-data/rnaseq/sea-ad
#   File     : SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad
#
#   ⚠  Raw data are subject to the Allen Institute data use agreement.
#      Users must download independently; see data/NOTE_restricted_data.md
#
# Paper reference
#   Methods — "SEA-AD snRNA-seq trajectory analysis"
#   Kim Y†, Heo W†, Park S†, et al. (2025). bioRxiv.
#   https://github.com/YoungOukKim/MCI-to-AD
#
# Usage
#   Rscript analysis/Fig1_SEA-AD/01_SEAAD_data_extraction.R
#
# Requirements
#   R >= 4.3.2  |  rhdf5, data.table
# ==============================================================================

suppressPackageStartupMessages({
  library(rhdf5)
  library(data.table)
})

set.seed(42)

# ── Paths (relative to repository root) ───────────────────────────────────────
# Set H5AD_PATH to the location of your downloaded SEA-AD h5ad file.
# All other paths are relative and require no editing.
H5AD_PATH <- Sys.getenv("H5AD_PATH",
  unset = "path/to/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad")
OUT_DIR   <- "output"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ── 1. Define target genes ────────────────────────────────────────────────────
# All genes used across Tables 1, 2, S1–S4 and module analyses
target_genes <- c(
  # Core cascade markers (Table 1, Table 2, Fig. 1)
  "NDUFS1", "PTGDS", "LCN2", "NGFR", "APOE",
  # NF-κB module (Table 2, Fig. S1)
  "PTGS2", "IL6ST", "NFKBIA",
  # Purinergic / Ca2+ module (Table 2, Fig. S1)
  "P2RY1", "P2RY12", "GJA1", "ITPR2",
  # Metabolic buffering module (Fig. S1)
  "HMOX1", "SOD2", "MT1E", "MT2A", "CLU", "SLC1A2",
  # Immediate early gene module (Fig. S1)
  "FOS", "JUN", "EGR1", "ATF3",
  # Microglia markers (Table 2)
  "TREM2", "C3",
  # Reactive astrocyte markers (Table S2)
  "CHI3L1", "SERPINA3", "GFAP",
  # Apoptosis / identity markers (Table S2)
  "BCL2", "BCL2L1", "MCL1", "BAX", "CASP3", "CASP9",
  "AQP4",
  # LCN2 receptor in neurons (Table S3)
  "SLC22A17",
  # SST interneuron marker (Table S5)
  "SST",
  # Ferroptosis protection (Fig. S1)
  "GPX4", "FTH1", "SLC7A11"
)

# ── 2. Map gene names to h5ad column indices ──────────────────────────────────
message(">>> [1/5] Mapping gene indices...")
gene_names_all <- h5read(H5AD_PATH, "var/_index")
gene_idx_map   <- match(target_genes, gene_names_all) - 1L   # 0-based
valid_mask     <- !is.na(gene_idx_map)

if (any(!valid_mask)) {
  message("  Genes not found in h5ad (skipped): ",
          paste(target_genes[!valid_mask], collapse = ", "))
}
target_genes <- target_genes[valid_mask]
gene_idx_map <- gene_idx_map[valid_mask]
message("  Extracting ", length(target_genes), " genes from ", H5AD_PATH)

# ── 3. Extract expression matrix (sparse CSC format) ─────────────────────────
message(">>> [2/5] Extracting expression from 1.3M cells...")
indptr  <- h5read(H5AD_PATH, "X/indptr", bit64conversion = "double")
n_cells <- length(indptr) - 1L
expr_mat <- matrix(0, nrow = n_cells, ncol = length(target_genes),
                   dimnames = list(NULL, target_genes))

block_size <- 200000L
for (start_cell in seq(1L, n_cells, by = block_size)) {
  end_cell  <- min(start_cell + block_size - 1L, n_cells)
  start_pos <- indptr[start_cell]
  end_pos   <- indptr[end_cell + 1L]

  if (end_pos > start_pos) {
    indices <- h5read(H5AD_PATH, "X/indices",
                      start = start_pos + 1L, count = end_pos - start_pos,
                      bit64conversion = "double")
    data    <- h5read(H5AD_PATH, "X/data",
                      start = start_pos + 1L, count = end_pos - start_pos)
    offset  <- 1L
    for (i in start_cell:end_cell) {
      nnz <- indptr[i + 1L] - indptr[i]
      if (nnz > 0L) {
        c_idx <- indices[offset:(offset + nnz - 1L)]
        c_dat <- data[offset:(offset + nnz - 1L)]
        for (j in seq_along(gene_idx_map)) {
          found <- which(c_idx == gene_idx_map[j])
          if (length(found) > 0L) expr_mat[i, j] <- c_dat[found[1L]]
        }
        offset <- offset + nnz
      }
    }
  }
  message(sprintf("  Progress: %5.1f%%  (%d / %d cells)",
                  end_cell / n_cells * 100, end_cell, n_cells))
}

# ── 4. Attach cell metadata ───────────────────────────────────────────────────
message(">>> [3/5] Attaching cell metadata...")
subclass_idx  <- h5read(H5AD_PATH, "obs/Subclass")
subclass_cats <- h5read(H5AD_PATH, "obs/__categories/Subclass")
cps_vec       <- as.numeric(
  h5read(H5AD_PATH, "obs/Continuous Pseudo-progression Score"))

df <- as.data.table(expr_mat)
df[, cell_type := subclass_cats[subclass_idx + 1L]]
df[, bin       := cps_vec]
df <- df[!is.na(bin)]

# ── 5. Classify cell types ────────────────────────────────────────────────────
message(">>> [4/5] Classifying cell types...")
df[, cell_type_label := "Other"]
df[grepl("^L[0-9]|IT|ET|CT|NP|L6b", cell_type), cell_type_label := "Exc_Neuron"]
df[cell_type == "Astrocyte",                      cell_type_label := "Astrocyte"]
df[grepl("Micro", cell_type),                     cell_type_label := "Microglia"]
df[, bin_rounded := round(bin, 1)]

# Report
message("\nCell type counts:")
print(df[, .N, by = cell_type_label])

# Sparse gene detection rates
for (g in c("LCN2", "NGFR", "APOE", "TREM2", "C3", "SLC22A17")) {
  if (g %in% colnames(df)) {
    nz <- sum(df[[g]] > 0)
    message(sprintf("  %-12s  %d non-zero cells  (%.2f%%)",
                    g, nz, nz / nrow(df) * 100))
  }
}

# ── 6. Save ──────────────────────────────────────────────────────────────────
message(">>> [5/5] Saving to ", file.path(OUT_DIR, "SEAAD_processed_data.RData"))
save(df, file = file.path(OUT_DIR, "SEAAD_processed_data.RData"))

# Verify bin × cell-type distribution (corresponds to Table S1 cell counts)
bin_dist <- df[cell_type_label %in% c("Exc_Neuron", "Astrocyte", "Microglia"),
               .N, by = .(bin_rounded, cell_type_label)]
message("\nBin × Cell type distribution:")
print(dcast(bin_dist, bin_rounded ~ cell_type_label, value.var = "N"))

message("\n>>> Script 01 complete.")
message("    Next: Rscript analysis/Fig1_SEA-AD/02_SEAAD_trajectory_analysis.R")
