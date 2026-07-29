# ==============================================================================
# 01_SEAAD_data_extraction.R
#
# Purpose : Extract gene expression data from the SEA-AD h5ad file (~1.38M nuclei)
# Output  : output/SEAAD_processed_data.RData  (input for Script 02)
#
# Data source
#   SEA-AD atlas — Gabitto et al. (2024), Nature Neuroscience
#   Download : https://portal.brain-map.org/atlases-and-data/rnaseq/sea-ad
#   File     : SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad
#
#   (!)  Raw data are subject to the Allen Institute data use agreement.
#        Users must download independently; see data/NOTE_restricted_data.md
#
# Paper reference
#   Methods — "SEA-AD snRNA-seq trajectory analysis"
#   Kim Y, Heo W, Park S, et al. (2026). Research Square. https://doi.org/10.21203/rs.3.rs-9499795
#   https://github.com/YoungOukKim/MCI-to-AD
#
# Usage
#   Rscript analysis/Fig1_SEA-AD/01_SEAAD_data_extraction.R
#
# Requirements
#   R >= 4.3.2  |  rhdf5, data.table
#
# SPEED NOTE
#   Section 3 now uses the fast block read: per cell, one vectorized match()
#   over the target-gene index vector, instead of a per-gene which() loop.
#   Output (expr_mat) is identical; all ~1.38M cells are still extracted because
#   downstream tables/figures use multiple cell types.
# ==============================================================================

suppressPackageStartupMessages({
  library(rhdf5)
  library(data.table)
})

set.seed(42)

# -- Paths (relative to repository root) ---------------------------------------
H5AD_PATH <- Sys.getenv("H5AD_PATH",
                        unset = "path/to/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad")
OUT_DIR   <- Sys.getenv("OUT_DIR", unset = "output")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# -- 1. Define target genes ----------------------------------------------------
target_genes <- c(
  "NDUFS1", "PTGDS", "LCN2", "NGFR", "APOE",
  "PTGS2", "IL6ST", "NFKBIA",
  "P2RY1", "P2RY12", "GJA1", "ITPR2",
  "HMOX1", "SOD2", "MT1E", "MT2A", "CLU", "SLC1A2",
  "FOS", "JUN", "EGR1", "ATF3",
  "TREM2", "C3",
  "CHI3L1", "SERPINA3", "GFAP",
  "BCL2", "BCL2L1", "MCL1", "BAX", "CASP3", "CASP9",
  "AQP4",
  "SLC22A17",
  "SST",
  "GPX4", "FTH1", "SLC7A11"
)

# -- 2. Map gene names to h5ad column indices ----------------------------------
message(">>> [1/5] Mapping gene indices...")
gene_names_all <- h5read(H5AD_PATH, "var/_index")
gene_idx_map   <- match(target_genes, gene_names_all) - 1L   # 0-based (== X/indices)
valid_mask     <- !is.na(gene_idx_map)
if (any(!valid_mask)) {
  message("  Genes not found in h5ad (skipped): ",
          paste(target_genes[!valid_mask], collapse = ", "))
}
target_genes <- target_genes[valid_mask]
gene_idx_map <- gene_idx_map[valid_mask]
message("  Extracting ", length(target_genes), " genes from ", H5AD_PATH)

# -- 3. Extract expression matrix (sparse CSC) — FAST block read ---------------
#   Per cell: a single match(gene_idx_map, cell_indices) maps all target genes
#   at once (vectorized), instead of one which() scan per gene.
message(">>> [2/5] Extracting expression from ~1.38M cells (fast block read)...")
indptr  <- h5read(H5AD_PATH, "X/indptr", bit64conversion = "double")
n_cells <- length(indptr) - 1L
expr_mat <- matrix(0, nrow = n_cells, ncol = length(gene_idx_map),
                   dimnames = list(NULL, target_genes))

block_size <- 200000L
for (s0 in seq(1L, n_cells, by = block_size)) {
  e0 <- min(s0 + block_size - 1L, n_cells)
  sp <- indptr[s0]
  cnt <- indptr[e0 + 1L] - sp
  if (cnt <= 0) next

  ci <- h5read(H5AD_PATH, "X/indices",
               start = sp + 1L, count = cnt, bit64conversion = "double")
  cd <- h5read(H5AD_PATH, "X/data",
               start = sp + 1L, count = cnt)

  for (i in s0:e0) {
    a <- indptr[i] - sp + 1L
    b <- indptr[i + 1L] - sp
    if (b < a) next
    h  <- match(gene_idx_map, ci[a:b])   # vectorized over genes
    ok <- !is.na(h)
    if (any(ok)) expr_mat[i, ok] <- cd[a:b][h[ok]]
  }
  message(sprintf("  Progress: %5.1f%%  (%d / %d cells)",
                  e0 / n_cells * 100, e0, n_cells))
}

# -- 4. Attach cell metadata ---------------------------------------------------
message(">>> [3/5] Attaching cell metadata...")
subclass_idx  <- h5read(H5AD_PATH, "obs/Subclass")
subclass_cats <- h5read(H5AD_PATH, "obs/__categories/Subclass")
cps_vec       <- as.numeric(
  h5read(H5AD_PATH, "obs/Continuous Pseudo-progression Score"))

# donor ID (required for donor-level statistics)
.read_obs_cat <- function(key) {
  v <- h5read(H5AD_PATH, paste0("obs/", key))
  cats <- tryCatch(h5read(H5AD_PATH, paste0("obs/__categories/", key)), error = function(e) NULL)
  if (!is.null(cats)) cats[as.integer(v) + 1L] else v
}
.obs_names <- { ll <- h5ls(H5AD_PATH); ll$name[ll$group == "/obs"] }
.donor_key <- .obs_names[grepl("donor", .obs_names, ignore.case = TRUE)][1]
donor_vec  <- if (!is.na(.donor_key)) as.character(.read_obs_cat(.donor_key)) else NA_character_
message("  donor field = ", ifelse(is.na(.donor_key), "<not found>", .donor_key))

df <- as.data.table(expr_mat)
df[, cell_type := subclass_cats[subclass_idx + 1L]]
df[, bin       := cps_vec]
df[, donor     := donor_vec]   # donor-level grouping
df <- df[!is.na(bin)]

# -- 5. Classify cell types ----------------------------------------------------
message(">>> [4/5] Classifying cell types...")
df[, cell_type_label := "Other"]
df[grepl("^L[0-9]|IT|ET|CT|NP|L6b", cell_type), cell_type_label := "Exc_Neuron"]
df[cell_type == "Astrocyte",                      cell_type_label := "Astrocyte"]
df[grepl("Micro", cell_type),                     cell_type_label := "Microglia"]
df[grepl("^Sst", cell_type, ignore.case = TRUE),  cell_type_label := "SST_Neuron"]   # SST interneurons (Table S4/S5)
df[, bin_rounded := round(bin, 1)]

message("\nCell type counts:")
print(df[, .N, by = cell_type_label])

for (g in c("LCN2", "NGFR", "APOE", "TREM2", "C3", "SLC22A17")) {
  if (g %in% colnames(df)) {
    nz <- sum(df[[g]] > 0)
    message(sprintf("  %-12s  %d non-zero cells  (%.2f%%)",
                    g, nz, nz / nrow(df) * 100))
  }
}

# -- 6. Save -------------------------------------------------------------------
message(">>> [5/5] Saving to ", file.path(OUT_DIR, "SEAAD_processed_data.RData"))
save(df, file = file.path(OUT_DIR, "SEAAD_processed_data.RData"))

bin_dist <- df[cell_type_label %in% c("Exc_Neuron", "Astrocyte", "Microglia"),
               .N, by = .(bin_rounded, cell_type_label)]
message("\nBin x Cell type distribution:")
print(dcast(bin_dist, bin_rounded ~ cell_type_label, value.var = "N"))

message("\n>>> Script 01 complete.")
message("    Next: Rscript analysis/Fig1_SEA-AD/02_SEAAD_trajectory_analysis.R")
