# ==============================================================================
# 02b_SEAAD_supplementary_tables.R
#
# Purpose : Generate supplementary Tables S2, S3, S5 from SEA-AD data
# Input   : output/SEAAD_processed_data.RData  (from Script 01)
#
# Outputs → Paper tables
#   output/TableS2_Spearman_Correlations.csv  →  Additional file 1: Table S2
#   output/TableS3_SLC22A17_Expression.csv    →  Additional file 1: Table S3
#   output/TableS5_1_Proportions.csv          →  Additional file 1: Table S5 (proportions)
#   output/TableS5_2_SST_Trajectories.csv     →  Additional file 1: Table S5 (SST+ markers)
#
# Methods
#   Table S2 : Spearman correlation coefficients (rho) for key astrocytic genes
#              against the continuous pseudo-progression score (CPS)
#   Table S3 : Constitutive expression rate of LCN2 receptor SLC22A17 in
#              excitatory neurons across pseudo-progression bins
#   Table S5 : Proportional changes of excitatory neurons and SST+ inhibitory
#              interneurons across pseudo-progression bins
#
# Paper reference
#   Kim Y†, Heo W†, Park S†, et al. (2025). bioRxiv.
#   https://github.com/YoungOukKim/MCI-to-AD
#
# Usage
#   Rscript analysis/Fig1_SEA-AD/02b_SEAAD_supplementary_tables.R
#
# Requirements
#   R >= 4.3.2  |  data.table
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

set.seed(42)

# ── Paths (relative to repository root) ───────────────────────────────────────
OUT_DIR <- "output"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ══════════════════════════════════════════════════════════════════════════════
# 0. Load data
# ══════════════════════════════════════════════════════════════════════════════
message(">>> Loading data...")
load(file.path(OUT_DIR, "SEAAD_processed_data.RData"))

# ══════════════════════════════════════════════════════════════════════════════
# 1. TABLE S2 — Spearman correlations of key astrocytic genes with CPS
# ══════════════════════════════════════════════════════════════════════════════
message(">>> Generating Table S2 (Spearman correlations with CPS)...")

astro_df <- df[cell_type_label == "Astrocyte" & !is.na(bin)]

s2_gene_list <- list(
  list(gene = "BCL2",    category = "Anti-apoptotic"),
  list(gene = "BCL2L1",  category = "Anti-apoptotic"),
  list(gene = "MCL1",    category = "Anti-apoptotic"),
  list(gene = "LCN2",    category = "Inflection marker"),
  list(gene = "PTGDS",   category = "Inflection marker"),
  list(gene = "AQP4",    category = "Astrocyte identity"),
  list(gene = "GFAP",    category = "Astrocyte identity"),
  list(gene = "SLC1A2",  category = "Astrocyte identity"),
  list(gene = "BAX",     category = "Pro-apoptotic"),
  list(gene = "CASP3",   category = "Pro-apoptotic"),
  list(gene = "CASP9",   category = "Pro-apoptotic")
)

# Warn about missing genes (need to be added to Script 01 target_genes)
available_mask <- sapply(s2_gene_list, function(x) x$gene %in% colnames(astro_df))
if (!all(available_mask)) {
  message("  Missing genes (add to target_genes in Script 01 and re-run):")
  message("    ", paste(sapply(s2_gene_list[!available_mask], `[[`, "gene"),
                        collapse = ", "))
}

s2_rows <- lapply(s2_gene_list[available_mask], function(item) {
  g  <- item$gene
  ct <- cor.test(astro_df[[g]], astro_df$bin, method = "spearman", exact = FALSE)
  data.table(
    Gene      = g,
    rho       = round(ct$estimate, 4),
    p.value   = signif(ct$p.value, 3),
    Direction = dplyr::case_when(
      ct$p.value > 0.05  ~ "Stable (ns)",
      ct$estimate > 0    ~ "Increase",
      TRUE               ~ "Decrease"
    ),
    Category  = item$category
  )
})

table_s2 <- rbindlist(s2_rows)
message("  Table S2 (n = ", nrow(astro_df), " astrocytes):")
print(table_s2)
write.csv(table_s2, file.path(OUT_DIR, "TableS2_Spearman_Correlations.csv"),
          row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# 2. TABLE S3 — SLC22A17 (LCN2 receptor) constitutive expression in neurons
# ══════════════════════════════════════════════════════════════════════════════
message(">>> Generating Table S3 (SLC22A17 expression in excitatory neurons)...")

if ("SLC22A17" %in% colnames(df)) {
  neuron_df <- df[
    cell_type_label == "Exc_Neuron" &
      !is.na(bin_rounded) & bin_rounded >= 0.1 & bin_rounded <= 0.9]

  table_s3 <- neuron_df[, .(
    n_cells         = .N,
    expressing_cells = sum(SLC22A17 > 0),
    expression_rate  = round(sum(SLC22A17 > 0) / .N * 100, 1)
  ), by = .(bin = bin_rounded)][order(bin)]

  message("  Table S3:")
  print(table_s3)
  message(sprintf("  Overall mean expression rate: %.1f%%",
                  mean(table_s3$expression_rate)))

  write.csv(table_s3, file.path(OUT_DIR, "TableS3_SLC22A17_Expression.csv"),
            row.names = FALSE)
} else {
  message("  SLC22A17 not found — add to target_genes in Script 01 and re-run")
}

# ══════════════════════════════════════════════════════════════════════════════
# 3. TABLE S5 — Neuronal subtype proportions and SST+ trajectories
# ══════════════════════════════════════════════════════════════════════════════
message(">>> Generating Table S5 (neuronal subtypes and SST+ trajectories)...")

# 3-a. Proportions ─────────────────────────────────────────────────────────────
all_neurons <- df[
  cell_type_label == "Exc_Neuron" |
    grepl("Sst|SST", cell_type, ignore.case = TRUE)]
all_neurons[, bin_r := round(bin, 1)]
all_neurons <- all_neurons[!is.na(bin_r) & bin_r >= 0.1 & bin_r <= 0.9]

all_neurons[, neuron_class := fcase(
  cell_type_label == "Exc_Neuron",                          "Excitatory_Neuron",
  grepl("Sst|SST", cell_type, ignore.case = TRUE),          "SST_Inhibitory",
  default = "Other_Neuron"
)]

# Denominator: all neuron subtypes per bin
all_subtypes <- df[
  cell_type_label == "Exc_Neuron" |
    grepl("^L[0-9]|IT|ET|CT|NP|L6b|Sst|SST|Pvalb|Vip|Lamp5|Sncg",
          cell_type, ignore.case = TRUE)]
all_subtypes[, bin_r := round(bin, 1)]
totals_per_bin <- all_subtypes[
  !is.na(bin_r) & bin_r >= 0.1 & bin_r <= 0.9,
  .(total_bin = .N), by = bin_r]

s5_counts <- all_neurons[
  neuron_class %in% c("Excitatory_Neuron", "SST_Inhibitory"),
  .N, by = .(bin = bin_r, cell_class = neuron_class)]

s5_1 <- merge(s5_counts, totals_per_bin, by.x = "bin", by.y = "bin_r")
s5_1[, proportion_pct := round(N / total_bin * 100, 2)]
setnames(s5_1, "N", "n_subtype")

message("  Table S5-1 (proportions):")
print(dcast(s5_1, bin ~ cell_class, value.var = "proportion_pct"))
write.csv(s5_1, file.path(OUT_DIR, "TableS5_1_Proportions.csv"), row.names = FALSE)

# 3-b. SST+ molecular trajectories ────────────────────────────────────────────
sst_cells <- df[
  grepl("Sst|SST", cell_type, ignore.case = TRUE) &
    !is.na(bin_rounded) & bin_rounded >= 0.1 & bin_rounded <= 0.9]

if (nrow(sst_cells) > 0L) {
  sst_markers <- intersect(
    c("PTGDS", "LCN2", "NGFR", "NDUFS1", "SST", "GPX4", "CASP3", "BAX", "BCL2"),
    colnames(sst_cells))

  s5_2_raw <- sst_cells[,
    c(lapply(.SD, mean, na.rm = TRUE), list(n_cells = .N)),
    by  = .(bin = bin_rounded),
    .SDcols = sst_markers][order(bin)]

  apply_3bin_ma <- function(x) {
    n  <- length(x)
    ma <- numeric(n)
    for (i in seq_len(n)) {
      idx   <- max(1L, i - 1L):min(n, i + 1L)
      ma[i] <- mean(x[idx], na.rm = TRUE)
    }
    ma
  }

  ma_names <- paste0(sst_markers, "_ma")
  s5_2 <- copy(s5_2_raw)
  s5_2[, (ma_names) := lapply(.SD, apply_3bin_ma), .SDcols = sst_markers]

  message("  Table S5-2 (SST+ trajectories):")
  print(s5_2[, c("bin", "n_cells", ma_names), with = FALSE])
  write.csv(s5_2, file.path(OUT_DIR, "TableS5_2_SST_Trajectories.csv"),
            row.names = FALSE)
} else {
  message("  No SST+ cells found — check cell_type annotations in h5ad")
}

# ══════════════════════════════════════════════════════════════════════════════
# 4. Summary
# ══════════════════════════════════════════════════════════════════════════════
message("\n", strrep("=", 65))
message("  Output files  →  Paper tables")
message(strrep("=", 65))
message("  TableS2_Spearman_Correlations.csv  →  Additional file 1: Table S2")
message("  TableS3_SLC22A17_Expression.csv    →  Additional file 1: Table S3")
message("  TableS5_1_Proportions.csv          →  Additional file 1: Table S5 (proportions)")
message("  TableS5_2_SST_Trajectories.csv     →  Additional file 1: Table S5 (SST+ markers)")
message(strrep("=", 65))
message(">>> Script 02b complete.")
