# ==============================================================================
# FigS1_module_trajectories.R
#
# Purpose : Additional file 1: Fig. S1
#           Module-level mechanistic dissection of the astrocyte-neuron cascade
#           in SEA-AD pseudo-progression
#
# Panels
#   A  Purinergic/Ca2+ module vs PTGDS  (Lag -1, r = -0.886)
#   B  NF-kB Priming module vs PTGDS    (Lag  0, r = -0.678)
#   C  IEG Stress module vs PTGDS       (Lag  0, concurrent)
#   D  Ferroptosis-Protective module vs NDUFS1 (r = -0.857)
#   E  Timeline schematic of module cascade within inflection window
#
# Matches Supp caption:
#   "Module-level trajectories of key astrocyte programs relative to PTGDS
#    across pseudo-progression bins (3-bin moving average, z-score normalized).
#    The inflection window (Bins 0.4-0.8) is highlighted. Lagged cross-
#    correlation values are indicated in boxed annotations."
#
# Input   : SEA-AD h5ad file (set H5AD_PATH env var)
# Output  : output/FigS1/FigS1.png  +  .tiff  (600 dpi, 14 x 11.5 in)
#
# Paper: Kim Y et al. (2025). bioRxiv. https://github.com/YoungOukKim/MCI-to-AD
# Usage: Rscript analysis/FigS/FigS1_module_trajectories.R
# Requirements: R >= 4.3.2 | rhdf5, dplyr, tidyr, ggplot2, patchwork
# ==============================================================================

suppressPackageStartupMessages({
  library(rhdf5)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})
set.seed(42)

# Paths (relative to repository root)
H5AD_PATH <- Sys.getenv("H5AD_PATH",
  unset = "path/to/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad")
OUT_DIR <- file.path("output", "FigS1")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Cell type classifier
classify_cell <- function(label) {
  label <- as.character(label)
  if (grepl("Sst", label, ignore.case = TRUE))          return("SST_Inhibitory")
  if (grepl("^L[0-9]|IT$|ET$|CT$|NP$|L6b", label))     return("Excitatory_Neuron")
  if (grepl("^Astro", label, ignore.case = TRUE))        return("Astrocyte")
  if (grepl("^Micro", label, ignore.case = TRUE))        return("Microglia")
  return("Other")
}

# ==============================================================================
# 1. Load SEA-AD h5ad
# ==============================================================================
message(">>> [1/4] Loading SEA-AD data...")
h5closeAll()

subclass_idx        <- h5read(H5AD_PATH, "/obs/Subclass")
subclass_categories <- h5read(H5AD_PATH, "/obs/__categories/Subclass")
cell_type_labels    <- subclass_categories[subclass_idx + 1]
bin_vec             <- as.numeric(
  h5read(H5AD_PATH, "/obs/Continuous Pseudo-progression Score"))
cell_class          <- sapply(cell_type_labels, classify_cell)

gene_names_all <- h5read(H5AD_PATH, "var/_index")
module_genes <- c(
  "P2RY1", "P2RY12", "GJA1", "ITPR2",            # Purinergic/Ca2+
  "PTGS2", "IL6ST", "NFKBIA",                      # NF-kB priming
  "FOS", "JUN", "EGR1", "ATF3",                    # IEG stress
  "HMOX1", "SOD2", "MT1E", "MT2A", "CLU", "SLC1A2", # Metabolic buffering
  "GPX4", "FTH1", "SLC7A11",                       # Ferroptosis protection
  "PTGDS", "LCN2", "NDUFS1"                        # Core cascade markers
)

gene_idx_map <- match(module_genes, gene_names_all) - 1L
valid_mask   <- !is.na(gene_idx_map)
module_genes <- module_genes[valid_mask]
gene_idx_map <- gene_idx_map[valid_mask]
message("  Extracting ", length(module_genes), " genes from ", length(bin_vec), " cells")

# ==============================================================================
# 2. Extract expression (block read)
# ==============================================================================
message(">>> [2/4] Extracting expression...")
indptr   <- h5read(H5AD_PATH, "X/indptr", bit64conversion = "double")
n_cells  <- length(indptr) - 1L
expr_mat <- matrix(0, nrow = n_cells, ncol = length(module_genes),
                   dimnames = list(NULL, module_genes))

block_size <- 200000L
for (start_cell in seq(1L, n_cells, by = block_size)) {
  end_cell  <- min(start_cell + block_size - 1L, n_cells)
  start_pos <- indptr[start_cell]
  end_pos   <- indptr[end_cell + 1L]
  count_pos <- end_pos - start_pos
  if (count_pos > 0L) {
    idx_block  <- h5read(H5AD_PATH, "X/indices",
                         start = start_pos + 1L, count = count_pos,
                         bit64conversion = "double")
    data_block <- h5read(H5AD_PATH, "X/data",
                         start = start_pos + 1L, count = count_pos)
    offset <- 1L
    for (i in start_cell:end_cell) {
      nnz <- indptr[i + 1L] - indptr[i]
      if (nnz > 0L) {
        ci <- idx_block[offset:(offset + nnz - 1L)]
        cd <- data_block[offset:(offset + nnz - 1L)]
        for (j in seq_along(gene_idx_map)) {
          m <- which(ci == gene_idx_map[j])
          if (length(m) > 0L) expr_mat[i, j] <- cd[m[1L]]
        }
        offset <- offset + nnz
      }
    }
  }
  message(sprintf("  Progress: %5.1f%%", end_cell / n_cells * 100))
}

# ==============================================================================
# 3. Bin-level summary and 3-bin moving average
# ==============================================================================
message(">>> [3/4] Computing trajectories...")
df_all <- data.frame(bin = round(bin_vec, 1), cell_class, expr_mat) %>%
  filter(!is.na(bin), bin >= 0.1, bin <= 0.9)

compute_ma3 <- function(x) {
  n <- length(x); ma <- numeric(n)
  for (i in seq_len(n)) { ma[i] <- mean(x[max(1,i-1):min(n,i+1)], na.rm=TRUE) }
  ma
}
zscore <- function(x) { m <- mean(x, na.rm=TRUE); s <- sd(x, na.rm=TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x))); (x - m) / s }

astro_bin <- df_all %>%
  filter(cell_class == "Astrocyte") %>%
  group_by(bin) %>%
  summarise(across(all_of(module_genes), ~ mean(.x, na.rm=TRUE)), .groups="drop")

neuron_bin <- df_all %>%
  filter(cell_class == "Excitatory_Neuron") %>%
  group_by(bin) %>%
  summarise(NDUFS1_neuron = mean(NDUFS1, na.rm=TRUE), .groups="drop")

astro_ma <- astro_bin %>%
  arrange(bin) %>%
  mutate(across(all_of(module_genes), compute_ma3, .names = "{.col}_ma"),
         PTGDS_ma = predict(loess(PTGDS ~ bin, data=., span=1.0)))

astro_ma <- astro_ma %>%
  mutate(
    Purinergic_Ca_ma       = rowMeans(pick(P2RY1_ma, P2RY12_ma, GJA1_ma, ITPR2_ma)),
    NFkB_Priming_ma        = rowMeans(pick(PTGS2_ma, IL6ST_ma, NFKBIA_ma)),
    IEG_Stress_ma          = rowMeans(pick(FOS_ma, JUN_ma, EGR1_ma, ATF3_ma)),
    Ferroptosis_Protect_ma = rowMeans(pick(GPX4_ma, FTH1_ma, SLC7A11_ma))
  ) %>%
  left_join(neuron_bin %>%
              arrange(bin) %>%
              mutate(NDUFS1_neuron_ma = compute_ma3(NDUFS1_neuron)),
            by = "bin")

# ==============================================================================
# 4. Build panels
# ==============================================================================
message(">>> [4/4] Building panels A-E...")

inflection_rect <- annotate("rect", xmin=0.4, xmax=0.8,
                             ymin=-Inf, ymax=Inf, fill="gray", alpha=0.12)
x_scale <- scale_x_continuous(breaks=seq(0.1,0.9,0.1))
theme_s1 <- theme_minimal(base_size=14) +
  theme(plot.title=element_blank(), legend.position="top",
        legend.title=element_blank(), axis.text=element_text(size=10),
        panel.grid.minor=element_blank(), legend.text=element_text(size=10))

plot_module_vs_ptgds <- function(mod_col, mod_name, color_mod, lag_text) {
  df <- astro_ma %>%
    select(bin, mod=all_of(mod_col), PTGDS=PTGDS_ma) %>%
    mutate(mod=zscore(mod), PTGDS=zscore(PTGDS)) %>%
    pivot_longer(c("mod","PTGDS"), names_to="Module", values_to="Expression") %>%
    mutate(Module=if_else(Module=="PTGDS", "PTGDS (Astrocyte)", mod_name))
  y_hi <- max(df$Expression, na.rm=TRUE)
  ggplot(df, aes(x=bin, y=Expression, color=Module)) +
    inflection_rect +
    geom_smooth(method="loess", span=0.9, se=FALSE, linewidth=1.5) +
    geom_point(size=1.8, alpha=0.5) +
    annotate("label", x=0.85, y=y_hi*0.92, label=lag_text,
             size=4.2, fill="white", alpha=0.85, hjust=1) +
    scale_color_manual(values=c("PTGDS (Astrocyte)"="#0072B2",
                                 setNames(color_mod, mod_name))) +
    x_scale +
    labs(x="Pseudo-progression Bin",
         y="3-bin MA Expression\n(z-score)") + theme_s1
}

# Panels A-C: Module vs PTGDS
# Table 3 values: A=Lag-1 r=-0.886, B=Lag0 r=-0.678, C=Lag0 concurrent
pA <- plot_module_vs_ptgds("Purinergic_Ca_ma",  "Purinergic/Ca\u00b2\u207a Module",
                            "#E74C3C", "Lag -1, r = -0.886")
pB <- plot_module_vs_ptgds("NFkB_Priming_ma",   "NF-\u03baB Priming Module",
                            "#ff7f0e", "Lag 0, r = -0.678")
pC <- plot_module_vs_ptgds("IEG_Stress_ma",     "IEG Stress Module",
                            "#2ca02c", "Lag 0 concurrent")

# Panel D: Ferroptosis vs NDUFS1 (r = -0.857, Table 3)
pD_df <- astro_ma %>%
  select(bin, Ferroptosis=Ferroptosis_Protect_ma, NDUFS1=NDUFS1_neuron_ma) %>%
  mutate(Ferroptosis=zscore(Ferroptosis), NDUFS1=zscore(NDUFS1)) %>%
  pivot_longer(c("Ferroptosis","NDUFS1"), names_to="Module", values_to="Expression") %>%
  mutate(Module=recode(Module,
    "Ferroptosis"="Ferroptosis-Protective Module",
    "NDUFS1"="NDUFS1 (Neuron)"))

pD <- ggplot(pD_df, aes(x=bin, y=Expression, color=Module)) +
  inflection_rect +
  geom_smooth(method="loess", span=0.9, se=FALSE, linewidth=1.5) +
  geom_point(size=1.8, alpha=0.5) +
  annotate("label", x=0.85, y=max(pD_df$Expression,na.rm=TRUE)*0.85,
           label="Inverse sync, r = -0.857", size=4.2, fill="white", alpha=0.85, hjust=1) +
  scale_color_manual(values=c(
    "Ferroptosis-Protective Module"="#9467bd", "NDUFS1 (Neuron)"="#2C3E50")) +
  x_scale +
  labs(x="Pseudo-progression Bin", y="3-bin MA Expression\n(z-score)") + theme_s1

# Panel E: Timeline schematic — Inflection Window
pE <- ggplot() +
  geom_segment(aes(x=0.1,xend=0.9,y=0,yend=0), color="black", linewidth=1.2) +
  annotate("rect", xmin=0.4,xmax=0.8,ymin=-0.3,ymax=0.3, fill="#d3d3d3", alpha=0.25) +
  annotate("text", x=0.6, y=0.45, label="Inflection Window (Bins 0.4-0.8)",
           size=5, fontface="bold") +
  geom_segment(aes(x=c(0.4,0.8),xend=c(0.4,0.8),y=-0.25,yend=0.25),
               color="black", linetype="dashed", linewidth=0.8) +
  geom_point(aes(x=c(0.2,0.35,0.48,0.6,0.75,0.85), y=0), size=5, color="black") +
  geom_text(aes(x=c(0.2,0.35,0.48,0.6,0.75,0.85),
                y=c(0.18,0.18,0.18,0.28,0.18,0.18),
                label=c("Purinergic/Ca\u00b2\u207a","NF-\u03baB Priming",
                        "IEG & Buffering","PTGDS Peak","LCN2 Surge",
                        "Ferroptosis \u2193")),
            size=4.2, fontface="bold", vjust=0) +
  annotate("label", x=0.6, y=0.28, label="PTGDS Peak (Bin ~0.6)",
           size=4.5, fill="#e6f3ff", alpha=0.9, color="darkblue") +
  theme_void() + theme(plot.margin=margin(t=30,b=10,l=5,r=5))

# Combine and save
fig_s1 <- (pA | pB) / (pC | pD) / pE +
  plot_layout(heights=c(1,1,0.85)) +
  plot_annotation(tag_levels='A') &
  theme(plot.tag=element_text(face="bold", size=28),
        plot.margin=margin(5,5,5,5))

out_png  <- file.path(OUT_DIR, "FigS1.png")
out_tiff <- file.path(OUT_DIR, "FigS1.tiff")
ggsave(out_png,  fig_s1, width=14, height=11.5, dpi=600, bg="white")
ggsave(out_tiff, fig_s1, width=14, height=11.5, dpi=600,
       device="tiff", compression="lzw", bg="white")

message(">>> Saved: ", out_png)
message(">>> Script FigS1 complete.")
