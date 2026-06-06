# ==============================================================================
# FigS2_module_trajectories.R
#
# Purpose : Generate Figure S2 — Module-level trajectory analysis
#           Cross-correlation of astrocytic regulatory modules with PTGDS
#           and neuronal NDUFS1 across SEA-AD pseudo-progression bins.
#
# Panels:
#   A = Purinergic/Ca²⁺ module vs PTGDS (z-score, Lag −1, r = −0.886)
#   B = NF-κB Priming module vs PTGDS   (z-score, Lag 0,  r = −0.678)
#   C = IEG Stress module vs PTGDS      (z-score, Lag 0 concurrent)
#   D = Ferroptosis-Protective vs NDUFS1 (z-score, r = −0.857)
#   E = Timeline schematic               (Inflection Window Bins 0.4–0.8)
#
# Outputs (-> output/FigS2/):
#   FigS2_final.png   (300 dpi)
#   FigS2_final.tiff  (300 dpi, LZW — journal submission)
#
# Paper : Kim Y et al. (2025). https://github.com/YoungOukKim/MCI-to-AD
# Usage : Rscript analysis/FigS/FigS2_module_trajectories.R
# Requirements : R >= 4.3.2 | rhdf5, dplyr, tidyr, ggplot2, patchwork
# ==============================================================================

suppressPackageStartupMessages({
  library(rhdf5)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)

# ------------------------------------------------------------------------------
# Paths (relative to repository root; override with environment variables)
# ------------------------------------------------------------------------------
H5AD_FILE <- Sys.getenv("SEAAD_H5AD",
               unset = "data/SEA-AD/SEAAD_MTG_RNAseq_final-nuclei.h5ad")
OUT_DIR   <- Sys.getenv("FIGS2_OUT",
               unset = "output/FigS2")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Cell type classifier
# ------------------------------------------------------------------------------
classify_cell <- function(label) {
  label <- as.character(label)
  if (grepl("Sst", label, ignore.case = TRUE))            return("SST_Inhibitory")
  if (grepl("^L[0-9]|IT$|ET$|CT$|NP$|L6b", label))      return("Excitatory_Neuron")
  if (grepl("^Astro", label, ignore.case = TRUE))         return("Astrocyte")
  if (grepl("^Micro", label, ignore.case = TRUE))         return("Microglia")
  return("Other")
}

# ------------------------------------------------------------------------------
# 1. Data extraction
# ------------------------------------------------------------------------------
message(">>> [1/4] Loading data from H5AD...")
h5closeAll()

subclass_idx        <- h5read(H5AD_FILE, "/obs/Subclass")
subclass_categories <- h5read(H5AD_FILE, "/obs/__categories/Subclass")
cell_type_labels    <- subclass_categories[subclass_idx + 1]
bin_vec             <- as.numeric(h5read(H5AD_FILE,
                         "/obs/Continuous Pseudo-progression Score"))
cell_class          <- sapply(cell_type_labels, classify_cell)

gene_names_all <- h5read(H5AD_FILE, "var/_index")

module_genes <- c(
  "P2RY1", "P2RY12", "GJA1", "ITPR2",           # Purinergic/Ca²⁺
  "PTGS2", "IL6ST", "NFKBIA",                    # NF-κB Priming
  "FOS", "JUN", "EGR1", "ATF3",                  # IEG Stress
  "HMOX1", "SOD2", "MT1E", "MT2A", "CLU", "SLC1A2",  # Metabolic Buffering
  "GPX4", "FTH1", "SLC7A11",                     # Ferroptosis-Protective
  "PTGDS", "LCN2", "NDUFS1"                      # Core markers
)

gene_idx_map <- match(module_genes, gene_names_all) - 1
valid_mask   <- !is.na(gene_idx_map)
module_genes <- module_genes[valid_mask]
gene_idx_map <- gene_idx_map[valid_mask]

indptr   <- h5read(H5AD_FILE, "X/indptr", bit64conversion = "double")
n_cells  <- length(indptr) - 1
expr_mat <- matrix(0, nrow = n_cells, ncol = length(module_genes))
colnames(expr_mat) <- module_genes

block_size <- 200000
for (start_cell in seq(1, n_cells, by = block_size)) {
  end_cell  <- min(start_cell + block_size - 1, n_cells)
  start_pos <- indptr[start_cell]
  end_pos   <- indptr[end_cell + 1]
  count_pos <- end_pos - start_pos
  if (count_pos > 0) {
    idx_block  <- h5read(H5AD_FILE, "X/indices",
                         start = start_pos + 1, count = count_pos,
                         bit64conversion = "double")
    data_block <- h5read(H5AD_FILE, "X/data",
                         start = start_pos + 1, count = count_pos)
    offset <- 1
    for (i in start_cell:end_cell) {
      cell_nnz <- indptr[i + 1] - indptr[i]
      if (cell_nnz > 0) {
        current_indices <- idx_block[offset:(offset + cell_nnz - 1)]
        current_data    <- data_block[offset:(offset + cell_nnz - 1)]
        for (j in seq_along(gene_idx_map)) {
          m <- which(current_indices == gene_idx_map[j])
          if (length(m) > 0) expr_mat[i, j] <- current_data[m[1]]
        }
        offset <- offset + cell_nnz
      }
    }
  }
  message("  Progress: ", round(end_cell / n_cells * 100, 1), "%")
}
h5closeAll()

# ------------------------------------------------------------------------------
# 2. Bin-level summary
# ------------------------------------------------------------------------------
message(">>> [2/4] Summarizing by bin...")

df_all <- data.frame(score = bin_vec, cell_class = cell_class, expr_mat) %>%
  filter(!is.na(score) & score >= 0.1) %>%
  mutate(bin = round(score, 1))

astro_bin <- df_all %>%
  filter(cell_class == "Astrocyte") %>%
  group_by(bin) %>%
  summarise(across(all_of(module_genes), \(x) mean(x, na.rm = TRUE)),
            n_cells = n(), .groups = "drop") %>%
  arrange(bin)

neuron_bin <- df_all %>%
  filter(cell_class == "Excitatory_Neuron") %>%
  group_by(bin) %>%
  summarise(NDUFS1_neuron = mean(NDUFS1, na.rm = TRUE),
            n_neurons = n(), .groups = "drop") %>%
  arrange(bin)

# ------------------------------------------------------------------------------
# 3. 3-bin moving average + module composites
# ------------------------------------------------------------------------------
compute_ma3 <- function(x) {
  n <- length(x); ma <- rep(NA, n)
  for (i in 1:n) ma[i] <- mean(x[max(1, i-1):min(n, i+1)], na.rm = TRUE)
  ma
}

astro_ma <- astro_bin %>%
  mutate(across(all_of(module_genes), compute_ma3, .names = "{col}_ma"))

# PTGDS: LOESS fit (consistent with Fig. 1)
ptgds_loess <- loess(PTGDS ~ bin, data = astro_bin, span = 1.0)
astro_ma <- astro_ma %>%
  mutate(PTGDS_ma = predict(ptgds_loess, newdata = data.frame(bin = bin)))

neuron_ma <- neuron_bin %>%
  mutate(NDUFS1_neuron_ma = compute_ma3(NDUFS1_neuron))

# Module composite scores
astro_ma <- astro_ma %>%
  mutate(
    Purinergic_Ca_ma       = rowMeans(pick(P2RY1_ma, P2RY12_ma, GJA1_ma, ITPR2_ma),
                                      na.rm = TRUE),
    NFkB_Priming_ma        = rowMeans(pick(PTGS2_ma, IL6ST_ma, NFKBIA_ma),
                                      na.rm = TRUE),
    IEG_Stress_ma          = rowMeans(pick(FOS_ma, JUN_ma, EGR1_ma, ATF3_ma),
                                      na.rm = TRUE),
    Ferroptosis_Protect_ma = rowMeans(pick(GPX4_ma, FTH1_ma, SLC7A11_ma),
                                      na.rm = TRUE)
  ) %>%
  left_join(neuron_ma[, c("bin", "NDUFS1_neuron_ma")], by = "bin")

# ------------------------------------------------------------------------------
# 4. Panel construction
# ------------------------------------------------------------------------------
message(">>> [3/4] Building panels...")

zscore <- function(x) {
  m <- mean(x, na.rm = TRUE); s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - m) / s
}

plot_module_vs_ptgds <- function(mod_col, mod_name, color_mod, lag_text,
                                  lag_pos_x = 0.85, lag_pos_y = 0.92) {
  df <- astro_ma %>%
    dplyr::select(bin, mod = all_of(mod_col), PTGDS = PTGDS_ma) %>%
    mutate(mod = zscore(mod), PTGDS = zscore(PTGDS)) %>%
    pivot_longer(c("mod", "PTGDS"), names_to = "Module", values_to = "Expression") %>%
    mutate(Module = ifelse(Module == "PTGDS", "PTGDS (Astro)", mod_name))

  color_vals        <- c("PTGDS (Astro)" = "#0072B2")
  color_vals[mod_name] <- color_mod
  y_range           <- range(df$Expression, na.rm = TRUE)

  ggplot(df, aes(x = bin, y = Expression, color = Module)) +
    annotate("rect", xmin = 0.4, xmax = 0.8, ymin = -Inf, ymax = Inf,
             fill = "gray", alpha = 0.12) +
    geom_smooth(method = "loess", span = 0.9, se = FALSE, linewidth = 1.5) +
    geom_point(size = 1.8, alpha = 0.5) +
    annotate("label", x = lag_pos_x, y = y_range[2] * lag_pos_y,
             label = lag_text, size = 4.2, fill = "white", alpha = 0.85, hjust = 1) +
    scale_color_manual(values = color_vals) +
    scale_x_continuous(breaks = seq(0.1, 0.9, 0.1)) +
    labs(x = "Pseudo-progression bin", y = "3-bin MA expression (z-score)") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "top", legend.title = element_blank(),
          axis.text = element_text(size = 10),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 10))
}

# Panel A
pA <- plot_module_vs_ptgds(
  "Purinergic_Ca_ma", "Purinergic/Ca\u00B2\u207A Module",
  "#E74C3C", "Lag \u22121, r = \u22120.886"
)

# Panel B
pB <- plot_module_vs_ptgds(
  "NFkB_Priming_ma", "NF-\u03BAB Priming Module",
  "#ff7f0e", "Lag 0, r = \u22120.678"
)

# Panel C
pC <- plot_module_vs_ptgds(
  "IEG_Stress_ma", "IEG Stress Module",
  "#2ca02c", "Lag 0 concurrent", lag_pos_y = 0.85
)

# Panel D — Ferroptosis-Protective vs Neuronal NDUFS1
df_d <- astro_ma %>%
  dplyr::select(bin,
                Ferroptosis = Ferroptosis_Protect_ma,
                NDUFS1      = NDUFS1_neuron_ma) %>%
  mutate(Ferroptosis = zscore(Ferroptosis), NDUFS1 = zscore(NDUFS1)) %>%
  pivot_longer(c("Ferroptosis", "NDUFS1"),
               names_to = "Module", values_to = "Expression") %>%
  mutate(Module = ifelse(Module == "Ferroptosis",
                         "Ferroptosis-Protective Module", "NDUFS1 (Neuron)"))

pD <- ggplot(df_d, aes(x = bin, y = Expression, color = Module)) +
  annotate("rect", xmin = 0.4, xmax = 0.8, ymin = -Inf, ymax = Inf,
           fill = "gray", alpha = 0.12) +
  geom_smooth(method = "loess", span = 0.9, se = FALSE, linewidth = 1.5) +
  geom_point(size = 1.8, alpha = 0.5) +
  annotate("label", x = 0.85, y = max(df_d$Expression, na.rm = TRUE) * 0.85,
           label = "Inverse sync, r = \u22120.857",
           size = 4.2, fill = "white", alpha = 0.85, hjust = 1) +
  scale_color_manual(values = c(
    "Ferroptosis-Protective Module" = "#9467bd",
    "NDUFS1 (Neuron)"              = "#2C3E50"
  )) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.1)) +
  labs(x = "Pseudo-progression bin", y = "3-bin MA expression (z-score)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top", legend.title = element_blank(),
        axis.text = element_text(size = 10),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 10))

# Panel E — Timeline schematic
pE <- ggplot() +
  geom_segment(aes(x = 0.1, xend = 0.9, y = 0, yend = 0),
               color = "black", linewidth = 1.2) +
  annotate("rect", xmin = 0.4, xmax = 0.8, ymin = -0.3, ymax = 0.3,
           fill = "#d3d3d3", alpha = 0.25) +
  annotate("text", x = 0.6, y = 0.45,
           label = "Inflection Window (Bins 0.4\u20130.8)",
           size = 5, fontface = "bold") +
  geom_segment(aes(x = 0.4, xend = 0.4, y = -0.25, yend = 0.25),
               color = "black", linetype = "dashed", linewidth = 0.8) +
  geom_segment(aes(x = 0.8, xend = 0.8, y = -0.25, yend = 0.25),
               color = "black", linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = 0.4, y = -0.25, label = "0.4", size = 4.5, vjust = 1.5) +
  annotate("text", x = 0.8, y = -0.25, label = "0.8", size = 4.5, vjust = 1.5) +
  geom_point(aes(x = c(0.2, 0.35, 0.48, 0.6, 0.75, 0.85), y = 0),
             size = 5, color = "black") +
  geom_text(aes(x = c(0.2, 0.35, 0.48, 0.6, 0.75, 0.85),
                y = c(0.18, 0.18, 0.18, 0.28, 0.18, 0.18),
                label = c("Purinergic/Ca\u00B2\u207A", "NF-\u03BAB Priming",
                          "IEG & Buffering", "PTGDS Peak",
                          "LCN2 Surge", "Ferroptosis \u2193")),
            size = 4.2, fontface = "bold", vjust = 0) +
  annotate("segment", x = 0.6, xend = 0.6, y = 0.05, yend = 0.20,
           color = "darkblue", linewidth = 1,
           arrow = arrow(length = unit(0.25, "cm"))) +
  annotate("label", x = 0.6, y = 0.28, label = "PTGDS Peak (Bin ~0.6)",
           size = 4.5, fill = "#e6f3ff", alpha = 0.9, color = "darkblue") +
  theme_void() +
  theme(plot.margin = margin(t = 30, b = 10, l = 5, r = 5))

# ------------------------------------------------------------------------------
# 5. Assemble and save
# ------------------------------------------------------------------------------
message(">>> [4/4] Assembling and saving...")

fig_s2 <- (pA | pB) / (pC | pD) / pE +
  plot_layout(heights = c(1, 1, 0.85)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag    = element_text(face = "bold", size = 28),
        plot.margin = margin(t = 5, b = 5, l = 5, r = 5))

W   <- 14
H   <- 11.5
DPI <- 300

out_png  <- file.path(OUT_DIR, "FigS2_final.png")
out_tiff <- file.path(OUT_DIR, "FigS2_final.tiff")

ggsave(out_png, fig_s2,
       width = W, height = H, dpi = DPI, bg = "white")
ggsave(out_tiff, fig_s2,
       width = W, height = H, dpi = DPI, bg = "white",
       device = "tiff", compression = "lzw")

message("\n", strrep("=", 55))
message("  Done! PNG : ", out_png)
message("  Done! TIFF: ", out_tiff)
message("  Size: ", W, " x ", H, " in @ ", DPI, " dpi")
message("  Panels: A=Purinergic/Ca | B=NF-kB | C=IEG | D=Ferroptosis | E=Timeline")
message(strrep("=", 55))
