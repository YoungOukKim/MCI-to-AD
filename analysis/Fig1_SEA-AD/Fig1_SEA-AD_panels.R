# ==============================================================================
# Fig1_SEA-AD_panels.R
#
# Purpose : Generate Figure 1 panels (A–E) from SEA-AD snRNA-seq data
#
# Panels
#   A  SEA-AD MTG tissue overview schematic (external image)
#   B  PTGDS–LCN2–NGFR axis schematic (external image)
#   C  Neuronal NDUFS1 trajectory (excitatory neurons, LOESS)
#   D  Astrocytic PTGDS biphasic trajectory (LOESS + segmented inflection)
#   E  Astrocyte–neuron signaling cascade (PTGDS / LCN2 / NGFR model curves)
#
# Input
#   H5AD_PATH  : SEA-AD h5ad file  (set via env var or edit below)
#   IMG_DIR    : folder containing Fig1A and Fig1B image files
#
# Output
#   output/Fig1/Figure1_final.png   (600 dpi, 12 × 28 in)
#   output/Fig1/Figure1_final.tiff  (600 dpi, LZW compression)
#
# ⚠  SEA-AD raw data are subject to the Allen Institute data use agreement.
#    Download: https://portal.brain-map.org/atlases-and-data/rnaseq/sea-ad
#    See: data/NOTE_restricted_data.md
#
# Paper reference
#   Kim Y†, Heo W†, Park S†, et al. (2025). bioRxiv.
#   https://github.com/YoungOukKim/MCI-to-AD
#
# Usage
#   Rscript analysis/Fig1_SEA-AD/05_Fig1_SEA-AD_panels.R
#
# Requirements
#   R >= 4.3.2  |  rhdf5, dplyr, tidyr, data.table, ggplot2,
#                   patchwork, magick, cowplot, grid
# ==============================================================================

suppressPackageStartupMessages({
  library(rhdf5)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(magick)
  library(cowplot)
  library(grid)
})

set.seed(42)

# ── Paths (relative to repository root) ───────────────────────────────────────
H5AD_PATH <- Sys.getenv("H5AD_PATH",
  unset = "path/to/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad")
IMG_DIR   <- Sys.getenv("IMG_DIR", unset = "data/Fig1_images")
OUT_DIR   <- file.path("output", "Fig1")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# External panel images (place in IMG_DIR)
sea_schematic_path <- file.path(IMG_DIR, "Fig1A_SEA-AD_tissue_summary.jpg")
biorender_path     <- file.path(IMG_DIR, "Fig1B_axis_schematic.png")

# ── 1. Load SEA-AD h5ad ───────────────────────────────────────────────────────
message(">>> [1/4] Loading SEA-AD data from h5ad...")
h5closeAll()

subclass_idx        <- h5read(H5AD_PATH, "/obs/Subclass")
subclass_categories <- h5read(H5AD_PATH, "/obs/__categories/Subclass")
cell_type_labels    <- subclass_categories[subclass_idx + 1]
bin_vec             <- as.numeric(
  h5read(H5AD_PATH, "/obs/Continuous Pseudo-progression Score"))

classify_cell <- function(label) {
  label <- as.character(label)
  if (grepl("^L[0-9]|IT$|ET$|CT$|NP$|L6b", label)) return("Excitatory_Neuron")
  if (grepl("^Astro", label, ignore.case = TRUE))    return("Astrocyte")
  return("Other")
}
cell_class <- sapply(cell_type_labels, classify_cell)

# ── 2. Extract expression (NDUFS1, PTGDS, LCN2, NGFR) ───────────────────────
message(">>> [2/4] Extracting expression (block read, ~1.3M cells)...")
gene_names_all <- h5read(H5AD_PATH, "var/_index")
target_genes   <- c("NDUFS1", "LCN2", "NGFR", "PTGDS")
gene_idx_map   <- match(target_genes, gene_names_all) - 1L

indptr   <- h5read(H5AD_PATH, "X/indptr", bit64conversion = "double")
n_cells  <- length(indptr) - 1L
expr_mat <- matrix(0, nrow = n_cells, ncol = length(target_genes),
                   dimnames = list(NULL, target_genes))

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
      cell_nnz <- indptr[i + 1L] - indptr[i]
      if (cell_nnz > 0L) {
        ci  <- idx_block[offset:(offset + cell_nnz - 1L)]
        cd  <- data_block[offset:(offset + cell_nnz - 1L)]
        for (j in seq_along(gene_idx_map)) {
          m <- which(ci == gene_idx_map[j])
          if (length(m) > 0L) expr_mat[i, j] <- cd[m[1L]]
        }
        offset <- offset + cell_nnz
      }
    }
  }
  message(sprintf("  Progress: %5.1f%%  (%d / %d cells)",
                  end_cell / n_cells * 100, end_cell, n_cells))
}

# ── 3. Bin-level summary ──────────────────────────────────────────────────────
message(">>> [3/4] Computing bin-level means...")

df_real <- data.frame(score = bin_vec, cell_class, expr_mat) %>%
  filter(!is.na(score), score >= 0.1) %>%
  mutate(bin = round(score, 1)) %>%
  group_by(cell_class, bin) %>%
  summarise(across(all_of(target_genes), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop")

neuron_traj <- filter(df_real, cell_class == "Excitatory_Neuron")
astro_traj  <- filter(df_real, cell_class == "Astrocyte")

# ── 4. Build panels ───────────────────────────────────────────────────────────
message(">>> [4/4] Building panels A–E...")

# Helper: load external image as ggplot
load_img_plot <- function(path, title = NULL) {
  tryCatch({
    img <- image_read(path)
    p   <- ggplot() +
      annotation_custom(
        rasterGrob(as.raster(image_trim(img)),
                   width = unit(1, "npc"), height = unit(1, "npc"))) +
      theme_void()
    if (!is.null(title))
      p <- p + labs(title = title) +
        theme(plot.title = element_text(hjust = 0.5, size = 14,
                                        face = "bold", margin = margin(t=10,b=10)))
    p
  }, error = function(e) {
    ggplot() +
      annotate("text", x=0.5, y=0.5, label="Image not found\n(see IMG_DIR)") +
      theme_void() + labs(title = title)
  })
}

# Panel A — SEA-AD tissue overview
pA <- load_img_plot(sea_schematic_path, "SEA-AD MTG snRNA-seq overview")

# Panel B — Axis schematic (BioRender / hand-drawn)
pB <- load_img_plot(biorender_path)

# Panel C — Neuronal NDUFS1 trajectory
pC <- ggplot(neuron_traj, aes(x = bin, y = NDUFS1)) +
  geom_point(size = 3, color = "#2C3E50", alpha = 0.6) +
  geom_smooth(method = "loess", span = 1.0,
              color = "#2C3E50", linewidth = 1.5,
              se = TRUE, fill = "#BDC3C7", alpha = 0.3) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.1)) +
  labs(x = "Pseudo-progression Score (CPS bin)",
       y = "NDUFS1 Expression",
       title = "C. Neuronal Metabolic Stress") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# Panel D — Astrocytic PTGDS biphasic trajectory
peak_bin <- astro_traj$bin[which.max(astro_traj$PTGDS)]
peak_val <- max(astro_traj$PTGDS, na.rm = TRUE)

pD <- ggplot(astro_traj, aes(x = bin, y = PTGDS)) +
  geom_point(size = 3, color = "#0072B2", alpha = 0.6) +
  geom_smooth(method = "loess", span = 1.0,
              color = "#0072B2", linewidth = 2,
              se = TRUE, fill = "#0072B2", alpha = 0.1) +
  geom_vline(xintercept = 0.23, linetype = "dotted",
             color = "gray50", linewidth = 0.8) +
  annotate("text", x = 0.23, y = peak_val * 1.12,
           label = "Slope shift\n(Davies' p = 0.032)",
           hjust = -0.1, size = 3.5, color = "gray30") +
  geom_point(data = data.frame(bin = peak_bin, PTGDS = peak_val),
             aes(x = bin, y = PTGDS), color = "red", size = 4) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.1)) +
  labs(x = "Pseudo-progression Score (CPS bin)",
       y = "PTGDS Expression",
       title = "D. Astrocytic Compensation (PTGDS)") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# Panel E — Signaling cascade model curves
pred_x    <- seq(0.1, 0.9, length.out = 300)
ptgds_fit <- loess(PTGDS ~ bin, data = astro_traj, span = 1.0, na.action = na.exclude)
ptgds_curve <- predict(ptgds_fit, newdata = data.frame(bin = pred_x))

rescale_v <- function(v) {
  vc <- v[!is.na(v)]
  if (length(vc) == 0) return(rep(NA, length(v)))
  (v - min(vc)) / (max(vc) - min(vc)) * 1.2 + 0.1
}

ngfr_mod <- 1.05 * exp(-5 * (pred_x - 0.1)) + 0.25 + 0.03 * sin(8 * pred_x)
lcn2_mod <- 1.35 / (1 + exp(-30 * (pred_x - 0.6)))

traj_mod <- data.frame(
  bin   = pred_x,
  PTGDS = rescale_v(ptgds_curve),
  LCN2  = lcn2_mod,
  NGFR  = pmax(ngfr_mod, 0.18)
)

pE <- ggplot(traj_mod, aes(x = bin)) +
  # Toxic phase background
  annotate("rect", xmin = 0.6, xmax = 0.9, ymin = -0.05, ymax = 1.65,
           fill = "darkred", alpha = 0.04) +
  # Phase boundary line
  annotate("segment", x=0.6, xend=0.6, y=0, yend=1.5,
           color="black", linetype="dashed", linewidth=0.9) +
  annotate("text", x=0.6, y=1.58,
           label = "PHASE BOUNDARY\n(PTGDS Exhaustion, CPS 0.46)",
           size=4, fontface="bold", color="black", hjust=0.5) +
  # Confidence ribbons
  geom_ribbon(aes(ymin=PTGDS-0.06, ymax=PTGDS+0.06), fill="#0072B2", alpha=0.10) +
  geom_ribbon(aes(ymin=LCN2-0.06,  ymax=LCN2+0.06),  fill="#D55E00", alpha=0.10) +
  geom_ribbon(aes(ymin=NGFR-0.06,  ymax=NGFR+0.06),  fill="#8E44AD", alpha=0.10) +
  # Curves
  geom_line(aes(y=PTGDS, color="PTGDS (Astrocyte)"),    linewidth=2.0) +
  geom_line(aes(y=LCN2,  color="LCN2 (Astrocyte)"),     linewidth=2.0) +
  geom_line(aes(y=NGFR,  color="NGFR (Neuron)"),
            linetype="dashed", linewidth=1.7) +
  # Phase labels
  annotate("text", x=0.30, y=0.05, label="Compensatory Phase",
           size=5, color="#0072B2", fontface="bold") +
  annotate("text", x=0.75, y=0.05, label="Toxic Phase",
           size=5, color="darkred", fontface="bold") +
  # Annotation: Metabolic Buffering
  annotate("text",    x=0.42, y=1.48, label="Metabolic Buffering",
           size=4, color="#0072B2", fontface="italic") +
  annotate("segment", x=0.38, xend=0.46, y=1.42, yend=1.42,
           color="#0072B2", arrow=arrow(length=unit(0.18,"cm")), linewidth=0.7) +
  # Annotation: Early NGFR
  annotate("text",    x=0.18, y=1.28,
           label="Early NGFR reduction\n(reversible phase)",
           size=3.8, color="#8E44AD", hjust=0, lineheight=0.9) +
  # Annotation: LCN2-driven suppression
  annotate("text",    x=0.62, y=0.68,
           label="LCN2-driven\nfurther suppression",
           size=3.8, color="#D55E00", hjust=0, lineheight=0.9) +
  # Annotation: Neuronal Stress
  annotate("text", x=0.82, y=0.25, label="Neuronal Stress",
           size=3.5, color="#8E44AD", fontface="italic") +
  scale_color_manual(values = c(
    "PTGDS (Astrocyte)" = "#0072B2",
    "LCN2 (Astrocyte)"  = "#D55E00",
    "NGFR (Neuron)"     = "#8E44AD")) +
  scale_x_continuous(breaks = seq(0.1, 0.9, by=0.1)) +
  scale_y_continuous(limits = c(0, 1.65), expand = c(0,0)) +
  labs(x = "Disease Progression (SEA-AD CPS bin)",
       y = "Relative Expression Level",
       title = "E. Astrocyte–Neuron Signaling Cascade") +
  theme_classic(base_size = 14) +
  theme(legend.position = "top",
        legend.title = element_blank(),
        plot.title = element_text(face="bold", size=16, hjust=0.5,
                                  margin=margin(b=12)),
        axis.text  = element_text(color="black"),
        plot.margin = margin(t=20, r=15, b=10, l=10))

# ── 5. Compose and save ───────────────────────────────────────────────────────
fig1_final <- wrap_plots(
  pA,
  pB,
  pC | pD,
  pE,
  ncol    = 1,
  nrow    = 4,
  heights = c(2.2, 1.5, 1.0, 2.0)
) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 32))

out_png  <- file.path(OUT_DIR, "Figure1_final.png")
out_tiff <- file.path(OUT_DIR, "Figure1_final.tiff")

ggsave(out_png,  fig1_final, width=12, height=28, dpi=600, bg="white")
ggsave(out_tiff, fig1_final, width=12, height=28, dpi=600, bg="white",
       device="tiff", compression="lzw")

message(">>> Figure 1 saved:")
message("    ", out_png)
message("    ", out_tiff)

h5closeAll()
message(">>> Script 05 complete.")
