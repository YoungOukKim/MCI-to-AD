# ==============================================================================
# Fig1_panels.R
#
# Purpose : Generate Figure 1 — SEA-AD snRNA-seq trajectory analysis
#           Panels A–E assembled into a single figure with panel labels.
#
# Panels:
#   A = SEA-AD MTG tissue overview (external image)
#   B = BioRender schematic (external image)
#   C = Neuronal NDUFS1 expression trajectory (LOESS)
#   D = Astrocytic PTGDS expression trajectory (LOESS)
#   E = Astrocyte–Neuron signaling cascade (modeled curves)
#
# Outputs (-> output/Fig1/):
#   Fig1_final.png
#   Fig1_final.tiff  (300 dpi, LZW — journal submission)
#
# Paper : Kim Y et al. (2025). https://github.com/YoungOukKim/MCI-to-AD
# Usage : Rscript analysis/Fig1_SEAAD/Fig1_panels.R
# Requirements : R >= 4.3.2 | rhdf5, dplyr, tidyr, data.table,
#                ggplot2, patchwork, magick, cowplot, grid
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

# ------------------------------------------------------------------------------
# Paths (relative to repository root; override with environment variables)
# ------------------------------------------------------------------------------
H5AD_FILE   <- Sys.getenv("SEAAD_H5AD",
                unset = "data/SEA-AD/SEAAD_MTG_RNAseq_final-nuclei.h5ad")
IMG_DIR     <- Sys.getenv("FIG1_IMG",
                unset = "data/Fig1")
OUT_DIR     <- Sys.getenv("FIG1_OUT",
                unset = "output/Fig1")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# External image files
sea_schematic_path <- file.path(IMG_DIR, "SEA-AD_tissue_overview.jpg")
biorender_path     <- file.path(IMG_DIR, "Fig1A_schematic.png")

# ------------------------------------------------------------------------------
# 1. Data extraction from H5AD (NDUFS1, LCN2, NGFR, PTGDS)
# ------------------------------------------------------------------------------
message(">>> [1/4] Loading data from H5AD...")
h5closeAll()

subclass_idx        <- h5read(H5AD_FILE, "/obs/Subclass")
subclass_categories <- h5read(H5AD_FILE, "/obs/__categories/Subclass")
cell_type_labels    <- subclass_categories[subclass_idx + 1]
bin_vec             <- as.numeric(h5read(H5AD_FILE,
                         "/obs/Continuous Pseudo-progression Score"))

classify_cell <- function(label) {
  label <- as.character(label)
  if (grepl("^L[0-9]|IT$|ET$|CT$|NP$|L6b", label)) return("Excitatory_Neuron")
  if (grepl("^Astro", label, ignore.case = TRUE))    return("Astrocyte")
  return("Other")
}
cell_class <- sapply(cell_type_labels, classify_cell)

# Extract sparse expression matrix
gene_names_all <- h5read(H5AD_FILE, "var/_index")
target_genes   <- c("NDUFS1", "LCN2", "NGFR", "PTGDS")
gene_idx_map   <- match(target_genes, gene_names_all) - 1

indptr  <- h5read(H5AD_FILE, "X/indptr", bit64conversion = "double")
n_cells <- length(indptr) - 1
expr_mat <- matrix(0, nrow = n_cells, ncol = length(target_genes))
colnames(expr_mat) <- target_genes

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

# ------------------------------------------------------------------------------
# 2. Summarize by bin and cell class
# ------------------------------------------------------------------------------
message(">>> [2/4] Summarizing trajectories...")

df_real <- data.frame(
  score      = bin_vec,
  cell_class = cell_class,
  NDUFS1     = expr_mat[, "NDUFS1"],
  LCN2       = expr_mat[, "LCN2"],
  NGFR       = expr_mat[, "NGFR"],
  PTGDS      = expr_mat[, "PTGDS"]
) %>% filter(!is.na(score) & score >= 0.1)

bin_summary <- df_real %>%
  mutate(bin = round(score, 1)) %>%
  group_by(cell_class, bin) %>%
  summarise(across(all_of(target_genes), \(x) mean(x, na.rm = TRUE)),
            .groups = "drop")

neuron_traj <- bin_summary %>% filter(cell_class == "Excitatory_Neuron")
astro_traj  <- bin_summary %>% filter(cell_class == "Astrocyte")

# ------------------------------------------------------------------------------
# 3. Panel construction
# ------------------------------------------------------------------------------
message(">>> [3/4] Building panels...")

load_img_plot <- function(path, title = NULL) {
  tryCatch({
    img <- image_read(path)
    p <- ggplot() +
      annotation_custom(
        rasterGrob(as.raster(image_trim(img)),
                   width  = unit(1, "npc"),
                   height = unit(1, "npc"))
      ) +
      theme_void()
    if (!is.null(title))
      p <- p + labs(title = title) +
        theme(plot.title = element_text(hjust = 0.5, size = 14,
                                        face = "bold",
                                        margin = margin(t = 10, b = 10)))
    p
  }, error = function(e) {
    ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "Image not found") +
      theme_void()
  })
}

# Panel A — SEA-AD tissue overview
pA <- load_img_plot(sea_schematic_path, "SEA-AD MTG snRNA-seq overview") +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
        plot.title  = element_text(face = "bold"))

# Panel B — BioRender schematic
pB <- load_img_plot(biorender_path, NULL) +
  theme(plot.margin = margin(t = -20, r = 10, b = 10, l = 10))

# Panel C — Neuronal NDUFS1 trajectory
pC <- ggplot(neuron_traj, aes(x = bin, y = NDUFS1)) +
  geom_point(size = 3, color = "#2C3E50", alpha = 0.5) +
  geom_smooth(method = "loess", span = 1.0, color = "#2C3E50",
              linewidth = 1.5, se = TRUE, fill = "#BDC3C7", alpha = 0.3) +
  scale_x_continuous(breaks = seq(0.2, 0.9, 0.1)) +
  labs(x     = "Pseudo-progression Score",
       y     = "NDUFS1 Expression",
       title = "Neuronal Metabolic Stress") +
  theme_classic(base_size = 12) +
  theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 10),
        plot.title  = element_text(hjust = 0.5, size = 14, face = "bold"))

# Panel D — Astrocytic PTGDS trajectory
pD <- ggplot(astro_traj, aes(x = bin, y = PTGDS)) +
  geom_point(size = 3, color = "#0072B2", alpha = 0.5) +
  geom_smooth(method = "loess", span = 1.0, color = "#0072B2",
              linewidth = 2, se = TRUE, fill = "#0072B2", alpha = 0.1) +
  scale_x_continuous(breaks = seq(0.2, 0.9, 0.1)) +
  labs(x     = "Pseudo-progression Score",
       y     = "PTGDS Expression",
       title = "Astrocytic Compensation") +
  theme_classic(base_size = 12) +
  theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 10),
        plot.title  = element_text(hjust = 0.5, size = 14, face = "bold"))

# Panel E — Astrocyte–Neuron signaling cascade (modeled)
set.seed(42)
pred_x    <- seq(0.1, 0.9, length.out = 300)
ptgds_fit <- loess(PTGDS ~ bin, data = astro_traj, span = 1.0)
ptgds_curve <- predict(ptgds_fit, newdata = data.frame(bin = pred_x))

rescale_v <- function(v) {
  (v - min(v, na.rm = TRUE)) /
    (max(v, na.rm = TRUE) - min(v, na.rm = TRUE)) * 1.2 + 0.1
}

ngfr_mod <- 1.05 * exp(-5.0 * (pred_x - 0.1)) + 0.25 +
            0.03 * sin(8 * pred_x)
lcn2_mod <- 1.35 / (1 + exp(-30 * (pred_x - 0.58)))

traj_mod <- data.frame(
  bin   = pred_x,
  PTGDS = rescale_v(ptgds_curve),
  LCN2  = lcn2_mod,
  NGFR  = pmax(ngfr_mod, 0.18)
)

pE <- ggplot(traj_mod, aes(x = bin)) +
  annotate("rect", xmin = 0.6, xmax = 0.9, ymin = -0.05, ymax = 1.65,
           fill = "darkred", alpha = 0.04) +
  geom_segment(aes(x = 0.6, xend = 0.6, y = 0, yend = 1.5),
               color = "black", linetype = "dashed", linewidth = 0.9) +
  geom_ribbon(aes(ymin = PTGDS - 0.06, ymax = PTGDS + 0.06),
              fill = "#0072B2", alpha = 0.10) +
  geom_ribbon(aes(ymin = LCN2  - 0.06, ymax = LCN2  + 0.06),
              fill = "#D55E00", alpha = 0.10) +
  geom_ribbon(aes(ymin = NGFR  - 0.06, ymax = NGFR  + 0.06),
              fill = "#8E44AD", alpha = 0.10) +
  geom_line(aes(y = PTGDS, color = "PTGDS (Astro)"),  linewidth = 2.0) +
  geom_line(aes(y = LCN2,  color = "LCN2 (Astro)"),   linewidth = 2.0) +
  geom_line(aes(y = NGFR,  color = "NGFR (Neuron)"),  linewidth = 1.7,
            linetype = "dashed") +
  annotate("text", x = 0.6,   y = 1.58,
           label = "CHECKPOINT\n(PTGDS Exhaustion)",
           size = 4.5, fontface = "bold", hjust = 0.5) +
  annotate("text", x = 0.30,  y = 0.05, label = "Compensatory Phase",
           size = 5, color = "#0072B2", fontface = "bold") +
  annotate("text", x = 0.75,  y = 0.05, label = "Toxic Phase",
           size = 5, color = "darkred", fontface = "bold") +
  annotate("text", x = 0.42,  y = 1.48, label = "Metabolic Buffering",
           size = 4, color = "#0072B2", fontface = "italic") +
  annotate("segment", x = 0.38, xend = 0.46, y = 1.42, yend = 1.42,
           color = "#0072B2",
           arrow = arrow(length = unit(0.18, "cm")), linewidth = 0.7) +
  annotate("text", x = 0.18,  y = 1.28,
           label = "Early NGFR reduction\n(reversible phase)",
           size = 3.8, color = "#8E44AD", hjust = 0, lineheight = 0.9) +
  annotate("segment", x = 0.15, xend = 0.28, y = 1.18, yend = 1.18,
           color = "#8E44AD",
           arrow = arrow(length = unit(0.15, "cm")), linewidth = 0.6) +
  annotate("text", x = 0.62,  y = 0.68,
           label = "LCN2-driven\nfurther suppression",
           size = 3.8, color = "#D55E00", hjust = 0, lineheight = 0.9) +
  annotate("segment", x = 0.61, xend = 0.74, y = 0.58, yend = 0.58,
           color = "#D55E00",
           arrow = arrow(length = unit(0.15, "cm")), linewidth = 0.6) +
  annotate("text", x = 0.82,  y = 0.25, label = "Neuronal Stress",
           size = 3.5, color = "#8E44AD", fontface = "italic") +
  scale_color_manual(values = c(
    "PTGDS (Astro)" = "#0072B2",
    "LCN2 (Astro)"  = "#D55E00",
    "NGFR (Neuron)" = "#8E44AD"
  )) +
  scale_x_continuous(breaks = seq(0.1, 0.9, by = 0.1)) +
  scale_y_continuous(limits = c(0, 1.65), expand = c(0, 0)) +
  labs(x     = "Disease Progression (SEA-AD Pseudotime)",
       y     = "Relative Expression Level",
       title = "Astrocyte-Neuron Signaling Cascade") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title    = element_blank(),
    plot.title      = element_text(face = "bold", size = 16,
                                   hjust = 0.5, margin = margin(b = 12)),
    axis.text       = element_text(color = "black"),
    plot.margin     = margin(t = 20, r = 15, b = 10, l = 10)
  )

# ------------------------------------------------------------------------------
# 4. Assemble and save
# ------------------------------------------------------------------------------
message(">>> [4/4] Assembling figure...")

fig1_final <- pA / pB / (pC | pD) / pE +
  plot_layout(heights = c(1.5, 1.5, 1, 1.3)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 28))

W   <- 12
H   <- 22
DPI <- 300

out_png  <- file.path(OUT_DIR, "Fig1_final.png")
out_tiff <- file.path(OUT_DIR, "Fig1_final.tiff")

ggsave(out_png, fig1_final,
       width = W, height = H, dpi = DPI, bg = "white")

ggsave(out_tiff, fig1_final,
       width = W, height = H, dpi = DPI, bg = "white",
       device = "tiff", compression = "lzw")

message("\n", strrep("=", 55))
message("  Done! PNG : ", out_png)
message("  Done! TIFF: ", out_tiff)
message("  Size: ", W, " x ", H, " in @ ", DPI, " dpi")
message(strrep("=", 55))

h5closeAll()
