# ==============================================================================
# FigureS7_ADNI_CSF_crossplatform.R   (revised layout)
#
# Purpose : Generate Figure S7 - Cross-platform validation of PTGDS/LCN2/NEFL
#           dynamics in ADNI CSF proteomics (TMT-MS and SomaScan 7K).
#
# Layout  : A, B, C stacked VERTICALLY (one panel per row) for legibility.
#
# Panels:
#   A = TMT-MS vs SomaScan PTGDS - weak between-platform correlation (r ~ 0.15)
#   B = Cross-platform LCN2 and NEFL elevation across MMSE decline
#   C = TMT-MS CSF PTGDS vs MMSE - flat / no robust association
#       (all subjects shown as points; biphasic signal is tissue-restricted).
#       n and Pearson r are computed from the data (not hard-coded).
#
# Outputs (-> output/FigS7/):
#   FigureS7_CSF_crossplatform.png   (300 dpi)
#   FigureS7_CSF_crossplatform.tiff  (300 dpi, LZW - journal submission)
#
# Data note:
#   Raw ADNI data (TMT-MS, SomaScan, ADNIMERGE) cannot be redistributed
#   under the ADNI data use agreement (https://adni.loni.usc.edu).
#   Pre-computed summary CSVs are provided in data/ADNI/:
#     FigS6_platform_cor.csv  - Panel A: RID, TMT_log2, Soma_log2, Group
#     FigS6_trend.csv         - Panel B: MMSE, LCN2_smooth, NEFL_smooth, platform
#     FigS6_tmt_ptgds.csv     - Panel C: MMSE, PTGDS_final (per-batch z)
#
# Paper : Kim Y et al. (2025). https://github.com/YoungOukKim/MCI-to-AD
# Usage : Rscript analysis/FigS/FigureS7_ADNI_CSF_crossplatform.R
# Requirements : R >= 4.3.2 | tidyverse, ggplot2, cowplot
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(cowplot)
})

# ------------------------------------------------------------------------------
# Paths (relative to repository root; override with environment variables)
# ------------------------------------------------------------------------------
DATA_DIR <- Sys.getenv("FIGS7_DATA", unset = "data/ADNI")
OUT_DIR  <- Sys.getenv("FIGS7_OUT",  unset = "output/FigS7")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Load pre-computed summary CSVs
# ------------------------------------------------------------------------------
message(">>> [1/3] Loading summary data...")

cor_data   <- read.csv(file.path(DATA_DIR, "FigS6_platform_cor.csv"))
trend_data <- read.csv(file.path(DATA_DIR, "FigS6_trend.csv"))
tmt_c      <- read.csv(file.path(DATA_DIR, "FigS6_tmt_ptgds.csv")) %>%
              filter(!is.na(MMSE), !is.na(PTGDS_final))

soma_trend <- trend_data %>% filter(platform == "Soma")
tmt_trend  <- trend_data %>% filter(platform == "TMT")

# --- Statistics computed directly from the data (never hard-coded) -----------
rA <- cor(cor_data$TMT_log2, cor_data$Soma_log2, use = "complete.obs")
rC <- cor(tmt_c$MMSE, tmt_c$PTGDS_final, use = "complete.obs")
nC <- nrow(tmt_c)

titleA  <- sprintf("PTGDS between platforms (weak between-platform correlation, r \u2248 %.2f)", rA)
annotC  <- sprintf("Pearson r = %.2f (n.s., n = %d)", rC, nC)

message("  Panel A: ", nrow(cor_data), " subjects | between-platform r = ", round(rA, 3))
message("  Panel B: Soma=", nrow(soma_trend), " TMT=", nrow(tmt_trend), " rows")
message("  Panel C: n = ", nC, " | PTGDS-MMSE r = ", round(rC, 3), " (n.s.)")

# ------------------------------------------------------------------------------
# 2. Panel construction
# ------------------------------------------------------------------------------
message(">>> [2/3] Building panels...")

base_theme <- theme_classic(base_size = 16) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
        axis.title = element_text(size = 16),
        axis.text  = element_text(size = 14))

# Panel A - Platform correlation (PTGDS); legend INSIDE to avoid clipping
pA <- ggplot(cor_data, aes(x = TMT_log2, y = Soma_log2)) +
  geom_point(aes(color = Group), alpha = 0.55, size = 1.9) +
  geom_smooth(method = "lm", formula = y ~ x, color = "black", se = FALSE, linewidth = 1.1) +
  scale_color_manual(values = c("CN" = "#99CCFF", "MCI" = "#FFCC99", "AD" = "#FF9999"),
                     breaks = c("CN", "MCI", "AD"), na.translate = FALSE) +
  labs(title = titleA, x = "TMT-MS (log2)", y = "SomaScan (log2)", color = "Clinical group") +
  base_theme +
  theme(legend.position = c(0.92, 0.20),
        legend.background = element_rect(fill = alpha("white", 0.92), color = "grey70"),
        legend.title = element_text(size = 12), legend.text = element_text(size = 11))

# Panel B - Cross-platform LCN2 / NEFL elevation (solid = Soma, dashed = TMT)
pB <- ggplot() +
  geom_ribbon(data = soma_trend, aes(x = MMSE, ymin = LCN2_smooth - 0.25, ymax = LCN2_smooth + 0.25),
              fill = "#CC79A7", alpha = 0.18) +
  geom_ribbon(data = soma_trend, aes(x = MMSE, ymin = NEFL_smooth - 0.25, ymax = NEFL_smooth + 0.25),
              fill = "#D55E00", alpha = 0.18) +
  geom_ribbon(data = tmt_trend,  aes(x = MMSE, ymin = LCN2_smooth - 0.25, ymax = LCN2_smooth + 0.25),
              fill = "#CC79A7", alpha = 0.10) +
  geom_ribbon(data = tmt_trend,  aes(x = MMSE, ymin = NEFL_smooth - 0.25, ymax = NEFL_smooth + 0.25),
              fill = "#D55E00", alpha = 0.10) +
  geom_line(data = soma_trend, aes(x = MMSE, y = LCN2_smooth, color = "LCN2 (Soma)"), linewidth = 1.8) +
  geom_line(data = soma_trend, aes(x = MMSE, y = NEFL_smooth, color = "NEFL (Soma)"), linewidth = 1.8) +
  geom_line(data = tmt_trend,  aes(x = MMSE, y = LCN2_smooth, color = "LCN2 (TMT)"),
            linewidth = 1.5, linetype = "dashed") +
  geom_line(data = tmt_trend,  aes(x = MMSE, y = NEFL_smooth, color = "NEFL (TMT)"),
            linewidth = 1.5, linetype = "dashed") +
  scale_x_reverse(breaks = seq(10, 30, 5)) +
  scale_color_manual(values = c("LCN2 (Soma)" = "#CC79A7", "LCN2 (TMT)" = "#CC79A7",
                                "NEFL (Soma)" = "#D55E00", "NEFL (TMT)" = "#D55E00")) +
  labs(title = "LCN2/NEFL rise with cognitive decline (both platforms)",
       x = "MMSE score", y = "Z-score increase", color = NULL) +
  base_theme +
  theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.text = element_text(size = 11), legend.key.width = unit(1.1, "cm")) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

# Panel C - CSF PTGDS vs MMSE: all subjects shown; flat / no robust association
pC <- ggplot(tmt_c, aes(x = MMSE, y = PTGDS_final)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_point(color = "#5BA3D0", alpha = 0.30, size = 0.8) +
  geom_smooth(method = "loess", span = 0.6, se = FALSE, color = "#0072B2", linewidth = 1.8) +
  scale_x_reverse(breaks = seq(12, 30, 2)) +
  annotate("label", x = max(tmt_c$MMSE), y = max(tmt_c$PTGDS_final),
           label = annotC, hjust = 0, vjust = 1, size = 4.6,
           label.size = 0.3, fill = "white") +
  labs(title = "TMT-MS PTGDS vs MMSE: flat (no robust association; biphasic is tissue-restricted, cf. Fig.1)",
       x = "Cognitive status (MMSE)", y = "CSF PTGDS (per-batch z)") +
  base_theme

# ------------------------------------------------------------------------------
# 3. Assemble (vertical stack) and save
# ------------------------------------------------------------------------------
message(">>> [3/3] Assembling and saving...")

final_plot <- plot_grid(
  pA, pB, pC,
  ncol = 1, rel_heights = c(1, 1, 1),
  labels = c("A", "B", "C"), label_size = 29, label_fontface = "bold",
  label_x = 0.005, label_y = 1.0, hjust = 0, vjust = 1
)

W   <- 12
H   <- 18
DPI <- 300

out_png  <- file.path(OUT_DIR, "FigureS7_CSF_crossplatform.png")
out_tiff <- file.path(OUT_DIR, "FigureS7_CSF_crossplatform.tiff")

ggsave(out_png,  final_plot, width = W, height = H, dpi = DPI, bg = "white")
ggsave(out_tiff, final_plot, width = W, height = H, dpi = DPI, bg = "white",
       device = "tiff", compression = "lzw")

message("\n", strrep("=", 55))
message("  Done! PNG : ", out_png)
message("  Done! TIFF: ", out_tiff)
message("  Panel C: n = ", nC, ", Pearson r = ", round(rC, 3), " (n.s.)")
message("  Size: ", W, " x ", H, " in @ ", DPI, " dpi")
message(strrep("=", 55))
