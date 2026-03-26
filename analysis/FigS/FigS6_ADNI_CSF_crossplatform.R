# ==============================================================================
# FigS6_ADNI_CSF_crossplatform.R
#
# Purpose : Generate Figure S6 — Cross-platform validation of PTGDS/LCN2/NEFL
#           dynamics in ADNI CSF proteomics (TMT-MS and SomaScan 7K).
#
# Panels:
#   A = Platform correlation: TMT-MS vs SomaScan PTGDS (scatter by clinical group)
#   B = Cross-platform LCN2 and NEFL elevation across MMSE decline
#   C = TMT-MS biphasic PTGDS dynamics across cognitive status
#
# Outputs (-> output/FigS6/):
#   FigS6_final.png   (300 dpi)
#   FigS6_final.tiff  (300 dpi, LZW — journal submission)
#
# Data note:
#   Raw ADNI data (TMT-MS, SomaScan, ADNIMERGE) cannot be redistributed
#   under the ADNI data use agreement (https://adni.loni.usc.edu).
#   Pre-computed summary CSVs are provided in data/ADNI/:
#     FigS6_platform_cor.csv  — Panel A: RID, TMT_log2, Soma_log2, Group
#     FigS6_trend.csv         — Panel B: MMSE, LCN2_smooth, NEFL_smooth, platform
#     FigS6_tmt_ptgds.csv     — Panel C: MMSE, PTGDS_final
#   To regenerate, uncomment the data preparation section at the end.
#
# Paper : Kim Y et al. (2025). https://github.com/YoungOukKim/MCI-to-AD
# Usage : Rscript analysis/FigS/FigS6_ADNI_CSF_crossplatform.R
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
DATA_DIR <- Sys.getenv("FIGS6_DATA", unset = "data/ADNI")
OUT_DIR  <- Sys.getenv("FIGS6_OUT",  unset = "output/FigS6")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Load pre-computed summary CSVs
# ------------------------------------------------------------------------------
message(">>> [1/3] Loading summary data...")

cor_data  <- read.csv(file.path(DATA_DIR, "FigS6_platform_cor.csv"))
trend_data <- read.csv(file.path(DATA_DIR, "FigS6_trend.csv"))
tmt_c     <- read.csv(file.path(DATA_DIR, "FigS6_tmt_ptgds.csv"))

soma_trend <- trend_data %>% filter(platform == "Soma")
tmt_trend  <- trend_data %>% filter(platform == "TMT")

message("  Panel A: ", nrow(cor_data), " subjects")
message("  Panel B: Soma=", nrow(soma_trend), " TMT=", nrow(tmt_trend), " rows")
message("  Panel C: ", nrow(tmt_c), " rows")

# ------------------------------------------------------------------------------
# 2. Panel construction
# ------------------------------------------------------------------------------
message(">>> [2/3] Building panels...")

# Panel A — Platform correlation (PTGDS)
pA <- ggplot(cor_data, aes(x = TMT_log2, y = Soma_log2)) +
  geom_point(aes(color = Group), alpha = 0.4, size = 1.8) +
  geom_smooth(method = "lm", color = "black") +
  scale_color_manual(values = c(
    "Healthy CN"    = "#99CCFF",
    "Early Peak"    = "#FFCC99",
    "Late Collapse" = "#FF9999"
  )) +
  labs(title    = "Platform correlation (PTGDS)",
       x        = "TMT-MS (log2)",
       y        = "SomaScan (log2)",
       color    = "Clinical group") +
  theme_classic(base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", hjust = 0.5))

# Panel B — Cross-platform LCN2 / NEFL elevation
pB <- ggplot() +
  # Soma ribbons
  geom_ribbon(data = soma_trend,
              aes(x = MMSE, ymin = LCN2_smooth - 0.25, ymax = LCN2_smooth + 0.25),
              fill = "#CC79A7", alpha = 0.18) +
  geom_ribbon(data = soma_trend,
              aes(x = MMSE, ymin = NEFL_smooth - 0.25, ymax = NEFL_smooth + 0.25),
              fill = "#D55E00", alpha = 0.18) +
  # TMT ribbons
  geom_ribbon(data = tmt_trend,
              aes(x = MMSE, ymin = LCN2_smooth - 0.25, ymax = LCN2_smooth + 0.25),
              fill = "#CC79A7", alpha = 0.10) +
  geom_ribbon(data = tmt_trend,
              aes(x = MMSE, ymin = NEFL_smooth - 0.25, ymax = NEFL_smooth + 0.25),
              fill = "#D55E00", alpha = 0.10) +
  # Lines: solid = Soma, dashed = TMT
  geom_line(data = soma_trend,
            aes(x = MMSE, y = LCN2_smooth, color = "LCN2 (Soma)"), linewidth = 1.8) +
  geom_line(data = soma_trend,
            aes(x = MMSE, y = NEFL_smooth, color = "NEFL (Soma)"), linewidth = 1.8) +
  geom_line(data = tmt_trend,
            aes(x = MMSE, y = LCN2_smooth, color = "LCN2 (TMT)"),
            linewidth = 1.5, linetype = "dashed") +
  geom_line(data = tmt_trend,
            aes(x = MMSE, y = NEFL_smooth, color = "NEFL (TMT)"),
            linewidth = 1.5, linetype = "dashed") +
  scale_x_reverse(breaks = seq(0, 30, 5)) +
  scale_color_manual(values = c(
    "LCN2 (Soma)" = "#CC79A7", "LCN2 (TMT)" = "#CC79A7",
    "NEFL (Soma)" = "#D55E00", "NEFL (TMT)" = "#D55E00"
  )) +
  coord_cartesian(ylim = c(-0.2, 2.6)) +
  labs(title = "Consistent elevation (Soma: \u2014, TMT: --)",
       x     = "MMSE score",
       y     = "Z-score increase",
       color = "Markers (platform)") +
  theme_classic(base_size = 14) +
  theme(legend.position  = "top",
        legend.key.width = unit(1.2, "cm"),
        legend.text      = element_text(size = 9),
        plot.title       = element_text(face = "bold", hjust = 0.5)) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

# Panel C — TMT-MS biphasic PTGDS dynamics
pC <- ggplot(tmt_c, aes(x = MMSE)) +
  geom_smooth(aes(y = PTGDS_final), color = "#0072B2",
              method = "loess", span = 0.6, linewidth = 2.0) +
  geom_vline(xintercept = 26, linetype = "dashed", color = "darkred") +
  scale_x_reverse(breaks = seq(18, 30, 2)) +
  labs(title = "TMT-MS validation: biphasic PTGDS dynamics",
       x     = "Cognitive status (MMSE)",
       y     = "Relative deviation") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# ------------------------------------------------------------------------------
# 3. Assemble and save
# ------------------------------------------------------------------------------
message(">>> [3/3] Assembling and saving...")

top_row <- plot_grid(pA, pB, ncol = 2, rel_widths = c(1.25, 1), scale = 0.95)

final_plot <- ggdraw() +
  draw_plot(top_row, x = 0,    y = 0.46, width = 1,   height = 0.52) +
  draw_plot(pC,      x = 0.05, y = 0.04, width = 0.9, height = 0.36) +
  draw_label("A", x = 0.02, y = 0.97, fontface = "bold", size = 36) +
  draw_label("B", x = 0.54, y = 0.97, fontface = "bold", size = 36) +
  draw_label("C", x = 0.02, y = 0.41, fontface = "bold", size = 36)

W   <- 18
H   <- 13
DPI <- 300

out_png  <- file.path(OUT_DIR, "FigS6_final.png")
out_tiff <- file.path(OUT_DIR, "FigS6_final.tiff")

ggsave(out_png,  final_plot,
       width = W, height = H, dpi = DPI, bg = "white")
ggsave(out_tiff, final_plot,
       width = W, height = H, dpi = DPI, bg = "white",
       device = "tiff", compression = "lzw")

message("\n", strrep("=", 55))
message("  Done! PNG : ", out_png)
message("  Done! TIFF: ", out_tiff)
message("  Size: ", W, " x ", H, " in @ ", DPI, " dpi")
message(strrep("=", 55))

# ==============================================================================
# (Optional) Regenerate summary CSVs from raw ADNI data
# Requires: EMORY_CSF_TMT_MS.csv, CruchagaLab SomaScan CSV, ADNIMERGE2 RData
# Uncomment and run once if raw data is available
# ==============================================================================
# message(">>> Regenerating summary CSVs from raw ADNI data...")
# library(tidyverse)
# tmt_raw  <- read.csv("data/ADNI/EMORY_CSF_TMT_MS.csv")
# soma_raw <- read.csv("data/ADNI/CruchagaLab_CSF_SOMAscan7k_postQC.csv")
# load("data/ADNI/ADNIMERGE2/DXSUM.rda")
# load("data/ADNI/ADNIMERGE2/MMSE.rda")
#
# mmse_clean <- MMSE %>%
#   mutate(RID = as.character(RID), MMSE = as.numeric(MMSCORE)) %>%
#   filter(!is.na(MMSE)) %>% distinct(RID, MMSE, .keep_all = TRUE)
#
# dx_ref <- DXSUM %>% mutate(RID = as.character(RID)) %>%
#   group_by(RID) %>% slice_tail(n=1) %>% ungroup()
#
# # Panel A
# cor_out <- inner_join(
#   tmt_raw  %>% mutate(RID=as.character(RID)) %>% select(RID, TMT=PTGDS_P41222),
#   soma_raw %>% mutate(RID=as.character(RID)) %>% select(RID, Soma=X10514.5),
#   by="RID") %>%
#   inner_join(dx_ref, by="RID") %>%
#   mutate(Group = case_when(
#     DIAGNOSIS=="CN"      ~ "Healthy CN",
#     DIAGNOSIS=="Dementia"~ "Late Collapse",
#     TRUE                 ~ "Early Peak"),
#     TMT_log2  = log2(TMT  + 1e-4),
#     Soma_log2 = log2(Soma + 1e-4))
# write.csv(cor_out %>% select(RID, TMT_log2, Soma_log2, Group),
#           file.path(DATA_DIR, "FigS6_platform_cor.csv"), row.names=FALSE)
#
# # Panel B (trend data) — run prepare_b_data() then bind rows
# # Panel C (tmt_c_data) — see original derivation
# # write.csv(..., file.path(DATA_DIR, "FigS6_trend.csv"), row.names=FALSE)
# # write.csv(..., file.path(DATA_DIR, "FigS6_tmt_ptgds.csv"), row.names=FALSE)
