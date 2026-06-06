# ==============================================================================
# FigS4_neuronal_subtypes.R
#
# Purpose : Generate Figure S4 — Differential vulnerability across
#           neuronal subtypes in SEA-AD pseudo-progression analysis.
#
# Panels:
#   A = Excitatory vs SST+ inhibitory neuron proportions per bin
#   B = NDUFS1 expression variance by subtype (metabolic heterogeneity)
#   C = SST+ BCL2 trajectory (apoptotic buffering)
#   D = SST+ BAX/BCL2 ratio trajectory (apoptotic stress index)
#
# Outputs (-> output/FigS4/):
#   FigS4_final.png   (300 dpi)
#   FigS4_final.tiff  (300 dpi, LZW — journal submission)
#
# Data note:
#   SEAAD_processed_data.RData contains cell-level data (~1.3M nuclei)
#   and cannot be redistributed (Allen Brain Atlas data use policy).
#   Pre-computed summary CSV files are provided in data/SEA-AD/:
#     FigS4_proportions.csv  — bin-level subtype proportions
#     FigS4_variance.csv     — NDUFS1 variance by subtype
#     FigS4_bcl2.csv         — SST+ BCL2 trajectory
#     FigS4_ratio.csv        — SST+ BAX/BCL2 ratio trajectory
#   To regenerate from raw data, uncomment write.csv() lines at the end.
#
# Paper : Kim Y et al. (2025). https://github.com/YoungOukKim/MCI-to-AD
# Usage : Rscript analysis/FigS/FigS4_neuronal_subtypes.R
# Requirements : R >= 4.3.2 | data.table, ggplot2, patchwork
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)

# ------------------------------------------------------------------------------
# Paths (relative to repository root; override with environment variables)
# ------------------------------------------------------------------------------
DATA_DIR <- Sys.getenv("FIGS4_DATA", unset = "data/SEA-AD")
OUT_DIR  <- Sys.getenv("FIGS4_OUT",  unset = "output/FigS4")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Load pre-computed summary CSVs
#    (generated from SEAAD_processed_data.RData — see data note above)
# ------------------------------------------------------------------------------
message(">>> [1/4] Loading summary data...")

df_prop  <- fread(file.path(DATA_DIR, "FigS4_proportions.csv"))
var_data <- fread(file.path(DATA_DIR, "FigS4_variance.csv"))
sst_bcl2 <- fread(file.path(DATA_DIR, "FigS4_bcl2.csv"))
sst_ratio <- fread(file.path(DATA_DIR, "FigS4_ratio.csv"))

# BCL2 decline annotation (Bin 0.1 → 0.4)
bcl2_bin01  <- sst_bcl2[bin_rounded == 0.1, BCL2_mean]
bcl2_bin04  <- sst_bcl2[bin_rounded == 0.4, BCL2_mean]
bcl2_decline <- (bcl2_bin01 - bcl2_bin04) / bcl2_bin01 * 100

message(sprintf("  BCL2 decline (Bin 0.1->0.4): %.1f%%", bcl2_decline))

# ------------------------------------------------------------------------------
# 2. Shared theme
# ------------------------------------------------------------------------------
theme_s4 <- theme_classic(base_size = 10) +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 12, face = "bold",
                                margin = margin(b = 5)),
    axis.title   = element_text(size = 10, face = "bold"),
    axis.text    = element_text(color = "black", size = 9),
    plot.margin  = margin(t = 5, r = 8, b = 5, l = 2),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text  = element_text(size = 9),
    legend.margin = margin(b = -8)
  )

# ------------------------------------------------------------------------------
# 3. Panel construction
# ------------------------------------------------------------------------------
message(">>> [2/4] Building panels...")

# Panel A — Subtype proportions per bin
pA <- ggplot(df_prop,
             aes(x = bin, y = proportion,
                 color = cell_class, group = cell_class)) +
  geom_smooth(method = "loess", span = 0.6, se = FALSE, linewidth = 1.2) +
  geom_point(size = 2, alpha = 0.6) +
  scale_color_manual(
    values = c("Excitatory_Neuron" = "#2C3E50",
               "SST_Inhibitory"   = "#8E44AD"),
    labels = c("Excitatory_Neuron" = "Excitatory",
               "SST_Inhibitory"   = "SST+")
  ) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) +
  labs(x = "Pseudo-progression bin",
       y = "Proportion (%)",
       title = "Neuron subtype proportions") +
  theme_s4

# Panel B — NDUFS1 variance by subtype
pB <- ggplot(var_data,
             aes(x = cell_class, y = mean_var, fill = cell_class)) +
  geom_bar(stat = "identity", width = 0.5, alpha = 0.8) +
  geom_errorbar(aes(ymin = min_var, ymax = max_var),
                width = 0.1, linewidth = 0.7) +
  scale_fill_manual(values = c("Excitatory"     = "#2C3E50",
                                "SST+ Inhibitory" = "#8E44AD")) +
  labs(x = "", y = "Expression variance",
       title = "NDUFS1 metabolic heterogeneity") +
  theme_s4 +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 9))

# Panel C — SST+ BCL2 trajectory
bcl2_label <- sprintf("-%.1f%% decline", bcl2_decline)

pC <- ggplot(sst_bcl2, aes(x = bin_rounded, y = BCL2_mean)) +
  geom_smooth(method = "loess", span = 1.0, se = FALSE,
              color = "#8E44AD", linewidth = 1.2) +
  geom_point(size = 2, color = "#8E44AD", alpha = 0.5) +
  annotate("segment",
           x = 0.15, xend = 0.35,
           y = bcl2_bin01 + 0.01, yend = bcl2_bin04 + 0.01,
           color = "red",
           arrow = arrow(length = unit(0.15, "cm")), linewidth = 0.8) +
  annotate("text", x = 0.25, y = bcl2_bin01 + 0.03,
           label = bcl2_label,
           color = "red", fontface = "bold", size = 3) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) +
  labs(x = "Pseudo-progression bin",
       y = "BCL2 expression",
       title = "SST+ apoptotic buffering") +
  theme_s4

# Panel D — SST+ BAX/BCL2 ratio
pD <- ggplot(sst_ratio, aes(x = bin_rounded, y = BAX_BCL2_ratio)) +
  geom_smooth(method = "loess", span = 1.0, se = FALSE,
              color = "darkred", linewidth = 1.2) +
  geom_point(size = 2, color = "darkred", alpha = 0.5) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) +
  labs(x = "Pseudo-progression bin",
       y = "BAX / BCL2 ratio",
       title = "SST+ apoptotic stress index") +
  theme_s4

# ------------------------------------------------------------------------------
# 4. Assemble and save
# ------------------------------------------------------------------------------
message(">>> [3/4] Assembling figure...")

fig_s4 <- (pA | pB) / (pC | pD) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag          = element_text(face = "bold", size = 20),
        plot.tag.position = c(0.03, 0.99))

W   <- 9
H   <- 8.5
DPI <- 300

out_png  <- file.path(OUT_DIR, "FigS4_final.png")
out_tiff <- file.path(OUT_DIR, "FigS4_final.tiff")

ggsave(out_png,  fig_s4,
       width = W, height = H, dpi = DPI, bg = "white")
ggsave(out_tiff, fig_s4,
       width = W, height = H, dpi = DPI, bg = "white",
       device = "tiff", compression = "lzw")

message("\n", strrep("=", 55))
message("  Done! PNG : ", out_png)
message("  Done! TIFF: ", out_tiff)
message("  Size: ", W, " x ", H, " in @ ", DPI, " dpi")
message(strrep("=", 55))

# ------------------------------------------------------------------------------
# 5. (Optional) Regenerate summary CSVs from raw RData
#    Uncomment and run once if SEAAD_processed_data.RData is available
# ------------------------------------------------------------------------------
# message(">>> [4/4] Regenerating summary CSVs from raw data...")
# load("data/SEA-AD/SEAAD_processed_data.RData")
#
# neurons <- df[cell_type_label %in% c("Exc_Neuron", "SST_Neuron") &
#               bin_rounded >= 0.1 & bin_rounded <= 0.9]
#
# # Proportions
# total_per_bin <- neurons[, .N, by = bin_rounded]; setnames(total_per_bin, "N", "total")
# sst_n <- neurons[cell_type_label == "SST_Neuron", .N, by = bin_rounded]; setnames(sst_n, "N", "SST_n")
# exc_n <- neurons[cell_type_label == "Exc_Neuron", .N, by = bin_rounded]; setnames(exc_n, "N", "Exc_n")
# prop_dt <- merge(merge(total_per_bin, sst_n, by = "bin_rounded", all.x = TRUE),
#                  exc_n, by = "bin_rounded", all.x = TRUE)
# prop_dt[is.na(SST_n), SST_n := 0]; prop_dt[is.na(Exc_n), Exc_n := 0]
# prop_dt[, `:=`(SST_pct = SST_n/total*100, Exc_pct = Exc_n/total*100)]
# df_prop_out <- rbind(
#   data.table(bin = prop_dt$bin_rounded, cell_class = "SST_Inhibitory",   proportion = prop_dt$SST_pct),
#   data.table(bin = prop_dt$bin_rounded, cell_class = "Excitatory_Neuron", proportion = prop_dt$Exc_pct)
# )
# write.csv(df_prop_out, file.path(DATA_DIR, "FigS4_proportions.csv"), row.names = FALSE)
#
# # NDUFS1 variance
# exc_var <- neurons[cell_type_label == "Exc_Neuron", .(NDUFS1_var = var(NDUFS1, na.rm=TRUE)), by=bin_rounded]
# sst_var <- neurons[cell_type_label == "SST_Neuron", .(NDUFS1_var = var(NDUFS1, na.rm=TRUE)), by=bin_rounded]
# var_out <- data.frame(
#   cell_class = c("Excitatory", "SST+ Inhibitory"),
#   mean_var   = c(mean(exc_var$NDUFS1_var), mean(sst_var$NDUFS1_var)),
#   min_var    = c(min(exc_var$NDUFS1_var),  min(sst_var$NDUFS1_var)),
#   max_var    = c(max(exc_var$NDUFS1_var),  max(sst_var$NDUFS1_var))
# )
# write.csv(var_out, file.path(DATA_DIR, "FigS4_variance.csv"), row.names = FALSE)
#
# # BCL2 trajectory
# sst_bcl2_out <- neurons[cell_type_label == "SST_Neuron",
#                         .(BCL2_mean = mean(BCL2, na.rm=TRUE)), by=bin_rounded][order(bin_rounded)]
# write.csv(sst_bcl2_out, file.path(DATA_DIR, "FigS4_bcl2.csv"), row.names = FALSE)
#
# # BAX/BCL2 ratio
# sst_ratio_out <- neurons[cell_type_label == "SST_Neuron",
#                          .(BAX_BCL2_ratio = mean(BAX, na.rm=TRUE)/mean(BCL2, na.rm=TRUE)),
#                          by=bin_rounded][order(bin_rounded)]
# write.csv(sst_ratio_out, file.path(DATA_DIR, "FigS4_ratio.csv"), row.names = FALSE)
#
# message("  CSVs saved to ", DATA_DIR)
