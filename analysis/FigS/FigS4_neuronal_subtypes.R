# ==============================================================================
# FigS4_neuronal_subtypes.R
#
# Purpose : Additional file 1: Fig. S4
#           Differential vulnerability and apoptotic buffering capacity
#           across neuronal subtypes
#
# Panels
#   A  Proportional changes of Excitatory vs SST+ neurons across CPS bins
#      (SST+: 11.24% -> 5.49% among classified subtypes)
#   B  NDUFS1 expression variance: SST+ higher than Excitatory (0.19-0.26 vs 0.12-0.16)
#   C  BCL2 trajectory in SST+ neurons: 6.8% decline (Bin 0.1->0.4)
#   D  BAX/BCL2 ratio in SST+ neurons: progressive increase
#
# Matches Supp caption:
#   "SST+ neurons exhibit a progressive decline (from 11.24% to 5.49% among
#    classified subtypes), while excitatory neuron proportion increases
#    correspondingly."
#
# NOTE: FigS4 code has incomplete sections (sst_bcl2, sst_ratio objects
#       not computed before use in Panels C/D). This script adds the missing
#       computation steps.
#
# Input   : output/SEAAD_processed_data.RData  (from 01_SEAAD_data_extraction.R)
# Output  : output/FigS4/FigS4.png + .tiff  (600 dpi, 9 x 8.5 in)
#
# Paper: Kim Y et al. (2025). bioRxiv. https://github.com/YoungOukKim/MCI-to-AD
# Usage: Rscript analysis/FigS/FigS4_neuronal_subtypes.R
# Requirements: R >= 4.3.2 | data.table, ggplot2, patchwork
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
set.seed(42)

# Paths (relative to repository root)
OUT_DIR <- file.path("output", "FigS4")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. Load data
# ==============================================================================
message(">>> Loading SEAAD_processed_data.RData...")
load(file.path("output", "SEAAD_processed_data.RData"))

# Ensure SST classification exists (add if missing)
if (!"SST_Neuron" %in% unique(df$cell_type_label)) {
  df[grepl("Sst|SST", cell_type, ignore.case = TRUE),
     cell_type_label := "SST_Neuron"]
}

neurons <- df[cell_type_label %in% c("Exc_Neuron", "SST_Neuron") &
                bin_rounded >= 0.1 & bin_rounded <= 0.9]

message(sprintf("  Neurons: %d total (Exc: %d, SST: %d)",
                nrow(neurons),
                neurons[cell_type_label == "Exc_Neuron", .N],
                neurons[cell_type_label == "SST_Neuron",  .N]))

# ==============================================================================
# 2. Panel A: Subtype proportions per bin
# ==============================================================================
total_per_bin <- neurons[, .(total = .N), by = bin_rounded]
sst_per_bin   <- neurons[cell_type_label == "SST_Neuron",  .(SST_n = .N), by = bin_rounded]
exc_per_bin   <- neurons[cell_type_label == "Exc_Neuron",  .(Exc_n = .N), by = bin_rounded]

prop_dt <- Reduce(function(a, b) merge(a, b, by="bin_rounded", all.x=TRUE),
                  list(total_per_bin, sst_per_bin, exc_per_bin))
prop_dt[is.na(SST_n), SST_n := 0L]
prop_dt[is.na(Exc_n), Exc_n := 0L]
prop_dt[, SST_pct := SST_n / total * 100]
prop_dt[, Exc_pct := Exc_n / total * 100]

message(sprintf("  SST%%: %.2f%% (Bin 0.1) -> %.2f%% (Bin 0.9)",
                prop_dt[bin_rounded == 0.1, SST_pct],
                prop_dt[bin_rounded == 0.9, SST_pct]))

df_prop <- rbind(
  data.table(bin=prop_dt$bin_rounded, cell_class="SST+ Inhibitory",
             proportion=prop_dt$SST_pct),
  data.table(bin=prop_dt$bin_rounded, cell_class="Excitatory Neuron",
             proportion=prop_dt$Exc_pct))

# ==============================================================================
# 3. Panel B: NDUFS1 variance per subtype
# ==============================================================================
exc_var <- neurons[cell_type_label == "Exc_Neuron",
                   .(NDUFS1_var = var(NDUFS1, na.rm=TRUE)), by = bin_rounded]
sst_var <- neurons[cell_type_label == "SST_Neuron",
                   .(NDUFS1_var = var(NDUFS1, na.rm=TRUE)), by = bin_rounded]

var_data <- data.frame(
  cell_class = c("Excitatory", "SST+ Inhibitory"),
  mean_var   = c(mean(exc_var$NDUFS1_var, na.rm=TRUE),
                 mean(sst_var$NDUFS1_var, na.rm=TRUE)),
  min_var    = c(min(exc_var$NDUFS1_var,  na.rm=TRUE),
                 min(sst_var$NDUFS1_var,  na.rm=TRUE)),
  max_var    = c(max(exc_var$NDUFS1_var,  na.rm=TRUE),
                 max(sst_var$NDUFS1_var,  na.rm=TRUE)))

message(sprintf("  NDUFS1 var  Exc: %.2f-%.2f | SST: %.2f-%.2f",
                var_data$min_var[1], var_data$max_var[1],
                var_data$min_var[2], var_data$max_var[2]))

# ==============================================================================
# 4. Panels C & D: SST+ BCL2 trajectory and BAX/BCL2 ratio
# ==============================================================================
sst_cells <- neurons[cell_type_label == "SST_Neuron"]

apoptosis_genes <- intersect(c("BCL2","BAX","CASP3"), colnames(sst_cells))
if (length(apoptosis_genes) == 0) {
  stop("BCL2/BAX not found in data. Add to target_genes in Script 01 and re-run.")
}

sst_bcl2  <- sst_cells[, .(BCL2_mean = mean(BCL2, na.rm=TRUE)),
                        by = .(bin_rounded)][order(bin_rounded)]
sst_ratio <- sst_cells[, .(BAX_BCL2_ratio = mean(BAX / (BCL2 + 1e-6), na.rm=TRUE)),
                        by = .(bin_rounded)][order(bin_rounded)]

bcl2_bin01  <- sst_bcl2[bin_rounded == 0.1, BCL2_mean]
bcl2_bin04  <- sst_bcl2[bin_rounded == 0.4, BCL2_mean]
bcl2_decline <- (bcl2_bin01 - bcl2_bin04) / bcl2_bin01 * 100

message(sprintf("  BCL2 decline (Bin 0.1->0.4): %.1f%%", bcl2_decline))

# ==============================================================================
# 5. Build panels
# ==============================================================================
theme_s4 <- theme_classic(base_size = 10) +
  theme(plot.title   = element_text(hjust=0.5, size=12, face="bold",
                                     margin=margin(b=5)),
        axis.title   = element_text(size=10, face="bold"),
        axis.text    = element_text(color="black", size=9),
        plot.margin  = margin(5,8,5,2),
        legend.position = "top",
        legend.title    = element_blank(),
        legend.text     = element_text(size=9),
        legend.margin   = margin(b=-8))

x_sc <- scale_x_continuous(breaks=seq(0.1,0.9,0.2))

# Panel A
pA <- ggplot(df_prop, aes(x=bin, y=proportion, color=cell_class, group=cell_class)) +
  geom_smooth(method="loess", span=0.6, se=FALSE, linewidth=1.2) +
  geom_point(size=2, alpha=0.6) +
  scale_color_manual(values=c("Excitatory Neuron"="#2C3E50",
                               "SST+ Inhibitory"="#8E44AD")) +
  x_sc +
  labs(x="Pseudo-progression", y="Proportion (%)",
       title="Neuronal Subtype Proportions") + theme_s4

# Panel B
pB <- ggplot(var_data, aes(x=cell_class, y=mean_var, fill=cell_class)) +
  geom_bar(stat="identity", width=0.5, alpha=0.8) +
  geom_errorbar(aes(ymin=min_var, ymax=max_var), width=0.1, linewidth=0.7) +
  scale_fill_manual(values=c("Excitatory"="#2C3E50",
                              "SST+ Inhibitory"="#8E44AD")) +
  labs(x="", y="NDUFS1 Expression Variance",
       title="Metabolic Heterogeneity (NDUFS1)") +
  theme_s4 + theme(legend.position="none",
                   axis.text.x=element_text(size=9))

# Panel C — BCL2 trajectory in SST+
bcl2_label <- sprintf("-%.1f%% Decline", bcl2_decline)
pC <- ggplot(sst_bcl2, aes(x=bin_rounded, y=BCL2_mean)) +
  geom_smooth(method="loess", span=1.0, se=FALSE,
              color="#8E44AD", linewidth=1.2) +
  geom_point(size=2, color="#8E44AD", alpha=0.5) +
  annotate("segment", x=0.15, xend=0.35,
           y=bcl2_bin01+0.01, yend=bcl2_bin04+0.01,
           color="red", arrow=arrow(length=unit(0.15,"cm")), linewidth=0.8) +
  annotate("text", x=0.25, y=bcl2_bin01+0.03,
           label=bcl2_label, color="red", fontface="bold", size=3) +
  x_sc +
  labs(x="Pseudo-progression", y="BCL2 Expression",
       title="SST+ Anti-apoptotic Buffering (BCL2)") + theme_s4

# Panel D — BAX/BCL2 ratio
pD <- ggplot(sst_ratio, aes(x=bin_rounded, y=BAX_BCL2_ratio)) +
  geom_smooth(method="loess", span=1.0, se=FALSE,
              color="darkred", linewidth=1.2) +
  geom_point(size=2, color="darkred", alpha=0.5) +
  x_sc +
  labs(x="Pseudo-progression", y="BAX / BCL2 Ratio",
       title="SST+ Apoptotic Stress Index") + theme_s4

# Combine
fig_s4 <- (pA | pB) / (pC | pD) +
  plot_annotation(tag_levels='A') &
  theme(plot.tag=element_text(face="bold", size=20))

out_png  <- file.path(OUT_DIR, "FigS4.png")
out_tiff <- file.path(OUT_DIR, "FigS4.tiff")
ggsave(out_png,  fig_s4, width=9, height=8.5, dpi=600, bg="white")
ggsave(out_tiff, fig_s4, width=9, height=8.5, dpi=600,
       device="tiff", compression="lzw", bg="white")

message(">>> Saved: ", out_png)
cat(strrep("=", 60), "\n")
cat(sprintf("  Panel A: SST%% %.2f%% -> %.2f%%\n",
            prop_dt[bin_rounded==0.1, SST_pct],
            prop_dt[bin_rounded==0.9, SST_pct]))
cat(sprintf("  Panel C: BCL2 decline = %.1f%%\n", bcl2_decline))
cat(strrep("=", 60), "\n")
message(">>> Script FigS4 complete.")
