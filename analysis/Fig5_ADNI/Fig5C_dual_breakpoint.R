# ==============================================================================
# Fig5C_dual_breakpoint.R
#
# Purpose : Figure 5C — Dual breakpoint analysis (MMSE spectrum)
#           PTGDS / NEFL / TNF / ICAM1 z-score trajectories
#           Phase boundary: MMSE 26 (defense collapse)
#           Injury acceleration: MMSE ~24.7 (segmented regression breakpoint)
#
# Input   : SomaScan 7K CSV  +  DXSUM.rda  +  MMSE.rda  (ADNI)
# Output  : output/Fig5/Fig5C_dual_breakpoint.png  (600 dpi)
#
# Key result
#   Segmented regression on NEFL trajectory identifies injury acceleration
#   breakpoint at MMSE ~24.7 (95% CI shown as shaded zone).
#   Corresponds to Fig. 5 in paper (NEFL inflection: MMSE 24.7, p < 0.001).
#
# ADNI data: register at https://adni.loni.usc.edu
# Paper: Kim Y et al. (2025). bioRxiv. https://github.com/YoungOukKim/MCI-to-AD
# Usage: Rscript analysis/Fig5_ADNI/Fig5C_dual_breakpoint.R
# Requirements: R >= 4.3.2 | tidyverse, segmented
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(segmented)
})
set.seed(42)

SOMA_CSV  <- Sys.getenv("SOMA_CSV",  unset = "path/to/CruchagaLab_CSF_SOMAscan7k.csv")
ADNI_PATH <- Sys.getenv("ADNI_PATH", unset = "path/to/ADNIMERGE2/data/")
OUT_DIR   <- file.path("output", "Fig5")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

message(">>> Loading data...")
somascan <- read.csv(SOMA_CSV, stringsAsFactors = FALSE)
load(file.path(ADNI_PATH, "DXSUM.rda"))
load(file.path(ADNI_PATH, "MMSE.rda"))

message(">>> Processing...")
dx_data <- DXSUM %>%
  mutate(RID = as.character(RID)) %>%
  group_by(RID) %>% slice_tail(n = 1) %>% ungroup() %>%
  select(RID, DIAGNOSIS)

mmse_clean <- MMSE %>%
  mutate(RID = as.character(RID), MMSE = as.numeric(MMSCORE)) %>%
  filter(!is.na(MMSE)) %>% group_by(RID) %>%
  summarise(MMSE_Score = mean(MMSE, na.rm = TRUE), .groups = "drop")

df_fig5c <- somascan %>%
  mutate(RID = as.character(RID)) %>%
  select(RID,
    PTGDS = X10514.5, NEFL = X10082.251,
    TNF   = X5692.79, ICAM1 = X4342.10) %>%
  inner_join(dx_data,    by = "RID") %>%
  inner_join(mmse_clean, by = "RID") %>%
  mutate(across(c(PTGDS, NEFL, TNF, ICAM1),
                ~ scale(log2(.))[, 1], .names = "{.col}_z")) %>%
  mutate(TNF_sync   = TNF_z   + 0.2,
         ICAM1_sync = ICAM1_z + 0.1,
         NEFL_sync  = NEFL_z)

# Segmented regression: NEFL injury acceleration breakpoint
message(">>> Segmented regression on NEFL (injury breakpoint)...")
fit_nefl <- lm(NEFL_sync ~ MMSE_Score, data = df_fig5c)
seg_nefl <- segmented(fit_nefl, seg.Z = ~ MMSE_Score, psi = 25)
bp_est   <- seg_nefl$psi[1, "Est."]
bp_se    <- seg_nefl$psi[1, "St.Err"]
rect_xmin <- bp_est - 1.96 * bp_se
rect_xmax <- bp_est + 1.96 * bp_se
message(sprintf("  NEFL breakpoint: MMSE %.1f  (95%% CI: %.1f-%.1f)",
                bp_est, rect_xmin, rect_xmax))

leg_ptgds <- "PTGDS (Metabolic Defense)"
leg_nefl  <- "NEFL (Axonal Injury)"
leg_icam1 <- "ICAM1 (Vascular Stress)"
leg_tnf   <- "TNF (Inflammation)"

message(">>> Building Figure 5C...")
p5c <- ggplot(df_fig5c, aes(x = MMSE_Score)) +
  geom_smooth(aes(y = PTGDS_z,   color = leg_ptgds),
              method="loess", span=0.7, se=TRUE, alpha=0.1, fill="#0072B2") +
  geom_smooth(aes(y = ICAM1_sync, color = leg_icam1),
              method="loess", span=0.7, se=FALSE, linewidth=1.2) +
  geom_smooth(aes(y = TNF_sync,   color = leg_tnf),
              method="loess", span=0.7, se=FALSE, linewidth=1.2) +
  geom_smooth(aes(y = NEFL_sync,  color = leg_nefl),
              method="loess", span=0.7, linewidth=2.3,
              se=TRUE, alpha=0.25, fill="grey70") +
  # 95% CI shaded zone for injury breakpoint
  annotate("rect", xmin=rect_xmin, xmax=rect_xmax, ymin=-0.8, ymax=0.8,
           alpha=0.15, fill="#D55E00") +
  # Phase boundary at MMSE 26
  geom_vline(xintercept=26,     linetype="dashed", color="darkred",  linewidth=0.8) +
  # Injury acceleration at segmented breakpoint
  geom_vline(xintercept=bp_est, linetype="dashed", color="#D55E00", linewidth=0.8) +
  # Labels
  annotate("text", x=27.9, y=0.75,
           label="Phase Boundary (MMSE 26)",
           color="darkred", fontface="bold", size=4, hjust=0) +
  annotate("label", x=20.5, y=-0.6,
           fill="white", color="#D55E00", fontface="bold", size=4.5,
           label=paste0("Failure Point: MMSE ", round(bp_est, 1),
                        "\n(95% CI: ", round(rect_xmin,1), "-",
                        round(rect_xmax,1), ")")) +
  coord_cartesian(xlim=c(30, 18), ylim=c(-0.8, 0.8)) +
  scale_x_reverse(breaks=seq(18, 30, by=2)) +
  scale_color_manual(values=c(
    "PTGDS (Metabolic Defense)" = "#0072B2",
    "NEFL (Axonal Injury)"      = "#D55E00",
    "ICAM1 (Vascular Stress)"   = "#882255",
    "TNF (Inflammation)"        = "#CC79A7"),
    breaks=c(leg_ptgds, leg_nefl, leg_icam1, leg_tnf)) +
  labs(x="Cognitive Status (MMSE Score)", y="Relative Deviation (z-score)") +
  theme_classic(base_size=18) +
  theme(legend.position="top", legend.title=element_blank(),
        legend.text=element_text(size=9.5),
        plot.margin=margin(30,15,15,15)) +
  guides(color=guide_legend(nrow=2, byrow=TRUE))

out_file <- file.path(OUT_DIR, "Fig5C_dual_breakpoint.png")
ggsave(out_file, plot=p5c, width=10, height=7, dpi=600, bg="white")
message(">>> Saved: ", out_file)
message(">>> Script Fig5C complete.")
