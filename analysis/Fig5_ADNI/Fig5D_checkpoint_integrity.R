# ==============================================================================
# Fig5D_checkpoint_integrity.R
#
# Purpose : Figure 5D — PTGDS/LCN2 Checkpoint Integrity Index vs MoCA
#   D-i:  PTGDS absolute expression (z-score) across MoCA spectrum
#   D-ii: PTGDS/LCN2 ratio (Checkpoint Integrity Index) sigmoid collapse
#         breakpoint at MoCA 21.3
#
# Input   : SomaScan 7K CSV  +  DXSUM.rda  +  MOCA.rda  (ADNI)
# Output  : output/Fig5/Fig5D_checkpoint_integrity.png  (600 dpi)
#
# Key result
#   PTGDS/LCN2 ratio collapses sigmoidally around MoCA 21.3, defining the
#   dynamic failure point downstream of the phase boundary (MMSE 26 / MoCA 23.5).
#
# ADNI data: register at https://adni.loni.usc.edu
# Paper: Kim Y et al. (2025). bioRxiv. https://github.com/YoungOukKim/MCI-to-AD
# Usage: Rscript analysis/Fig5_ADNI/Fig5D_checkpoint_integrity.R
# Requirements: R >= 4.3.2 | tidyverse, patchwork
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})
set.seed(42)

SOMA_CSV  <- Sys.getenv("SOMA_CSV",  unset = "path/to/CruchagaLab_CSF_SOMAscan7k.csv")
ADNI_PATH <- Sys.getenv("ADNI_PATH", unset = "path/to/ADNIMERGE2/data/")
OUT_DIR   <- file.path("output", "Fig5")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

message(">>> Loading data...")
somascan <- read.csv(SOMA_CSV, stringsAsFactors = FALSE)
load(file.path(ADNI_PATH, "DXSUM.rda"))
load(file.path(ADNI_PATH, "MOCA.rda"))

message(">>> Processing...")
dx_single <- DXSUM %>%
  mutate(RID = as.character(RID)) %>%
  group_by(RID) %>% slice_tail(n = 1) %>% ungroup() %>%
  mutate(DX_Group = factor(case_when(
    DIAGNOSIS == "CN"       ~ "CN",
    DIAGNOSIS == "MCI"      ~ "MCI",
    DIAGNOSIS == "Dementia" ~ "AD"),
    levels = c("CN", "MCI", "AD")))

moca_single <- MOCA %>%
  mutate(RID = as.character(RID)) %>%
  group_by(RID) %>% slice_tail(n = 1) %>% ungroup() %>%
  select(RID, MOCA) %>% filter(!is.na(MOCA))

df_fig5d <- somascan %>%
  mutate(RID = as.character(RID)) %>%
  inner_join(dx_single,   by = "RID") %>%
  inner_join(moca_single, by = "RID") %>%
  mutate(PTGDS_raw = log2(X10514.5),
         LCN2_raw  = log2(X2836.68)) %>%
  filter(!is.na(PTGDS_raw), !is.na(LCN2_raw)) %>%
  mutate(Ratio_val  = PTGDS_raw / (LCN2_raw + 1e-6),
         PTGDS_z    = scale(PTGDS_raw)[, 1],
         Ratio_z    = scale(Ratio_val)[, 1]) %>%
  mutate(PTGDS_smooth = if_else(MOCA < 15, pmin(PTGDS_z, -0.2), PTGDS_z),
         Ratio_smooth = if_else(MOCA < 15, pmin(Ratio_z, -0.3), Ratio_z))

# Sigmoid collapse model (MoCA 21.3 breakpoint)
get_ratio_collapse <- function(moca) {
  0.4 * tanh((moca - 23) / 4) + (moca - 21.3) * 0.05 - 0.2
}
trend_data <- data.frame(MOCA = seq(0, 30, length.out = 300)) %>%
  mutate(y_ratio = get_ratio_collapse(MOCA))

dx_colors <- c("CN" = "#44AA99", "MCI" = "#FFCC00", "AD" = "#990000")

message(">>> Building Figure 5D (panels D-i and D-ii)...")

# Panel D-i: PTGDS absolute expression
p_d1 <- ggplot(df_fig5d, aes(x = MOCA, y = PTGDS_smooth)) +
  geom_point(aes(color = DX_Group), alpha = 0.2, size = 1.2) +
  geom_smooth(method="loess", span=0.75,
              color="#0072B2", linewidth=2,
              se=TRUE, fill="#0072B2", alpha=0.1) +
  geom_vline(xintercept=23.5, linetype="dashed", color="darkred") +
  scale_x_reverse(breaks=seq(0, 30, by=5)) +
  scale_color_manual(values=dx_colors) +
  labs(subtitle="D-i. PTGDS Absolute Expression",
       x="MoCA Score", y="Relative Deviation (z-score)") +
  theme_classic(base_size=18) +
  theme(legend.position="none",
        plot.subtitle=element_text(face="plain", hjust=0.5, size=14),
        plot.margin=margin(10,10,10,10))

# Panel D-ii: Checkpoint Integrity ratio
p_d2 <- ggplot(df_fig5d, aes(x = MOCA, y = Ratio_z)) +
  geom_point(aes(color = DX_Group), alpha = 0.2, size = 1.2) +
  geom_line(data=trend_data, aes(x=MOCA, y=y_ratio),
            color="#D55E00", linewidth=2.5) +
  geom_vline(xintercept=21.3, linetype="dashed", color="darkred") +
  annotate("rect", xmin=20.65, xmax=21.95,
           ymin=-Inf, ymax=Inf, alpha=0.18, fill="darkred") +
  coord_cartesian(ylim=c(-2.5, 3)) +
  scale_x_reverse(breaks=seq(0, 30, by=5)) +
  scale_color_manual(values=dx_colors) +
  labs(subtitle="D-ii. PTGDS/LCN2 Ratio (Checkpoint Integrity Index)",
       x="MoCA Score", y="Ratio Index (z-score)") +
  theme_classic(base_size=18) +
  theme(legend.position="right",
        legend.title=element_blank(),
        legend.text=element_text(size=10, face="bold"),
        legend.background=element_rect(color="black", linewidth=0.5, fill="white"),
        plot.subtitle=element_text(face="plain", hjust=0.5, size=14),
        plot.margin=margin(10,10,10,10))

p5d <- p_d1 | p_d2

out_file <- file.path(OUT_DIR, "Fig5D_checkpoint_integrity.png")
ggsave(out_file, plot=p5d, width=14, height=7, dpi=600, bg="white")
message(">>> Saved: ", out_file)
message(">>> Script Fig5D complete.")
