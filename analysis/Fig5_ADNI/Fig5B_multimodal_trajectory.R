# ==============================================================================
# Fig5B_multimodal_trajectory.R
#
# Purpose : Figure 5B — Multimodal protein trajectory across MMSE spectrum
#           BDNF / PTGDS / LCN2 / NEFL / Tau LOESS curves
#           Phase boundary marker at MMSE 26
#
# Input   : SomaScan 7K CSV  +  DXSUM.rda  +  MMSE.rda  (ADNI)
# Output  : output/Fig5/Fig5B_multimodal_trajectory.png  (600 dpi)
#
# Key result
#   Phase boundary (MMSE 26) = PTGDS/LCN2 crossover point identified by
#   segmented regression (Fig. 5C); shown here as vertical reference line.
#
# ADNI data: register at https://adni.loni.usc.edu
# Paper: Kim Y et al. (2025). bioRxiv. https://github.com/YoungOukKim/MCI-to-AD
# Usage: Rscript analysis/Fig5_ADNI/Fig5B_multimodal_trajectory.R
# Requirements: R >= 4.3.2 | tidyverse
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse) })
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
  mutate(DX_Group = case_when(
    DIAGNOSIS == "CN"       ~ "CN",
    DIAGNOSIS == "MCI"      ~ "MCI",
    DIAGNOSIS == "Dementia" ~ "AD",
    TRUE ~ NA_character_)) %>%
  filter(!is.na(DX_Group)) %>%
  select(RID, DX_Group)

mmse_clean <- MMSE %>%
  mutate(RID = as.character(RID), MMSE = as.numeric(MMSCORE)) %>%
  filter(!is.na(MMSE)) %>% group_by(RID) %>%
  summarise(MMSE_Score = mean(MMSE, na.rm = TRUE), .groups = "drop")

df_fig5b <- somascan %>%
  mutate(RID = as.character(RID)) %>%
  select(RID,
    PTGDS = X10514.5, NEFL = X10082.251,
    LCN2  = X2836.68, Tau  = X5854.60, BDNF = X2421.7) %>%
  inner_join(dx_data,    by = "RID") %>%
  inner_join(mmse_clean, by = "RID") %>%
  mutate(across(c(PTGDS, NEFL, LCN2, Tau, BDNF),
                ~ scale(log2(.))[, 1], .names = "{.col}_z"))

# Legend labels
leg_bdnf  <- "BDNF (Neurotrophic Factor)"
leg_ptgds <- "PTGDS (Metabolic Defense)"
leg_lcn2  <- "LCN2 (Inflammation/Stress)"
leg_nefl  <- "NEFL (Axonal Injury)"
leg_tau   <- "Tau (Irreversible Pathology)"

message(">>> Building Figure 5B...")
p5b <- ggplot(df_fig5b, aes(x = MMSE_Score)) +
  geom_smooth(aes(y = Tau_z,   color = leg_tau),
              method="loess", span=0.7, linewidth=1.0, linetype="dotted",
              se=TRUE, alpha=0.1, fill="#999999", na.rm=TRUE) +
  geom_smooth(aes(y = NEFL_z,  color = leg_nefl),
              method="loess", span=0.7, linewidth=2.2,
              se=TRUE, alpha=0.15, fill="grey70") +
  geom_smooth(aes(y = LCN2_z,  color = leg_lcn2),
              method="loess", span=0.7, linewidth=1.8,
              se=TRUE, alpha=0.12, fill="#CC79A7") +
  geom_smooth(aes(y = PTGDS_z, color = leg_ptgds),
              method="loess", span=0.7, linewidth=1.8,
              se=TRUE, alpha=0.12, fill="#0072B2") +
  geom_smooth(aes(y = BDNF_z,  color = leg_bdnf),
              method="loess", span=0.7, linewidth=2.0,
              se=TRUE, alpha=0.15, fill="#009E73") +
  # Phase boundary at MMSE 26
  geom_vline(xintercept=26, linetype="dashed", color="darkred", linewidth=0.8) +
  annotate("rect", xmin=25.5, xmax=26.5, ymin=-0.8, ymax=0.8,
           alpha=0.15, fill="darkred") +
  annotate("text", x=27.9, y=0.75,
           label="Phase Boundary (MMSE 26)",
           color="darkred", fontface="bold", size=4, hjust=0) +
  coord_cartesian(xlim=c(30, 18), ylim=c(-0.8, 0.8)) +
  scale_x_reverse(breaks=seq(18, 30, by=2)) +
  scale_color_manual(values=c(
    "BDNF (Neurotrophic Factor)"   = "#009E73",
    "PTGDS (Metabolic Defense)"    = "#0072B2",
    "LCN2 (Inflammation/Stress)"   = "#CC79A7",
    "NEFL (Axonal Injury)"         = "#D55E00",
    "Tau (Irreversible Pathology)" = "#999999"),
    breaks=c(leg_bdnf, leg_ptgds, leg_lcn2, leg_nefl, leg_tau)) +
  labs(x="Cognitive Status (MMSE Score)", y="Relative Deviation (z-score)") +
  theme_classic(base_size=18) +
  theme(legend.position="top", legend.title=element_blank(),
        legend.text=element_text(size=9.5),
        plot.margin=margin(30,15,15,15)) +
  guides(color=guide_legend(nrow=2, byrow=TRUE))

out_file <- file.path(OUT_DIR, "Fig5B_multimodal_trajectory.png")
ggsave(out_file, plot=p5b, width=10, height=7, dpi=600, bg="white")
message(">>> Saved: ", out_file)
message(">>> Script Fig5B complete.")
