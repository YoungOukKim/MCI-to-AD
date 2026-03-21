# ==============================================================================
# Fig5A_PTGDS_clinical_stages.R
#
# Purpose : Figure 5A — PTGDS expression across clinical stages (boxplot)
#           Healthy CN  /  Early Peak (Compensatory)  /  Late-stage Collapse
#
# Input   : SomaScan 7K CSV  +  DXSUM.rda  +  MMSE.rda  (ADNI)
# Output  : output/Fig5/Fig5A_PTGDS_clinical_stages.png  (600 dpi)
#
# Clinical stage classification
#   Healthy CN          : DIAGNOSIS == "CN"  &  MMSE >= 29.5
#   Early Peak (Comp.)  : MMSE 26-28  (pre-boundary compensatory peak)
#   Late-stage Collapse : DIAGNOSIS == "Dementia"  OR  MMSE < 24
#
# ADNI data: register at https://adni.loni.usc.edu
# Paper: Kim Y et al. (2025). bioRxiv. https://github.com/YoungOukKim/MCI-to-AD
# Usage: Rscript analysis/Fig5_ADNI/Fig5A_PTGDS_clinical_stages.R
# Requirements: R >= 4.3.2 | tidyverse, ggplot2
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})
set.seed(42)

# Paths (relative to repository root)
SOMA_CSV  <- Sys.getenv("SOMA_CSV",  unset = "path/to/CruchagaLab_CSF_SOMAscan7k.csv")
ADNI_PATH <- Sys.getenv("ADNI_PATH", unset = "path/to/ADNIMERGE2/data/")
OUT_DIR   <- file.path("output", "Fig5")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Load data
message(">>> Loading data...")
somascan <- read.csv(SOMA_CSV, stringsAsFactors = FALSE)
load(file.path(ADNI_PATH, "DXSUM.rda"))
load(file.path(ADNI_PATH, "MMSE.rda"))

# Merge and classify
message(">>> Processing clinical stages...")
df_fig5a <- somascan %>%
  mutate(RID = as.character(RID)) %>%
  select(RID, PTGDS_Raw = X10514.5) %>%
  inner_join(
    DXSUM %>% mutate(RID = as.character(RID)) %>%
      group_by(RID) %>% slice_tail(n = 1) %>% ungroup(),
    by = "RID") %>%
  inner_join(
    MMSE %>% mutate(RID = as.character(RID), MMSE = as.numeric(MMSCORE)) %>%
      filter(!is.na(MMSE)) %>% group_by(RID) %>%
      summarise(MMSE_Score = mean(MMSE, na.rm = TRUE), .groups = "drop"),
    by = "RID") %>%
  mutate(
    PTGDS_z = scale(log2(PTGDS_Raw))[, 1],
    Stage = case_when(
      DIAGNOSIS == "CN"       & MMSE_Score >= 29.5 ~ "Healthy CN",
      MMSE_Score >= 26        & MMSE_Score  < 29   ~ "Early Peak\n(Compensatory)",
      DIAGNOSIS == "Dementia" | MMSE_Score  < 24   ~ "Late-stage\nCollapse",
      TRUE ~ NA_character_)) %>%
  filter(!is.na(PTGDS_z), !is.na(Stage)) %>%
  mutate(Stage = factor(Stage,
    levels = c("Healthy CN", "Early Peak\n(Compensatory)", "Late-stage\nCollapse")))

message(sprintf("  N = %d", nrow(df_fig5a)))
print(table(df_fig5a$Stage))

# Build plot
message(">>> Building Figure 5A...")
p5a <- ggplot(df_fig5a, aes(x = Stage, y = PTGDS_z, fill = Stage)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.5,
               lwd = 0.8, color = "#333333", fatten = 0) +
  stat_summary(fun = mean, geom = "errorbar",
               aes(ymax = after_stat(y), ymin = after_stat(y)),
               width = 0.5, lwd = 1.0, color = "#333333") +
  geom_jitter(width = 0.12, alpha = 0.12, size = 1.1, color = "grey40") +
  stat_summary(fun = mean, geom = "line",  aes(group = 1),
               linewidth = 1.2, color = "firebrick") +
  stat_summary(fun = mean, geom = "point", size = 4.5,
               color = "firebrick", shape = 18) +
  scale_fill_manual(values = c(
    "Healthy CN"                  = "#99CCFF",
    "Early Peak\n(Compensatory)"  = "#FFCC99",
    "Late-stage\nCollapse"        = "#FF9999")) +
  coord_cartesian(ylim = c(-1.5, 1.5)) +
  labs(y = "Relative Deviation (z-score)", x = NULL) +
  theme_classic(base_size = 18) +
  theme(
    legend.position  = "none",
    axis.text.x      = element_text(face = "bold", size = 13, color = "black"),
    axis.title.y     = element_text(face = "plain", size = 15),
    plot.margin      = margin(25, 10, 10, 10))

# Save
out_file <- file.path(OUT_DIR, "Fig5A_PTGDS_clinical_stages.png")
ggsave(out_file, plot = p5a, width = 7, height = 6, dpi = 600, bg = "white")
message(">>> Saved: ", out_file)
message(">>> Script Fig5A complete.")
