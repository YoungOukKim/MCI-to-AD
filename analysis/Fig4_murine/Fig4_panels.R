# ==============================================================================
# Fig4_panels.R
#
# Purpose : Generate individual Fig. 4 panels (A–F), NO panel labels.
#           Labels are added separately in Fig4_combine.R
#
# Outputs (all -> OUT_DIR)
#   Fig4A_ymaze.png
#   Fig4B_PAT.png
#   Fig4C_inflammatory.png      (LCN2, TNF-alpha, IL-6, IL-1beta)
#   Fig4E_trophic.png           (NGFR, BDNF)
#   Fig4F_oxidative.png         (Keap1, APOE, ABCA1)
#
# Panel mapping:
#   A = Y-maze spontaneous alternation (%)
#   B = Passive Avoidance Test — retention latency (sec)
#   C = Astrocytic LCN2 + pro-inflammatory cytokines (qPCR)
#   D = (immunohistochemistry images — supplied separately)
#   E = Neuronal NGFR and BDNF (qPCR)
#   F = Keap1, APOE, ABCA1 (qPCR)
#
# Data source : data/murine/murine_behavioral_molecular_data.xlsx
# Paper       : Kim Y et al. (2025). https://github.com/YoungOukKim/MCI-to-AD
# Usage       : Rscript analysis/Fig4_murine/Fig4_panels.R
# Requirements: R >= 4.3.2 | ggplot2, dplyr, readxl, tidyr, ggpubr
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readxl)
  library(tidyr)
})

# ------------------------------------------------------------------------------
# Paths  (override with environment variables or edit defaults)
# All defaults are relative to the repository root.
# ------------------------------------------------------------------------------
DATA_FILE <- Sys.getenv("MURINE_DATA",
  unset = "data/murine/murine_behavioral_molecular_data.xlsx")
OUT_DIR   <- Sys.getenv("FIG4_OUT", unset = "output/Fig4")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Shared colour palette  (7 groups)
# ------------------------------------------------------------------------------
GROUPS <- c("Sham",
            "A\u03b2 (20\u03bcM)",
            "BXP-101 50mg/kg",
            "BXP-101 100mg/kg",
            "BXP-101 200mg/kg",
            "BXP-101 400mg/kg",
            "Donepezil 5mg/kg")

COL7 <- c(
  "Sham"              = "#4472C4",
  "A\u03b2 (20\u03bcM)"  = "#ED7D31",
  "BXP-101 50mg/kg"   = "#C5E0B4",
  "BXP-101 100mg/kg"  = "#92D050",
  "BXP-101 200mg/kg"  = "#70AD47",
  "BXP-101 400mg/kg"  = "#375623",
  "Donepezil 5mg/kg"  = "#A9A9A9"
)

# PCR groups (4 groups only)
PCR_GROUPS <- c("Sham","A\u03b2 (20\u03bcM)",
                "BXP-101 200mg/kg","BXP-101 400mg/kg")
COL4 <- COL7[PCR_GROUPS]

# ------------------------------------------------------------------------------
# Shared theme
# ------------------------------------------------------------------------------
BASE <- theme_classic(base_size = 16) +
  theme(
    axis.text.x     = element_text(size = 13, angle = 35, hjust = 1),
    axis.text.y     = element_text(size = 13),
    axis.title      = element_text(size = 15),
    legend.position = "none",
    plot.margin     = margin(10, 10, 10, 10)
  )

DPI   <- 600
WHITE <- "white"

# ==============================================================================
# PANEL A — Y-maze spontaneous alternation
# ==============================================================================
message(">>> Panel A ...")

ym <- read_excel(DATA_FILE, sheet = "Ab_Y-maze_summary") %>%
  rename(Group = 1, n = 2, Mean = 3, SEM = 4, p_Sham = 5) %>%
  filter(!is.na(Mean)) %>%
  mutate(
    Group = factor(trimws(Group), levels = GROUPS),
    Mean  = as.numeric(Mean),
    SEM   = as.numeric(SEM)
  )

# significance brackets (from spreadsheet)
sig_ab <- data.frame(
  x1 = c(1, 4, 7), x2 = c(2, 2, 2),
  y  = c(42, 45, 48),
  lab = c("*", "**", "**")
)

pA <- ggplot(ym, aes(x = Group, y = Mean, fill = Group)) +
  geom_bar(stat = "identity", width = 0.65, alpha = 0.85,
           color = "gray30", linewidth = 0.4) +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM),
                width = 0.25, linewidth = 0.7) +
  # significance: * Sham vs Ab
  annotate("segment", x=1, xend=2, y=42, yend=42,
           color="black", linewidth=0.6) +
  annotate("text", x=1.5, y=43.2, label="*", size=7, color="red") +
  # ** Ab vs BXP-101 100
  annotate("segment", x=2, xend=4, y=45, yend=45,
           color="black", linewidth=0.6) +
  annotate("text", x=3, y=46.2, label="**", size=6, color="blue") +
  # ** Ab vs DNZ
  annotate("segment", x=2, xend=7, y=49, yend=49,
           color="black", linewidth=0.6) +
  annotate("text", x=4.5, y=50.2, label="**", size=6, color="blue") +
  scale_fill_manual(values = COL7) +
  coord_cartesian(ylim = c(0, 55)) +
  labs(x = NULL, y = "Spontaneous Alternation (%)") +
  BASE

ggsave(file.path(OUT_DIR, "Fig4A_ymaze.png"), pA,
       width = 8, height = 6.5, dpi = DPI, bg = WHITE)
message("  Saved Fig4A_ymaze.png")

# ==============================================================================
# PANEL B — Passive Avoidance Test (Retention latency)
# ==============================================================================
message(">>> Panel B ...")

pat <- read_excel(DATA_FILE, sheet = "Ab_PAT_summary", skip = 1) %>%
  select(Group = 1, n = 2,
         Acq_mean = 3, Acq_SEM = 4,
         Ret_mean = 5, Ret_SEM = 6) %>%
  filter(!is.na(Ret_mean)) %>%
  mutate(
    Group    = factor(trimws(Group), levels = GROUPS),
    Ret_mean = as.numeric(Ret_mean),
    Ret_SEM  = as.numeric(Ret_SEM)
  )

pB <- ggplot(pat, aes(x = Group, y = Ret_mean, fill = Group)) +
  geom_bar(stat = "identity", width = 0.65, alpha = 0.85,
           color = "gray30", linewidth = 0.4) +
  geom_errorbar(aes(ymin = Ret_mean - Ret_SEM, ymax = Ret_mean + Ret_SEM),
                width = 0.25, linewidth = 0.7) +
  # * Sham vs Ab
  annotate("segment", x=1, xend=2, y=270, yend=270,
           color="black", linewidth=0.6) +
  annotate("text", x=1.5, y=277, label="*", size=7, color="red") +
  # significance: BXP-101 400 vs Ab
  annotate("segment", x=2, xend=6, y=290, yend=290,
           color="black", linewidth=0.6) +
  annotate("text", x=4, y=297, label="*", size=7, color="blue") +
  scale_fill_manual(values = COL7) +
  coord_cartesian(ylim = c(0, 310)) +
  labs(x = NULL, y = "Retention Latency (sec)") +
  BASE

ggsave(file.path(OUT_DIR, "Fig4B_PAT.png"), pB,
       width = 8, height = 6.5, dpi = DPI, bg = WHITE)
message("  Saved Fig4B_PAT.png")

# ==============================================================================
# PANELS C, E, F — PCR molecular markers
# ==============================================================================
message(">>> Panels C + E + F ...")

pcr_raw <- read_excel(DATA_FILE, sheet = "Ab_PCR_molecular") %>%
  rename(Marker=1, Panel=2, Group=3, Rep1=4, Rep2=5, Rep3=6, Mean=7, SEM=8) %>%
  filter(!is.na(Mean), !is.na(Group)) %>%
  mutate(
    Group  = factor(trimws(Group), levels = PCR_GROUPS),
    Mean   = as.numeric(Mean),
    SEM    = as.numeric(SEM),
    Marker = trimws(Marker)
  ) %>%
  filter(!is.na(Group))

# Fill down Marker column (merged cells read as NA after first row)
pcr_raw$Marker[pcr_raw$Marker == ""] <- NA
pcr_raw <- pcr_raw %>%
  tidyr::fill(Marker, .direction = "down") %>%
  tidyr::fill(Panel,  .direction = "down")

# Helper: multi-marker faceted barplot
make_pcr_plot <- function(markers, ylab, width_in = 10) {
  df <- pcr_raw %>%
    filter(Marker %in% markers) %>%
    mutate(Marker = factor(Marker, levels = markers))

  ggplot(df, aes(x = Group, y = Mean, fill = Group)) +
    geom_bar(stat = "identity", width = 0.6, alpha = 0.85,
             color = "gray30", linewidth = 0.4) +
    geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM),
                  width = 0.22, linewidth = 0.7) +
    geom_hline(yintercept = 1, linetype = "dashed",
               color = "gray50", linewidth = 0.5) +
    facet_wrap(~Marker, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = COL4) +
    labs(x = NULL, y = ylab) +
    BASE +
    theme(
      strip.text  = element_text(face = "bold.italic", size = 14),
      axis.text.x = element_text(size = 11, angle = 40, hjust = 1)
    )
}

# Panel C — inflammatory markers
pC <- make_pcr_plot(
  c("LCN2","TNF-\u03b1","IL-6","IL-1\u03b2"),
  "Relative expression (vs Sham)",
  width_in = 13
)
ggsave(file.path(OUT_DIR, "Fig4C_inflammatory.png"), pC,
       width = 13, height = 5.5, dpi = DPI, bg = WHITE)
message("  Saved Fig4C_inflammatory.png")

# Panel E — trophic markers
pE <- make_pcr_plot(
  c("NGFR","BDNF"),
  "Relative expression (vs Sham)",
  width_in = 8
)
ggsave(file.path(OUT_DIR, "Fig4E_trophic.png"), pE,
       width = 8, height = 5.5, dpi = DPI, bg = WHITE)
message("  Saved Fig4E_trophic.png")

# Panel F — oxidative / lipid markers
pF <- make_pcr_plot(
  c("Keap1","APOE","ABCA1"),
  "Relative expression (vs Sham)",
  width_in = 11
)
ggsave(file.path(OUT_DIR, "Fig4F_oxidative.png"), pF,
       width = 11, height = 5.5, dpi = DPI, bg = WHITE)
message("  Saved Fig4F_oxidative.png")

message("\n", strrep("=", 60))
message("  All panels saved to: ", OUT_DIR)
message("  Next step: run Fig4_combine.R to assemble with labels.")
message(strrep("=", 60))
