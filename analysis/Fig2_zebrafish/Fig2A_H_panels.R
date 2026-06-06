# ==============================================================================
# Fig2A_H_panels.R
#
# Purpose : Generate individual Fig. 2 panels (A-H), NO panel labels.
#           Labels are added separately in Fig2_combine.R
#
# Outputs (all -> OUT_DIR)
#   Fig2A_timeline.png
#   Fig2B_avoidance.png
#   Fig2C_BLBP_quant.png
#   Fig2D_Nestin_quant.png
#   Fig2E_WIF_images.png        (3x3 grid assembled from 9 raw images)
#   Fig2F_qPCR.png
#   Fig2G_ELISA.png
#   Fig2H_schematic.png
#
# Panel mapping:
#   A = Experimental timeline schematic
#   B = Red ball avoidance (14 dpf)
#   C = BLBP quantification (WIF, reactive gliosis)
#   D = Nestin quantification (WIF, neural progenitor)
#   E = Representative WIF images (3x3 grid)
#   F = ptgdsb.1/2 longitudinal qPCR (14 + 21 dpf)
#   G = ELISA: TNF-alpha + IL-6 (14 dpf)
#   H = Molecular recovery schematic
#      (NF-kB suppression -> ptgdsb capacity -> TNF-a/IL-6 down -> NGFR re-expression)
#
# Paper : Kim Y et al. (2025). https://github.com/YoungOukKim/MCI-to-AD
# Usage : Rscript analysis/Fig2_zebrafish/Fig2A_H_panels.R
# Requirements : R >= 4.3.2 | ggplot2, dplyr, tidyr, magick
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(magick)
})

# ------------------------------------------------------------------------------
# Paths (relative to repository root; override with environment variables)
# ------------------------------------------------------------------------------
DATA_DIR <- Sys.getenv("ZEBRAFISH_DATA", unset = "data/zebrafish")
IMG_DIR  <- Sys.getenv("WIF_IMG_DIR",    unset = "data/zebrafish/WIF")
OUT_DIR  <- Sys.getenv("FIG2_OUT",       unset = "output/Fig2")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Shared colour palette
# ------------------------------------------------------------------------------
COL <- c(
  Control   = "#4472C4",
  MCI       = "#ED7D31",
  Donepezil = "#A9A9A9",
  BXP_0.3   = "#C5E0B4",
  BXP_0.4   = "#70AD47",
  BXP_0.6   = "#375623"
)

LV6 <- c(
  "Control",
  "MCI",
  "MCI+Donepezil\n(10 \u03bcM)",
  "MCI+BXP-101\n(0.3 \u03bcg/ml)",
  "MCI+BXP-101\n(0.4 \u03bcg/ml)",
  "MCI+BXP-101\n(0.6 \u03bcg/ml)"
)

LV3 <- c("Control", "MCI", "MCI+BXP-101\n(0.4 \u03bcg/ml)")

COL6 <- setNames(
  unname(COL[c("Control","MCI","Donepezil","BXP_0.3","BXP_0.4","BXP_0.6")]),
  LV6
)

COL3 <- c(
  "Control"                        = unname(COL["Control"]),
  "MCI"                            = unname(COL["MCI"]),
  "MCI+BXP-101\n(0.4 \u03bcg/ml)" = unname(COL["BXP_0.4"])
)

# ------------------------------------------------------------------------------
# Shared ggplot2 theme
# ------------------------------------------------------------------------------
BASE <- theme_classic(base_size = 18) +
  theme(
    axis.text.x     = element_text(size = 14, angle = 30, hjust = 1),
    axis.text.y     = element_text(size = 14),
    axis.title      = element_text(size = 16),
    plot.title      = element_text(size = 17, face = "bold", hjust = 0.5),
    legend.position = "none",
    plot.margin     = margin(10, 10, 10, 10)
  )

DPI   <- 600
WHITE <- "white"

# ==============================================================================
# PANEL A — Experimental timeline
# ==============================================================================
message(">>> Panel A ...")

pA <- ggplot() +
  annotate("segment", x = 0, xend = 23, y = 0, yend = 0,
           color = "black", linewidth = 1.2,
           arrow = arrow(length = unit(0.25, "cm"), ends = "last")) +
  geom_vline(xintercept = c(2, 4, 6, 14, 15, 21),
             linetype = "dashed", color = "gray60", linewidth = 0.5) +
  annotate("text",
           x = c(0, 2, 4, 6, 14, 15, 21), y = -0.25,
           label = c("0","2","4","6","14","15","21"),
           size = 5.5, color = "black") +
  annotate("text", x = 23, y = -0.25, label = "dpf",
           size = 5.5, fontface = "italic") +
  annotate("rect", xmin = 2,  xmax = 14, ymin = 1.65, ymax = 1.95,
           fill = "#E41A1C", alpha = 0.75) +
  annotate("rect", xmin = 4,  xmax = 14, ymin = 1.25, ymax = 1.55,
           fill = "#FF7F00", alpha = 0.75) +
  annotate("rect", xmin = 6,  xmax = 14, ymin = 0.85, ymax = 1.15,
           fill = COL["BXP_0.4"], alpha = 0.80) +
  annotate("rect", xmin = 15, xmax = 21, ymin = 0.85, ymax = 1.15,
           fill = "#377EB8", alpha = 0.80) +
  annotate("text", x = 8,  y = 1.80,
           label = "LPS (2\u201314 dpf)",
           size = 4.5, color = "white", fontface = "bold") +
  annotate("text", x = 9,  y = 1.40,
           label = "D-galactose + HFD (4\u201314 dpf)",
           size = 4.5, color = "white", fontface = "bold") +
  annotate("text", x = 10, y = 1.00,
           label = "BXP-101 early (6\u201314 dpf)",
           size = 4.5, color = "white", fontface = "bold") +
  annotate("text", x = 18, y = 1.00,
           label = "BXP-101 late (15\u201321 dpf)",
           size = 4.5, color = "white", fontface = "bold") +
  annotate("point", x = 14, y = 0.45, size = 4, color = "#984EA3") +
  annotate("point", x = 21, y = 0.45, size = 4, color = "#984EA3") +
  annotate("text", x = 14, y = 0.62,
           label = "Behavior \u00b7 qPCR \u00b7 ELISA \u00b7 WIF",
           size = 4.5, color = "#984EA3") +
  annotate("text", x = 21, y = 0.62,
           label = "qPCR", size = 4.5, color = "#984EA3") +
  scale_y_continuous(limits = c(-0.4, 2.2)) +
  scale_x_continuous(limits = c(-0.5, 24)) +
  theme_void() +
  theme(plot.margin = margin(10, 10, 10, 10))

ggsave(file.path(OUT_DIR, "Fig2A_timeline.png"), pA,
       width = 10, height = 3, dpi = DPI, bg = WHITE)
message("  Saved Fig2A_timeline.png")

# ==============================================================================
# PANEL B — Red ball avoidance
# ==============================================================================
message(">>> Panel B ...")

beh <- read.csv(file.path(DATA_DIR, "behavior_14dpf.csv"),
                stringsAsFactors = FALSE)

beh_avg <- beh %>%
  filter(Well == "Average") %>%
  select(-Well) %>%
  pivot_longer(-Replicate, names_to = "Group", values_to = "Pct") %>%
  mutate(
    Group = dplyr::recode(Group,
      "Control"            = "Control",
      "MCI"                = "MCI",
      "MCI_Donepezil_10uM" = "MCI+Donepezil\n(10 \u03bcM)",
      "MCI_BXP101_0.3ugml" = "MCI+BXP-101\n(0.3 \u03bcg/ml)",
      "MCI_BXP101_0.4ugml" = "MCI+BXP-101\n(0.4 \u03bcg/ml)",
      "MCI_BXP101_0.6ugml" = "MCI+BXP-101\n(0.6 \u03bcg/ml)"),
    Group = factor(Group, levels = LV6)
  )

beh_sum <- beh_avg %>%
  group_by(Group) %>%
  summarise(Mean = mean(Pct), SD = sd(Pct), .groups = "drop")

pB <- ggplot(beh_sum, aes(x = Group, y = Mean, fill = Group)) +
  geom_bar(stat = "identity", width = 0.65, alpha = 0.85,
           color = "gray30", linewidth = 0.4) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0.25, linewidth = 0.7) +
  geom_jitter(data = beh_avg, aes(x = Group, y = Pct, fill = Group),
              width = 0.12, size = 2, alpha = 0.6,
              shape = 21, color = "gray30") +
  annotate("segment", x = 1, xend = 2, y = 83, yend = 83,
           color = "black", linewidth = 0.6) +
  annotate("text", x = 1.5, y = 84.5, label = "*",
           size = 7, color = "red") +
  annotate("segment", x = 2, xend = 6, y = 88, yend = 88,
           color = "black", linewidth = 0.6) +
  annotate("text", x = 4, y = 89.5, label = "#",
           size = 7, color = "blue") +
  scale_fill_manual(values = COL6) +
  coord_cartesian(ylim = c(0, 95)) +
  labs(x = NULL, y = "% larvae in non-stimuli area") +
  BASE +
  theme(plot.margin = margin(8, 8, 20, 8))

ggsave(file.path(OUT_DIR, "Fig2B_avoidance.png"), pB,
       width = 7, height = 7, dpi = DPI, bg = WHITE)
message("  Saved Fig2B_avoidance.png")

# ==============================================================================
# PANELS C & D — BLBP and Nestin quantification (WIF)
# ==============================================================================
message(">>> Panels C + D ...")

wif_long <- read.csv(file.path(DATA_DIR, "WIF_14dpf.csv"),
                     stringsAsFactors = FALSE) %>%
  filter(Larva != "Average") %>%
  mutate(Larva = as.numeric(Larva)) %>%
  pivot_longer(c(Control, MCI, MCI_BXP101_0.4ugml),
               names_to = "Group", values_to = "Intensity") %>%
  mutate(
    Group = dplyr::recode(Group,
      "Control"            = "Control",
      "MCI"                = "MCI",
      "MCI_BXP101_0.4ugml" = "MCI+BXP-101\n(0.4 \u03bcg/ml)"),
    Group = factor(Group, levels = c(
      "Control", "MCI", "MCI+BXP-101\n(0.4 \u03bcg/ml)"))
  )

wif_sum <- wif_long %>%
  group_by(Marker, Group) %>%
  summarise(Mean = mean(Intensity), SD = sd(Intensity), .groups = "drop")

make_wif_bar <- function(marker, ylab, sigs) {
  d_sum <- wif_sum %>% filter(Marker == marker)
  d_raw <- wif_long %>% filter(Marker == marker)
  p <- ggplot(d_sum, aes(x = Group, y = Mean, fill = Group)) +
    geom_bar(stat = "identity", width = 0.6, alpha = 0.85,
             color = "gray30", linewidth = 0.4) +
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                  width = 0.2, linewidth = 0.7) +
    geom_jitter(data = d_raw, aes(x = Group, y = Intensity, fill = Group),
                width = 0.1, size = 2, alpha = 0.6,
                shape = 21, color = "gray30") +
    scale_fill_manual(values = COL3) +
    labs(x = NULL, y = ylab) +
    BASE
  for (s in sigs) {
    p <- p +
      annotate("segment",
               x = s$x1, xend = s$x2, y = s$y, yend = s$y,
               color = "black", linewidth = 0.6) +
      annotate("text",
               x = (s$x1 + s$x2) / 2, y = s$y + s$off,
               label = s$lab, size = 7,
               color = ifelse(grepl("#", s$lab), "blue", "red"))
  }
  p
}

pC <- make_wif_bar(
  "BLBP_TRITC",
  "Mean fluorescence intensity (a.u.)\nBLBP (Reactive Gliosis)",
  list(
    list(x1 = 1, x2 = 2, y = 27,  off = 0.8, lab = "**"),
    list(x1 = 2, x2 = 3, y = 30,  off = 0.8, lab = "###")
  )
)

pD <- make_wif_bar(
  "Nestin_GFP",
  "Mean fluorescence intensity (a.u.)\nNestin (Neural Progenitor)",
  list(
    list(x1 = 1, x2 = 2, y = 13.5, off = 0.3, lab = "*"),
    list(x1 = 1, x2 = 3, y = 14.5, off = 0.3, lab = "***")
  )
)

ggsave(file.path(OUT_DIR, "Fig2C_BLBP_quant.png"), pC,
       width = 4, height = 6, dpi = DPI, bg = WHITE)
ggsave(file.path(OUT_DIR, "Fig2D_Nestin_quant.png"), pD,
       width = 4, height = 6, dpi = DPI, bg = WHITE)
message("  Saved Fig2C_BLBP_quant.png + Fig2D_Nestin_quant.png")

# ==============================================================================
# PANEL E — Representative WIF images (3x3 grid)
# ==============================================================================
message(">>> Panel E ...")

load_wif <- function(fname, scalebar = TRUE) {
  path <- file.path(IMG_DIR, fname)
  if (!file.exists(path)) {
    img <- image_blank(800, 620, color = "black")
    img <- image_annotate(img, paste0("Missing:\n", fname),
                          color = "white", size = 20, gravity = "Center")
    return(img)
  }
  img <- image_read(path)
  img <- image_trim(img)
  img <- image_scale(img, "800x620!")
  if (scalebar) {
    bar_w <- 160L
    bar_h <- 5L
    bar   <- image_blank(bar_w, bar_h, color = "white")
    img   <- image_composite(img, bar,
                              operator = "Over",
                              offset   = paste0("+", 800L - bar_w - 20L,
                                                "+", 620L - 30L))
    img   <- image_annotate(img, "1000 \u03bcm",
                             color = "white", size = 24,
                             gravity = "SouthEast", location = "+10+8")
  }
  img
}

add_row_label <- function(img, label) {
  image_annotate(img, label,
                 color = "white", size = 48, font = "bold",
                 gravity = "NorthWest", location = "+10+10")
}

imgs <- list(
  add_row_label(load_wif("Vehicle_MERGE.png"),   "Control"),
  load_wif("Vehicle_RFP_BLBP.png"),
  load_wif("Vehicle_GFP_Nestin.png"),
  add_row_label(load_wif("MCI_MERGE.png"),        "MCI"),
  load_wif("MCI_RFP_BLBP.png"),
  load_wif("MCI_GFP_Nestin.png"),
  add_row_label(load_wif("BXP_101_MERGE.png"),    "MCI+BXP-101"),
  load_wif("BXP_101_RFP_BLBP.png"),
  load_wif("BXP_101_GFP_Nestin.png")
)

row_e1 <- image_append(c(imgs[[1]], imgs[[2]], imgs[[3]]))
row_e2 <- image_append(c(imgs[[4]], imgs[[5]], imgs[[6]]))
row_e3 <- image_append(c(imgs[[7]], imgs[[8]], imgs[[9]]))
grid_e <- image_append(c(row_e1, row_e2, row_e3), stack = TRUE)

hdr <- image_blank(2400, 130, color = "white")
hdr <- image_annotate(hdr, "MERGE",
                       color = "black", size = 65, font = "bold",
                       location = "+260+15")
hdr <- image_annotate(hdr, "BLBP (TRITC)",
                       color = "black", size = 65, font = "bold",
                       location = "+1030+15")
hdr <- image_annotate(hdr, "Nestin (GFP)",
                       color = "black", size = 65, font = "bold",
                       location = "+1790+15")

pE_img <- image_append(c(hdr, grid_e), stack = TRUE)
image_write(pE_img, file.path(OUT_DIR, "Fig2E_WIF_images.png"))
message("  Saved Fig2E_WIF_images.png")

# ==============================================================================
# PANEL F — ptgdsb.1/2 longitudinal qPCR
# ==============================================================================
message(">>> Panel F ...")

q14 <- read.csv(file.path(DATA_DIR, "qPCR_14dpf.csv"),
                stringsAsFactors = FALSE) %>%
  filter(Replicate == "Average",
         Gene %in% c("ptgdsb.1", "ptgdsb.2")) %>%
  pivot_longer(c(Control, MCI, MCI_Donepezil_10uM,
                 MCI_BXP101_0.3ugml, MCI_BXP101_0.4ugml, MCI_BXP101_0.6ugml),
               names_to = "Group", values_to = "Expr") %>%
  mutate(dpf = "14dpf")

q21 <- read.csv(file.path(DATA_DIR, "qPCR_21dpf.csv"),
                stringsAsFactors = FALSE) %>%
  filter(Replicate == "Average",
         Gene %in% c("ptgdsb.1", "ptgdsb.2")) %>%
  pivot_longer(c(Control, MCI, MCI_BXP101_0.4ugml),
               names_to = "Group", values_to = "Expr") %>%
  mutate(dpf = "21dpf")

qPCR_long <- bind_rows(q14, q21) %>%
  filter(Group %in% c("Control", "MCI", "MCI_BXP101_0.4ugml")) %>%
  mutate(
    Group = dplyr::recode(Group,
      "Control"            = "Control",
      "MCI"                = "MCI",
      "MCI_BXP101_0.4ugml" = "MCI+BXP-101\n(0.4 \u03bcg/ml)"),
    Group = factor(Group, levels = c(
      "Control", "MCI", "MCI+BXP-101\n(0.4 \u03bcg/ml)")),
    dpf  = factor(dpf,  levels = c("14dpf", "21dpf")),
    Gene = factor(Gene, levels = c("ptgdsb.1", "ptgdsb.2"))
  )

pF <- ggplot(qPCR_long,
             aes(x = dpf, y = Expr, color = Group, group = Group)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray60") +
  facet_wrap(~Gene, ncol = 2) +
  scale_color_manual(values = COL3) +
  scale_x_discrete(labels = c("14dpf" = "14 dpf", "21dpf" = "21 dpf")) +
  labs(x = "Time point (dpf)",
       y = "Relative expression\n(vs Control)",
       color = NULL) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.text     = element_text(size = 14),
    strip.text      = element_text(face = "bold.italic", size = 16),
    plot.margin     = margin(8, 8, 8, 8)
  )

ggsave(file.path(OUT_DIR, "Fig2F_qPCR.png"), pF,
       width = 7, height = 5.5, dpi = DPI, bg = WHITE)
message("  Saved Fig2F_qPCR.png")

# ==============================================================================
# PANEL G — ELISA: TNF-alpha + IL-6
# ==============================================================================
message(">>> Panel G ...")

elisa <- read.csv(file.path(DATA_DIR, "ELISA_14dpf.csv"),
                  stringsAsFactors = FALSE) %>%
  pivot_longer(c(Control, MCI, MCI_BXP101_0.3ugml,
                 MCI_BXP101_0.4ugml, MCI_BXP101_0.6ugml),
               names_to = "Group", values_to = "Conc") %>%
  mutate(
    Group = dplyr::recode(Group,
      "Control"            = "Control",
      "MCI"                = "MCI",
      "MCI_BXP101_0.3ugml" = "BXP-101\n0.3",
      "MCI_BXP101_0.4ugml" = "BXP-101\n0.4",
      "MCI_BXP101_0.6ugml" = "BXP-101\n0.6"),
    Group = factor(Group, levels = c(
      "Control", "MCI", "BXP-101\n0.3", "BXP-101\n0.4", "BXP-101\n0.6")),
    Marker = factor(Marker,
      levels = c("TNF_alpha", "IL_6"),
      labels = c("TNF-\u03b1 (pg/ml)", "IL-6 (pg/ml)"))
  )

# Note: ELISA values are single measurements per group (n=1);
# error bars are not applicable.
col5 <- setNames(
  unname(COL[c("Control","MCI","BXP_0.3","BXP_0.4","BXP_0.6")]),
  c("Control","MCI","BXP-101\n0.3","BXP-101\n0.4","BXP-101\n0.6")
)

pG <- ggplot(elisa, aes(x = Group, y = Conc, fill = Group)) +
  geom_bar(stat = "identity", width = 0.65, alpha = 0.85,
           color = "gray30", linewidth = 0.4) +
  facet_wrap(~Marker, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = col5) +
  labs(x = NULL, y = "Concentration (pg/ml)",
       caption = "BXP-101 concentrations in \u03bcg/ml") +
  BASE +
  theme(
    strip.text   = element_text(face = "bold", size = 16),
    axis.text.x  = element_text(size = 13, angle = 45, hjust = 1),
    plot.caption = element_text(size = 11, hjust = 0.5,
                                color = "gray40", margin = margin(t = 4)),
    plot.margin  = margin(8, 16, 8, 8)
  )

ggsave(file.path(OUT_DIR, "Fig2G_ELISA.png"), pG,
       width = 9, height = 6.5, dpi = DPI, bg = WHITE)
message("  Saved Fig2G_ELISA.png")

# ==============================================================================
# PANEL H — Molecular recovery schematic
# ==============================================================================
message(">>> Panel H ...")

nodes <- data.frame(
  x     = 1:5,
  y     = rep(1, 5),
  label = c(
    "NF-\u03baB\nSuppression",
    "ptgdsb\nCapacity\nPreserved",
    "TNF-\u03b1/IL-6\n\u2193",
    "NGFR\nRe-expression",
    "Cognitive\nRecovery"
  ),
  col = c("#E41A1C","#4DAF4A","#FF7F00","#377EB8","#984EA3"),
  stringsAsFactors = FALSE
)

pH <- ggplot(nodes, aes(x = x, y = y)) +
  annotate("segment",
           x    = c(1.40, 2.40, 3.40, 4.40),
           xend = c(1.60, 2.60, 3.60, 4.60),
           y = 1, yend = 1,
           arrow = arrow(length = unit(0.50, "cm")),
           color = "gray30", linewidth = 1.8) +
  geom_point(aes(color = col), size = 48, alpha = 0.14) +
  geom_point(aes(color = col), size = 48, shape = 1, stroke = 2.0) +
  geom_text(aes(label = label, color = col),
            size = 8.0, fontface = "bold", lineheight = 0.95) +
  scale_color_identity() +
  scale_x_continuous(limits = c(0.3, 5.7)) +
  scale_y_continuous(limits = c(0.62, 1.38)) +
  theme_void() +
  theme(plot.margin = margin(t = 2, r = 10, b = 2, l = 10))

ggsave(file.path(OUT_DIR, "Fig2H_schematic.png"), pH,
       width = 15, height = 3.8, dpi = DPI, bg = WHITE)
message("  Saved Fig2H_schematic.png")

message("\n", strrep("=", 60))
message("  All panels saved to: ", OUT_DIR)
message("  Next step: run Fig2_combine.R to assemble with labels.")
message(strrep("=", 60))
