# ==============================================================================
# Fig2_combine.R  (v8)
#
# Fix: E/G label overlap — labels placed at EXACT row boundaries
#      E row label moved INSIDE image (white text on MERGE column)
#      G label x position corrected for 1:1.2 column split
#
# Layout (15 inches wide, 600 dpi):
#   Row 1  A              full width       h = 4.0
#   Row 2  B | C          1 : 1            h = 6.0
#   Row 3  D | E          1 : 1.5          h = 6.0
#   Row 4  F | G          1 : 1.2          h = 6.0
#   Row 5  H | I          1 : 2            h = 5.0
#   Total                                  h = 27.0
# ==============================================================================

suppressPackageStartupMessages({
  library(cowplot)
  library(magick)
})

IMG_DIR <- Sys.getenv("FIG2_OUT", unset = "D:/work/paper-results/Fig2")
OUT_DIR <- IMG_DIR
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

load_panel <- function(fname, scale = 0.97) {
  path <- file.path(IMG_DIR, fname)
  if (!file.exists(path))
    stop("Panel not found: ", path, "\n  Run Fig2A_I_panels.R first.")
  ggdraw() + draw_image(image_read(path), scale = scale)
}

message(">>> Loading panels...")
pA <- load_panel("Fig2A_timeline.png")
pB <- load_panel("Fig2B_avoidance.png")
pC <- load_panel("Fig2C_BLBP_quant.png")
pD <- load_panel("Fig2D_Nestin_quant.png")
pE <- load_panel("Fig2E_WIF_images.png")
pF <- load_panel("Fig2F_qPCR.png")
pG <- load_panel("Fig2G_ELISA.png")
pH <- load_panel("Fig2H_correlation.png")
pI <- load_panel("Fig2I_schematic.png")

h <- c(4.0, 6.0, 6.0, 6.0, 5.0)
H_TOTAL <- sum(h)   # 27.0
h_rel   <- h / H_TOTAL

message(">>> Assembling rows...")
pad <- theme(plot.margin = margin(t = 10, r = 35, b = 10, l = 35))

row1 <- pA + pad
row2 <- plot_grid(pB, pC, ncol = 2, rel_widths = c(1,   1  )) + pad
row3 <- plot_grid(pD, pE, ncol = 2, rel_widths = c(1,   1.5)) + pad
row4 <- plot_grid(pF, pG, ncol = 2, rel_widths = c(1,   1.2)) + pad
row5 <- plot_grid(pH, pI, ncol = 2, rel_widths = c(1,   2  )) + pad

fig2_base <- plot_grid(
  row1, row2, row3, row4, row5,
  ncol = 1, rel_heights = h_rel
)

# ── Label positions ────────────────────────────────────────────────────────────
row_top <- c(
  1.000,
  1.000 - h_rel[1],
  1.000 - sum(h_rel[1:2]),
  1.000 - sum(h_rel[1:3]),
  1.000 - sum(h_rel[1:4])
)
Y_OFF <- 0.010
y_lbl <- row_top - Y_OFF

# Left column x
X_L <- 0.02

# Right column x — calculated from rel_widths per row
# Row2 (1:1):   right starts at 0.5  → X_R2 = 0.50
# Row3 (1:1.5): right starts at 1/2.5 = 0.40 of row → X_R3 = 0.40
# Row4 (1:1.2): right starts at 1/2.2 = 0.455 → X_R4 = 0.455
# Row5 (1:2):   I at 1/3 → x_I = 0.34
X_R2 <- 0.50
X_R3 <- 0.40    # D|E split: D=1, E=1.5 → E starts at 40% of figure width
X_R4 <- 0.455   # F|G split: F=1, G=1.2 → G starts at 45.5%
x_I  <- 0.34    # H|I split: H=1, I=2 → I starts at 33%

LSIZE <- 32

message(">>> Adding labels...")
fig2_labeled <- ggdraw(fig2_base) +
  draw_label("A", x = X_L,  y = y_lbl[1], fontface = "bold", size = LSIZE) +
  draw_label("B", x = X_L,  y = y_lbl[2], fontface = "bold", size = LSIZE) +
  draw_label("C", x = X_R2, y = y_lbl[2], fontface = "bold", size = LSIZE) +
  draw_label("D", x = X_L,  y = y_lbl[3], fontface = "bold", size = LSIZE) +
  draw_label("E", x = X_R3, y = y_lbl[3], fontface = "bold", size = LSIZE) +
  draw_label("F", x = X_L,  y = y_lbl[4], fontface = "bold", size = LSIZE) +
  draw_label("G", x = X_R4, y = y_lbl[4], fontface = "bold", size = LSIZE) +
  draw_label("H", x = X_L,  y = y_lbl[5], fontface = "bold", size = LSIZE) +
  draw_label("I", x = x_I,  y = y_lbl[5], fontface = "bold", size = LSIZE)

# ── Save ───────────────────────────────────────────────────────────────────────
W   <- 15
DPI <- 600
tmp_png  <- file.path(OUT_DIR, "Fig2_temp.png")
out_png  <- file.path(OUT_DIR, "Fig2_final.png")
out_tiff <- file.path(OUT_DIR, "Fig2_final.tiff")

message(">>> Saving (", W, " x ", H_TOTAL, " in @ ", DPI, " dpi)...")
ggsave(tmp_png, fig2_labeled,
       width = W, height = H_TOTAL, dpi = DPI, bg = "white")

img  <- image_read(tmp_png)
info <- image_info(img)
img  <- image_crop(img,
         geometry = paste0(info$width, "x", info$height - 240, "+0+0"))
image_write(img, out_png)
if (file.exists(tmp_png)) file.remove(tmp_png)

ggsave(out_tiff, fig2_labeled,
       width = W, height = H_TOTAL, dpi = DPI, bg = "white",
       device = "tiff", compression = "lzw")

message("\n", strrep("=", 60))
message("  Done! PNG : ", out_png)
message("  Done! TIFF: ", out_tiff)
message("  Size: ", W, " x ", H_TOTAL, " in @ ", DPI, " dpi")
message(strrep("=", 60))
