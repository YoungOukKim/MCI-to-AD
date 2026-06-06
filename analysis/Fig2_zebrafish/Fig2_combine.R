# ==============================================================================
# Fig2_combine.R
#
# Assembles individual Fig. 2 panels (A-H) into a single figure with labels.
#
# Layout (15 inches wide, 300 dpi):
#   Row 1  A              full width       h = 4.0
#   Row 2  B | C          1 : 1            h = 6.0
#   Row 3  D | E          1 : 1.5          h = 6.0
#   Row 4  F | G          1 : 1.2          h = 6.0
#   Row 5  H (full width)                  h = 4.0
#   Total                                  h = 26.0
#
# Run after Fig2A_H_panels.R has been executed.
#
# Paper : Kim Y et al. (2025). https://github.com/YoungOukKim/MCI-to-AD
# Usage : Rscript analysis/Fig2_zebrafish/Fig2_combine.R
# Requirements : R >= 4.3.2 | cowplot, magick
# ==============================================================================

suppressPackageStartupMessages({
  library(cowplot)
  library(magick)
})

# ------------------------------------------------------------------------------
# Paths (relative to repository root; override with environment variables)
# ------------------------------------------------------------------------------
IMG_DIR <- Sys.getenv("FIG2_OUT", unset = "output/Fig2")
OUT_DIR <- IMG_DIR
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

load_panel <- function(fname, scale = 0.97) {
  path <- file.path(IMG_DIR, fname)
  if (!file.exists(path))
    stop("Panel not found: ", path, "\n  Run Fig2A_H_panels.R first.")
  ggdraw() + draw_image(image_read(path), scale = scale)
}

# ------------------------------------------------------------------------------
# Load panels
# ------------------------------------------------------------------------------
message(">>> Loading panels...")
pA <- load_panel("Fig2A_timeline.png")
pB <- load_panel("Fig2B_avoidance.png")
pC <- load_panel("Fig2C_BLBP_quant.png")
pD <- load_panel("Fig2D_Nestin_quant.png")
pE <- load_panel("Fig2E_WIF_images.png")
pF <- load_panel("Fig2F_qPCR.png")
pG <- load_panel("Fig2G_ELISA.png")

# Row heights (inches)
h       <- c(4.0, 6.0, 6.0, 6.0, 4.0)
H_TOTAL <- sum(h)   # 26.0
h_rel   <- h / H_TOTAL

# ------------------------------------------------------------------------------
# Assemble rows
# ------------------------------------------------------------------------------
message(">>> Assembling rows...")
pad <- theme(plot.margin = margin(t = 10, r = 35, b = 10, l = 35))

row1 <- pA + pad
row2 <- plot_grid(pB, pC, ncol = 2, rel_widths = c(1,   1  )) + pad
row3 <- plot_grid(pD, pE, ncol = 2, rel_widths = c(1,   1.5)) + pad
row4 <- plot_grid(pF, pG, ncol = 2, rel_widths = c(1,   1.2)) + pad

# Row 5: Panel H schematic placed at top of its row (y=0.65 pushes image up)
row5 <- ggdraw() +
  draw_image(magick::image_read(file.path(IMG_DIR, "Fig2H_schematic.png")),
             x = 0, y = 0.65, width = 1, height = 1.0)

fig2_base <- plot_grid(
  row1, row2, row3, row4, row5,
  ncol = 1, rel_heights = h_rel
)

# ------------------------------------------------------------------------------
# Panel label positions
# ------------------------------------------------------------------------------
row_top <- c(
  1.000,
  1.000 - h_rel[1],
  1.000 - sum(h_rel[1:2]),
  1.000 - sum(h_rel[1:3]),
  1.000 - sum(h_rel[1:4])
)
Y_OFF    <- 0.010
y_lbl    <- row_top - Y_OFF
y_lbl[5] <- row_top[5] + 0.025   # H label aligned with schematic top

X_L  <- 0.02    # left column labels
X_R2 <- 0.50    # Row2: B|C  (1:1)
X_R3 <- 0.40    # Row3: D|E  (1:1.5)
X_R4 <- 0.455   # Row4: F|G  (1:1.2)

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
  draw_label("H", x = X_L,  y = y_lbl[5], fontface = "bold", size = LSIZE)

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------
W   <- 15
DPI <- 300
tmp_png  <- file.path(OUT_DIR, "Fig2_temp.png")
out_png  <- file.path(OUT_DIR, "Fig2_final.png")
out_tiff <- file.path(OUT_DIR, "Fig2_final.tiff")

message(">>> Saving (", W, " x ", H_TOTAL, " in @ ", DPI, " dpi)...")
ggsave(tmp_png, fig2_labeled,
       width = W, height = H_TOTAL, dpi = DPI, bg = "white")

# Crop bottom whitespace below Panel H schematic
# Row5 h=4.0 in, y=0.65 -> schematic occupies top ~1.4 in -> crop ~1100 px
img  <- image_read(tmp_png)
info <- image_info(img)
img  <- image_crop(img,
         geometry = paste0(info$width, "x", info$height - 1100, "+0+0"))
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
