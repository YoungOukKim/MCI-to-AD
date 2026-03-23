# ==============================================================================
# Fig4_combine.R
#
# Purpose : Assemble individual Fig. 4 panels into final figure with labels.
#           Panel label size/position follows Fig. 6 convention (size=30, 15in).
#
# Layout (15 inches wide, 600 dpi):
#
#   Row 1  A  |  B               1 : 1          h = 6.0
#   Row 2  C  (full width)                      h = 5.5
#   Row 3  D  |  E               1.2 : 1        h = 5.5
#   Row 4  F  (full width)                      h = 5.5
#   Total                                       h = 22.5
#
# Note: Panel D (immunohistochemistry images) must be supplied as
#       Fig4D_IHC.png in OUT_DIR. A placeholder is shown if missing.
#
# Paper : Kim Y et al. (2025). https://github.com/YoungOukKim/MCI-to-AD
# Usage : Rscript analysis/Fig4_murine/Fig4_combine.R
# Requirements: R >= 4.3.2 | cowplot, magick
# ==============================================================================

suppressPackageStartupMessages({
  library(cowplot)
  library(magick)
})

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
OUT_DIR <- Sys.getenv("FIG4_OUT", unset = "output/Fig4")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Helper: load panel image as ggdraw object
# ------------------------------------------------------------------------------
load_panel <- function(fname, scale = 0.97) {
  path <- file.path(OUT_DIR, fname)
  if (!file.exists(path)) {
    # placeholder
    img <- magick::image_blank(800, 600, color = "white")
    img <- magick::image_annotate(img,
      paste0("Missing:\n", fname, "\n(run Fig4_panels.R first)"),
      color = "gray50", size = 28, gravity = "Center")
    magick::image_write(img, file.path(OUT_DIR, fname))
  }
  ggdraw() + draw_image(magick::image_read(path), scale = scale)
}

# Panel D placeholder if IHC image not yet available
make_placeholder <- function(label) {
  img <- magick::image_blank(900, 700, color = "#F5F5F5")
  img <- magick::image_annotate(img,
    paste0("Panel ", label, "\n(IHC images — supplied separately)"),
    color = "gray60", size = 36, gravity = "Center", font = "bold")
  magick::image_write(img, file.path(OUT_DIR, paste0("Fig4", label, "_placeholder.png")))
  ggdraw() + draw_image(img, scale = 0.97)
}

message(">>> Loading panels...")
pA <- load_panel("Fig4A_ymaze.png")
pB <- load_panel("Fig4B_PAT.png")
pC <- load_panel("Fig4C_inflammatory.png")
pD <- if (file.exists(file.path(OUT_DIR, "Fig4D_IHC.png")))
        load_panel("Fig4D_IHC.png") else make_placeholder("D")
pE <- load_panel("Fig4E_trophic.png")
pF <- load_panel("Fig4F_oxidative.png")

# ------------------------------------------------------------------------------
# Row heights
# ------------------------------------------------------------------------------
h     <- c(6.0, 5.5, 5.5, 5.5)   # rows 1–4
H_TOT <- sum(h)                    # 22.5 inches
h_rel <- h / H_TOT

# ------------------------------------------------------------------------------
# Assemble rows
# ------------------------------------------------------------------------------
message(">>> Assembling rows...")
pad <- theme(plot.margin = margin(t = 8, r = 30, b = 8, l = 30))

row1 <- plot_grid(pA, pB, ncol = 2, rel_widths = c(1, 1)) + pad
row2 <- pC + pad
row3 <- plot_grid(pD, pE, ncol = 2, rel_widths = c(1.2, 1)) + pad
row4 <- pF + pad

fig4_base <- plot_grid(row1, row2, row3, row4,
                       ncol = 1, rel_heights = h_rel)

# ------------------------------------------------------------------------------
# Panel label positions
# ------------------------------------------------------------------------------
row_top <- c(
  1.000,
  1.000 - h_rel[1],
  1.000 - sum(h_rel[1:2]),
  1.000 - sum(h_rel[1:3])
)
Y_OFF  <- 0.010
y_lbl  <- row_top - Y_OFF

X_L    <- 0.02   # left column
X_R    <- 0.51   # right column (row1: 1:1 split)
# Row3: D=1.2, E=1 -> E starts at 1.2/2.2 = 0.545
x_E    <- 1.2 / 2.2

LSIZE  <- 30     # matches Fig. 6 convention

message(">>> Adding panel labels...")
fig4_labeled <- ggdraw(fig4_base) +
  draw_label("A", x = X_L,  y = y_lbl[1], fontface = "bold", size = LSIZE) +
  draw_label("B", x = X_R,  y = y_lbl[1], fontface = "bold", size = LSIZE) +
  draw_label("C", x = X_L,  y = y_lbl[2], fontface = "bold", size = LSIZE) +
  draw_label("D", x = X_L,  y = y_lbl[3], fontface = "bold", size = LSIZE) +
  draw_label("E", x = x_E,  y = y_lbl[3], fontface = "bold", size = LSIZE) +
  draw_label("F", x = X_L,  y = y_lbl[4], fontface = "bold", size = LSIZE)

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------
W   <- 15
DPI <- 600
tmp_png  <- file.path(OUT_DIR, "Fig4_temp.png")
out_png  <- file.path(OUT_DIR, "Fig4_final.png")
out_tiff <- file.path(OUT_DIR, "Fig4_final.tiff")

message(">>> Saving (", W, " x ", H_TOT, " in @ ", DPI, " dpi)...")
ggsave(tmp_png, fig4_labeled,
       width = W, height = H_TOT, dpi = DPI, bg = "white")

img  <- magick::image_read(tmp_png)
info <- magick::image_info(img)
img  <- magick::image_crop(img,
         geometry = paste0(info$width, "x", info$height - 240, "+0+0"))
magick::image_write(img, out_png)
if (file.exists(tmp_png)) file.remove(tmp_png)

ggsave(out_tiff, fig4_labeled,
       width = W, height = H_TOT, dpi = DPI, bg = "white",
       device = "tiff", compression = "lzw")

message("\n", strrep("=", 60))
message("  Done! PNG : ", out_png)
message("  Done! TIFF: ", out_tiff)
message("  Size : ", W, " x ", H_TOT, " in @ ", DPI, " dpi")
message("  Note : Place IHC images as Fig4D_IHC.png to replace placeholder.")
message(strrep("=", 60))
