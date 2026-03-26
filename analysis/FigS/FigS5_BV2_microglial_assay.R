# ==============================================================================
# FigS5_BV2_microglial_assay.R
#
# Purpose : Generate Figure S5 — BV-2 microglial assay results
#           Cell viability, nitric oxide production, and NF-kB activation
#           in response to BXP-101 treatment.
#
# Panels:
#   A-C = Row 1 (3 panels)
#   D-F = Row 2 (3 panels)
#   G   = Row 3 (full width)
#   H   = Row 4 (full width)
#
# Outputs (-> output/FigS5/):
#   FigS5_final.png   (300 dpi)
#   FigS5_final.tiff  (300 dpi, LZW — journal submission)
#
# Paper : Kim Y et al. (2025). https://github.com/YoungOukKim/MCI-to-AD
# Usage : Rscript analysis/FigS/FigS5_BV2_microglial_assay.R
# Requirements : R >= 4.3.2 | magick, cowplot, ggplot2, grid
# ==============================================================================

suppressPackageStartupMessages({
  library(magick)
  library(cowplot)
  library(ggplot2)
  library(grid)
})

# ------------------------------------------------------------------------------
# Paths (relative to repository root; override with environment variables)
# ------------------------------------------------------------------------------
IMG_DIR <- Sys.getenv("FIGS5_IMG", unset = "data/FigS5")
OUT_DIR <- Sys.getenv("FIGS5_OUT", unset = "output/FigS5")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Load and trim images
# ------------------------------------------------------------------------------
message(">>> [1/3] Loading images...")

load_img <- function(fname) {
  path <- file.path(IMG_DIR, fname)
  if (!file.exists(path)) {
    warning("Image not found: ", path)
    return(ggdraw() + draw_label(paste("Missing:", fname), color = "red", size = 9))
  }
  image_trim(image_read(path))
}

img_a <- load_img("FigS5A.png"); img_b <- load_img("FigS5B.png")
img_c <- load_img("FigS5C.png"); img_d <- load_img("FigS5D.png")
img_e <- load_img("FigS5E.png"); img_f <- load_img("FigS5F.png")
img_g <- load_img("FigS5G.png"); img_h <- load_img("FigS5H.png")

# ------------------------------------------------------------------------------
# 2. Assemble layout
# ------------------------------------------------------------------------------
message(">>> [2/3] Assembling layout...")

row1 <- plot_grid(
  ggdraw() + draw_image(img_a, scale = 0.95),
  ggdraw() + draw_image(img_b, scale = 0.95),
  ggdraw() + draw_image(img_c, scale = 0.95),
  ncol = 3
)
row2 <- plot_grid(
  ggdraw() + draw_image(img_d, scale = 0.95),
  ggdraw() + draw_image(img_e, scale = 0.95),
  ggdraw() + draw_image(img_f, scale = 0.95),
  ncol = 3
)
row3 <- ggdraw() + draw_image(img_g, scale = 0.98)
row4 <- ggdraw() + draw_image(img_h, scale = 0.98)

final_plot <- plot_grid(row1, row2, row3, row4,
                        ncol = 1, rel_heights = c(1, 1, 1.5, 1.5))

# Panel labels
label_size <- 24
y_row1 <- 0.99; y_row2 <- 0.795; y_row3 <- 0.595; y_row4 <- 0.295

final_labeled <- ggdraw(final_plot) +
  draw_label("A", x = 0.01, y = y_row1, fontface = "bold", size = label_size, hjust = 0, vjust = 1) +
  draw_label("B", x = 0.34, y = y_row1, fontface = "bold", size = label_size, hjust = 0, vjust = 1) +
  draw_label("C", x = 0.67, y = y_row1, fontface = "bold", size = label_size, hjust = 0, vjust = 1) +
  draw_label("D", x = 0.01, y = y_row2, fontface = "bold", size = label_size, hjust = 0, vjust = 1) +
  draw_label("E", x = 0.34, y = y_row2, fontface = "bold", size = label_size, hjust = 0, vjust = 1) +
  draw_label("F", x = 0.67, y = y_row2, fontface = "bold", size = label_size, hjust = 0, vjust = 1) +
  draw_label("G", x = 0.01, y = y_row3, fontface = "bold", size = label_size, hjust = 0, vjust = 1) +
  draw_label("H", x = 0.01, y = y_row4, fontface = "bold", size = label_size, hjust = 0, vjust = 1) +
  # Dashed annotation boxes
  draw_grob(rectGrob(x = 0.505, y = 0.89, width = 0.32, height = 0.17,
                     gp = gpar(col = "#BFBFBF", lwd = 1.2, lty = 2, fill = NA))) +
  draw_grob(rectGrob(x = 0.505, y = 0.69, width = 0.32, height = 0.18,
                     gp = gpar(col = "#BFBFBF", lwd = 1.2, lty = 2, fill = NA)))

# ------------------------------------------------------------------------------
# 3. Save (300 dpi, LZW TIFF) with bottom whitespace crop
# ------------------------------------------------------------------------------
message(">>> [3/3] Saving...")

W   <- 12
H   <- 16
DPI <- 300

tmp_tiff <- file.path(OUT_DIR, "temp_figS5.tiff")
out_tiff <- file.path(OUT_DIR, "FigS5_final.tiff")
out_png  <- file.path(OUT_DIR, "FigS5_final.png")

ggsave(tmp_tiff, final_labeled,
       width = W, height = H, dpi = DPI, bg = "white",
       device = "tiff", compression = "lzw")

# Crop bottom whitespace (scan from bottom for last non-white row)
final_img <- image_read(tmp_tiff)
info      <- image_info(final_img)
px        <- as.integer(final_img[[1]])
arr       <- array(px, dim = c(4, info$width, info$height))
last_row  <- info$height
for (row in seq(info$height, 1, by = -1)) {
  if (any(arr[1,,row] < 250) || any(arr[2,,row] < 250) || any(arr[3,,row] < 250)) {
    last_row <- row; break
  }
}
crop_height <- min(last_row - 60, info$height)
final_img   <- image_crop(image_read(tmp_tiff),
                           geometry = paste0(info$width, "x", crop_height, "+0+0"))
image_write(final_img, out_tiff, format = "tiff", compression = "lzw")
image_write(final_img, out_png)
if (file.exists(tmp_tiff)) file.remove(tmp_tiff)

message("\n", strrep("=", 55))
message("  Done! PNG : ", out_png)
message("  Done! TIFF: ", out_tiff)
message("  Size: ", W, " x ", round(crop_height/DPI, 2), " in @ ", DPI, " dpi")
message(strrep("=", 55))
