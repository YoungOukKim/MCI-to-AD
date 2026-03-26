# ==============================================================================
# FigS7_compound_characterization.R
#
# Purpose : Generate Figure S7 — HPLC chromatograms and chemical structures
#           of BXP-101 active compounds (Atractylodin, Wedelolactone, Honokiol).
#
# Panels:
#   A = Atractylodin: HPLC chromatogram + chemical structure
#   B = Wedelolactone: HPLC chromatogram + chemical structure
#   C = Honokiol: HPLC chromatogram + chemical structure
#
# Outputs (-> output/FigS7/):
#   FigS7_final.png   (300 dpi)
#   FigS7_final.tiff  (300 dpi, LZW — journal submission)
#
# Paper : Kim Y et al. (2025). https://github.com/YoungOukKim/MCI-to-AD
# Usage : Rscript analysis/FigS/FigS7_compound_characterization.R
# Requirements : R >= 4.3.2 | ggplot2, cowplot, magick
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(magick)
})

# ------------------------------------------------------------------------------
# Paths (relative to repository root; override with environment variables)
# ------------------------------------------------------------------------------
IMG_DIR <- Sys.getenv("FIGS7_IMG", unset = "data/FigS7")
OUT_DIR <- Sys.getenv("FIGS7_OUT", unset = "output/FigS7")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Load and standardize images
# ------------------------------------------------------------------------------
message(">>> [1/3] Loading images...")

load_img <- function(fname, rotate_90 = FALSE) {
  path <- file.path(IMG_DIR, fname)
  if (!file.exists(path)) {
    warning("Image not found: ", path)
    return(image_blank(2000, 1500, color = "white") %>%
           image_annotate(paste("Missing:", fname), color = "red",
                          size = 40, gravity = "Center"))
  }
  img <- image_read(path) %>% image_trim()
  if (rotate_90) img <- img %>% image_rotate(90) %>% image_trim()
  image_scale(img, "2000x1500!")
}

img_a_p <- load_img("FigS7A_atractylodin_peak.png")
img_a_c <- load_img("FigS7A_atractylodin_chem.png")
img_b_p <- load_img("FigS7B_wedelolactone_peak.png")
img_b_c <- load_img("FigS7B_wedelolactone_chem.png")
img_c_p <- load_img("FigS7C_honokiol_peak.png")
img_c_c <- load_img("FigS7C_honokiol_chem.png", rotate_90 = TRUE)

# ------------------------------------------------------------------------------
# 2. Assemble layout
# ------------------------------------------------------------------------------
message(">>> [2/3] Assembling layout...")

draw_peak <- function(img) {
  ggdraw() + draw_image(img, x = 0, y = 0, width = 1, height = 1, scale = 1.0)
}

draw_chem <- function(img) {
  ggdraw() + draw_image(img, x = -0.05, y = 0.12, width = 1, height = 1, scale = 0.75)
}

make_row <- function(p_img, c_img) {
  plot_grid(draw_peak(p_img), draw_chem(c_img),
            ncol = 2, rel_widths = c(1, 0.9)) +
    theme(plot.margin = margin(t = 20, r = 20, b = 20, l = 50))
}

row1 <- make_row(img_a_p, img_a_c)
row2 <- make_row(img_b_p, img_b_c)
row3 <- make_row(img_c_p, img_c_c)

final_plot <- plot_grid(row1, row2, row3, ncol = 1)

# Panel labels
final_labeled <- ggdraw(final_plot) +
  draw_label("A", x = 0.02, y = 0.985, fontface = "bold", size = 28) +
  draw_label("B", x = 0.02, y = 0.655, fontface = "bold", size = 28) +
  draw_label("C", x = 0.02, y = 0.325, fontface = "bold", size = 28)

# ------------------------------------------------------------------------------
# 3. Save (300 dpi, LZW TIFF) with right/bottom whitespace crop
# ------------------------------------------------------------------------------
message(">>> [3/3] Saving...")

W   <- 14
H   <- 18
DPI <- 300

tmp_tiff <- file.path(OUT_DIR, "temp_figS7.tiff")
out_tiff <- file.path(OUT_DIR, "FigS7_final.tiff")
out_png  <- file.path(OUT_DIR, "FigS7_final.png")

ggsave(tmp_tiff, final_labeled,
       width = W, height = H, dpi = DPI, bg = "white",
       device = "tiff", compression = "lzw")

# Crop: remove right excess and bottom whitespace
final_img <- image_read(tmp_tiff)
info      <- image_info(final_img)

# Scan bottom for last non-white row
px      <- as.integer(final_img[[1]])
arr     <- array(px, dim = c(4, info$width, info$height))
last_row <- info$height
for (row in seq(info$height, 1, by = -1)) {
  if (any(arr[1,,row] < 250) || any(arr[2,,row] < 250) || any(arr[3,,row] < 250)) {
    last_row <- row; break
  }
}
crop_h <- min(last_row + 30, info$height)
crop_w <- info$width - 150   # trim right margin (adjust if needed)

final_img <- image_crop(image_read(tmp_tiff),
                         geometry = paste0(crop_w, "x", crop_h, "+0+0"))
image_write(final_img, out_tiff, format = "tiff", compression = "lzw")
image_write(final_img, out_png)
if (file.exists(tmp_tiff)) file.remove(tmp_tiff)

message("\n", strrep("=", 55))
message("  Done! PNG : ", out_png)
message("  Done! TIFF: ", out_tiff)
message("  Size: ", round(crop_w/DPI, 2), " x ", round(crop_h/DPI, 2),
        " in @ ", DPI, " dpi")
message(strrep("=", 55))
