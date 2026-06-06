# ==============================================================================
# FigS3_PTGDS_conservation.R
#
# Purpose : Generate Figure S3 — Cross-species conservation of PTGDS
#           NF-kB binding motif and structural homology analysis.
#
# Panels:
#   A = Structural homology overview (external image)
#   B = Sequence alignment / conservation (external image)
#   C = NF-kB (RELA) binding affinity scores across species (bar chart)
#
# Outputs (-> output/FigS3/):
#   FigS3_final.png   (300 dpi)
#   FigS3_final.tiff  (300 dpi, LZW — journal submission)
#
# Paper : Kim Y et al. (2025). https://github.com/YoungOukKim/MCI-to-AD
# Usage : Rscript analysis/FigS/FigS3_PTGDS_conservation.R
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
IMG_DIR <- Sys.getenv("FIGS3_IMG",
              unset = "data/FigS3")
OUT_DIR <- Sys.getenv("FIGS3_OUT",
              unset = "output/FigS3")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Load external images
# ------------------------------------------------------------------------------
message(">>> [1/3] Loading images...")

img_a <- image_trim(image_read(file.path(IMG_DIR, "FigS3A.png")))
img_b <- image_trim(image_read(file.path(IMG_DIR, "FigS3B.png")))

# ------------------------------------------------------------------------------
# 2. Panel C — NF-kB (RELA) binding affinity scores
# ------------------------------------------------------------------------------
message(">>> [2/3] Building Panel C...")

motif_data <- data.frame(
  Gene    = c("Human PTGDS", "Zebrafish ptgdsb.1", "Zebrafish ptgdsb.2"),
  Species = c("Mammalian", "Danio rerio", "Danio rerio"),
  Max_Score = c(0.8870, 0.8539, 0.7980)
)
motif_data$Species <- factor(motif_data$Species,
                              levels = c("Mammalian", "Danio rerio"))
motif_data$Gene    <- factor(motif_data$Gene,
                              levels = rev(motif_data$Gene))

pC <- ggplot(motif_data, aes(x = Gene, y = Max_Score, fill = Species)) +
  geom_bar(stat = "identity", width = 0.5, alpha = 0.9) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "#E63946") +
  scale_fill_manual(values = c("Mammalian"  = "#1D3557",
                                "Danio rerio" = "#457B9D")) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1.1), expand = c(0, 0)) +
  theme_classic(base_size = 18) +
  labs(y = "NF-kB (RELA) Binding Affinity Score",
       x = "", fill = "Species") +
  theme(legend.position = "right")

# ------------------------------------------------------------------------------
# 3. Assemble layout
# ------------------------------------------------------------------------------
message(">>> [3/3] Assembling and saving...")

row_bottom <- plot_grid(
  ggdraw() + draw_image(img_b, scale = 0.95),
  pC,
  ncol = 2, rel_widths = c(0.8, 2.2)
)

label_size <- 36

final_plot <- ggdraw() +
  draw_plot(ggdraw() + draw_image(img_a, scale = 1.0),
            x = 0.02, y = 0.76, width = 0.96, height = 0.22) +
  draw_plot(row_bottom,
            x = 0, y = 0, width = 1, height = 0.72) +
  draw_label("A", x = 0.02, y = 0.97, fontface = "bold", size = label_size) +
  draw_label("B", x = 0.02, y = 0.73, fontface = "bold", size = label_size) +
  draw_label("C", x = 0.30, y = 0.73, fontface = "bold", size = label_size)

# ------------------------------------------------------------------------------
# 4. Save (300 dpi, LZW TIFF)
# ------------------------------------------------------------------------------
W   <- 18
H   <- 10
DPI <- 300

out_png  <- file.path(OUT_DIR, "FigS3_final.png")
out_tiff <- file.path(OUT_DIR, "FigS3_final.tiff")

ggsave(out_png,  final_plot,
       width = W, height = H, dpi = DPI, bg = "white")
ggsave(out_tiff, final_plot,
       width = W, height = H, dpi = DPI, bg = "white",
       device = "tiff", compression = "lzw")

message("\n", strrep("=", 55))
message("  Done! PNG : ", out_png)
message("  Done! TIFF: ", out_tiff)
message("  Size: ", W, " x ", H, " in @ ", DPI, " dpi")
message(strrep("=", 55))
