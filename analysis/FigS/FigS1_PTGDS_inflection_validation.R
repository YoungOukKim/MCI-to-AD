# ==============================================================================
# FigS1_PTGDS_inflection_validation.R
#
# Purpose : Generate Figure S1 — Statistical validation of biphasic
#           PTGDS trajectory and inflection point identification.
#
# Panels:
#   A = Raw bin-level PTGDS expression with SEM error bars
#   B = 3-bin moving average trajectory
#   C = LOESS regression (consistent with Fig. 1D style)
#   D = Segmented regression with breakpoint CI and Davies' test
#
# Outputs (-> output/FigS1/):
#   FigS1_final.png   (300 dpi)
#   FigS1_final.tiff  (300 dpi, LZW — journal submission)
#
# Paper : Kim Y et al. (2025). https://github.com/YoungOukKim/MCI-to-AD
# Usage : Rscript analysis/FigS/FigS1_PTGDS_inflection_validation.R
# Requirements : R >= 4.3.2 | ggplot2, dplyr, segmented, data.table,
#                cowplot, patchwork
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(segmented)
  library(data.table)
  library(cowplot)
  library(patchwork)
})

set.seed(42)

# ------------------------------------------------------------------------------
# Paths (relative to repository root; override with environment variables)
# ------------------------------------------------------------------------------
# NOTE: SEAAD_processed_data.RData contains cell-level data (~1.3M nuclei)
#       and cannot be redistributed due to Allen Brain Atlas data use policy.
#       Instead, use the pre-computed bin-level summary provided in this repo:
#         data/SEA-AD/FigS1_bin_means.csv
#       To regenerate from raw data, run: analysis/01_SEAAD_data_processing.R
BIN_CSV <- Sys.getenv("FIGS1_BIN_CSV",
              unset = "data/SEA-AD/FigS1_bin_means.csv")
OUT_DIR <- Sys.getenv("FIGS1_OUT",
              unset = "output/FigS1")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Data loading
# ------------------------------------------------------------------------------
# Uses pre-computed bin-level summary (9 rows x bin).
# Columns: bin, PTGDS_mean, PTGDS_sd, n_cells, PTGDS_se, PTGDS_ma
message(">>> [1/4] Loading bin-level summary...")
bin_means <- data.table::as.data.table(read.csv(BIN_CSV))

# Verify required columns
required_cols <- c("bin", "PTGDS_mean", "PTGDS_sd", "n_cells", "PTGDS_se", "PTGDS_ma")
missing_cols  <- setdiff(required_cols, names(bin_means))
if (length(missing_cols) > 0)
  stop("Missing columns in BIN_CSV: ", paste(missing_cols, collapse = ", "))

message("  Loaded ", nrow(bin_means), " bins from: ", BIN_CSV)

# ------------------------------------------------------------------------------
# 2. Segmented regression
# ------------------------------------------------------------------------------
message(">>> [2/4] Fitting segmented regression...")

lm_fit  <- lm(PTGDS_mean ~ bin, data = bin_means)
seg_fit <- segmented(lm_fit, seg.Z = ~ bin, psi = 0.5)

# Breakpoint estimate
bp <- summary(seg_fit)$psi[, "Est."]

# 95% CI via confint.segmented() — consistent with Table 2
ci_result <- tryCatch({
  ci_raw <- confint.segmented(seg_fit)
  ci_mat <- if (is.list(ci_raw) && !is.matrix(ci_raw)) ci_raw[[1]] else ci_raw
  if (is.matrix(ci_mat)) list(lo = ci_mat[1, 2], hi = ci_mat[1, 3])
  else                   list(lo = ci_mat[2],     hi = ci_mat[3])
}, error = function(e) {
  message("  confint.segmented() failed, using ±1.96*SE: ", e$message)
  bp_se <- summary(seg_fit)$psi[, "St.Err"]
  list(lo = bp - 1.96 * bp_se, hi = bp + 1.96 * bp_se)
})
ci_lo_bp <- ci_result$lo
ci_hi_bp <- ci_result$hi

# Davies' test
davies_p <- davies.test(lm_fit, ~ bin)$p.value

# Pre/post slopes
sl         <- slope(seg_fit)
slope_pre  <- round(sl$bin[1, "Est."], 2)
slope_post <- round(sl$bin[2, "Est."], 2)

# Fine-grid prediction for smooth segmented line
x_fine          <- data.frame(bin = seq(0.1, 0.9, by = 0.005))
x_fine$seg_pred <- predict(seg_fit, newdata = x_fine)

message(sprintf("  Breakpoint = %.3f (95%% CI: %.2f-%.2f), Davies p = %.4f",
                bp, round(ci_lo_bp, 2), round(ci_hi_bp, 2), davies_p))
message(sprintf("  Slopes: pre = %+.2f, post = %+.2f", slope_pre, slope_post))

# ------------------------------------------------------------------------------
# 3. Shared plot settings
# ------------------------------------------------------------------------------
y_lo    <- min(bin_means$PTGDS_mean) - 0.08
y_hi    <- max(bin_means$PTGDS_mean) + 0.15
x_scale <- scale_x_continuous(breaks = seq(0.1, 0.9, 0.1))

theme_panel <- theme_classic(base_size = 14) +
  theme(
    plot.margin = margin(t = 10, r = 15, b = 10, l = 15),
    axis.title  = element_text(size = 12),
    axis.text   = element_text(size = 11, color = "black")
  )

# ------------------------------------------------------------------------------
# 4. Panel construction
# ------------------------------------------------------------------------------
message(">>> [3/4] Building panels...")

# Panel A — Raw expression with SEM error bars
pA <- ggplot(bin_means, aes(x = bin, y = PTGDS_mean)) +
  geom_pointrange(aes(ymin = PTGDS_mean - PTGDS_se,
                      ymax = PTGDS_mean + PTGDS_se),
                  size = 0.6, color = "#0072B2", linewidth = 0.5) +
  geom_line(color = "#0072B2", linewidth = 0.8) +
  x_scale + coord_cartesian(ylim = c(y_lo, y_hi)) +
  theme_panel +
  labs(x = "Pseudo-progression bin",
       y = "Mean PTGDS expression\n(log-normalized)")

# Panel B — 3-bin moving average
pB <- ggplot(bin_means, aes(x = bin)) +
  geom_point(aes(y = PTGDS_mean), size = 2.5, alpha = 0.35, color = "gray55") +
  geom_line(aes(y = PTGDS_ma),  linewidth = 1.3, color = "#0072B2") +
  geom_point(aes(y = PTGDS_ma), size = 3,        color = "#0072B2") +
  x_scale + coord_cartesian(ylim = c(y_lo, y_hi)) +
  theme_panel +
  labs(x = "Pseudo-progression bin",
       y = "PTGDS expression")

# Panel C — LOESS regression (consistent with Fig. 1D: span=1.0, se ribbon)
pC <- ggplot(bin_means, aes(x = bin, y = PTGDS_mean)) +
  geom_point(size = 3, color = "#0072B2", alpha = 0.6) +
  geom_smooth(method = "loess", span = 1.0, se = TRUE,
              color = "#0072B2", fill = "#0072B2",
              linewidth = 2, alpha = 0.1) +
  x_scale + coord_cartesian(ylim = c(y_lo, y_hi)) +
  theme_panel +
  labs(x = "Pseudo-progression bin",
       y = "PTGDS expression")

# Panel D — Segmented regression with breakpoint CI
bp_label <- sprintf(
  "Breakpoint = %.3f\n(95%% CI: %.2f\u2013%.2f)\nDavies' p = %.4f\nSlope: %+.2f \u2192 %+.2f",
  bp, round(ci_lo_bp, 2), round(ci_hi_bp, 2), davies_p, slope_pre, slope_post
)

pD <- ggplot() +
  annotate("rect", xmin = ci_lo_bp, xmax = ci_hi_bp,
           ymin = -Inf, ymax = Inf, alpha = 0.12, fill = "#D6604D") +
  geom_vline(xintercept = bp, linetype = "dashed",
             color = "#D6604D", linewidth = 0.8) +
  geom_smooth(data = bin_means, aes(x = bin, y = PTGDS_mean),
              method = "loess", span = 1.0, se = FALSE,
              color = "gray70", linewidth = 1.0) +
  geom_line(data = x_fine, aes(x = bin, y = seg_pred),
            linewidth = 1.5, color = "#B2182B") +
  geom_point(data = bin_means, aes(x = bin, y = PTGDS_mean),
             size = 3.5, color = "#0072B2") +
  annotate("label",
           x = 0.92, y = y_lo + 0.08,
           label = bp_label,
           hjust = 1, vjust = 0, size = 3.3, fontface = "italic",
           fill = "white", label.size = 0.3, alpha = 0.95) +
  x_scale + coord_cartesian(ylim = c(y_lo, y_hi)) +
  theme_panel +
  labs(x = "Pseudo-progression bin",
       y = "PTGDS expression")

# ------------------------------------------------------------------------------
# 5. Assemble and save
# ------------------------------------------------------------------------------
message(">>> [4/4] Assembling and saving...")

fig_s1 <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag    = element_text(face = "bold", size = 28),
        plot.margin = margin(t = 5, b = 5, l = 5, r = 5))

W   <- 14
H   <- 10
DPI <- 300

out_png  <- file.path(OUT_DIR, "FigS1_final.png")
out_tiff <- file.path(OUT_DIR, "FigS1_final.tiff")

ggsave(out_png,  fig_s1,
       width = W, height = H, dpi = DPI, bg = "white")
ggsave(out_tiff, fig_s1,
       width = W, height = H, dpi = DPI, bg = "white",
       device = "tiff", compression = "lzw")

message("\n", strrep("=", 55))
message("  Done! PNG : ", out_png)
message("  Done! TIFF: ", out_tiff)
message("  Size: ", W, " x ", H, " in @ ", DPI, " dpi")
message(sprintf("  Breakpoint: %.3f (95%% CI: %.2f-%.2f), Davies p=%.4f",
                bp, round(ci_lo_bp, 2), round(ci_hi_bp, 2), davies_p))
message(strrep("=", 55))

# Save bin-level summary for repository (run once to generate data/SEA-AD/FigS1_bin_means.csv)
# Uncomment below if regenerating from raw SEAAD_processed_data.RData:
# write.csv(bin_means, file.path("data/SEA-AD", "FigS1_bin_means.csv"), row.names = FALSE)
# message("  Saved bin-level summary to data/SEA-AD/FigS1_bin_means.csv")
