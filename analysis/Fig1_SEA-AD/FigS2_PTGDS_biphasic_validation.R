# ==============================================================================
# 04_FigS2_PTGDS_biphasic_validation.R
#
# Purpose : Statistical validation of the PTGDS biphasic trajectory and
#           inflection point (Additional file 1: Fig. S2)
# Input   : output/SEAAD_processed_data.RData  (from Script 01)
#
# Outputs → Paper figure panels
#   output/FigS2/FigS2A_raw.png           →  Fig. S2A (unsmoothed expression)
#   output/FigS2/FigS2B_moving_average.png →  Fig. S2B (3-bin MA)
#   output/FigS2/FigS2C_loess.png         →  Fig. S2C (LOESS regression)
#   output/FigS2/FigS2D_segmented.png     →  Fig. S2D (segmented regression)
#   output/FigS2/FigS2_combined.png/pdf   →  Combined 2×2 panel
#
# Key statistical results (reproduced from paper)
#   Inflection breakpoint : Bin ≈ 0.23  (95% CI: 0.13–0.33)
#   Davies' test          : p = 0.032
#   Pre-breakpoint slope  : +4.50  (compensatory rise)
#   Post-breakpoint slope : −0.52  (gradual deceleration)
#   Raw peak              : Bin 0.5 (mean = 2.02)
#
# Paper reference
#   Results: "segmented regression on unsmoothed bin-level means (n = 9 bins)
#             identified a significant early inflection at Bin 0.23
#             (95% CI: 0.13–0.33; Davies' p = 0.032)"
#   Kim Y†, Heo W†, Park S†, et al. (2025). bioRxiv.
#   https://github.com/YoungOukKim/MCI-to-AD
#
# Usage
#   Rscript analysis/Fig1_SEA-AD/04_FigS2_PTGDS_biphasic_validation.R
#
# Requirements
#   R >= 4.3.2  |  ggplot2, dplyr, segmented, data.table
#   Optional   : patchwork  (for combined 2×2 panel)
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(segmented)
  library(data.table)
})

set.seed(42)

# ── Paths (relative to repository root) ───────────────────────────────────────
OUT_DIR <- file.path("output", "FigS2")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Shared ggplot2 theme
theme_paper <- function(base_size = 14) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      plot.title   = element_text(face = "bold", size = base_size - 1),
      axis.text    = element_text(size = base_size - 2),
      axis.title   = element_text(size = base_size - 1)
    )
}

# ══════════════════════════════════════════════════════════════════════════════
# 1. Load and prepare astrocyte PTGDS bin means
# ══════════════════════════════════════════════════════════════════════════════
message(">>> [1/5] Loading data and computing bin means...")

load(file.path("output", "SEAAD_processed_data.RData"))

astro_df <- df[
  cell_type_label == "Astrocyte" &
    !is.na(bin_rounded) & bin_rounded >= 0.1 & bin_rounded <= 0.9]

bin_means <- astro_df[, .(
  PTGDS_mean = mean(PTGDS, na.rm = TRUE),
  n_cells    = .N
), by = .(bin = bin_rounded)][order(bin)]

message("  Astrocyte PTGDS bin means (unsmoothed):")
print(bin_means)

# 3-bin centered moving average
bin_means[, PTGDS_ma := {
  x  <- PTGDS_mean
  n  <- length(x)
  ma <- numeric(n)
  for (i in seq_len(n)) {
    idx   <- max(1L, i - 1L):min(n, i + 1L)
    ma[i] <- mean(x[idx])
  }
  ma
}]

# ══════════════════════════════════════════════════════════════════════════════
# 2. Panel A — Raw unsmoothed expression
# ══════════════════════════════════════════════════════════════════════════════
message(">>> [2/5] Panel A: raw expression...")

p_A <- ggplot(bin_means, aes(x = bin, y = PTGDS_mean)) +
  geom_point(size = 3.5, color = "#2166AC") +
  geom_line( color = "#2166AC", linewidth = 0.8) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.1),
                     labels = paste0("Bin ", seq(0.1, 0.9, 0.1))) +
  theme_paper() +
  labs(title = "A. Raw PTGDS expression (unsmoothed)",
       x = "Pseudo-progression score (CPS bin)",
       y = "Mean PTGDS expression\n(log-normalized)")

# ══════════════════════════════════════════════════════════════════════════════
# 3. Panel B — 3-bin centered moving average
# ══════════════════════════════════════════════════════════════════════════════
message(">>> [3/5] Panel B: 3-bin moving average...")

p_B <- ggplot(bin_means, aes(x = bin)) +
  geom_point(aes(y = PTGDS_mean), size = 2.5, alpha = 0.35, color = "gray50") +
  geom_line( aes(y = PTGDS_ma),   linewidth = 1.3,  color = "#2166AC") +
  geom_point(aes(y = PTGDS_ma),   size = 3.0,       color = "#2166AC") +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.1)) +
  theme_paper() +
  labs(title = "B. 3-bin centered moving average",
       x = "Pseudo-progression score (CPS bin)",
       y = "PTGDS expression")

# ══════════════════════════════════════════════════════════════════════════════
# 4. Panel C — LOESS regression (span = 0.6)
# ══════════════════════════════════════════════════════════════════════════════
message(">>> [4/5] Panel C: LOESS regression (span = 0.6)...")

p_C <- ggplot(bin_means, aes(x = bin, y = PTGDS_mean)) +
  geom_point(size = 3.5, color = "#2166AC") +
  geom_smooth(method = "loess", span = 0.6, se = TRUE,
              color = "#B2182B", fill = "#FDDBC7", linewidth = 1.2) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.1)) +
  theme_paper() +
  labs(title = "C. LOESS regression (span = 0.6)",
       x = "Pseudo-progression score (CPS bin)",
       y = "PTGDS expression")

# ══════════════════════════════════════════════════════════════════════════════
# 5. Panel D — Segmented regression with Davies' test
# ══════════════════════════════════════════════════════════════════════════════
message(">>> [5/5] Panel D: segmented regression + Davies' test...")

lm_fit  <- lm(PTGDS_mean ~ bin, data = bin_means)
seg_fit <- segmented(lm_fit, seg.Z = ~ bin, psi = 0.5)

bp       <- summary(seg_fit)$psi[, "Est."]
bp_se    <- summary(seg_fit)$psi[, "St.Err"]
ci_lo    <- bp - 1.96 * bp_se
ci_hi    <- bp + 1.96 * bp_se
davies_p <- davies.test(lm_fit, ~ bin)$p.value
slopes   <- slope(seg_fit)

message(sprintf("  Breakpoint: %.3f  (95%% CI: %.3f–%.3f)",   bp, ci_lo, ci_hi))
message(sprintf("  Davies' test p-value: %.4f",                davies_p))
message(sprintf("  Pre-breakpoint slope:  %.4f",  slopes$bin[1, "Est."]))
message(sprintf("  Post-breakpoint slope: %.4f",  slopes$bin[2, "Est."]))

bin_means[, seg_pred := predict(seg_fit)]

p_D <- ggplot(bin_means, aes(x = bin, y = PTGDS_mean)) +
  annotate("rect", xmin = ci_lo, xmax = ci_hi,
           ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "#D6604D") +
  geom_point(size = 3.5, color = "#2166AC") +
  geom_line( aes(y = seg_pred), linewidth = 1.3, color = "#B2182B") +
  geom_vline(xintercept = bp, linetype = "dashed",
             color = "#D6604D", linewidth = 0.9) +
  annotate("text",
           x = bp + 0.045,
           y = max(bin_means$PTGDS_mean) * 0.985,
           label = sprintf("BP = %.3f\n(%.3f\u2013%.3f)\nDavies p = %.4f",
                           bp, ci_lo, ci_hi, davies_p),
           hjust = 0, size = 3.4, fontface = "italic") +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.1)) +
  theme_paper() +
  labs(title = "D. Segmented regression \u2014 inflection point",
       x = "Pseudo-progression score (CPS bin)",
       y = "PTGDS expression")

# ══════════════════════════════════════════════════════════════════════════════
# 6. Save individual panels (300 dpi PNG)
# ══════════════════════════════════════════════════════════════════════════════
for (panel in list(list(p_A,"FigS2A_raw"),
                   list(p_B,"FigS2B_moving_average"),
                   list(p_C,"FigS2C_loess"),
                   list(p_D,"FigS2D_segmented"))) {
  ggsave(file.path(OUT_DIR, paste0(panel[[2]], ".png")), panel[[1]],
         width = 5, height = 4, dpi = 300, bg = "white")
}
message("  Individual panels saved to ", OUT_DIR)

# ══════════════════════════════════════════════════════════════════════════════
# 7. Combined 2×2 panel (patchwork, if available)
# ══════════════════════════════════════════════════════════════════════════════
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  combined <- (p_A | p_B) / (p_C | p_D)
  ggsave(file.path(OUT_DIR, "FigS2_combined.png"), combined,
         width = 11, height = 9, dpi = 300, bg = "white")
  ggsave(file.path(OUT_DIR, "FigS2_combined.pdf"), combined,
         width = 11, height = 9, bg = "white")
  message("  Combined panel saved.")
} else {
  message("  patchwork not installed — skipping combined panel.")
  message("  Install with: install.packages('patchwork')")
}

# ══════════════════════════════════════════════════════════════════════════════
# 8. Manuscript-ready summary statistics
# ══════════════════════════════════════════════════════════════════════════════
cat("\n", strrep("=", 65), "\n")
cat("  Fig. S2 — Results for manuscript\n")
cat(strrep("=", 65), "\n")
cat(sprintf("  Breakpoint estimate : Bin %.3f\n",           bp))
cat(sprintf("  95%% CI             : %.3f \u2013 %.3f\n",   ci_lo, ci_hi))
cat(sprintf("  Davies' test p      : %.4f\n",               davies_p))
cat(sprintf("  Pre-BP slope        : %.4f  (rising)\n",     slopes$bin[1, "Est."]))
cat(sprintf("  Post-BP slope       : %.4f  (declining)\n",  slopes$bin[2, "Est."]))
cat(sprintf("  Raw peak            : Bin %.1f  (%.4f)\n",
            bin_means$bin[which.max(bin_means$PTGDS_mean)],
            max(bin_means$PTGDS_mean)))
cat(sprintf("  MA  peak            : Bin %.1f  (%.4f)\n",
            bin_means$bin[which.max(bin_means$PTGDS_ma)],
            max(bin_means$PTGDS_ma)))
cat(strrep("=", 65), "\n")

cat("\nFig. S2 legend (fill breakpoint values from output above):\n")
cat("Statistical validation of the PTGDS biphasic trajectory. ")
cat("A. Unsmoothed mean PTGDS expression across 9 pseudo-progression bins ")
cat("(1.3M nuclei, 84 donors; SEA-AD MTG). ")
cat("B. 3-bin centered moving average used for cross-correlation analyses (Table 2). ")
cat("C. LOESS regression (span = 0.6) confirming biphasic curvature. ")
cat(sprintf("D. Segmented regression identifying an inflection at Bin %.3f ", bp))
cat(sprintf("(95%% CI: %.3f\u2013%.3f; Davies' test p = %.4f).\n", ci_lo, ci_hi, davies_p))

message("\n>>> Script 04 complete.")
