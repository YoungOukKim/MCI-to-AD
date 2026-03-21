# ==============================================================================
# FigS2_PTGDS_inflection_validation.R
#
# Purpose : Additional file 1: Fig. S2
#           Statistical validation of the biphasic PTGDS trajectory and
#           inflection point
#
# Panels
#   A  Unsmoothed bin-level mean PTGDS expression (± SEM)
#   B  3-bin centered moving average trajectory
#   C  LOESS regression (span = 1.0, Fig. 1D style)
#   D  Segmented regression with breakpoint CI and Davies' test
#
# Key results (matches Supp caption)
#   Breakpoint: Bin 0.23 (95% CI: 0.13-0.33)
#   Davies' test p = 0.032
#   Pre-breakpoint slope: +4.50 (rising)
#   Post-breakpoint slope: -0.52 (declining)
#
# Input   : output/SEAAD_processed_data.RData  (from 01_SEAAD_data_extraction.R)
# Output  : output/FigS2/FigS2.png + .tiff  (600 dpi, 14 x 10 in)
#
# Paper: Kim Y et al. (2025). bioRxiv. https://github.com/YoungOukKim/MCI-to-AD
# Usage: Rscript analysis/FigS/FigS2_PTGDS_inflection_validation.R
# Requirements: R >= 4.3.2 | ggplot2, dplyr, segmented, data.table, patchwork
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(segmented)
  library(data.table)
  library(patchwork)
})
set.seed(42)

# Paths (relative to repository root)
OUT_DIR <- file.path("output", "FigS2")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. Load pre-processed data (from 01_SEAAD_data_extraction.R)
# ==============================================================================
message(">>> Loading SEAAD_processed_data.RData...")
load(file.path("output", "SEAAD_processed_data.RData"))

astro_df <- df[cell_type_label == "Astrocyte" &
                 !is.na(bin_rounded) & bin_rounded >= 0.1 & bin_rounded <= 0.9]

bin_means <- astro_df[, .(
  PTGDS_mean = mean(PTGDS, na.rm = TRUE),
  PTGDS_sd   = sd(PTGDS,   na.rm = TRUE),
  n_cells    = .N
), by = .(bin = bin_rounded)][order(bin)]

bin_means[, PTGDS_se := PTGDS_sd / sqrt(n_cells)]
bin_means[, PTGDS_ma := {
  x <- PTGDS_mean; n <- length(x); ma <- numeric(n)
  for (i in seq_len(n)) {
    ma[i] <- mean(x[max(1L, i-1L):min(n, i+1L)], na.rm = TRUE)
  }
  ma
}]

message(sprintf("  Astrocytes: %d cells across %d bins",
                sum(bin_means$n_cells), nrow(bin_means)))

# ==============================================================================
# 2. Segmented regression
#    CI via confint.segmented() — matches Table 2 in paper (0.13-0.33)
# ==============================================================================
lm_fit  <- lm(PTGDS_mean ~ bin, data = bin_means)
seg_fit <- segmented(lm_fit, seg.Z = ~ bin, psi = 0.5)

bp <- summary(seg_fit)$psi[, "Est."]

ci_result <- tryCatch({
  ci_raw <- confint.segmented(seg_fit)
  ci_mat <- if (is.list(ci_raw) && !is.matrix(ci_raw)) ci_raw[[1]] else ci_raw
  list(lo = ci_mat[1, 2], hi = ci_mat[1, 3])
}, error = function(e) {
  message("  confint.segmented() fallback to +/-1.96*SE: ", e$message)
  bp_se <- summary(seg_fit)$psi[, "St.Err"]
  list(lo = bp - 1.96 * bp_se, hi = bp + 1.96 * bp_se)
})
ci_lo <- ci_result$lo
ci_hi <- ci_result$hi

davies_p  <- davies.test(lm_fit, ~ bin)$p.value
sl        <- slope(seg_fit)
slope_pre  <- round(sl$bin[1, "Est."], 2)
slope_post <- round(sl$bin[2, "Est."], 2)

x_fine <- data.frame(bin = seq(0.1, 0.9, by = 0.005))
x_fine$seg_pred <- predict(seg_fit, newdata = x_fine)

message(sprintf("  Breakpoint: %.3f (95%% CI: %.2f-%.2f)", bp, ci_lo, ci_hi))
message(sprintf("  Davies' p = %.4f  |  Slopes: %+.2f / %+.2f",
                davies_p, slope_pre, slope_post))

# ==============================================================================
# 3. Build panels
# ==============================================================================
y_lo <- min(bin_means$PTGDS_mean) - 0.08
y_hi <- max(bin_means$PTGDS_mean) + 0.15
x_sc <- scale_x_continuous(breaks = seq(0.1, 0.9, 0.1))
ylims <- coord_cartesian(ylim = c(y_lo, y_hi))

theme_s2 <- theme_classic(base_size = 14) +
  theme(plot.title  = element_blank(),
        plot.margin = margin(10, 15, 10, 15),
        axis.title  = element_text(size = 12),
        axis.text   = element_text(size = 11, color = "black"))

# Panel A — Raw expression with SEM error bars
pA <- ggplot(bin_means, aes(x = bin, y = PTGDS_mean)) +
  geom_pointrange(aes(ymin = PTGDS_mean - PTGDS_se,
                      ymax = PTGDS_mean + PTGDS_se),
                  size = 0.6, color = "#0072B2", linewidth = 0.5) +
  geom_line(color = "#0072B2", linewidth = 0.8) +
  x_sc + ylims + theme_s2 +
  labs(x = "Pseudo-progression bin",
       y = "Mean PTGDS expression\n(log-normalized)")

# Panel B — 3-bin moving average
pB <- ggplot(bin_means, aes(x = bin)) +
  geom_point(aes(y = PTGDS_mean), size = 2.5, alpha = 0.35, color = "gray55") +
  geom_line(aes(y = PTGDS_ma), linewidth = 1.3, color = "#0072B2") +
  geom_point(aes(y = PTGDS_ma), size = 3, color = "#0072B2") +
  x_sc + ylims + theme_s2 +
  labs(x = "Pseudo-progression bin", y = "PTGDS expression")

# Panel C — LOESS regression (Fig. 1D style: span=1.0, #0072B2)
pC <- ggplot(bin_means, aes(x = bin, y = PTGDS_mean)) +
  geom_point(size = 3, color = "#0072B2", alpha = 0.6) +
  geom_smooth(method = "loess", span = 1.0, se = TRUE,
              color = "#0072B2", fill = "#0072B2",
              linewidth = 2, alpha = 0.1) +
  x_sc + ylims + theme_s2 +
  labs(x = "Pseudo-progression bin", y = "PTGDS expression")

# Panel D — Segmented regression + LOESS underlay
bp_label <- sprintf(
  "Breakpoint = %.3f\n(95%% CI: %.2f-%.2f)\nDavies' p = %.4f\nSlope: %+.2f -> %+.2f",
  bp, round(ci_lo, 2), round(ci_hi, 2), davies_p, slope_pre, slope_post)

pD <- ggplot() +
  annotate("rect", xmin=ci_lo, xmax=ci_hi,
           ymin=-Inf, ymax=Inf, alpha=0.12, fill="#D6604D") +
  geom_vline(xintercept=bp, linetype="dashed", color="#D6604D", linewidth=0.8) +
  geom_smooth(data=bin_means, aes(x=bin, y=PTGDS_mean),
              method="loess", span=1.0, se=FALSE,
              color="gray70", linewidth=1.0) +
  geom_line(data=x_fine, aes(x=bin, y=seg_pred),
            linewidth=1.5, color="#B2182B") +
  geom_point(data=bin_means, aes(x=bin, y=PTGDS_mean),
             size=3.5, color="#0072B2") +
  annotate("label", x=0.92, y=y_lo+0.08,
           label=bp_label, hjust=1, vjust=0,
           size=3.3, fontface="italic",
           fill="white", label.size=0.3, alpha=0.95) +
  x_sc + ylims + theme_s2 +
  labs(x="Pseudo-progression bin", y="PTGDS expression")

# Combine
fig_s2 <- (pA | pB) / (pC | pD) +
  plot_layout(heights=c(1,1)) +
  plot_annotation(tag_levels='A') &
  theme(plot.tag=element_text(face="bold", size=28),
        plot.margin=margin(5,5,5,5))

out_png  <- file.path(OUT_DIR, "FigS2.png")
out_tiff <- file.path(OUT_DIR, "FigS2.tiff")
ggsave(out_png,  fig_s2, width=14, height=10, dpi=600, bg="white")
ggsave(out_tiff, fig_s2, width=14, height=10, dpi=600,
       device="tiff", compression="lzw", bg="white")

message(">>> Saved: ", out_png)
cat(strrep("=", 60), "\n")
cat(sprintf("  Fig. S2 -- Manuscript values\n"))
cat(sprintf("  Breakpoint: Bin %.3f (95%% CI: %.2f-%.2f)\n", bp, ci_lo, ci_hi))
cat(sprintf("  Davies' p = %.4f\n", davies_p))
cat(sprintf("  Pre-slope: %+.2f  |  Post-slope: %+.2f\n", slope_pre, slope_post))
cat(sprintf("  Raw peak: Bin %.1f (%.4f)\n",
            bin_means$bin[which.max(bin_means$PTGDS_mean)],
            max(bin_means$PTGDS_mean)))
cat(strrep("=", 60), "\n")
message(">>> Script FigS2 complete.")
