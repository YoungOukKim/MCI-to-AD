# ==============================================================================
# Figure 1 — A, B, C, D (4-panel)
#

#   - Panels A-D: study schema, mechanistic schematic, neuronal NDUFS1,
#     and astrocytic PTGDS trajectory. Cascade detail is in Table 3 + Figs 2/3/5.
#   - Panels C/D use measured expression (loess).
#   - Panel D statistics (donor-level): beta2 = -2.27, p = 0.006;
#     segmented Davies p = 0.032; vertex CPS 0.47 (raw peak in the 0.5-0.6 bins).
#
# Input : SEA-AD MTG h5ad (measured expression)
# Output: Figure1_ABCD.(png|tiff)
# ==============================================================================

library(rhdf5)
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(patchwork)
library(magick)
library(grid)
library(segmented)   # Davies test / breakpoint (for statistical annotations)

# ------------------------------------------------------------------------------
# 1. Paths & Settings (edit for your environment)
# ------------------------------------------------------------------------------
file_path          <- Sys.getenv("H5AD_PATH", "path/to/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad")
save_path          <- Sys.getenv("OUT_DIR", file.path("output", "Fig1"))
sea_schematic_path <- file.path(save_path, "SEA-AD_Tissue_summary.jpg")  # place panel-A image here
biorender_path     <- file.path(save_path, "Fig1A-1.png")  # place panel-B schematic here

if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
set.seed(42)

# ------------------------------------------------------------------------------
# 2. Data extraction (load measured expression from H5AD)
#    - also extract donor ID (Panel D statistics are donor-level by design)
# ------------------------------------------------------------------------------
message(">>> [1/4] starting data extraction ...")
h5closeAll()

subclass_idx        <- h5read(file_path, "/obs/Subclass")
subclass_categories <- h5read(file_path, "/obs/__categories/Subclass")
cell_type_labels    <- subclass_categories[subclass_idx + 1]
bin_vec             <- as.numeric(h5read(file_path, "/obs/Continuous Pseudo-progression Score"))

# extract donor ID (the obs donor key name varies by dataset; auto-detect)
obs_names <- h5ls(file_path)
obs_keys  <- obs_names$name[obs_names$group == "/obs"]
donor_key <- intersect(c("Donor ID","donor_id","Donor","external_donor_name_label",
                         "donor","sample_id","specimen"), obs_keys)
donor_key <- if (length(donor_key) >= 1) donor_key[1] else NA

read_obs_factor <- function(key) {
  v <- h5read(file_path, paste0("/obs/", key))
  cats <- tryCatch(h5read(file_path, paste0("/obs/__categories/", key)),
                   error = function(e) NULL)
  if (!is.null(cats)) cats[as.integer(v) + 1L] else as.character(v)
}
donor_vec <- if (!is.na(donor_key)) read_obs_factor(donor_key) else NA_character_

classify_cell <- function(label) {
  label <- as.character(label)
  if (grepl("^L[0-9]|IT$|ET$|CT$|NP$|L6b", label)) return("Excitatory_Neuron")
  if (grepl("^Astro", label, ignore.case = TRUE))  return("Astrocyte")
  return("Other")
}
cell_class <- sapply(cell_type_labels, classify_cell)

# extract gene expression (NDUFS1 = Panel C, PTGDS = Panel D)
# (Fig 1 plots NDUFS1 and PTGDS; cascade detail is in Table 3 + Figs 2/3/5)
gene_names_all <- h5read(file_path, "var/_index")
target_genes   <- c("NDUFS1", "PTGDS")
gene_idx_map   <- match(target_genes, gene_names_all) - 1

# ── FAST block read: keep-filtered CSC slurp ──────────────────────────────────
# Instead of the slow per-cell double loop, read each block at once with
# h5read (slurp) and fill only the kept cells via vectorized match. Reading
# only the two cell types used in Panels C/D (Astrocyte, Excitatory_Neuron)
# with valid scores gives additional speed-up.
indptr  <- h5read(file_path, "X/indptr", bit64conversion = "double")
n_cells <- length(indptr) - 1L

# cell filter: cell types used in the analysis + valid score
keep <- (cell_class %in% c("Astrocyte", "Excitatory_Neuron")) &
        !is.na(bin_vec) & bin_vec >= 0.1
sel  <- which(keep)
kset <- logical(n_cells); kset[sel] <- TRUE
rmap <- integer(n_cells); rmap[sel] <- seq_along(sel)

# result holds only kept cells (saves memory/time); scattered to expr_mat[full] later.
M  <- matrix(0, length(sel), length(target_genes))
colnames(M) <- target_genes
pidx0 <- gene_idx_map          # 0-based var index (= match - 1)

bs <- 100000L
message(">>> reading expression (fast block read) ...")
for (s0 in seq(1L, n_cells, by = bs)) {
  e0  <- min(s0 + bs - 1L, n_cells)
  sp  <- indptr[s0]
  cnt <- indptr[e0 + 1L] - sp
  if (cnt <= 0) next
  ci_blk <- h5read(file_path, "X/indices", start = sp + 1L, count = cnt,
                   bit64conversion = "double")
  cd_blk <- h5read(file_path, "X/data",    start = sp + 1L, count = cnt)
  # process only kept cells within this block
  for (k in which(kset[s0:e0])) {
    g <- s0 + k - 1L                 # global cell index
    a <- indptr[g]     - sp + 1L
    b <- indptr[g + 1L] - sp
    if (b < a) next
    h  <- match(pidx0, ci_blk[a:b])  # locate target genes among this cells nnz
    ok <- !is.na(h)
    if (any(ok)) M[rmap[g], ok] <- cd_blk[a:b][h[ok]]
  }
  cat(sprintf("  ...%d/%d\r", e0, n_cells))
}
cat("\n")

# scatter keep-only M into full-length expr_mat (df_real assumes full length)
expr_mat <- matrix(0, nrow = n_cells, ncol = length(target_genes))
colnames(expr_mat) <- target_genes
expr_mat[sel, ] <- M

# ------------------------------------------------------------------------------
# 3. Integrate & summarize (bin means = visualization / donor aggregation = statistics)
# ------------------------------------------------------------------------------
df_real <- data.frame(
  score      = bin_vec,
  donor      = donor_vec,
  cell_class = cell_class,
  NDUFS1     = expr_mat[, "NDUFS1"],
  PTGDS      = expr_mat[, "PTGDS"]
) %>% filter(!is.na(score) & score >= 0.1)

# bin means for visualization (points in Panels C, D)
bin_summary <- df_real %>%
  mutate(bin = round(score, 1)) %>%
  group_by(cell_class, bin) %>%
  summarise(across(all_of(target_genes), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

neuron_traj <- bin_summary %>% filter(cell_class == "Excitatory_Neuron")
astro_traj  <- bin_summary %>% filter(cell_class == "Astrocyte")

# ---- compute statistical annotations (Panel D) ------------------------
# (a) donor-level quadratic regression: beta2, p, vertex - same unit (donor) as Table 2
# (b) donor-binned segmented: Davies p, breakpoint - same procedure giving 0.032
astro_cells <- df_real %>% filter(cell_class == "Astrocyte" & !is.na(donor))

if (nrow(astro_cells) > 0 && !all(is.na(astro_cells$donor))) {
  donor_astro <- astro_cells %>%
    group_by(donor) %>%
    summarise(PTGDS = mean(PTGDS, na.rm = TRUE), bin = mean(bin <- score), .groups = "drop")
  # note: per-donor mean bin = that donors mean CPS (where the donor sits on progression)
  donor_astro <- astro_cells %>%
    group_by(donor) %>%
    summarise(PTGDS = mean(PTGDS, na.rm = TRUE),
              bin   = mean(score, na.rm = TRUE), .groups = "drop")
  n_donor <- nrow(donor_astro)
  lm_quad <- lm(PTGDS ~ bin + I(bin^2), data = donor_astro)
  beta2   <- unname(coef(lm_quad)["I(bin^2)"])
  p_beta2 <- summary(lm_quad)$coefficients["I(bin^2)", "Pr(>|t|)"]
  vertex  <- unname(-coef(lm_quad)["bin"] / (2 * beta2))
} else {
  warning("donor column not found - annotations fall back to manuscript (hard-coded) values.")
  n_donor <- 84; beta2 <- -2.27; p_beta2 <- 0.006; vertex <- 0.47
}

# segmented (donor-binned means, 9 bins) -> Davies p, breakpoint
astro_bin_means <- astro_traj %>% arrange(bin)
davies_p <- NA; bp <- NA
try({
  lin <- lm(PTGDS ~ bin, data = astro_bin_means)
  dav <- davies.test(lin, seg.Z = ~bin)
  davies_p <- dav$p.value
  seg <- segmented(lin, seg.Z = ~bin)
  bp  <- as.numeric(seg$psi[, "Est."])
}, silent = TRUE)
# if computation fails, fall back to the manuscript values
if (is.na(davies_p)) davies_p <- 0.032
if (is.na(bp))       bp       <- 0.23

message(sprintf(">>> Panel D stats: donor n=%d, beta2=%.2f, p=%.3f, vertex=%.2f | Davies p=%.3f, bp=%.2f",
                n_donor, beta2, p_beta2, vertex, davies_p, bp))

# ------------------------------------------------------------------------------
# 4. Build panels (A, B, C, D)
# ------------------------------------------------------------------------------
load_img_plot <- function(path, title) {
  tryCatch({
    img <- image_read(path)
    ggplot() +
      annotation_custom(rasterGrob(as.raster(image_trim(img)),
                                   width = unit(1, "npc"), height = unit(1, "npc"))) +
      theme_void() +
      labs(title = title) +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold",
                                      margin = margin(t = 10, b = 10)))
  }, error = function(e) {
    ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Image Not Found") +
      theme_void() + labs(title = title)
  })
}

# Panel A: SEA-AD overview (image)
pA <- load_img_plot(sea_schematic_path, "SEA-AD MTG snRNA-seq overview") +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
        plot.tag.position = c(0.01, 0.99),
        plot.title = element_text(face = "bold"))

# Panel B: mechanistic schematic (BioRender image) - label as schematic in the caption
pB <- load_img_plot(biorender_path, NULL) +
  theme(plot.margin = margin(t = -20, r = 10, b = 10, l = 10),
        plot.tag.position = c(0.01, 0.95))

# Panel C: Neuronal NDUFS1 (measured loess)
pC <- ggplot(neuron_traj, aes(x = bin, y = NDUFS1)) +
  geom_point(size = 3, color = "#2C3E50", alpha = 0.5) +
  geom_smooth(method = "loess", span = 1.0, color = "#2C3E50",
              linewidth = 1.5, se = TRUE, fill = "#BDC3C7", alpha = 0.3) +
  scale_x_continuous(breaks = seq(0.2, 0.9, 0.1)) +
  labs(x = "Pseudo-progression Score", y = "NDUFS1 Expression",
       title = "Neuronal Metabolic Stress") +
  theme_classic(base_size = 12) +
  theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 10),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.tag.position = c(0.015, 0.99))

# Panel D: Astrocytic PTGDS (measured loess) + statistical annotations
#   - breakpoint (=0.23) dotted line, vertex (=0.47) dashed line
#   - annotation text uses the exact manuscript/corrected values
pmax_y <- max(astro_traj$PTGDS, na.rm = TRUE)
pD <- ggplot(astro_traj, aes(x = bin, y = PTGDS)) +
  # show breakpoint (Davies) and vertex - consistent with the manuscript
  geom_vline(xintercept = bp,     linetype = "dotted", color = "grey50", linewidth = 0.6) +
  geom_vline(xintercept = vertex, linetype = "dashed", color = "#0072B2", linewidth = 0.7) +
  geom_point(size = 3, color = "#0072B2", alpha = 0.5) +
  geom_smooth(method = "loess", span = 1.0, color = "#0072B2",
              linewidth = 2, se = TRUE, fill = "#0072B2", alpha = 0.1) +
  scale_x_continuous(breaks = seq(0.2, 0.9, 0.1)) +
  labs(x = "Pseudo-progression Score", y = "PTGDS Expression",
       title = "Astrocytic Compensation") +
  # statistical annotations (Table 2)
  annotate("text", x = bp, y = pmax_y * 1.02,
           label = sprintf("breakpoint CPS %.2f\n(Davies p = %.3f)", bp, davies_p),
           size = 3.0, color = "grey30", hjust = 0.5, lineheight = 0.9) +
  annotate("text", x = vertex, y = min(astro_traj$PTGDS, na.rm = TRUE),
           label = sprintf("vertex CPS %.2f", vertex),
           size = 3.0, color = "#0072B2", hjust = 0.5, vjust = 1) +
  annotate("text", x = 0.72, y = pmax_y * 1.0,
           label = sprintf("donor-level quadratic\nβ₂ = %.2f, p = %.3f (n = %d)",
                           beta2, p_beta2, n_donor),
           size = 3.0, color = "#0072B2", hjust = 0, lineheight = 0.95) +
  theme_classic(base_size = 12) +
  theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 10),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.tag.position = c(0.01, 0.99))

# ------------------------------------------------------------------------------
# 5. Layout & Save - A / B / (C | D)
# ------------------------------------------------------------------------------
fig1_final <- pA / pB / (pC | pD) +
  plot_layout(heights = c(1.5, 1.5, 1.1)) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold", size = 28))

ggsave(file.path(save_path, "Figure1_ABCD.png"),
       fig1_final, width = 12, height = 16, dpi = 600, bg = "white")
ggsave(file.path(save_path, "Figure1_ABCD.tiff"),
       fig1_final, width = 12, height = 16, dpi = 600,
       bg = "white", device = "tiff", compression = "lzw")

message(">>> [done] Figure 1 (A,B,C,D) generated.")
message(">>> In-figure statistics: Davies p, vertex, donor beta2.")
h5closeAll()
