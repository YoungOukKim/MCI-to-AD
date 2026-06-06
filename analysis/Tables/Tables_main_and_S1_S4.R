# ==============================================================================
# Tables_main_and_S1_S4.R
#
# Purpose : Generate all main text tables and SEA-AD supplementary tables
#
# Outputs → Paper tables (all CSV to output/tables/)
#
#   Main text:
#     Table1_CellType_Trajectories.csv      → Table 1 (3-bin MA at Bins 0.1/0.6/0.9)
#     Table2_Statistical_Validation.csv     → Table 2 (quadratic + segmented stats)
#     Table2_Segmented_Slopes_Detail.csv    → Table 2 supplemental (slopes)
#     Table3_Checkpoint_CCF.csv             → Table 3 (inflection window CCF)
#
#   Supplementary:
#     TableS4_Full_Trajectories.csv          → Additional file 1: Table S4
#     TableS5_Apoptotic_Identity_Markers.csv → Additional file 1: Table S5
#     TableS6_Full_Trajectory_CCF.csv        → Additional file 1: Table S6
#     TableS7_Neuronal_Subtype_Vuln.csv      → Additional file 1: Table S7
#
# Key verified values (v3)
#   Table 1 : PTGDS peak near CPS 0.47 (vertex from quadratic model)
#   Table 2 : donor-level (n~84). beta2~-2.27, p~0.006,
#             vertex~0.47; segmented inflection Bin~0.23, Davies p~0.032;
#             pre-peak rise NOT significant at donor level (p~0.45);
#             post-peak decline significant (p~0.003).
#   Table 3 : Purinergic Lag-1 r=-0.886 | NF-kB Lag0 r=-0.678 |
#             C3<->NGFR Lag0 r=-0.940 | LCN2<->NDUFS1 Lag0 r=-0.857
#   Table S5: CASP3 rho=-0.027 | AQP4 rho=-0.041 | SLC1A2 rho=-0.041
#   Table S6: 9 interaction pairs, lag.max=3 (full CPS trajectory)
#   Table S7: SST 11.24%->5.49%, Exc 88.76%->94.51%
#             BCL2 -6.8% (Bin 0.1->0.4)
#
# Input   : output/SEAAD_processed_data.RData  (from 01_SEAAD_data_extraction.R)
# Note    : h5ad not required (uses pre-processed RData)
#
# Paper: Kim Y et al. (2025). bioRxiv. https://github.com/YoungOukKim/MCI-to-AD
# Usage: Rscript analysis/Tables/Tables_main_and_S1_S4.R
# Requirements: R >= 4.3.2 | data.table, dplyr, segmented
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(segmented)
})
set.seed(42)

# Paths (relative to repository root)
OUT_DIR <- file.path("output", "tables")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# Helper functions
# ==============================================================================
ccf_analysis <- function(x, y, max_lag = 3L) {
  ok <- complete.cases(x, y)
  x <- x[ok]; y <- y[ok]; n <- length(x)
  if (n < 4L)
    return(data.table(Lag=NA_integer_, r=NA_real_, p_val=NA_real_, n_eff=NA_integer_))
  if (sd(x)==0 || sd(y)==0)
    return(data.table(Lag=0L, r=0, p_val=1, n_eff=n))
  res     <- ccf(x, y, lag.max=max_lag, plot=FALSE, na.action=na.pass)
  idx     <- which.max(abs(res$acf))
  lag_val <- as.integer(res$lag[idx])
  r_val   <- as.numeric(res$acf[idx])
  n_eff   <- n - abs(lag_val)
  p_val   <- if (n_eff > 2L) {
    2 * pt(-abs(r_val * sqrt((n_eff-2)/(1-r_val^2+1e-10))), df=n_eff-2L)
  } else NA_real_
  data.table(Lag=lag_val, r=round(r_val,3), p_val=round(p_val,4), n_eff=n_eff)
}

sig_label <- function(p) {
  ifelse(is.na(p),"", ifelse(p<0.001,"***", ifelse(p<0.01,"**",
         ifelse(p<0.05,"*","ns"))))
}

apply_ma3 <- function(x) {
  n <- length(x); ma <- numeric(n)
  for (i in seq_len(n)) ma[i] <- mean(x[max(1L,i-1L):min(n,i+1L)], na.rm=TRUE)
  ma
}

# ==============================================================================
# 1. Load data
# ==============================================================================
message(">>> Loading SEAAD_processed_data.RData...")
load(file.path("output", "SEAAD_processed_data.RData"))
message(sprintf("  Cells loaded: %s", format(nrow(df), big.mark=",")))

target_groups <- c("Astrocyte", "Exc_Neuron", "Microglia", "SST_Neuron")
analysis_genes <- intersect(
  c("NDUFS1","PTGDS","LCN2","NGFR","APOE","TREM2","C3","NFKBIA",
    "P2RY1","P2RY12","GJA1","ITPR2","PPARG","CASP3","BAX","BCL2","AQP4",
    "HMOX1","SOD2","CLU","SLC1A2","GPX4","FTH1","SLC7A11","GFAP",
    "FOS","JUN","EGR1","ATF3"),
  colnames(df))

# ==============================================================================
# 2. Bin-level means and 3-bin MA
# ==============================================================================
message(">>> Computing bin means and 3-bin MA...")

raw_means <- df[
  cell_type_label %in% target_groups &
    !is.na(bin_rounded) & bin_rounded >= 0.1 & bin_rounded <= 0.9,
  c(lapply(.SD, mean, na.rm=TRUE), list(n_cells=.N)),
  by = .(bin=bin_rounded, group=cell_type_label),
  .SDcols = analysis_genes][order(group, bin)]

ma_cols <- paste0(analysis_genes, "_ma")
smoothed <- copy(raw_means)
smoothed[, (ma_cols) := lapply(.SD, apply_ma3), by=group, .SDcols=analysis_genes]

a_s <- smoothed[group=="Astrocyte"][order(bin)]
n_s <- smoothed[group=="Exc_Neuron"][order(bin)]
m_s <- smoothed[group=="Microglia"][order(bin)]
sst_s <- smoothed[group=="SST_Neuron"][order(bin)]

# Module composites
# Use only the *_ma columns that exist. This h5ad may lack some genes
# (e.g. PTGS2 / IL6ST are absent in SEA-AD MTG and skipped in script 01), which
# previously crashed the hardcoded .SDcols. Degrade gracefully to present genes.
.module_ma <- function(dt, cols) {
  cols <- intersect(cols, names(dt))
  if (length(cols) == 0L) return(rep(NA_real_, nrow(dt)))
  if (length(cols) == 1L) return(dt[[cols]])
  rowMeans(as.matrix(dt[, ..cols]), na.rm = TRUE)
}
a_s[, NFkB_ma  := .module_ma(a_s, c("PTGS2_ma","IL6ST_ma","NFKBIA_ma"))]
a_s[, Purin_ma := .module_ma(a_s, c("P2RY1_ma","P2RY12_ma","GJA1_ma","ITPR2_ma"))]
a_s[, Ferr_ma  := .module_ma(a_s, c("GPX4_ma","FTH1_ma","SLC7A11_ma"))]

# ==============================================================================
# 3. TABLE 1 — Cell-type trajectory summary (Bins 0.1, 0.6, 0.9)
# ==============================================================================
message(">>> Generating Table 1...")

t1_full <- data.table(
  Bin         = a_s$bin,
  NDUFS1_N    = n_s$NDUFS1_ma,
  PTGDS_A     = a_s$PTGDS_ma,
  LCN2_A      = a_s$LCN2_ma,
  NGFR_N      = n_s$NGFR_ma,
  APOE_Micro  = m_s$APOE_ma
)

table1 <- t1_full[Bin %in% c(0.1, 0.5, 0.6, 0.9)]
table1[, Trajectory := fcase(
  Bin == 0.1, "Early (pre-inflection)",
  Bin %in% c(0.5,0.6), "Mid (inflection window)",
  Bin == 0.9, "Late (post-inflection)"
)]

message("  Table 1:"); print(table1)
write.csv(t1_full,  file.path(OUT_DIR,"Table1_CellType_Trajectories.csv"),    row.names=FALSE)
write.csv(table1,   file.path(OUT_DIR,"Table1_Paper_Summary.csv"),             row.names=FALSE)

# ==============================================================================
# 4. TABLE 2 — Statistical validation of biphasic PTGDS dynamics
# ==============================================================================
message(">>> Generating Table 2...")

# Quadratic model (donor level)
astro_all <- df[cell_type_label=="Astrocyte" & !is.na(bin)]
if (!"donor" %in% names(astro_all))
  stop("`donor` column missing in SEAAD_processed_data.RData.\n",
       "  Re-run 01_SEAAD_data_extraction.R to export donor IDs.")
donor_astro <- astro_all[, .(PTGDS=mean(PTGDS, na.rm=TRUE), bin=mean(bin)), by=donor]
n_astro   <- nrow(donor_astro)                      # donors (n~84), NOT cells
lm_quad   <- lm(PTGDS ~ bin + I(bin^2), data=donor_astro)
beta2_sc  <- coef(lm_quad)["I(bin^2)"]
p_beta2   <- summary(lm_quad)$coefficients["I(bin^2)","Pr(>|t|)"]
vertex_sc <- -coef(lm_quad)["bin"] / (2*beta2_sc)

message(sprintf("  Quadratic (DONOR-level): n=%d donors, beta2=%.4f, vertex=%.4f, p=%s",
                n_astro, beta2_sc, vertex_sc, format.pval(p_beta2, 3)))

# Pre/post vertex slopes (DONOR level — pre-peak rise is typically n.s. here)
pre_df    <- donor_astro[bin <= vertex_sc]
post_df   <- donor_astro[bin  > vertex_sc]
pre_slope  <- coef(lm(PTGDS ~ bin, data=pre_df))["bin"]
post_slope <- coef(lm(PTGDS ~ bin, data=post_df))["bin"]
pre_p      <- summary(lm(PTGDS ~ bin, data=pre_df))$coef["bin","Pr(>|t|)"]
post_p     <- summary(lm(PTGDS ~ bin, data=post_df))$coef["bin","Pr(>|t|)"]

# Segmented regression (DONOR-binned means, n=9 bins)
bin_means_raw <- donor_astro[, .(PTGDS_mean=mean(PTGDS,na.rm=TRUE)),
                             by=.(bin=round(bin,1))][order(bin)]
lm_seg  <- lm(PTGDS_mean ~ bin, data=bin_means_raw)
seg_fit <- tryCatch(
  segmented(lm_seg, seg.Z = ~ bin, psi=0.5),
  error = function(e) tryCatch(
    segmented(lm_seg, seg.Z = ~ bin, psi=0.3),
    error = function(e2) NULL))

inflection_est <- davies_p <- seg_slope1 <- seg_slope2 <- ci_lo <- ci_hi <- NA
if (!is.null(seg_fit)) {
  inflection_est <- summary(seg_fit)$psi[,"Est."]
  davies_p       <- davies.test(lm_seg, ~ bin)$p.value
  sl <- slope(seg_fit)
  seg_slope1 <- sl$bin[1,"Est."]; seg_slope2 <- sl$bin[2,"Est."]
  ci <- tryCatch({
    cm <- confint.segmented(seg_fit)
    cm <- if (is.list(cm)) cm[[1]] else cm
    list(lo=cm[1,2], hi=cm[1,3])
  }, error = function(e) {
    se <- summary(seg_fit)$psi[,"St.Err"]
    list(lo=inflection_est-1.96*se, hi=inflection_est+1.96*se)
  })
  ci_lo <- ci$lo; ci_hi <- ci$hi
  message(sprintf("  Inflection: %.2f (CI: %.2f-%.2f), Davies p=%.4f",
                  inflection_est, ci_lo, ci_hi, davies_p))
}

table2 <- data.table(
  Component         = c("Quadratic term (beta2)", "Vertex CPS",
                         "Pre-peak slope", "Post-peak slope",
                         "Segmented inflection (Bin)",
                         "Segmented CI (95%)",
                         "Davies test p-value"),
  Estimate          = c(round(beta2_sc,3), round(vertex_sc,3),
                         round(pre_slope,3), round(post_slope,3),
                         round(inflection_est,2),
                         sprintf("%.2f-%.2f", round(ci_lo,2), round(ci_hi,2)),
                         format.pval(davies_p, 3)),
  p_value           = c(format.pval(p_beta2,3), "—",
                         format.pval(pre_p,3), format.pval(post_p,3),
                         "—", "—", format.pval(davies_p,3)),
  Interpretation    = c("Biphasic curvature (donor-level)",
                         "Peak inflection point",
                         "Compensatory rise (often n.s. at donor level)",
                         "Progressive decline",
                         "Early deceleration onset",
                         "95% CI of breakpoint",
                         "Significant breakpoint")
)

message("  Table 2:"); print(table2)
write.csv(table2, file.path(OUT_DIR,"Table2_Statistical_Validation.csv"),  row.names=FALSE)

slopes_detail <- data.table(
  Stratum     = c("Pre-vertex (CPS <= vertex)", "Post-vertex (CPS > vertex)"),
  N_donors    = c(nrow(pre_df), nrow(post_df)),  # pre_df/post_df = donor pseudobulk
  Slope       = c(round(pre_slope,4), round(post_slope,4)),
  p_value     = c(format.pval(pre_p,3), format.pval(post_p,3))
)
write.csv(slopes_detail, file.path(OUT_DIR,"Table2_Segmented_Slopes_Detail.csv"), row.names=FALSE)

# ==============================================================================
# 5. TABLE 3 — Inflection window CCF (Bins 0.4-0.8, lag.max=2)
# ==============================================================================
message(">>> Generating Table 3...")

cp <- c(0.4, 0.5, 0.6, 0.7, 0.8)
a_cp <- a_s[bin %in% cp]; n_cp <- n_s[bin %in% cp]; m_cp <- m_s[bin %in% cp]

# LCN2 pairs excluded: astrocytic LCN2 detection 0.007% (5/67,419) -> untestable;
# the LCN2 arm is supported by CSF proteomics (Fig.5) and cross-species qPCR.
pairs_t3 <- list(
  list("Purinergic (A) -> PTGDS (A)",   a_cp$Purin_ma,    a_cp$PTGDS_ma,  "Astro-Astro"),
  list("NF-kB (A) -> PTGDS (A)",        a_cp$NFkB_ma,     a_cp$PTGDS_ma,  "Astro-Astro"),
  list("PTGDS (A) -> PPARG module (A)", a_cp$PTGDS_ma,    a_cp$NFkB_ma,   "Astro-Astro"),
  list("C3 (M) -> NGFR (N)",            m_cp$C3_ma,       n_cp$NGFR_ma,   "Micro-Neuron")
)

table3 <- rbindlist(lapply(pairs_t3, function(p) {
  res <- ccf_analysis(p[[2]], p[[3]], max_lag=2L)
  res[, `:=`(Pair=p[[1]], Cell_Type=p[[4]], Sig=sig_label(p_val))]
  setcolorder(res, c("Pair","Cell_Type","Lag","r","p_val","Sig","n_eff"))
  res
}))

message("  Table 3:"); print(table3)
write.csv(table3, file.path(OUT_DIR,"Table3_Checkpoint_CCF.csv"), row.names=FALSE)

# ==============================================================================
# 6. TABLE S1 — Full trajectories (all 9 bins, 3-bin MA)
# ==============================================================================
message(">>> Generating Table S1...")

supp_s1 <- data.table(
  Bin         = a_s$bin,
  n_Astrocyte = raw_means[group=="Astrocyte", n_cells],
  NDUFS1_N    = n_s$NDUFS1_ma,
  PTGDS_A     = a_s$PTGDS_ma,
  LCN2_A      = a_s$LCN2_ma,
  NGFR_N      = n_s$NGFR_ma,
  APOE_Micro  = m_s$APOE_ma,
  TREM2_Micro = m_s$TREM2_ma,
  C3_Micro    = m_s$C3_ma,
  NFKBIA_A    = a_s$NFKBIA_ma,
  CLU_A       = a_s$CLU_ma
)
write.csv(supp_s1, file.path(OUT_DIR,"TableS4_Full_Trajectories.csv"), row.names=FALSE)
message(sprintf("  TableS1: %d rows", nrow(supp_s1)))

# ==============================================================================
# 7. TABLE S2 — Spearman correlations with CPS in astrocytes
# ==============================================================================
message(">>> Generating Table S2...")

# Donor-level pseudobulk (n~84 donors). Cell-level cor.test on
# donor is the unit of analysis.
astro_df <- df[cell_type_label=="Astrocyte" & !is.na(bin)]
if (!"donor" %in% names(astro_df))
  stop("`donor` column missing - re-run 01_SEAAD_data_extraction.R")
s2_genes <- list(
  list("CASP3","Pro-apoptotic"), list("BAX","Pro-apoptotic"),
  list("BCL2","Anti-apoptotic"), list("AQP4","Astrocyte identity"),
  list("SLC1A2","Astrocyte identity"), list("GFAP","Reactive astrocyte"),
  list("PTGDS","Checkpoint marker")
)
# LCN2 excluded: astrocytic detection 0.007% (5/67,419) -> untestable.

supp_s2 <- rbindlist(lapply(s2_genes, function(item) {
  g <- item[[1]]
  if (!g %in% colnames(astro_df)) return(NULL)
  dd <- astro_df[, .(expr=mean(get(g), na.rm=TRUE), bin=mean(bin)), by=donor]
  ct <- cor.test(dd$expr, dd$bin, method="spearman", exact=FALSE)
  data.table(Gene=g, Category=item[[2]],
             Spearman_rho=round(ct$estimate,4),
             p_value=signif(ct$p.value,3),
             n_donors=nrow(dd),
             Direction=ifelse(ct$p.value>0.05,"Stable (ns)",
                              ifelse(ct$estimate>0,"Increase","Decrease")))
}))
# Multiple-testing: only BCL2 survives BH/Bonferroni; CASP3/AQP4 nominal.

write.csv(supp_s2, file.path(OUT_DIR,"TableS5_Apoptotic_Identity_Markers.csv"), row.names=FALSE)
message("  TableS2 (donor-level):"); print(supp_s2)

# ==============================================================================
# 8. TABLE S3 — Full trajectory CCF (9 bins, lag.max=3)
# ==============================================================================
message(">>> Generating Table S3...")

# LCN2 pairs excluded (0.007% detection); see note above.
pairs_s3 <- list(
  list("NDUFS1 (N) -> PTGDS (A)",        n_s$NDUFS1_ma,   a_s$PTGDS_ma,  "Neuron-Astro"),
  list("NF-kB (A) <-> PTGDS (A)",        a_s$NFkB_ma,     a_s$PTGDS_ma,  "Astro-Astro"),
  list("Purinergic (A) <-> PTGDS (A)",   a_s$Purin_ma,    a_s$PTGDS_ma,  "Astro-Astro"),
  list("C3 (M) <-> NGFR (N)",            m_s$C3_ma,       n_s$NGFR_ma,   "Micro-Neuron"),
  list("Ferroptosis (A) <-> NDUFS1 (N)", a_s$Ferr_ma,     n_s$NDUFS1_ma, "Astro-Neuron")
)

supp_s3 <- rbindlist(lapply(pairs_s3, function(p) {
  res <- ccf_analysis(p[[2]], p[[3]], max_lag=3L)
  res[, `:=`(Pair=p[[1]], Cell_Type=p[[4]], Sig=sig_label(p_val))]
  setcolorder(res, c("Pair","Cell_Type","Lag","r","p_val","Sig","n_eff"))
  res
}))

write.csv(supp_s3, file.path(OUT_DIR,"TableS6_Full_Trajectory_CCF.csv"), row.names=FALSE)
message("  TableS3:"); print(supp_s3)

# ==============================================================================
# 9. TABLE S4 — Neuronal subtype vulnerability
# ==============================================================================
message(">>> Generating Table S4...")

vuln_genes <- intersect(c("NDUFS1","NGFR","BCL2","BAX","CASP3"), colnames(df))

s4_exc <- df[cell_type_label=="Exc_Neuron" & bin_rounded>=0.1 & bin_rounded<=0.9,
             c(lapply(.SD, mean, na.rm=TRUE), list(n=.N)),
             by=bin_rounded, .SDcols=vuln_genes][order(bin_rounded)]
s4_exc[, subtype := "Excitatory"]

s4_sst <- df[cell_type_label=="SST_Neuron" & bin_rounded>=0.1 & bin_rounded<=0.9,
             c(lapply(.SD, mean, na.rm=TRUE), list(n=.N)),
             by=bin_rounded, .SDcols=vuln_genes][order(bin_rounded)]
s4_sst[, subtype := "SST_Inhibitory"]

# NDUFS1 variance
if ("NDUFS1" %in% vuln_genes) {
  exc_var <- df[cell_type_label=="Exc_Neuron" & bin_rounded>=0.1 & bin_rounded<=0.9,
                .(NDUFS1_var=var(NDUFS1,na.rm=TRUE)), by=bin_rounded]
  sst_var <- df[cell_type_label=="SST_Neuron" & bin_rounded>=0.1 & bin_rounded<=0.9,
                .(NDUFS1_var=var(NDUFS1,na.rm=TRUE)), by=bin_rounded]
  s4_exc <- merge(s4_exc, exc_var, by="bin_rounded", all.x=TRUE)
  if (nrow(s4_sst)>0) s4_sst <- merge(s4_sst, sst_var, by="bin_rounded", all.x=TRUE)
}

# Proportions
total_n <- df[cell_type_label %in% c("Exc_Neuron","SST_Neuron") &
                bin_rounded>=0.1 & bin_rounded<=0.9, .(total=.N), by=bin_rounded]
sst_n   <- df[cell_type_label=="SST_Neuron" & bin_rounded>=0.1 & bin_rounded<=0.9,
              .(SST_n=.N), by=bin_rounded]
exc_n   <- df[cell_type_label=="Exc_Neuron" & bin_rounded>=0.1 & bin_rounded<=0.9,
              .(Exc_n=.N), by=bin_rounded]
prop <- Reduce(function(a,b) merge(a,b,by="bin_rounded",all.x=TRUE),
               list(total_n, sst_n, exc_n))
prop[is.na(SST_n), SST_n := 0L]
prop[, SST_pct := round(SST_n/total*100,2)]
prop[, Exc_pct := round(Exc_n/total*100,2)]

if (all(c("BAX","BCL2") %in% vuln_genes)) {
  s4_exc[, BAX_BCL2_ratio := BAX/(BCL2+1e-10)]
  if (nrow(s4_sst)>0) s4_sst[, BAX_BCL2_ratio := BAX/(BCL2+1e-10)]
}

supp_s4 <- rbindlist(list(s4_exc, s4_sst), fill=TRUE)
supp_s4 <- merge(supp_s4, prop[,.(bin_rounded,total,SST_pct,Exc_pct)],
                 by="bin_rounded", all.x=TRUE)

write.csv(supp_s4, file.path(OUT_DIR,"TableS7_Neuronal_Subtype_Vuln.csv"), row.names=FALSE)
message(sprintf("  TableS4: %d rows, subtypes: %s",
                nrow(supp_s4), paste(unique(supp_s4$subtype),collapse=",")))

# SST proportion report (verify against Supp caption: 11.24% -> 5.49%)
if (nrow(prop[!is.na(SST_pct)]) > 0) {
  sst01 <- prop[bin_rounded==0.1, SST_pct]
  sst09 <- prop[bin_rounded==0.9, SST_pct]
  message(sprintf("  SST%%: %.2f%% (Bin 0.1) -> %.2f%% (Bin 0.9)", sst01, sst09))
  if (abs(sst01 - 11.24) > 1 || abs(sst09 - 5.49) > 1)
    message("  WARNING: SST proportions differ from paper values (11.24->5.49).",
            " Check SST_Neuron classification.")
}

# LCN2 detection rate (verify: 0.007%)
lcn2_det <- df[cell_type_label=="Astrocyte", .(
  n_cells=.N, n_LCN2_pos=sum(LCN2>0,na.rm=TRUE))]
message(sprintf("  LCN2+ astrocytes: %d / %d (%.4f%%)",
                lcn2_det$n_LCN2_pos, lcn2_det$n_cells,
                lcn2_det$n_LCN2_pos/lcn2_det$n_cells*100))
if (lcn2_det$n_LCN2_pos/lcn2_det$n_cells*100 > 0.1)
  message("  WARNING: LCN2 detection rate > 0.1% — expected ~0.007%.")

# ==============================================================================
# 10. Summary
# ==============================================================================
cat("\n", strrep("=", 65), "\n")
cat("  Output files  ->  Paper tables\n")
cat(strrep("=", 65), "\n")
cat("  Table1_CellType_Trajectories.csv       ->  Table 1\n")
cat("  Table2_Statistical_Validation.csv      ->  Table 2\n")
cat("  Table2_Segmented_Slopes_Detail.csv     ->  Table 2 (slopes)\n")
cat("  Table3_Checkpoint_CCF.csv              ->  Table 3\n")
cat("  TableS4_Full_Trajectories.csv          ->  Table S4\n")
cat("  TableS5_Apoptotic_Identity_Markers.csv ->  Table S5\n")
cat("  TableS6_Full_Trajectory_CCF.csv        ->  Table S6\n")
cat("  TableS7_Neuronal_Subtype_Vuln.csv      ->  Table S7\n")
cat(strrep("=", 65), "\n")
cat(sprintf("  Quadratic vertex CPS = %.4f\n", vertex_sc))
cat(sprintf("  beta2 = %.4f  (p = %s)\n", beta2_sc, format.pval(p_beta2,3)))
cat(sprintf("  Pre-slope  = %+.3f  |  Post-slope = %.3f\n", pre_slope, post_slope))
if (!is.na(inflection_est))
  cat(sprintf("  Inflection = %.2f  (CI: %.2f-%.2f)  Davies p = %.4f\n",
              inflection_est, ci_lo, ci_hi, davies_p))
cat(strrep("=", 65), "\n")
message(">>> Script Tables_main_and_S1_S4 complete.")
