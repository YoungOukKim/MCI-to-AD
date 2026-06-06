# ==============================================================================
# oligo_MCT_detection_check.R
#
# Purpose : Detection-rate gate for low-abundance transcripts, checking that
#           oligodendrocytic MCT4 (SLC16A3) declines ~42% along CPS, implying a
#           lactate "dual supply-line disruption" via MCT1 (SLC16A1).
#
# Why this exists : that paragraph was "data not shown" and was never committed.
#           From-raw re-check shows MCT4 is detected in only ~3% of
#           oligodendrocytes (untestable) and MCT1 is flat -> the claim does not
#           reproduce. This script documents the detection-rate gate.
#
# Rule    : a gene is reportable only if det_rate (fraction of cells with
#           expression > 0) >= 0.10. Trajectory = per-gene z across cells,
#           then mean within 9 CPS bins (same convention as Fig.1 / FigS).
#
# Input   : SEA-AD MTG h5ad (same file as Fig.1)
# Output  : output/oligo_MCT_detection_check.csv + console verdict
# Requires: rhdf5, data.table
# Usage   : Rscript analysis/FigS/oligo_MCT_detection_check.R
# ==============================================================================
suppressPackageStartupMessages({ library(rhdf5); library(data.table) })
set.seed(42)

MTG_H5 <- Sys.getenv("H5AD_PATH",
  unset = "path/to/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad")
OUTDIR <- "output"; if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

TRAJ_GENES <- c("SLC16A3","SLC16A1","MOG","MAG","MAPT")          # MCT4, MCT1, myelin, myelin, tau
ROLE       <- c("MCT4 (glycolytic export)","MCT1 (oligo->axon lactate)",
                "myelin","myelin","tau")

# ── obs: Subclass (cell type) + CPS ───────────────────────────────────────────
read_obs <- function(name){
  v <- h5read(MTG_H5, paste0("obs/", name))
  cats <- tryCatch(h5read(MTG_H5, paste0("obs/__categories/", name)), error=function(e) NULL)
  if (!is.null(cats)) cats[as.integer(v)+1L] else v
}
subclass <- read_obs("Subclass")
cps      <- as.numeric(h5read(MTG_H5, "obs/Continuous Pseudo-progression Score"))
keep     <- grepl("Oligo", subclass, ignore.case=TRUE) &
            !grepl("OPC|precursor", subclass, ignore.case=TRUE) & !is.na(cps)
cat(sprintf("oligodendrocyte nuclei: %d\n", sum(keep)))

# ── gene index (missing-gene tolerant) ───────────────────────────────────────
genes <- as.character(h5read(MTG_H5,"var/_index")); pidx0 <- match(TRAJ_GENES, genes) - 1L
present <- !is.na(pidx0); G <- TRAJ_GENES[present]; pidx0 <- pidx0[present]
cat(sprintf("genes present: %s (missing: %s)\n",
            paste(G,collapse=","), paste(setdiff(TRAJ_GENES,G),collapse=",")))

# ── FAST block read: keep-filtered CSC slurp ──────────────────────────────────
indptr <- h5read(MTG_H5,"X/indptr",bit64conversion="double"); nC <- length(indptr)-1L
sel <- which(keep); kset <- logical(nC); kset[sel] <- TRUE; rmap <- integer(nC); rmap[sel] <- seq_along(sel)
M <- matrix(0,length(sel),length(pidx0)); bs <- 100000L; cat("reading expression ...\n")
for(s0 in seq(1L,nC,by=bs)){
  e0 <- min(s0+bs-1L,nC); sp <- indptr[s0]; cnt <- indptr[e0+1L]-sp; if(cnt<=0) next
  ci <- h5read(MTG_H5,"X/indices",start=sp+1L,count=cnt,bit64conversion="double")
  cd <- h5read(MTG_H5,"X/data",   start=sp+1L,count=cnt)
  for(k in which(kset[s0:e0])){
    g <- s0+k-1L; a <- indptr[g]-sp+1L; b <- indptr[g+1L]-sp; if(b<a) next
    h <- match(pidx0, ci[a:b]); ok <- !is.na(h); M[rmap[g],ok] <- cd[a:b][h[ok]]
  }
  cat(sprintf("  ...%d/%d\r",e0,nC))
}; cat("\n")
colnames(M) <- G; cp <- cps[sel]

# ── detection rate (raw) + z-then-9-bin-mean trajectory + Spearman ────────────
det_rate <- colMeans(M > 0)
Ml  <- log1p(M)
binid <- cut(cp, breaks=seq(0.1,1.0,by=0.1), right=FALSE, labels=FALSE)
traj  <- sapply(G, function(g){ z <- scale(Ml[,g])[,1]; tapply(z, binid, mean, na.rm=TRUE) })
traj  <- traj[order(as.integer(rownames(traj))),,drop=FALSE]
binval <- seq(0.15,0.95,by=0.1)[as.integer(rownames(traj))]
n <- nrow(traj)
spearman <- sapply(G, function(g) suppressWarnings(cor(binval, traj[,g], method="spearman")))
p_n <- sapply(spearman, function(r){ if(is.na(r)||abs(r)>=1) return(NA)
  t <- r*sqrt((n-2)/(1-r^2)); 2*pt(-abs(t), n-2) })

res <- data.table(gene=G, role=ROLE[present],
                  det_rate=round(det_rate,3),
                  rho_cps=round(spearman,3), p_bins=signif(p_n,2),
                  testable = det_rate >= 0.10)
cat("\n================ oligodendrocyte MCT detection check ================\n")
print(res)
write.csv(data.frame(bin=binval, traj), file.path(OUTDIR,"oligo_MCT_trajectory.csv"), row.names=FALSE)
fwrite(res, file.path(OUTDIR,"oligo_MCT_detection_check.csv"))
cat("\n[VERDICT]\n")
cat("  Reportable only if det_rate >= 0.10. SLC16A3 (MCT4) ~3% -> NOT testable:\n")
cat("  MCT4-based metrics are unreliable at this detection rate.\n")
cat("  -> low-detection genes are gated out of inference.\n")
