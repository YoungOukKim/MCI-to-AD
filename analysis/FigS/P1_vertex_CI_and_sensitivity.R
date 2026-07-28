# ==============================================================================
# P1_vertex_CI_and_sensitivity.R
#
# Purpose : Robustness of the astrocytic PTGDS inflection reported in Table 2.
#           Addresses four questions not covered by the original analysis:
#
#   Q1  Bootstrap 95% confidence interval for the quadratic vertex (CPS 0.47).
#       The original analysis reports the vertex as a point estimate only.
#   Q2  Bin 0.1 contains a single donor (CPS 0.150). Does the curvature survive
#       removal of that donor, and of every other donor in turn?
#   Q3  Does the curvature survive truncation of the lowest-CPS donors, and are
#       the two limbs of the trajectory equally supported?
#   Q4  Donor counts per CPS bin, which are markedly uneven and were not
#       tabulated in the original submission.
#
# Pre-specified decision rules (fixed before results were read):
#   - Anchor gate: the full-cohort fit must reproduce the published values
#     (84 donors, 67,419 astrocytes, beta2 = -2.27, p = 0.006, vertex = 0.47).
#     If it does not, all downstream output is advisory only.
#   - Q2 verdict: if any leave-one-donor-out refit yields p(beta2) >= 0.05, that
#     donor is named and the "biphasic" description is reconsidered.
#   - Vertex CI: donor resampling, 10,000 replicates, percentile method. The
#     proportion of resampled vertices falling outside the observed CPS range is
#     reported alongside; if that proportion is large the vertex is not
#     determined by the data and no interval should be quoted.
#   - X is already log-normalised; no further transformation is applied.
#
# Method  : Matches analysis/Tables/Tables_main_and_S1_S4.R lines 169-174.
#           Donor-level PTGDS = mean log-normalised expression across that
#           donor's astrocytes; lm(PTGDS ~ cps + I(cps^2)); vertex = -b1/(2*b2).
#           Note that in the original scripts the variable named `bin` holds the
#           continuous CPS value (01_SEAAD_data_extraction.R line 128,
#           df[, bin := cps_vec]); the rounded value is `bin_rounded`.
#
# Inputs (first available source is used):
#   [1] data/SEA-AD/P1_donor_level_PTGDS.csv   in this repository; default path,
#                                              runs in seconds, fully reproducible
#   [2] output/SEAAD_processed_data.RData      local only; cell-level data, not
#                                              redistributable under the Allen
#                                              Brain Atlas data use policy
#   [3] SEA-AD h5ad file                       local only; reads astrocyte rows
#                                              only, set P1_H5AD to enable
#
# Outputs (-> output/P1_vertex/):
#   P1_vertex_CI_and_sensitivity.txt
#   P1_leave_one_donor_out.csv
#   P1_donor_level_PTGDS.csv        written only when derived from [2] or [3]
#
# Paper : Kim Y et al. (2026). https://github.com/YoungOukKim/MCI-to-AD
# Usage : Rscript analysis/FigS/P1_vertex_CI_and_sensitivity.R
# Requirements : R >= 4.3.2 | base stats only for source [1];
#                rhdf5 additionally for source [3]
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(42)

# ------------------------------------------------------------------------------
# Paths (relative to repository root; override with environment variables)
# ------------------------------------------------------------------------------
DONOR_CSV <- Sys.getenv("P1_DONOR_CSV", unset = "data/SEA-AD/P1_donor_level_PTGDS.csv")
RDATA     <- Sys.getenv("P1_RDATA",     unset = "output/SEAAD_processed_data.RData")
H5AD      <- Sys.getenv("P1_H5AD",      unset = "")
OUT_DIR   <- Sys.getenv("P1_OUT",       unset = "output/P1_vertex")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

LOG <- file.path(OUT_DIR, "P1_vertex_CI_and_sensitivity.txt")

CPS_OBS   <- "Continuous Pseudo-progression Score"
SUBCL_OBS <- "Subclass"
GENE      <- "PTGDS"
TARGET    <- "Astrocyte"
NBOOT     <- 10000L
GAP       <- 5L      # source [3]: astrocytes separated by <= GAP cells share a read
BLOCK     <- 2000L   # source [3]: fallback block size for a full sequential scan

REF_N_DONORS <- 84L
REF_N_CELLS  <- 67419L
REF_BETA2    <- -2.27
REF_VERTEX   <- 0.47

sink(LOG, split = TRUE)
cat("P1 robustness: vertex confidence interval and single-donor sensitivity\n")
cat("Run:", format(Sys.time()), "\n\n")

DF <- NULL

# ------------------------------------------------------------------------------
# Source [1] - donor-level summary deposited in this repository
# ------------------------------------------------------------------------------
cat("[1] Donor-level summary CSV\n")
if (file.exists(DONOR_CSV)) {
  d0 <- read.csv(DONOR_CSV)
  need <- c("donor", "cps", "ptgds")
  if (all(need %in% names(d0))) {
    if (!"ncell" %in% names(d0)) d0$ncell <- NA_real_
    DF <- d0[, c("donor", "cps", "ptgds", "ncell")]
    cat("    Found:", DONOR_CSV, "-", nrow(DF), "donors\n")
    cat("    Source [1] used; no further input required.\n")
  } else {
    cat("    Missing columns:", paste(setdiff(need, names(d0)), collapse = ", "), "\n")
  }
} else {
  cat("    Not found:", DONOR_CSV, "\n")
}

# ------------------------------------------------------------------------------
# Source [2] - cell-level RData produced by 01_SEAAD_data_extraction.R
# ------------------------------------------------------------------------------
if (is.null(DF)) {
  cat("\n[2] Cell-level RData\n")
  if (file.exists(RDATA)) {
    ev <- new.env(); load(RDATA, envir = ev)
    df <- get(ls(ev)[1], envir = ev)
    need <- c(GENE, "donor", "bin", "cell_type_label")
    if (all(need %in% names(df))) {
      dd <- as.data.frame(df)[, need]
      dd <- dd[dd$cell_type_label == TARGET & !is.na(dd$bin) & !is.na(dd$donor), ]
      cat("    Found:", RDATA, "-", format(nrow(dd), big.mark = ","),
          TARGET, "cells,", length(unique(dd$donor)), "donors\n")
      pb <- tapply(dd[[GENE]], dd$donor, mean)
      cp <- tapply(dd$bin,     dd$donor, function(z) z[1])
      nc <- tapply(rep(1L, nrow(dd)), dd$donor, sum)
      DF <- data.frame(donor = names(pb), cps = as.numeric(cp),
                       ptgds = as.numeric(pb), ncell = as.numeric(nc))
    } else {
      cat("    Missing columns:", paste(setdiff(need, names(df)), collapse = ", "), "\n")
    }
  } else {
    cat("    Not found:", RDATA, "\n")
    cat("    This file holds cell-level data and is not redistributed under the\n")
    cat("    Allen Brain Atlas data use policy. Regenerate it with\n")
    cat("    analysis/Fig1_SEA-AD/01_SEAAD_data_extraction.R.\n")
  }
}

# ------------------------------------------------------------------------------
# Source [3] - direct read from the SEA-AD h5ad, astrocyte rows only
# ------------------------------------------------------------------------------
if (is.null(DF)) {
  cat("\n[3] Direct h5ad read\n")
  if (!nzchar(H5AD) || !file.exists(H5AD)) {
    cat("    No h5ad path set. Set P1_H5AD to enable this source.\n")
    cat("    No usable input found; stopping.\n")
    sink(); quit(save = "no", status = 1)
  }
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    cat("    Package rhdf5 is required for this source; stopping.\n")
    sink(); quit(save = "no", status = 1)
  }
  suppressMessages(library(rhdf5))

  read_obs_cat <- function(key) {
    p  <- paste0("obs/", key)
    cp <- paste0("obs/__categories/", key)
    tryCatch({
      codes <- as.integer(h5read(H5AD, p))
      cats  <- as.character(h5read(H5AD, cp))
      cats[codes + 1L]
    }, error = function(e) as.character(h5read(H5AD, p)))
  }

  subcl <- read_obs_cat(SUBCL_OBS)
  obs_names <- { ll <- h5ls(H5AD); ll$name[ll$group == "/obs"] }
  donor_key <- obs_names[grepl("donor", obs_names, ignore.case = TRUE)][1]
  cat("    donor field =", ifelse(is.na(donor_key), "<not found>", donor_key), "\n")
  donor <- read_obs_cat(donor_key)
  cps   <- as.numeric(h5read(H5AD, paste0("obs/", CPS_OBS)))
  genes <- as.character(h5read(H5AD, "var/_index"))
  ncell <- length(subcl); ngene <- length(genes)

  cat("\n    Subclass labels, top 10. Do not hard-code these without checking.\n")
  tb <- sort(table(subcl), decreasing = TRUE)
  for (i in seq_len(min(10, length(tb)))) {
    cat(sprintf("      %-22s %d\n", names(tb)[i], tb[i]))
  }
  if (!(TARGET %in% names(tb))) {
    cat("\n    Label '", TARGET, "' not present; stopping.\n", sep = "")
    sink(); quit(save = "no", status = 1)
  }

  indptr <- as.numeric(h5read(H5AD, "X/indptr", bit64conversion = "double"))
  nnz <- indptr[length(indptr)]
  cat(sprintf("\n    nnz %s | length(indptr)-1 %s | cells %s | density %.4f\n",
              format(nnz, big.mark = ","), format(length(indptr) - 1L, big.mark = ","),
              format(ncell, big.mark = ","), nnz / (as.numeric(ncell) * ngene)))
  if ((length(indptr) - 1L) != ncell) {
    cat("    Matrix is not CSR; this script handles CSR only. Stopping.\n")
    sink(); quit(save = "no", status = 1)
  }
  cat("    CSR confirmed.\n")

  gi <- match(GENE, genes)
  if (is.na(gi)) {
    cat("    ", GENE, " absent from var; stopping.\n", sep = "")
    sink(); quit(save = "no", status = 1)
  }
  cat(sprintf("    %s var row index %d\n", GENE, gi))

  keep <- !is.na(cps) & !is.na(subcl) & (subcl == TARGET)
  keep[is.na(keep)] <- FALSE
  ac <- which(keep)
  cat(sprintf("\n    %s cells, %d donors\n",
              format(length(ac), big.mark = ","), length(unique(donor[ac]))))
  cat(sprintf("    Published counts %s / %d -> %s\n",
              format(REF_N_CELLS, big.mark = ","), REF_N_DONORS,
              ifelse(length(ac) == REF_N_CELLS &&
                     length(unique(donor[ac])) == REF_N_DONORS, "match", "MISMATCH")))

  # Only one gene is required, so read the astrocyte rows rather than the whole
  # matrix. Astrocytes lying close together are merged into one contiguous read.
  cutp <- c(0L, which(diff(ac) > GAP), length(ac))
  ngrp <- length(cutp) - 1L
  span <- sum(vapply(seq_len(ngrp),
                     function(k) ac[cutp[k + 1]] - ac[cutp[k] + 1L] + 1L, numeric(1)))
  cat(sprintf("\n    %s read spans covering %s cells (%.1f%% of the matrix, about %.1fx less I/O)\n",
              format(ngrp, big.mark = ","), format(span, big.mark = ","),
              100 * span / ncell, ncell / span))

  expr <- numeric(ncell)
  fid <- H5Fopen(H5AD, flags = "H5F_ACC_RDONLY")
  t0 <- Sys.time()
  res <- try({
    for (k in seq_len(ngrp)) {
      i1 <- ac[cutp[k] + 1L]; i2 <- ac[cutp[k + 1]]
      a0 <- indptr[i1]; a1 <- indptr[i2 + 1L]
      cnt <- a1 - a0
      if (cnt > 0) {
        ci  <- as.integer(h5read(fid, "X/indices", start = a0 + 1, count = cnt)) + 1L
        sel <- which(ci == gi)
        if (length(sel)) {
          xv <- as.numeric(h5read(fid, "X/data", start = a0 + 1, count = cnt))
          cell_of <- rep(i1:i2, diff(indptr[i1:(i2 + 1L)]))
          expr[cell_of[sel]] <- xv[sel]
        }
      }
      if (k %% 5000L == 0L) {
        el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
        cat(sprintf("    ... span %s / %s (%.1f min elapsed, about %.0f min total)\n",
                    format(k, big.mark = ","), format(ngrp, big.mark = ","),
                    el, el * ngrp / k))
      }
    }
    TRUE
  }, silent = TRUE)
  H5Fclose(fid)

  if (inherits(res, "try-error")) {
    cat("    Span-based read failed:", conditionMessage(attr(res, "condition")), "\n")
    cat("    Falling back to a full sequential scan, which is substantially slower.\n")
    expr <- numeric(ncell); nb <- 0L; t0 <- Sys.time()
    for (b in seq.int(1L, ncell, by = BLOCK)) {
      e2 <- min(b + BLOCK - 1L, ncell)
      if (!any(keep[b:e2])) next
      a0 <- indptr[b]; a1 <- indptr[e2 + 1L]
      if (a1 <= a0) next
      ci  <- as.integer(h5read(H5AD, "X/indices", start = a0 + 1, count = a1 - a0)) + 1L
      sel <- which(ci == gi)
      if (length(sel)) {
        xv <- as.numeric(h5read(H5AD, "X/data", start = a0 + 1, count = a1 - a0))
        cell_of <- rep(b:e2, diff(indptr[b:(e2 + 1L)]))
        expr[cell_of[sel]] <- xv[sel]
      }
      nb <- nb + 1L
      if (nb %% 100L == 0L) {
        el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
        cat(sprintf("    ... block %d / about %d (%.1f min elapsed, about %.0f min total)\n",
                    nb, ceiling(ncell / BLOCK), el, el * ceiling(ncell / BLOCK) / nb))
      }
    }
  }

  v  <- expr[ac]
  nz <- v[v > 0]
  if (!length(nz)) {
    cat("    All values are zero; check the reader or the gene index. Stopping.\n")
    sink(); quit(save = "no", status = 1)
  }
  cat(sprintf("    Detection rate %.4f | non-zero range %.3f to %.3f, no log1p applied\n",
              mean(v > 0), min(nz), max(nz)))

  d1 <- donor[ac]; c1 <- cps[ac]
  pb <- tapply(v,  d1, mean)
  cp <- tapply(c1, d1, function(z) z[1])
  nc <- tapply(rep(1L, length(d1)), d1, sum)
  DF <- data.frame(donor = names(pb), cps = as.numeric(cp),
                   ptgds = as.numeric(pb), ncell = as.numeric(nc))
}

# ------------------------------------------------------------------------------
# Donor-level table
# ------------------------------------------------------------------------------
DF <- DF[!is.na(DF$cps) & !is.na(DF$ptgds), ]
DF <- DF[order(DF$cps), ]
cat(sprintf("\n[4] Donor-level table: n = %d | CPS %.3f to %.3f\n",
            nrow(DF), min(DF$cps), max(DF$cps)))

cat("\n[5] Q4. Donors per CPS bin, not tabulated in the original submission\n")
bn <- round(DF$cps, 1)
for (k in sort(unique(bn))) {
  idx <- which(bn == k)
  cells <- if (all(is.na(DF$ncell[idx]))) "NA" else format(sum(DF$ncell[idx]), big.mark = ",")
  cat(sprintf("      bin %.1f : %2d donors | %8s cells | PTGDS %.4f\n",
              k, length(idx), cells, mean(DF$ptgds[idx])))
}

fit_v <- function(d) {
  m  <- lm(ptgds ~ cps + I(cps^2), data = d)
  s  <- summary(m)$coefficients
  b1 <- coef(m)[["cps"]]; b2 <- coef(m)[["I(cps^2)"]]
  list(b1 = b1, b2 = b2, p2 = s["I(cps^2)", "Pr(>|t|)"], vertex = -b1 / (2 * b2))
}

# ------------------------------------------------------------------------------
# Anchor gate
# ------------------------------------------------------------------------------
cat("\n[6] Anchor gate: full-cohort quadratic fit\n")
f0 <- fit_v(DF)
cat(sprintf("    beta1 %+.4f | beta2 %+.4f | p(beta2) %.5f | vertex %.4f\n",
            f0$b1, f0$b2, f0$p2, f0$vertex))
cat(sprintf("    Published: beta2 %.2f, p 0.006, vertex %.2f\n", REF_BETA2, REF_VERTEX))
gate <- (f0$b2 < 0) && (abs(f0$vertex - REF_VERTEX) < 0.08) && (f0$p2 < 0.05)
cat(sprintf("    -> %s\n", ifelse(gate,
      "Reproduced. Proceeding.",
      "NOT reproduced. Treat all output below as advisory only.")))

# ------------------------------------------------------------------------------
# Q1. Vertex bootstrap confidence interval
# ------------------------------------------------------------------------------
cat(sprintf("\n[7] Q1. Vertex bootstrap 95%% CI, %s donor resamples\n",
            format(NBOOT, big.mark = ",")))
vb <- rep(NA_real_, NBOOT); b2b <- rep(NA_real_, NBOOT)
for (i in seq_len(NBOOT)) {
  s <- DF[sample.int(nrow(DF), replace = TRUE), ]
  f <- try(fit_v(s), silent = TRUE)
  if (!inherits(f, "try-error")) { vb[i] <- f$vertex; b2b[i] <- f$b2 }
}
vb2 <- vb[is.finite(vb)]
inr <- vb2[vb2 >= min(DF$cps) & vb2 <= max(DF$cps)]
cat(sprintf("    Valid %s | inside the observed CPS range %s (%.1f%%)\n",
            format(length(vb2), big.mark = ","), format(length(inr), big.mark = ","),
            100 * length(inr) / max(1, length(vb2))))
if (length(inr) < 100) {
  cat("    Too few resampled vertices fall inside the observed range. The vertex is\n")
  cat("    not determined by the data and no confidence interval should be quoted.\n")
} else {
  ci <- quantile(inr, c(0.025, 0.975))
  cat(sprintf("    Vertex 95%% CI = [%.3f, %.3f], point estimate %.3f\n",
              ci[1], ci[2], f0$vertex))
}
cat(sprintf("    Proportion of resamples with beta2 < 0 = %.1f%%\n",
            100 * mean(b2b[is.finite(b2b)] < 0)))

# ------------------------------------------------------------------------------
# Q2. Leave-one-donor-out
# ------------------------------------------------------------------------------
cat("\n[8] Q2. Leave-one-donor-out refits\n")
lo <- do.call(rbind, lapply(seq_len(nrow(DF)), function(i) {
  f <- fit_v(DF[-i, ])
  data.frame(donor = DF$donor[i], cps = DF$cps[i],
             b2 = f$b2, p2 = f$p2, vertex = f$vertex)
}))
cat(sprintf("    beta2 %+.3f to %+.3f | p %.4f to %.4f | vertex %.3f to %.3f\n",
            min(lo$b2), max(lo$b2), min(lo$p2), max(lo$p2),
            min(lo$vertex), max(lo$vertex)))
bad <- lo[lo$p2 >= 0.05, ]
if (nrow(bad) > 0) {
  cat("    Curvature loses significance when these donors are removed:\n")
  for (i in seq_len(nrow(bad))) {
    cat(sprintf("       donor %s (CPS %.3f) -> beta2 %+.3f, p %.4f\n",
                bad$donor[i], bad$cps[i], bad$b2[i], bad$p2[i]))
  }
  cat("    The biphasic description depends on these donors; reconsider it.\n")
} else {
  cat("    Curvature remains significant in every refit, so it does not depend on\n")
  cat("    any single donor, including the single donor occupying bin 0.1.\n")
}
cat("\n    Five most influential donors, largest p for beta2\n")
top <- lo[order(-lo$p2), ][seq_len(min(5, nrow(lo))), ]
for (i in seq_len(nrow(top))) {
  cat(sprintf("      donor %-10s CPS %.3f -> beta2 %+.3f, p %.4f, vertex %.3f\n",
              top$donor[i], top$cps[i], top$b2[i], top$p2[i], top$vertex[i]))
}

# ------------------------------------------------------------------------------
# Q3. Low-CPS truncation and limb slopes
# ------------------------------------------------------------------------------
cat("\n[9] Q3. Low-CPS truncation\n")
for (cut in c(0.00, 0.16, 0.20, 0.25)) {
  sub <- DF[DF$cps >= cut, ]
  if (nrow(sub) < 20) {
    cat(sprintf("    CPS >= %.2f : n = %d, too few donors\n", cut, nrow(sub)))
  } else {
    f <- fit_v(sub)
    cat(sprintf("    CPS >= %.2f : n = %2d | beta2 %+.3f (p %.4f) | vertex %.3f\n",
                cut, nrow(sub), f$b2, f$p2, f$vertex))
  }
}

cat("\n[10] Q3b. Separate linear slopes either side of the vertex\n")
seg_slope <- function(lo_bound, hi_bound, label) {
  s <- DF[DF$cps > lo_bound & DF$cps <= hi_bound, ]
  if (nrow(s) < 5) {
    cat(sprintf("    %-38s n = %d, too few donors\n", label, nrow(s)))
  } else {
    m <- summary(lm(ptgds ~ cps, data = s))$coefficients
    cat(sprintf("    %-38s n = %2d | slope %+.4f | p %.4f\n",
                label, nrow(s), m["cps", "Estimate"], m["cps", "Pr(>|t|)"]))
  }
}
seg_slope(0.00, REF_VERTEX, "Ascending limb, all donors")
seg_slope(0.16, REF_VERTEX, "Ascending limb, lowest donor removed")
seg_slope(0.20, REF_VERTEX, "Ascending limb, CPS 0.20 to 0.47")
seg_slope(REF_VERTEX, 1.00, "Descending limb, reference")

# ------------------------------------------------------------------------------
# Model comparison
# ------------------------------------------------------------------------------
cat("\n[11] Model comparison by AIC\n")
nn     <- nrow(DF)
rss_of <- function(m) sum(residuals(m)^2)
aic_of <- function(rss, k) nn * log(rss / nn) + 2 * (k + 1)
m_lin  <- lm(ptgds ~ cps, data = DF)
m_qua  <- lm(ptgds ~ cps + I(cps^2), data = DF)
best <- NULL
for (kk in seq(0.20, 0.75, by = 0.005)) {
  mh <- lm(ptgds ~ pmax(0, cps - kk), data = DF)
  r  <- rss_of(mh)
  if (is.null(best) || r < best$rss) {
    best <- list(k = kk, rss = r, slope = unname(coef(mh)[2]))
  }
}
a_lin <- aic_of(rss_of(m_lin), 2)
a_qua <- aic_of(rss_of(m_qua), 3)
a_hin <- aic_of(best$rss, 3)
cat(sprintf("    Linear                  AIC %8.2f\n", a_lin))
cat(sprintf("    Quadratic, published    AIC %8.2f\n", a_qua))
cat(sprintf("    Hinge, plateau-decline  AIC %8.2f | breakpoint %.3f, slope %+.3f\n",
            a_hin, best$k, best$slope))
cat(sprintf("    Quadratic minus hinge = %+.2f; negative favours the quadratic\n",
            a_qua - a_hin))

# ------------------------------------------------------------------------------
# Write outputs
# ------------------------------------------------------------------------------
if (!file.exists(DONOR_CSV)) {
  write.csv(DF, file.path(OUT_DIR, "P1_donor_level_PTGDS.csv"), row.names = FALSE)
}
write.csv(lo, file.path(OUT_DIR, "P1_leave_one_donor_out.csv"), row.names = FALSE)

cat("\n[12] Outputs written to", OUT_DIR, "\n")
cat("\n[13] Notes\n")
cat("    - The vertex is the ratio -b1/(2*b2) and becomes unstable as b2 approaches\n")
cat("      zero. Read the interval together with the proportion of resampled\n")
cat("      vertices falling inside the observed range, reported in section [7].\n")
cat("    - All inference is at donor level; no cell-level testing is performed.\n")
cat("    - X is treated as already log-normalised and is not transformed further.\n")
sink()
