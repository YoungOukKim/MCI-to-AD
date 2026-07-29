#!/usr/bin/env Rscript
# P1_age_cps.R -- Donor-level relationship between chronological age and CPS in SEA-AD MTG.
#
# Reads ONLY the /obs metadata from the .h5ad (not the expression matrix), so it is fast
# and low-memory even for ~1.38M nuclei. Produces the real donor-level Spearman rho between
# age-at-death and CPS. Nothing is fabricated.
#
# REQUIREMENTS:  install.packages("hdf5r")
# RUN:           Rscript analysis/FigS/P1_age_cps.R      (run from the repository root)
#                Override the defaults with the H5AD_PATH / OUTDIR environment variables,
#                or pass the h5ad path as the first argument.
# OUTPUT (OUTDIR): age_cps_donor.csv, age_cps_scatter.png

# ---- paths (relative to the repository root; override via H5AD_PATH / OUTDIR) ----
# The SEA-AD h5ad is not redistributed here. Obtain it from the Allen Brain Map portal
# (see data/NOTE_restricted_data.md) and place it under data/SEA-AD/, or point
# H5AD_PATH at wherever you keep it.
MTG_H5 <- Sys.getenv("H5AD_PATH",
                     unset = file.path("data", "SEA-AD",
                                       "SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad"))
OUTDIR <- Sys.getenv("OUTDIR", unset = file.path("output", "FigS8"))

suppressWarnings(suppressMessages({
  if (!requireNamespace("hdf5r", quietly = TRUE))
    stop("Package 'hdf5r' is required. Install with: install.packages('hdf5r')")
  library(hdf5r)
}))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && nzchar(args[[1]])) MTG_H5 <- args[[1]]

if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
csv_out <- file.path(OUTDIR, "age_cps_donor.csv")
png_out <- file.path(OUTDIR, "age_cps_scatter.png")

cat(sprintf("[load]   %s\n", MTG_H5))
cat(sprintf("[outdir] %s\n", OUTDIR))
if (!file.exists(MTG_H5)) stop(sprintf("h5ad not found: %s", MTG_H5))

f <- H5File$new(MTG_H5, mode = "r")
close_f <- function() try(f$close_all(), silent = TRUE)

obs <- f[["obs"]]

# --- robust listing of obs columns: try AnnData 'column-order' attr, then names(), then ls() ---
get_members <- function(grp) {
  m <- tryCatch(as.character(hdf5r::h5attributes(grp)[["column-order"]]),
                error = function(e) NULL)
  if (is.null(m) || !length(m))
    m <- tryCatch(names(grp), error = function(e) NULL)
  if (is.null(m) || !length(m))
    m <- tryCatch(grp$ls()$name, error = function(e) NULL)
  if (is.null(m)) stop("could not list /obs columns")
  setdiff(m, c("_index", "__categories"))
}
members <- get_members(obs)
cat(sprintf("[obs] %d columns; first few: %s\n",
            length(members), paste(utils::head(members, 8), collapse = ", ")))

# --- read one obs column, decoding AnnData categoricals (group: categories + codes) ---
read_obs_col <- function(obs, name) {
  item <- obs[[name]]
  if (inherits(item, "H5Group")) {                 # categorical
    cats  <- item[["categories"]]$read()
    codes <- item[["codes"]]$read()
    out <- rep(NA_character_, length(codes))
    ok  <- codes >= 0
    out[ok] <- cats[codes[ok] + 1L]                # codes are 0-based; -1 = NA
    try(item$close(), silent = TRUE)
    return(out)
  } else {                                         # plain dataset (numeric / string)
    v <- item$read(); try(item$close(), silent = TRUE); return(v)
  }
}

find_col <- function(cols, candidates) {
  for (cand in candidates) { hit <- cols[tolower(cols) == tolower(cand)]; if (length(hit)) return(hit[1]) }
  for (cand in candidates) { hit <- cols[grepl(tolower(cand), tolower(cols), fixed = TRUE)]; if (length(hit)) return(hit[1]) }
  NA_character_
}

age_col   <- find_col(members, c("Age at Death","Age_at_Death","age_at_death","AgeDeath","Age","age"))
cps_col   <- find_col(members, c("Continuous Pseudo-progression Score","Continuous Pseudo progression Score",
                                 "pseudo-progression","pseudoprogression","CPS","cps","pseudo_progression_score"))
donor_col <- find_col(members, c("Donor ID","donor_id","DonorID","donor","external_donor_name",
                                 "specimen","individualID","Donor"))
cat(sprintf("[cols] age=%s  cps=%s  donor=%s\n",
            ifelse(is.na(age_col),"<none>",age_col),
            ifelse(is.na(cps_col),"<none>",cps_col),
            ifelse(is.na(donor_col),"<none>",donor_col)))

if (is.na(age_col) || is.na(cps_col)) {
  cat("\n[!] Could not auto-detect age/CPS columns. Available columns:\n")
  for (c in members) cat("    ", c, "\n")
  cat("\nEdit the candidate vectors near the top, then re-run.\n")
  close_f(); quit(status = 0)
}

age_raw <- read_obs_col(obs, age_col)
cps_raw <- as.numeric(read_obs_col(obs, cps_col))
donor   <- if (!is.na(donor_col)) read_obs_col(obs, donor_col) else seq_along(cps_raw)
if (is.na(donor_col))
  cat("[!] No donor column found -- falling back to per-nucleus (NOT donor-level).\n")

extract_num <- function(x) vapply(as.character(x), function(s) {
  mm <- regmatches(s, regexpr("[0-9]+(\\.[0-9]+)?", s))
  if (length(mm) == 0) NA_real_ else as.numeric(mm)
}, numeric(1), USE.NAMES = FALSE)
censored <- grepl("[+>]", as.character(age_raw))
age      <- extract_num(age_raw)

df <- data.frame(donor = as.character(donor), age = age, cps = cps_raw,
                 censored = censored, stringsAsFactors = FALSE)
df <- df[is.finite(df$age) & is.finite(df$cps), ]

agg_age  <- tapply(df$age,  df$donor, median, na.rm = TRUE)
agg_cps  <- tapply(df$cps,  df$donor, median, na.rm = TRUE)
agg_cens <- tapply(df$censored, df$donor, mean, na.rm = TRUE)
g <- data.frame(donor = names(agg_age), age = as.numeric(agg_age),
                cps = as.numeric(agg_cps), censored_frac = as.numeric(agg_cens))
g <- g[is.finite(g$age) & is.finite(g$cps), ]
n <- nrow(g)
cens_donor <- mean(g$censored_frac > 0)

sp <- suppressWarnings(cor.test(g$age, g$cps, method = "spearman"))
pe <- cor.test(g$age, g$cps, method = "pearson")
rho <- unname(sp$estimate); p_s <- sp$p.value
r   <- unname(pe$estimate); p_p <- pe$p.value

write.csv(g, csv_out, row.names = FALSE)

cat("\n========== DONOR-LEVEL RESULT ==========\n")
cat(sprintf("n donors            : %d\n", n))
cat(sprintf("age censored (>=1)  : %.1f%% of donors had a '90+'-type capped age\n", 100 * cens_donor))
cat(sprintf("Spearman rho        : %.3f   (p = %.3g)\n", rho, p_s))
cat(sprintf("Pearson  r          : %.3f   (p = %.3g)\n", r, p_p))
cat("========================================\n")
cat("\nSuggested supplementary sentence (choose weak/moderate per the rho printed above):\n")
cat(sprintf('  "Across SEA-AD donors (n = %d), donor age at death correlated only [weakly/moderately]\n   with CPS (Spearman rho = %.2f, p = %.2g), confirming that CPS indexes neuropathological\n   progression rather than chronological age."\n', n, rho, p_s))

ok <- tryCatch({
  png(png_out, width = 1100, height = 1050, res = 260)
  par(mar = c(4.2, 4.4, 2.2, 1))
  plot(g$age, g$cps, pch = 21, bg = "#4F81BD", col = "#333333", cex = 1.0, lwd = 0.5,
       xlab = "Donor age at death (years)", ylab = "CPS (donor-level)",
       main = "Age vs CPS (SEA-AD donors)")
  lab <- sprintf("Spearman rho = %.2f\np = %.2g\nn = %d%s",
                 rho, p_s, n, ifelse(cens_donor > 0, "\n(age censored at 90+)", ""))
  usr <- par("usr")
  text(usr[1] + 0.03 * (usr[2] - usr[1]), usr[4] - 0.03 * (usr[4] - usr[3]),
       lab, adj = c(0, 1), cex = 0.8)
  dev.off(); TRUE
}, error = function(e) { message("[plot skipped] ", conditionMessage(e)); FALSE })

cat(sprintf("\n[saved] %s\n", csv_out))
if (isTRUE(ok)) cat(sprintf("[saved] %s\n", png_out))
close_f()
