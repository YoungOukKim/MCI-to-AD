# Fig. 5 — Clinical Translation: PTGDS/LCN2 Phase Boundary in ADNI CSF

Clinical validation of the astrocytic PTGDS phase boundary using longitudinal
ADNI CSF proteomics (n=735), demonstrating alignment with MMSE 26.0 and
conversion risk prediction.

## Scripts

| File | Order | Description |
|------|-------|-------------|
| `03_ADNI_clinical_analysis.R` | 1st | Data loading, ComBat batch correction, LOESS smoothing, covariate adjustment |
| `Fig5A_PTGDS_clinical_stages.R` | 2nd | Fig. 5A — Biphasic CSF PTGDS across CN/EMCI/LMCI/AD |
| `Fig5B_multimodal_trajectory.R` | 3rd | Fig. 5B–C — Segmented regression: PTGDS/LCN2 crossover and NEFL inflection |
| `Fig5C_dual_breakpoint.R` | 4th | Fig. 5D — MoCA-based concordant breakpoint |
| `Fig5D_checkpoint_integrity.R` | 5th | Fig. 5E–F — ROC analysis and Kaplan–Meier survival curve |

## Panels

| Panel | Description |
|-------|-------------|
| A | Biphasic CSF PTGDS: EMCI +23% (p=0.008), LMCI +41% (p<0.001), AD −16% vs LMCI (p=0.02) |
| B | Segmented regression: PTGDS/LCN2 crossover at MMSE 26.0 (95% CI: 24.8–27.2) |
| C | NEFL inflection at MMSE 24.7 (95% CI: 23.1–26.3; p<0.001) |
| D | Cross-platform validation — TMT-MS and SomaScan concordance; MoCA breakpoint 21.3 |
| E | ROC analysis: AUC=0.743 (95% CI: 0.698–0.788) for MCI vs AD discrimination |
| F | Kaplan–Meier: PTGDS/LCN2 index predicts conversion risk (HR 3.2, 95% CI: 2.3–4.4) |

## Data source

**ADNI (Alzheimer's Disease Neuroimaging Initiative)**
- Access: https://adni.loni.usc.edu (application required)
- Reference: Petersen et al. (2010) Neurology doi:10.1212/WNL.0b013e3181cb3e25

### Proteomics platforms

| Platform | n | Markers | Dataset |
|----------|---|---------|---------|
| TMT-MS (Emory cohort) | 415 | PTGDS absolute quantification | ADNI proteomics |
| SomaScan 7K | 735 | PTGDS, LCN2, NEFL, Tau | SeqIds below |

SomaScan SeqIds: PTGDS X10514-5 · LCN2 X2836-68 · NEFL X10082-251 · Tau X5854-60

### Clinical variables
- Source: ADNIMERGE (MMSE, MoCA, diagnosis, age, sex, education, APOE ε4)
- CSF Aβ42, total Tau, pTau181: Roche Elecsys (UPENNBIOMK dataset)

## Key results

| Metric | Value |
|--------|-------|
| PTGDS/LCN2 crossover | MMSE 26.0 (95% CI: 24.8–27.2) |
| NEFL inflection | MMSE 24.7 (95% CI: 23.1–26.3) |
| MoCA breakpoint | 21.3 |
| MMSE window | 1.3 points (26.0 → 24.7) |
| AUC (MCI vs AD) | 0.743 (95% CI: 0.698–0.788) |
| Hazard ratio | HR 3.2 (95% CI: 2.3–4.4) |
| PTGDS/LCN2 vs MMSE | r = 0.34, p < 0.001 |
| PTGDS/LCN2 vs NEFL | β = −0.37, p < 10⁻²⁸ |

## Preprocessing pipeline

```
Raw SomaScan/TMT-MS data
  → ComBat batch correction (Johnson et al. 2007)
  → LOESS smoothing
  → Segmented regression (Davies' test)
  → Covariate adjustment: age, sex, education, APOE ε4
  → ROC (DeLong 95% CI) + Cox proportional hazards
      stratified by 75th percentile risk score
```

## Usage

```r
# ADNI data must be downloaded separately (login required)
Sys.setenv(ADNI_DATA = "data/ADNI")
Sys.setenv(FIG5_OUT  = "output/Fig5")

# Run in order
source("analysis/Fig5_ADNI/03_ADNI_clinical_analysis.R")
source("analysis/Fig5_ADNI/Fig5A_PTGDS_clinical_stages.R")
source("analysis/Fig5_ADNI/Fig5B_multimodal_trajectory.R")
source("analysis/Fig5_ADNI/Fig5C_dual_breakpoint.R")
source("analysis/Fig5_ADNI/Fig5D_checkpoint_integrity.R")
```

## R package requirements

```r
install.packages(c("ggplot2","dplyr","tidyr","readr",
                   "segmented","survival","survminer",
                   "pROC","lme4","sva"))
```

## ADNI data access note

> ADNI data are available to qualified researchers at https://adni.loni.usc.edu.
> Raw data files are not included in this repository per ADNI data use agreement.
> Place downloaded files in `data/ADNI/` before running scripts.

## Citation

ADNI:
> Petersen RC, Aisen PS, Beckett LA, et al. (2010).
> Alzheimer's Disease Neuroimaging Initiative (ADNI).
> *Neurology*. doi:10.1212/WNL.0b013e3181cb3e25

SomaScan proteomics:
> Johnson ECB, Dammer EB, Duong DM, et al. (2020).
> Large-scale proteomic analysis of Alzheimer's disease brain.
> *Nature Medicine*. doi:10.1038/s41591-020-0815-6

All analyses reproducible with `set.seed(42)` in R v4.3.2.
