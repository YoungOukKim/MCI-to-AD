# Fig. 3 — ADNI CSF Concordance for Downstream Effectors

CSF-level analysis of the astrocytic PTGDS axis using longitudinal ADNI
proteomics (n = 735) across two orthogonal platforms (TMT-MS, SomaScan).
The analysis evaluates which arm of the axis reaches the CSF compartment and
distinguishes the tissue-restricted upstream marker (PTGDS) from downstream
injury effectors (NEFL, LCN2).

> **Note on figure numbering:** the CSF concordance analysis is reported as
> **Figure 3** in the manuscript. Earlier internal drafts placed this content
> at Fig. 5; the script and outputs here use the current Fig. 3 numbering.

---

## Scripts

| File | Order | Description |
|------|-------|-------------|
| `03_ADNI_clinical_analysis.R` | 1st | Data loading, per-batch normalization, covariate adjustment, CSF marker–cognition associations against a NEFL positive control |
| `Fig5_CSF_concordance.R` | 2nd | Figure 3 panels: NEFL / LCN2 / PTGDS CSF associations with cognitive decline |

---

## Key results

| Metric | Value | Notes |
|--------|-------|-------|
| CSF NEFL vs cognitive decline | r = −0.35 | Strong; positive control for the injury arm |
| CSF LCN2 vs cognitive decline | r = −0.14 (TMT-MS) | Weak; clearest on TMT-MS |
| CSF PTGDS (upstream) | between-platform r ≈ 0.15 | Tissue-restricted; no robust cognitive association |

**Interpretation.** Upstream PTGDS is tissue-restricted and does not behave as a
standalone fluid biomarker; the CSF data instead provide mechanistic
concordance for the **downstream** effectors (NEFL, LCN2). The PTGDS/LCN2 ratio
does not act as an independent predictor in cross-validated analysis (ratio-only
CV-AUC ≈ chance), and accordingly **no ROC/hazard-ratio prediction claim is made
from the CSF data**.

---

## Data source

**ADNI (Alzheimer's Disease Neuroimaging Initiative)**
- Access: https://adni.loni.usc.edu (application + data use agreement required)
- Reference: Petersen et al. (2010) *Neurology* doi:10.1212/WNL.0b013e3181cb3e25

### Proteomics platforms

| Platform | n | Markers | Dataset |
|----------|---|---------|---------|
| TMT-MS (Emory cohort) | 415 | PTGDS, NEFL, LCN2 | ADNI proteomics |
| SomaScan 7K | 735 | PTGDS, LCN2, NEFL, Tau | ADNI proteomics |

### Clinical variables
- Source: ADNIMERGE (MMSE, MoCA, diagnosis, age, sex, education, APOE ε4)
- CSF Aβ42, total Tau, pTau181: Roche Elecsys (UPENNBIOMK dataset)

---

## Preprocessing pipeline

```
Raw SomaScan / TMT-MS data
  → per-batch normalization (z-scoring within platform/batch)
  → covariate adjustment: age, sex, education, APOE ε4
  → marker–cognition association vs NEFL positive control
  → incremental-value check (cross-validated; ratio tested as predictor)
```

---

## Usage

```r
# ADNI data must be downloaded separately (login required)
Sys.setenv(ADNI_DATA = "data/ADNI")
Sys.setenv(FIG3_OUT  = "output/Fig3")

# Run in order
source("analysis/Fig5_ADNI/03_ADNI_clinical_analysis.R")
source("analysis/Fig5_ADNI/Fig5_CSF_concordance.R")
```

## R package requirements

```r
install.packages(c("ggplot2","dplyr","tidyr","readr"))
```

## ADNI data access note

> ADNI data are available to qualified researchers at https://adni.loni.usc.edu.
> Raw data files are not included in this repository per the ADNI data use
> agreement. The CSVs under `data/ADNI/` are de-identified, group-level summary
> statistics used to render the figure; they are not subject-level data.

## Citation

> Petersen RC, Aisen PS, Beckett LA, et al. (2010).
> Alzheimer's Disease Neuroimaging Initiative (ADNI).
> *Neurology*. doi:10.1212/WNL.0b013e3181cb3e25

All analyses use `set.seed(42)` in R v4.3.2.
