# Fig. 3 — ADNI CSF Concordance for Downstream Effectors

CSF-level analysis of the astrocytic PTGDS axis using longitudinal ADNI CSF
proteomics across two orthogonal platforms (Emory TMT-MS and SomaScan 7K).
The analysis evaluates which arm of the axis reaches the CSF compartment and
distinguishes the tissue-restricted upstream marker (PTGDS) from downstream
injury effectors (NEFL, LCN2). TMT-MS marker–cognition analyses are reported at
**n = 1,104** in the manuscript.

> **Note on figure numbering:** the CSF concordance analysis is reported as
> **Figure 3** in the manuscript. Earlier internal drafts placed this content
> at Fig. 5; the folder, scripts, and outputs now use the current Fig. 3
> numbering.

---

## Scripts

| File | Order | Description |
|------|-------|-------------|
| `03_ADNI_clinical_analysis.R` | 1st | Data loading, per-batch normalization, covariate adjustment, CSF marker–cognition associations against a NEFL positive control; reports associations only (no prediction claim) |
| `Fig3_CSF_concordance.R` | 2nd | Figure 3 panels: NEFL / LCN2 / PTGDS CSF associations with cognitive decline |

---

## Key results

| Metric | Value | Notes |
|--------|-------|-------|
| CSF NEFL vs cognition | r = −0.35 (TMT-MS) | Strong; injury-arm positive control (NEFL vs pathology axis r = +0.48 [TMT], clears absolute r > 0.3) |
| CSF LCN2 vs cognition | r = −0.14 (TMT-MS, p = 2×10⁻⁶) | Weak; clearest on TMT-MS where the NEFL control is satisfied |
| CSF LCN2 (SomaScan) | r = −0.10 | Direction-consistent but supportive only (SomaScan NEFL control did not clear threshold) |
| CSF PTGDS (upstream) | TMT r = −0.02, n.s. (n = 1,104); cross-platform r ≈ 0.15 | Tissue-restricted; no robust cognitive association |

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

| Platform | n | Markers | Role |
|----------|---|---------|------|
| Emory TMT-MS | 1,104 | PTGDS, LCN2, NEFL | Primary (NEFL positive control satisfied) |
| SomaScan 7K | computed at runtime by `03_ADNI_clinical_analysis.R`; not separately reported in the manuscript | PTGDS, LCN2, NEFL, Tau | Secondary / supportive |

SeqIds (SomaScan): PTGDS X10514-5, LCN2 X2836-68, NEFL X10082-251.

### Clinical variables
- Source: ADNIMERGE (MMSE, MoCA, diagnosis, age, sex, education, APOE ε4)
- CSF Aβ42, total Tau, pTau181: Roche Elecsys (UPENNBIOMK dataset)

---

## Preprocessing pipeline

```
Raw SomaScan / TMT-MS data
  → log-transform + per-batch z-normalization (batch = TMT plex / assay run)
  → covariate adjustment: age, sex, education, APOE ε4
  → marker–cognition association vs NEFL positive control
  → incremental-value check (cross-validated; ratio tested as predictor, null)
```

---

## Usage

```r
# ADNI raw data must be downloaded separately (login required); see access note.
Sys.setenv(ADNI_PATH = "path/to/ADNIMERGE2/data/")
Sys.setenv(SOMA_CSV  = "path/to/CruchagaLab_CSF_SOMAscan7k.csv")

# Run in order
source("analysis/Fig3_ADNI/03_ADNI_clinical_analysis.R")  # -> output/TableS8_Demographics.csv, output/Fig3_CSF_associations.csv
source("analysis/Fig3_ADNI/Fig3_CSF_concordance.R")       # -> output/Fig3/Fig3_CSF_concordance.png (.tif)
```

> `03_ADNI_clinical_analysis.R` reads `ADNI_PATH` and `SOMA_CSV` and writes to
> `output/`; `Fig3_CSF_concordance.R` writes to `output/Fig3/` (override either
> with the `OUT_DIR` environment variable).

## R package requirements

```r
install.packages(c("tidyverse","survival","lubridate","ggplot2"))
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
