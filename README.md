# Multi-faceted research platform for MCI-to-AD progression: integrating data from zebrafish, mouse models, and human patient cohorts

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R version](https://img.shields.io/badge/R-4.3.2-276DC3)](https://www.r-project.org/)

---

## Authors

**Youngouk Kim**¹†\*, **Woomyung Heo**¹†, **Sejin Park**¹†, and colleagues

¹ BioXP Research Institute, Donghae, Gangwon-do, Republic of Korea
† Co-first authors  · \* Corresponding author

---

## Abstract

Alzheimer’s disease shows prolonged prodromal stability before accelerating decline, yet molecular markers resolving this heterogeneity are limited. Using pseudo-progression analysis of 1.3 million SEA-AD single nuclei (84 donors), we identify a reproducible biphasic astrocytic trajectory anchored to prostaglandin D2 synthase (PTGDS), with a statistically resolved donor-level inflection (quadratic β₂ = −2.27, p = 0.006; vertex CPS 0.47). The same directional change is independently reproduced in external brain proteomics (ROSMAP and Banner; AD versus control p = 3.4 × 10⁻³), and the biphasic pattern reconciles previously conflicting CSF reports as stage-dependent. In ADNI CSF, downstream NEFL tracks cognitive decline strongly and LCN2 weakly, whereas PTGDS itself is tissue-restricted and not a stand-alone predictor. We propose, but do not establish, that post-inflection PTGDS attenuation accompanies LCN2-linked inflammation and NGFR suppression. These data position astrocytic PTGDS as a candidate stage marker, not a causal driver, of the compensatory-to-vulnerable shift in the aging brain.

---

## Repository Structure

```
MCI-to-AD/
├── README.md
│
├── analysis/
│   ├── Fig1_SEA-AD/
│   │   ├── 01_SEAAD_data_extraction.R       # h5ad -> donor-aware processed data
│   │   └── Figure1_panels.R                 # Fig. 1 (A-D): trajectories + donor-level stats
│   ├── Tables/
│   │   └── Tables_main_and_S1_S4.R          # Table 2 (donor-level); Tables S4-S7
│   ├── Fig5_ADNI/
│   │   ├── 03_ADNI_clinical_analysis.R      # CSF associations, demographics
│   │   └── Fig5_CSF_concordance.R           # NEFL/LCN2/PTGDS CSF panels (Fig. 3)
│   └── FigS/
│       └── FigureS7_ADNI_CSF_crossplatform.R   # Fig. S7 panels (between-platform r ~ 0.15)
│
├── data/ADNI/                          # de-identified summary CSVs (not raw ADNI)
│   ├── FigS6_platform_cor.csv
│   ├── FigS6_trend.csv
│   └── FigS6_tmt_ptgds.csv
│
└── figures/                            # PNG + TIFF (300/600 dpi)
    ├── Figure1_ABCD.png
    ├── Fig3_CSF_concordance.png
    └── FigureS7_CSF_crossplatform.png
```

---

## Data Availability

| Dataset | Access | Link |
| --- | --- | --- |
| SEA-AD snRNA-seq (1.3M nuclei) | Public | Allen Brain Cell Atlas (human MTG 10x SEA-AD) |
| ADNI CSF proteomics (n = 1104) | Registration required | adni.loni.usc.edu |
| Zebrafish qPCR / ELISA / Behavior | This repository | `data/zebrafish/` |
| Murine qPCR / Behavior | Available upon request | Corresponding author |

> ⚠️ ADNI and SEA-AD raw data are subject to their data-use agreements and **cannot
> be redistributed**. The CSVs under `data/ADNI/` are de-identified, group-level
> summary statistics used to render Fig. S7; they are not subject-level data.

---

## Key Findings

| Finding | Value | Method |
| --- | --- | --- |
| PTGDS inflection (vertex) | CPS 0.47 | Quadratic model, donor-level: β₂ = −2.27, p = 0.006 (n = 84) |
| Segmented breakpoint | CPS 0.23 (95% CI 0.13–0.33) | Davies' test, donor-level p = 0.032 |
| External bulk-protein replication | PTGDS AD vs control p = 3.4 × 10⁻³ | ROSMAP/Banner/Emory (Johnson et al. 2022) |
| CSF downstream arm | NEFL r = −0.35; LCN2 r = −0.14 (TMT) | per-batch z; NEFL positive control |
| CSF upstream PTGDS | tissue-restricted (between-platform r ≈ 0.15) | TMT-MS / SomaScan |
| Zebrafish avoidance deficit | −22.7% (p = 0.020) | Red-ball avoidance, 14 dpf |
| BXP-101 rescue | +18.6 pp (p = 0.016) | BXP-101 0.6 μg/ml vs MCI |
| IL-6 suppression (protein) | −75% vs MCI | ELISA, 14 dpf |

---

## Reproducibility

- R version 4.3.2
- Key packages: `segmented`, `lme4`, `ggplot2`, `pROC`, `survival`, `rhdf5`, `data.table`, `patchwork`

```r
# Fig. 1 — SEA-AD trajectories + donor-level statistics
source("analysis/Fig1_SEA-AD/Figure1_panels.R")

# Fig. 3 — ADNI CSF concordance
source("analysis/Fig5_ADNI/Fig5_CSF_concordance.R")

# Fig. S7 — cross-platform CSF
source("analysis/FigS/FigureS7_ADNI_CSF_crossplatform.R")
```

Script paths use relative `path/to/...` placeholders or `Sys.getenv()` environment
variables; set these to your local data locations before running.

---

## Figure & Table Index

| Figure | Content | Script / asset |
|---|---|---|
| Fig. 1 | SEA-AD snRNA PTGDS biphasic trajectory & inflection (A–D) | `Fig1_SEA-AD/Figure1_panels.R` → `figures/Figure1_ABCD.png` |
| Fig. 2 | Cross-species validation (zebrafish, murine) | (external) |
| Fig. 3 | ADNI CSF concordance | `Fig5_ADNI/Fig5_CSF_concordance.R` → `figures/Fig3_CSF_concordance.png` |
| Fig. 4 | Network pharmacology / docking | (external) |
| Fig. 5 | BXP-101 murine rescue | (external) |
| Fig. S7 | Cross-platform CSF validation | `FigS/FigureS7_ADNI_CSF_crossplatform.R` → `figures/FigureS7_CSF_crossplatform.png` |

| Supp. table | Content | Script → CSV |
| --- | --- | --- |
| Table S4 | Bin-resolved cell-type trajectories | `Tables_main_and_S1_S4.R` → `TableS4_Full_Trajectories.csv` |
| Table S5 | Apoptotic / identity markers (donor-level Spearman) | `Tables_main_and_S1_S4.R` → `TableS5_Apoptotic_Identity_Markers.csv` |
| Table S6 | Lagged cross-correlation (full CPS trajectory) | `Tables_main_and_S1_S4.R` → `TableS6_Full_Trajectory_CCF.csv` |
| Table S7 | Neuronal subtype vulnerability | `Tables_main_and_S1_S4.R` → `TableS7_Neuronal_Subtype_Vuln.csv` |
| Table S8 | ADNI participant demographics | `03_ADNI_clinical_analysis.R` → `TableS8_Demographics.csv` |
| Table S9 | Independent bulk-proteomic validation | external (ROSMAP / Banner / Emory) |

Tables S1–S3 (QC & Braak mapping; antibody & primer resources; mouse behavior) are
curated, not script-generated.

---

## License

Code: [MIT](https://opensource.org/licenses/MIT). Zebrafish data: CC BY 4.0.

## Citation

If you use this code or data, please cite:

```
Kim Y, Heo W, Park SJ, et al. (2025).
A biphasic astrocytic PTGDS trajectory marks a metabolic vulnerability stage in prodromal Alzheimer’s disease.
Manuscript under review.
```

## Contact

**Youngouk Kim** — BioXP Research Institute, Donghae, Gangwon-do, Republic of Korea
