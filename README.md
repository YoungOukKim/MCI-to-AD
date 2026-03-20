# Multi-faceted research platform for MCI-to-AD progression: integrating deep learning models with data from zebrafish, mouse models, and human patient cohorts

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![bioRxiv](https://img.shields.io/badge/bioRxiv-preprint-blue)](https://www.biorxiv.org/[TO_BE_INSERTED])
[![R version](https://img.shields.io/badge/R-4.3.2-276DC3)](https://www.r-project.org/)

---

## Authors

**Youngouk Kim**¹†*, **Woomyung Heo**¹†, **Sejin Park**¹†, and colleagues

¹ BioXP Research Institute, Donghae, Gangwon-do, Republic of Korea  
† Co-first authors  
\* Corresponding author: [email — to be inserted]

---

## Abstract

We identify a lipid-metabolic exhaustion event that defines a discrete phase boundary
separating compensatory resilience from inflammatory vulnerability during the
MCI-to-AD transition. Integrating SEA-AD single-cell pseudo-progression
(1.3M nuclei, 84 donors), a reversible zebrafish MCI model, murine validation,
and longitudinal ADNI CSF proteomics (n = 735), we show that astrocytic PTGDS
undergoes a statistically bounded biphasic trajectory — compensatory induction
followed by abrupt exhaustion — at CPS 0.46 (Davies' p = 0.032).
Metabolic collapse precedes inflammatory amplification, repositioning LCN2-mediated
neuroinflammation as a downstream consequence rather than a causal driver.
This phase boundary corresponds clinically to MMSE 26.0 (95% CI: 24.8–27.2),
defining a quantifiable therapeutic window.

> **Preprint:** [bioRxiv link — to be inserted upon posting]

---

## Repository Structure

```
ptgds-phase-boundary/
│
├── README.md
│
├── data/
│   ├── zebrafish/
│   │   ├── behavior_14dpf.csv          # Red ball avoidance + locomotor (14 dpf)
│   │   ├── qPCR_14dpf.csv              # Gene expression, 6 groups, 14 dpf
│   │   ├── qPCR_21dpf.csv              # Gene expression, 3 groups, 21 dpf
│   │   └── ELISA_14dpf.csv             # TNF-α, IL-6 protein (pg/ml)
│   │
│   └── NOTE_restricted_data.md         # ADNI & SEA-AD access instructions
│
├── analysis/
│   ├── Fig1_SEA-AD/
│   │   ├── segmented_regression.R      # Davies' test, breakpoint estimation
│   │   ├── quadratic_model.R           # Biphasic curvature, vertex CPS 0.46
│   │   └── lagged_crosscorrelation.R   # PTGDS–LCN2–NGFR temporal ordering
│   │
│   ├── Fig2_zebrafish/
│   │   └── Fig2_panels.R               # Behavior, qPCR, ELISA, correlation
│   │
│   ├── Fig3_network_pharmacology/
│   │   └── network_pharmacology.R      # Target prediction, PPI, enrichment
│   │
│   ├── Fig4_murine/
│   │   └── murine_validation.R         # LCN2, NGFR, inflammatory markers
│   │
│   └── Fig5_ADNI/
│       ├── ADNI_proteomics.R           # LOESS, segmented regression, MMSE
│       ├── ROC_analysis.R              # AUC, DeLong CI
│       └── survival_analysis.R         # Cox regression, Kaplan–Meier
│
├── figures/
│   └── (high-resolution PNG/TIFF — added upon publication)
│
└── environment/
    ├── session_info.txt                # R sessionInfo() output
    └── packages.R                      # Package installation script
```

---

## Data Availability

| Dataset | Access | Link |
|---------|--------|------|
| SEA-AD snRNA-seq (1.3M nuclei) | Public | [Allen Brain Cell Atlas](https://portal.brain-map.org/atlases-and-data/rnaseq/human-mtg-10x_sea-ad) |
| ADNI CSF proteomics (n=735) | Registration required | [adni.loni.usc.edu](https://adni.loni.usc.edu) |
| Zebrafish qPCR / ELISA / Behavior | This repository | `data/zebrafish/` |
| Murine qPCR / Behavior | Available upon request | Contact corresponding author |

> ⚠️ **Note:** ADNI and SEA-AD raw data are subject to their respective data use agreements and **cannot be redistributed**. Analysis code is provided for reproducibility; users must obtain data independently through the above links.

---

## Reproducibility

### Requirements

- R version 4.3.2
- Key packages: `segmented`, `lme4`, `ggplot2`, `pROC`, `survival`, `CellChat`

### Installation

```r
source("environment/packages.R")
```

### Running analyses

```r
# Fig. 1 — SEA-AD segmented regression
source("analysis/Fig1_SEA-AD/segmented_regression.R")

# Fig. 2 — Zebrafish panels
source("analysis/Fig2_zebrafish/Fig2_panels.R")

# Fig. 5 — ADNI proteomics
source("analysis/Fig5_ADNI/ADNI_proteomics.R")
```

---

## Key Findings

| Finding | Value | Method |
|---------|-------|--------|
| PTGDS inflection point | CPS 0.46 | Quadratic model (β₂ = −2.07, p < 2.2e−16) |
| Segmented regression breakpoint | Bin 0.23 (95% CI: 0.13–0.33) | Davies' test (p = 0.032) |
| Clinical phase boundary | MMSE 26.0 (95% CI: 24.8–27.2) | Segmented regression, ADNI |
| Zebrafish avoidance deficit | −22.7% (p = 0.020) | Red ball avoidance, 14 dpf |
| BXP-101 rescue | +18.6 pp (p = 0.016) | BXP-101 0.6 μg/ml vs MCI |
| IL-6 suppression (protein) | −75% vs MCI | ELISA, 14 dpf |
| MCI→AD conversion risk | HR 3.2 (95% CI: 2.3–4.4) | Cox regression, ADNI |
| MCI vs AD discrimination | AUC 0.743 (95% CI: 0.698–0.788) | ROC, PTGDS/LCN2 ratio |

---

## Citation

If you use this code or data, please cite:

```
Kim Y*, Heo W†, Park S†, et al. (2025).
A lipid-metabolic exhaustion event defines a phase boundary structuring
the reversible-to-irreversible transition in MCI-to-AD progression.
bioRxiv. https://doi.org/[TO_BE_INSERTED]
```

---

## License

Code in this repository is released under the [MIT License](LICENSE).  
Zebrafish experimental data (`data/zebrafish/`) are released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

---

## Contact

**Youngouk Kim** (Corresponding author)  
BioXP Research Institute, Donghae, Gangwon-do, Republic of Korea  
sola16@gmail.com

---

*This work was conducted at BioXP Research Institute, Donghae, Gangwon-do, Republic of Korea.*  
*Zebrafish experiments were performed at Zefit Co., Ltd. (Seoul, Republic of Korea).*
