# A biphasic astrocytic PTGDS trajectory marks a metabolic vulnerability stage in prodromal Alzheimer's disease

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Preprint: Research Square](https://img.shields.io/badge/Preprint-Research%20Square-blue)](https://doi.org/10.21203/rs.3.rs-9499795)
[![R version](https://img.shields.io/badge/R-4.3.2-276DC3)](https://www.r-project.org/)

Analysis code and openly shareable data for the study of a reproducible **astrocytic PTGDS inflection (CPS 0.47)** that marks a metabolic vulnerability stage during the MCI-to-AD continuum.

---

## Authors

**YoungOuk Kim**†*, **WooMyung Heo**†, **Se Jin Park**†, YoungChul Kim, Ye Eun Cho, Ye-Won Lee, JungYeon Kim

† Co-first authors
* Corresponding author: YoungOuk Kim — yo.kim@bioxp.biz (BioXP Research Institute)

> Full author affiliations are listed in the manuscript.

---

## Abstract

Alzheimer's disease shows prolonged prodromal stability before accelerating decline, yet molecular markers resolving this heterogeneity are limited. Using pseudo-progression analysis of 1.3 million SEA-AD single nuclei (84 donors), we identify a reproducible biphasic astrocytic trajectory anchored to prostaglandin D2 synthase (PTGDS), with a statistically resolved donor-level inflection (quadratic beta2 = -2.27, p = 0.006; vertex CPS 0.47). The same directional change is independently reproduced in external brain proteomics (ROSMAP and Banner; AD versus control p = 3.4 x 10^-3), and the biphasic pattern reconciles previously conflicting CSF reports as stage-dependent. In ADNI CSF, downstream NEFL tracks cognitive decline strongly and LCN2 weakly, whereas PTGDS itself is tissue-restricted and not a stand-alone predictor. We propose, but do not establish, that post-inflection PTGDS attenuation accompanies LCN2-linked inflammation and NGFR suppression. These data position astrocytic PTGDS as a candidate stage marker, not a causal driver, of the compensatory-to-vulnerable shift in the aging brain — the astrocytic PTGDS inflection (CPS 0.47).

> **Preprint:** Research Square — https://doi.org/10.21203/rs.3.rs-9499795 (not peer reviewed; under journal consideration)

---

## Repository structure

> Adjust this section to match the files actually committed (run `git ls-files`).

```
MCI-to-AD/
├── README.md
├── LICENSE
├── .gitignore
│
├── data/
│   └── SEA-AD/
│       └── age_cps_donor.csv          # Donor-level age vs CPS (n = 84; Supplementary Fig. S8)
│
├── analysis/
│   └── FigS/
│       └── P1_age_cps.R               # Donor-level age–CPS correlation + scatter
│
└── figures/
    ├── Figure1_ABCD.png               # SEA-AD trajectories; inflection vertex CPS 0.47
    ├── Fig3_CSF_concordance.png       # ADNI CSF: NEFL / LCN2 / PTGDS vs MMSE
    └── age_cps_scatter.png            # Age vs CPS (Supplementary Fig. S8)
```

Manuscript figure mapping: **Fig. 1** SEA-AD single-nucleus trajectories · **Fig. 2** zebrafish MCI model · **Fig. 3** ADNI CSF concordance · **Fig. 4** network pharmacology and in silico docking · **Fig. 5** murine validation.

---

## Data availability

| Dataset | Access | Link |
| --- | --- | --- |
| SEA-AD snRNA-seq (1.3M nuclei) | Public | https://portal.brain-map.org |
| ADNI CSF proteomics | Registration + DUA | https://adni.loni.usc.edu |
| ROSMAP / Banner proteomic statistics | Published (Johnson et al.) via AD Knowledge Portal | https://adknowledgeportal.synapse.org |

> ADNI and SEA-AD raw data are subject to their respective data use agreements and **cannot be redistributed**. ROSMAP/Banner values used here were taken from the published study of Johnson et al. (Supplementary Tables 2 and 17); no raw spectra were reanalysed. Analysis code is provided for reproducibility; users must obtain restricted data independently.

---

## Key findings

| Finding | Value | Method |
| --- | --- | --- |
| Astrocytic PTGDS inflection (donor-level) | vertex CPS 0.47 (beta2 = -2.27, p = 0.006; n = 84) | Donor-level quadratic model |
| Segmented breakpoint | CPS 0.23 (95% CI 0.13–0.33) | Davies' test (p = 0.032) |
| External proteomic reproduction | AD vs control p = 3.4 x 10^-3 | ROSMAP/Banner consensus TMT-MS (Johnson et al.) |
| CSF NEFL (positive control) | r = -0.35, p < 10^-30 | ADNI TMT-MS vs MMSE |
| CSF LCN2 | TMT r = -0.14 (n = 1104); SomaScan r = -0.10 | ADNI proteomics vs MMSE |
| CSF PTGDS | r = -0.02 (n.s.; n = 1104) — tissue-restricted | ADNI proteomics vs MMSE |
| MCI-to-AD conversion | **not predicted** (PTGDS HR 1.04, p = 0.59; PTGDS/LCN2 ratio HR 0.93, p = 0.37) | Pre-registered Cox regression |
| Age independence of CPS | Spearman rho = -0.14 (p = 0.19, n = 84) | Donor-level correlation |
| Clinical alignment | MMSE 24–27 range | Cross-scale alignment |
| Zebrafish avoidance deficit | -22.7% (p = 0.020) | Red-ball avoidance, 14 dpf |
| Cross-species PTGDS conservation | 85.6% amino-acid identity | zebrafish ptgdsb vs human PTGDS |

> BXP-101 is used as a broad, multi-component NF-kB **pharmacological probe** to test whether the pre-inflection window is modifiable; it is not a PTGDS-specific intervention and does not establish PTGDS causality. See the competing-interests statement.

---

## Reproducibility

* R version 4.3.2
* Key packages: `hdf5r`, `segmented`, `ggplot2`, `cowplot`

```r
# Donor-level age vs CPS (Supplementary Fig. S8)
source("analysis/FigS/P1_age_cps.R")
```

---

## Citation

If you use this code or data, please cite:

```
Kim Y, Heo W, Park SJ, et al. (2026).
A biphasic astrocytic PTGDS trajectory marks a metabolic vulnerability stage in prodromal Alzheimer's disease.
Research Square. https://doi.org/10.21203/rs.3.rs-9499795
```

---

## Competing interests

Y.O. Kim, W.M. Heo, and Y.C. Kim are affiliated with BioXP, Inc. (BioXP Research Institute), the developer of BXP-101, which was used in this study as a pharmacological tool. W.M. Heo is the Chief Executive Officer of BioXP, Inc., and Y.O. Kim is Director of its Research Institute. These authors therefore declare a competing financial interest. The remaining authors (S.J. Park, Y.E. Cho, Y.-W. Lee, J.Y. Kim) declare no competing interests.

---

## License

Code is released under the [MIT License](LICENSE). Openly shared data files are released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

---

## Contact

**YoungOuk Kim** (Corresponding author) — BioXP Research Institute, Republic of Korea — yo.kim@bioxp.biz
