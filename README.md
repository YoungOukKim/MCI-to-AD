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

Alzheimer's disease shows prolonged prodromal stability before accelerating decline, yet molecular markers resolving this heterogeneity are limited. Using pseudo-progression analysis of 1.24 million SEA-AD single nuclei retained after quality control (84 donors), we identify a reproducible biphasic astrocytic trajectory anchored to prostaglandin D2 synthase (PTGDS), with a statistically resolved donor-level inflection (quadratic beta2 = -2.27, p = 0.006; vertex CPS 0.47). The same directional change is independently reproduced in external brain proteomics (ROSMAP and Banner; AD versus control p = 3.4 x 10^-3), and the biphasic pattern reconciles previously conflicting CSF reports as stage-dependent. In ADNI CSF, downstream NEFL tracks cognitive decline strongly and LCN2 weakly, whereas PTGDS itself is tissue-restricted and not a stand-alone predictor. We propose, but do not establish, that post-inflection PTGDS attenuation accompanies LCN2-linked inflammation and NGFR suppression. These data position astrocytic PTGDS as a candidate stage marker, not a causal driver, of the compensatory-to-vulnerable shift in the aging brain — the astrocytic PTGDS inflection (CPS 0.47).

> **Preprint:** Research Square — https://doi.org/10.21203/rs.3.rs-9499795 (not peer reviewed). The manuscript is currently under review at *Alzheimer's Research & Therapy*.

---

## Repository structure

```
MCI-to-AD/
├── README.md
├── LICENSE                              # MIT (code); shared data CC BY 4.0
├── .gitignore
│
├── analysis/
│   ├── Fig1_SEA-AD/                     # Fig. 1 — SEA-AD snRNA-seq PTGDS trajectories
│   │   ├── 01_SEAAD_data_extraction.R
│   │   ├── 02_SEAAD_trajectory_analysis.R          # main Tables 1 and 3, plus supplementary tables
│   │   ├── Figure1_panels.R
│   │   └── README.md
│   ├── Fig2_zebrafish/                  # Fig. 2 — reversible zebrafish MCI model
│   │   ├── Fig2A_H_panels.R
│   │   ├── Fig2_combine.R
│   │   └── README.md
│   ├── Fig3_ADNI/                       # Fig. 3 — ADNI CSF concordance (NEFL / LCN2 / PTGDS)
│   │   ├── 03_ADNI_clinical_analysis.R
│   │   ├── Fig3_CSF_concordance.R
│   │   └── README.md
│   ├── Fig4_network_pharmacology/       # Fig. 4 — network pharmacology and in silico docking
│   │   ├── Fig4_network_pharmacology.R
│   │   └── README.md
│   ├── Fig5_murine/                     # Fig. 5 — murine validation
│   │   ├── Fig5_panels.R
│   │   ├── Fig5_combine.R
│   │   └── README.md
│   ├── FigS/                            # supplementary figures S1–S8
│   │   ├── FigS1_PTGDS_inflection_validation.R    # Supp. Fig. S1
│   │   ├── FigS2_module_trajectories.R            # Supp. Fig. S2
│   │   ├── FigS3_compound_characterization.R      # Supp. Fig. S8 (BXP-101 HPLC)
│   │   ├── FigS4_BV2_microglial_assay.R           # Supp. Fig. S5 (BV-2 in vitro)
│   │   ├── FigS5_PTGDS_conservation.R             # Supp. Fig. S4 (conservation)
│   │   ├── FigS6_neuronal_subtypes.R              # Supp. Fig. S3 (neuronal subtypes)
│   │   ├── FigS7_ADNI_CSF_crossplatform.R         # Supp. Fig. S6 (CSF cross-platform)
│   │   ├── P1_age_cps.R                           # Add. file 1, Fig. S8 (age vs CPS)
│   │   ├── P1_vertex_CI_and_sensitivity.R         # Add. file 1, Table S5 (vertex CI, LODO, donor counts)
│   │   └── README.md
│   └── Tables/
│       ├── TableS1_QC_Braak_mapping.R
│       └── Tables_main_and_S1_S4.R
│
├── data/
│   ├── SEA-AD/                          # derived donor-level CSVs (see mapping tables below)
│   ├── ADNI/                            # de-identified group-level summary CSVs (cross-platform)
│   ├── murine/                          # murine_behavioral_molecular_data.xlsx
│   ├── zebrafish/                       # raw 14/21 dpf behavior, qPCR, ELISA, WIF (+ README)
│   ├── FigS3/  FigS4/  FigS5/           # panel images (repo FigS3=Supp S8 HPLC, FigS4=Supp S5 BV-2, FigS5=Supp S4 conservation)
│   └── NOTE_restricted_data.md          # SEA-AD / ADNI access terms (not redistributable)
│
├── assets/                              # repository cover images
├── environment/
│   └── packages.R                       # R package manifest
├── figures/                             # assembled manuscript figure outputs
└── output/                              # script-generated tables/figures (created on run)
```

> **Main figures:** **Fig. 1** SEA-AD trajectories · **Fig. 2** zebrafish MCI model · **Fig. 3** ADNI CSF concordance · **Fig. 4** network pharmacology / docking · **Fig. 5** murine validation. Main-figure folders match the manuscript numbering.

> **Supplementary numbering.** The manuscript supplementary items are delivered as **Additional file 1** (figures S1–S8, tables S1–S10) and **Additional file 2** (machine-readable copy of Table S10). Most script names now match the final numbering; the three exceptions are listed below. Files are intentionally **not renamed** to avoid breaking hard-coded data paths and reproducibility.

### Supplementary figure -> repository script

| Additional file 1 | Content | Repository script |
| --- | --- | --- |
| **Fig. S1** | Statistical validation of the biphasic PTGDS trajectory and inflection point | `analysis/FigS/FigS1_PTGDS_inflection_validation.R` |
| **Fig. S2** | Module-level mechanistic dissection of the astrocyte-neuron cascade | `analysis/FigS/FigS2_module_trajectories.R` |
| **Fig. S3** | Chemical standardization of BXP-101 (HPLC-DAD) | `analysis/FigS/FigS3_compound_characterization.R` |
| **Fig. S4** | In vitro BV-2 microglial validation of BXP-101 | `analysis/FigS/FigS4_BV2_microglial_assay.R` |
| **Fig. S5** | Differential vulnerability and apoptotic buffering across neuronal subtypes | `analysis/FigS/FigS6_neuronal_subtypes.R` *(name mismatch)* |
| **Fig. S6** | Evolutionary conservation and structural homology of PTGDS | `analysis/FigS/FigS5_PTGDS_conservation.R` *(name mismatch)* |
| **Fig. S7** | Cross-platform CSF biomarker comparison (ADNI; TMT-MS vs SomaScan) | `analysis/FigS/FigS7_ADNI_CSF_crossplatform.R` |
| **Fig. S8** | Donor age at death is not correlated with CPS | `analysis/FigS/P1_age_cps.R` *(name mismatch)* |

### Supplementary table -> repository source

> Script-internal `TableS#` labels do **not** match the final numbering; use the **Additional file 1** column.

| Additional file 1 | Content | Repository source (internal label) |
| --- | --- | --- |
| **Table S1** | QC metrics and bin-to-Braak stage mapping | `analysis/Tables/TableS1_QC_Braak_mapping.R` |
| **Table S2** | Antibody and primer resources | manuscript supplementary (not script-generated) |
| **Table S3** | Behavioural test results, amyloid-beta murine model | `data/murine/murine_behavioral_molecular_data.xlsx` |
| **Table S4** | Bin-resolved cell-type expression trajectories | `analysis/Tables/Tables_main_and_S1_S4.R` -> `TableS4_Full_Trajectories.csv` |
| **Table S5** | Donor counts per CPS bin and robustness of the inflection | `analysis/FigS/P1_vertex_CI_and_sensitivity.R`; data `data/SEA-AD/FigS1_bin_means_with_donors.csv` |
| **Table S6** | Spearman correlations of apoptotic, identity and reactive markers with CPS | `Tables_main_and_S1_S4.R` -> `TableS5_Apoptotic_Identity_Markers.csv` |
| **Table S7** | Lagged cross-correlation across the full CPS trajectory | `Tables_main_and_S1_S4.R` -> `TableS6_Full_Trajectory_CCF.csv` |
| **Table S8** | Neuronal subtype expression and vulnerability across CPS bins | `Tables_main_and_S1_S4.R` -> `TableS7_Neuronal_Subtype_Vuln.csv` |
| **Table S9** | ADNI participant demographics | `analysis/Fig3_ADNI/03_ADNI_clinical_analysis.R` -> `TableS8_Demographics.csv` |
| **Table S10** | Independent brain bulk proteomic validation (ROSMAP / Banner) | standalone workbook, delivered as **Additional file 2** |

---

## Data availability

| Dataset | Access | Link |
| --- | --- | --- |
| SEA-AD MTG snRNA-seq (~1.38M nuclei in the released file; 1.24M analysed after QC) | Public | https://portal.brain-map.org/explore/seattle-alzheimers-disease |
| ADNI CSF proteomics | Registration + DUA | https://adni.loni.usc.edu |
| ROSMAP / Banner proteomic statistics | Published (Johnson et al.) via AD Knowledge Portal | https://adknowledgeportal.synapse.org |

> ADNI and SEA-AD raw data are subject to their respective data use agreements and **cannot be redistributed**. ROSMAP/Banner values used here were taken from the published study of Johnson et al. (Supplementary Tables 2 and 17); no raw spectra were reanalysed. Analysis code is provided for reproducibility; users must obtain restricted data independently.

---

## Key findings

| Finding | Value | Method |
| --- | --- | --- |
| Astrocytic PTGDS inflection (donor-level) | vertex CPS 0.47 (beta2 = -2.27, p = 0.006; bootstrap 95% CI 0.32-0.53; n = 84) | Donor-level quadratic model |
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
* Key packages: `hdf5r`, `segmented`, `ggplot2`, `cowplot` (see `environment/packages.R`)

```r
# Donor-level age vs CPS (Supplementary Fig. S7)
source("analysis/FigS/P1_age_cps.R")
```

Each `analysis/` subfolder has its own README describing inputs, outputs, and run order.

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
