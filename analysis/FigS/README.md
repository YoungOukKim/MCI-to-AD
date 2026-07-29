# Supplementary Figures (Additional file 1: Fig. S1–S8)

Scripts in this folder generate the supplementary figures delivered as
**Additional file 1** of the manuscript. Each script reads openly shareable
derived or summary data from `data/` and writes to `output/FigS#/`.

> **Numbering.** Supplementary figures were renumbered during manuscript
> revision so that they appear in order of first citation in the text.
> Most script names now match the final numbering. Three do not
> (`FigS5_*`, `FigS6_*`, `P1_age_cps.R`) and are flagged below. Scripts are
> intentionally **not renamed** to avoid breaking hard-coded data paths.
> Always use the **Additional file 1** column.

| Script | Additional file 1 | Description | Data |
| --- | :---: | --- | --- |
| `FigS1_PTGDS_inflection_validation.R` | **Fig. S1** | Statistical validation of the biphasic PTGDS trajectory and inflection point (raw, 3-bin MA, LOESS, segmented) | `data/SEA-AD/FigS1_bin_means.csv` |
| `FigS2_module_trajectories.R` | **Fig. S2** | Module-level mechanistic dissection of the astrocyte-neuron cascade in SEA-AD pseudo-progression | `data/SEA-AD/` |
| `FigS3_compound_characterization.R` | **Fig. S3** | Chemical standardization of BXP-101 (HPLC-DAD chromatograms, structures) | `data/FigS3/` |
| `FigS4_BV2_microglial_assay.R` | **Fig. S4** | In vitro BV-2 microglial validation of BXP-101 (non-cytotoxic anti-inflammatory synergy) | `data/FigS4/` |
| `FigS6_neuronal_subtypes.R` | **Fig. S5** *(name mismatch)* | Differential neuronal-subtype vulnerability and apoptotic buffering | `data/SEA-AD/FigS6_*.csv` |
| `FigS5_PTGDS_conservation.R` | **Fig. S6** *(name mismatch)* | PTGDS evolutionary conservation and structural homology (MSA, structure, NF-kB RELA) | `data/FigS5/` |
| `FigS7_ADNI_CSF_crossplatform.R` | **Fig. S7** | Cross-platform CSF biomarker comparison (ADNI; TMT-MS vs SomaScan) | `data/ADNI/` |
| `P1_age_cps.R` | **Fig. S8** *(name mismatch)* | Donor age at death vs pseudo-progression (CPS) — not correlated | `data/SEA-AD/age_cps_donor.csv` |

## Supplementary table generated in this folder

| Script | Additional file 1 | Description | Data |
| --- | :---: | --- | --- |
| `P1_vertex_CI_and_sensitivity.R` | **Table S5** | Bootstrap 95% CI for the quadratic vertex, leave-one-donor-out refits, low-CPS truncation, and donor counts per CPS bin | `data/SEA-AD/P1_donor_level_PTGDS.csv`, `data/SEA-AD/P1_leave_one_donor_out.csv`, `data/SEA-AD/FigS1_bin_means_with_donors.csv` |

This script also produces the bootstrap confidence interval reported in main
Figure 1 and Table 2. It runs an anchor gate first: the full-cohort fit must
reproduce the published values (84 donors, 67,419 astrocytes, beta2 = -2.27,
p = 0.006, vertex = 0.47) before any downstream output is interpreted.

## Key results referenced

- **Fig. S1:** inflection breakpoint Bin ~ 0.23 (95% CI 0.13–0.33; Davies' p = 0.032)
- **Fig. S6:** ptgdsb.1 vs human PTGDS 85.6% identity; RMSD < 1.3 A; NF-kB RELA motifs conserved
- **Fig. S8:** Spearman rho = -0.14 (p = 0.19, n = 84) — CPS is not an age proxy
- **Table S5:** vertex CPS 0.47, bootstrap 95% CI 0.32–0.53; curvature retained in all 84 leave-one-donor-out refits

## Usage

```r
# example: Additional file 1, Fig. S8 (donor-level age vs CPS)
source("analysis/FigS/P1_age_cps.R")

# example: Additional file 1, Table S5 (vertex CI and sensitivity)
source("analysis/FigS/P1_vertex_CI_and_sensitivity.R")
```

Each script reads its data path from an environment variable (default shown in
the script header) and writes to `output/FigS#/`.

> **Note:** `FigS7_ADNI_CSF_crossplatform.R` (Fig. S7) reads its cross-platform
> CSV inputs from `data/ADNI/FigS6_*.csv` — a legacy data-file prefix that does
> not match the figure number. The script is correct; only the input CSV names
> are mislabeled and can be normalized in a later pass.
