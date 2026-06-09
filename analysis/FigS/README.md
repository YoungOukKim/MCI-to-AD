# Supplementary Figures (S1–S8)

Scripts in this folder are named by their **final manuscript supplementary figure
number**. Each reads openly shareable derived/summary data from `data/` and
writes to `output/FigS#/`.

| Script | Supp. Fig. | Description | Data |
|--------|:---:|-------------|------|
| `FigS1_PTGDS_inflection_validation.R` | **S1** | Statistical validation of the biphasic PTGDS trajectory and inflection point (raw, 3-bin MA, LOESS, segmented) | `data/SEA-AD/FigS1_bin_means.csv` |
| `FigS2_module_trajectories.R` | **S2** | Module-level astrocyte–neuron cascade in SEA-AD pseudo-progression | `data/SEA-AD/` |
| `FigS3_compound_characterization.R` | **S3** | Chemical standardization of BXP-101 (HPLC-DAD chromatograms, structures) | `data/FigS3/` |
| `FigS4_BV2_microglial_assay.R` | **S4** | In vitro BV-2 microglia validation of BXP-101 (non-cytotoxic anti-inflammatory synergy) | `data/FigS4/` |
| `FigS5_PTGDS_conservation.R` | **S5** | PTGDS evolutionary conservation and structural homology (MSA, structure, NF-κB RELA) | `data/FigS5/` |
| `FigS6_neuronal_subtypes.R` | **S6** | Differential neuronal-subtype vulnerability and apoptotic buffering | `data/SEA-AD/FigS6_*.csv` |
| `FigS7_ADNI_CSF_crossplatform.R` | **S7** | Cross-platform CSF biomarker comparison (ADNI; TMT-MS vs SomaScan) | `data/ADNI/` |
| `P1_age_cps.R` | **S8** | Donor age at death vs pseudo-progression (CPS) — not correlated | `data/SEA-AD/age_cps_donor.csv` |

## Key results referenced

- **S1:** inflection breakpoint Bin ≈ 0.23 (95% CI 0.13–0.33; Davies' p = 0.032)
- **S5:** ptgdsb.1 vs human PTGDS 85.6% identity; RMSD < 1.3 Å; NF-κB RELA motifs conserved
- **S8:** Spearman ρ = −0.14 (p = 0.19, n = 84) — CPS is not an age proxy

## Usage

```r
# example: Supplementary Fig. S8 (donor-level age vs CPS)
source("analysis/FigS/P1_age_cps.R")
```

Each script reads its data path from an environment variable (default shown in
the script header) and writes to `output/FigS#/`.

> **Note:** `FigS7_ADNI_CSF_crossplatform.R` (Supp. Fig. S7) currently reads its
> cross-platform CSV inputs from `data/ADNI/FigS6_*.csv` — a legacy data-file
> prefix that does not match the figure number. The script is correct; only the
> input CSV names are mislabeled and can be normalized in a later pass.
