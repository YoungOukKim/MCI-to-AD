# Supplementary Figures — Analysis Scripts

R scripts for generating Supplementary Figures.

---

## Scripts

| Script | Figure | Description |
|--------|--------|-------------|
| `FigS1_PTGDS_inflection_validation.R` | Fig. S1 | Statistical validation of the biphasic PTGDS trajectory (segmented regression, LOESS, Davies' test) |
| `FigS2_module_trajectories.R`         | Fig. S2 | Module-level trajectory analysis (Purinergic/Ca²⁺, NF-κB, IEG, Ferroptosis vs PTGDS/NDUFS1) |
| `FigS3_PTGDS_conservation.R`          | Fig. S5 | Cross-species PTGDS conservation and NF-κB binding affinity (human, zebrafish) |
| `FigS4_neuronal_subtypes.R`           | Fig. S6 | Differential vulnerability across neuronal subtypes (excitatory vs SST+) |
| `FigS5_BV2_microglial_assay.R`        | Fig. S4 | BV-2 microglial assay (cell viability, NO production, NF-κB activation) |
| `FigS7_compound_characterization.R`   | Fig. S3 | HPLC chromatograms and chemical structures (atractylodin, wedelolactone, honokiol) |
| `FigureS7_ADNI_CSF_crossplatform.R`   | Fig. S7 | Cross-platform CSF comparison in ADNI (between-platform r ≈ 0.15 for PTGDS) |

> **Note on figure numbering:** script and data-file prefixes (`FigS3`, `FigS4`,
> `FigS5`, `FigS7`) are **legacy identifiers from earlier drafts** and do NOT match
> the published Supplementary Figure numbers. Authoritative mapping (script → paper):
> `FigS3_PTGDS_conservation` → **S5**; `FigS4_neuronal_subtypes` → **S6**;
> `FigS5_BV2_microglial_assay` → **S4**; `FigS7_compound_characterization` → **S3**;
> `FigureS7_ADNI_CSF_crossplatform` → **S7**. `FigS1`/`FigS2` are unchanged.

---

## Data Notes

### Fig. S1 — SEA-AD bin-level summary
`SEAAD_processed_data.RData` (~1.3M nuclei) cannot be redistributed under the
[Allen Brain Atlas data use policy](https://portal.brain-map.org/explore/policies).
A pre-computed summary is provided in `data/SEA-AD/FigS1_bin_means.csv` (9 rows).
To regenerate: download the H5AD from the
[Allen Brain Cell Atlas](https://portal.brain-map.org/atlases-and-data/rnaseq/human-mtg-10x_sea-ad),
run `analysis/Fig1_SEA-AD/01_SEAAD_data_extraction.R`, then enable the
`write.csv()` step in `FigS1_PTGDS_inflection_validation.R`.

### Fig. S2 — SEA-AD H5AD direct
Reads directly from the SEA-AD H5AD file (set via the `SEAAD_H5AD` env variable).

### Fig. S5 (PTGDS conservation) — External images
Place `FigS3A.png` and `FigS3B.png` in `data/FigS3/` (legacy `FigS3` prefix).

### Fig. S6 (neuronal subtypes) — SEA-AD summary CSVs
Pre-computed CSVs in `data/SEA-AD/`: `FigS4_proportions.csv`,
`FigS4_variance.csv`, `FigS4_bcl2.csv`, `FigS4_ratio.csv`.

### Fig. S4 (BV-2 microglial assay) — External images
Place `FigS5A.png` through `FigS5H.png` in `data/FigS5/` (legacy `FigS5` prefix).

### Fig. S7 (CSF cross-platform) — ADNI summary CSVs
Raw ADNI data cannot be redistributed under the
[ADNI data use agreement](https://adni.loni.usc.edu).
Pre-computed, de-identified group-level CSVs in `data/ADNI/`:
`FigS6_platform_cor.csv`, `FigS6_trend.csv`, `FigS6_tmt_ptgds.csv`.
(The `FigS6_` filename prefix is retained only as a stable data-file identifier.)

### Fig. S3 (compound characterization / HPLC) — External images
Place 6 images in `data/FigS7/`:
`FigS7A_atractylodin_peak.png`, `FigS7A_atractylodin_chem.png`,
`FigS7B_wedelolactone_peak.png`, `FigS7B_wedelolactone_chem.png`,
`FigS7C_honokiol_peak.png`, `FigS7C_honokiol_chem.png`.

---

## Usage

```r
Rscript analysis/FigS/FigS1_PTGDS_inflection_validation.R
Rscript analysis/FigS/FigS2_module_trajectories.R
Rscript analysis/FigS/FigS3_PTGDS_conservation.R
Rscript analysis/FigS/FigS4_neuronal_subtypes.R
Rscript analysis/FigS/FigS5_BV2_microglial_assay.R
Rscript analysis/FigS/FigS7_compound_characterization.R
Rscript analysis/FigS/FigureS7_ADNI_CSF_crossplatform.R
```

Output paths can be overridden via environment variables (see each script header).
