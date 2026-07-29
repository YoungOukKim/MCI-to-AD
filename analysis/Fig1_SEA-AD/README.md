# Fig. 1 — SEA-AD snRNA-seq PTGDS Inflection

Generates Figure 1 panels A–D from the SEA-AD snRNA-seq dataset, with
donor-level statistical validation of the astrocytic PTGDS biphasic trajectory.

## Scripts

| File | Description |
|------|-------------|
| `01_SEAAD_data_extraction.R` | Extracts donor-aware processed data from the SEA-AD H5AD |
| `02_SEAAD_trajectory_analysis.R` | Builds main Tables 1 and 3 and the SEA-AD supplementary trajectory tables |
| `Figure1_panels.R` | Builds Figure 1 (A–D): trajectories + donor-level statistics |

## Usage
```r
Rscript analysis/Fig1_SEA-AD/01_SEAAD_data_extraction.R
Rscript analysis/Fig1_SEA-AD/02_SEAAD_trajectory_analysis.R
Rscript analysis/Fig1_SEA-AD/Figure1_panels.R
```

## Input
| File | Path | Description |
|------|------|-------------|
| SEA-AD H5AD | `data/SEA-AD/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad` | ~1.38M nuclei in file; 1.24M analysed after QC; 84 donors |

> SEA-AD raw data are subject to the Allen Brain Atlas data use policy and are
> not redistributed here. See `data/NOTE_restricted_data.md`. A pre-computed
> bin-level summary is provided in `data/SEA-AD/FigS1_bin_means.csv`.

## Output
- `figures/Figure1_ABCD.png` — Figure 1, panels A–D
- `output/Fig1/` — intermediate files

## Key statistics (donor-level, n = 84)
- Quadratic vertex: **CPS 0.47** (β₂ = −2.27, p = 0.006; bootstrap 95% CI 0.32–0.53)
- Segmented breakpoint: **Bin 0.23** (95% CI 0.13–0.33; Davies' p = 0.032)
- Pre-peak rise not significant at donor level; post-peak decline significant
- Curvature retained in all 84 leave-one-donor-out refits and after truncating the lowest-CPS donors
  (see `analysis/FigS/P1_vertex_CI_and_sensitivity.R`)
