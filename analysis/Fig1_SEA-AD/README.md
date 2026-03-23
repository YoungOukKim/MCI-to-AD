# Fig. 1 — SEA-AD Single-cell Trajectories

Single-nucleus trajectories in SEA-AD reveal an astrocytic PTGDS phase boundary
structuring the MCI-to-AD transition.

## Scripts

| File | Order | Description |
|------|-------|-------------|
| `01_SEAAD_data_extraction.R` | 1st | Download and preprocess SEA-AD snRNA-seq data (MTG atlas) |
| `02_SEAAD_trajectory_analysis.R` | 2nd | CPS binning, trajectory smoothing, segmented regression, CCF |
| `02b_SEAAD_supplementary_tables.R` | 3rd | Generate Additional file 1: Tables S1–S4, S7 |
| `Fig1_SEA-AD_panels.R` | 4th | Generate individual Fig. 1 panels A–E (no labels) |
| `FigS2_PTGDS_biphasic_validation.R` | 5th | Generate Additional file 1: Fig. S2 (change-point validation) |

## Panels

| Panel | Description |
|-------|-------------|
| A | Study schema — SEA-AD (1.3M nuclei, n=84), zebrafish model, ADNI CSF |
| B | Mechanistic framework — PTGDS–LCN2–NGFR axis |
| C | Neuronal NDUFS1 decline from Bin 0.1; concurrent PTGDS upregulation |
| D | Biphasic PTGDS trajectory — slope change at Bin 0.23 (Davies' p=0.032), peak at CPS 0.46 |
| E | Quadratic model curvature (β₂=−2.07, p<2.2e−16); vertex at CPS 0.46 |

## Data source

**SEA-AD Middle Temporal Gyrus snRNA-seq atlas**
- Allen Brain Cell Atlas: https://portal.brain-map.org/atlasAndData/bulk-data
- Version: MTG snRNA-seq (2024-02-13)
- 1,240,908 nuclei (post-QC) from 84 donors
- CPS (Continuous Pseudo-progression Score) provided by SEA-AD consortium

```r
# Data access — Allen Brain Cell Atlas API
# Reference: Gabitto et al. (2024) Nat Neurosci doi:10.1038/s41593-024-01774-5
```

## Key results

| Metric | Value |
|--------|-------|
| Segmented regression inflection | Bin 0.23 (95% CI: 0.13–0.33) |
| Davies' test | p = 0.032 |
| Quadratic curvature β₂ | −2.07 (p < 2.2e−16) |
| PTGDS peak (vertex) | CPS 0.46 |
| Pre-peak slope | +0.244 (p = 0.000345) |
| Post-peak slope | −1.213 (p < 2.2e−16) |
| Astrocytes analyzed | n = 67,419 |
| Total nuclei | 1,240,908 (84 donors) |

## QC parameters

- Minimum genes per nucleus: ≥ 500
- Maximum mitochondrial read fraction: < 20%
- Doublet removal: Scrublet (Wolock et al. 2019)
- CPS bins: 0.1–0.9 (rounded to 1 decimal place)
- Trajectory smoothing: 3-bin moving average
- CCF lag.max = 3; inflection window = Bins 0.4–0.8

## Usage

```r
# Set paths (defaults are relative to repo root)
Sys.setenv(SEAAD_DATA = "data/SEA-AD")
Sys.setenv(FIG1_OUT   = "output/Fig1")

# Run in order
source("analysis/Fig1_SEA-AD/01_SEAAD_data_extraction.R")
source("analysis/Fig1_SEA-AD/02_SEAAD_trajectory_analysis.R")
source("analysis/Fig1_SEA-AD/02b_SEAAD_supplementary_tables.R")
source("analysis/Fig1_SEA-AD/Fig1_SEA-AD_panels.R")
source("analysis/Fig1_SEA-AD/FigS2_PTGDS_biphasic_validation.R")
```

## R package requirements

```r
install.packages(c("ggplot2","dplyr","tidyr","readr",
                   "segmented","strucchange","lme4"))
# Bioconductor
BiocManager::install(c("Seurat","SingleCellExperiment"))
```

## Citation

SEA-AD consortium data:
> Gabitto MI, Travaglini KJ, Rachleff VM, et al. (2024).
> Integrated multimodal cell atlas of Alzheimer's disease.
> *Nature Neuroscience*. doi:10.1038/s41593-024-01774-5

All analyses reproducible with `set.seed(42)` in R v4.3.2.
