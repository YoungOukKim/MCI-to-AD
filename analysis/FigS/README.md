# Supplementary Figures — Analysis Scripts

R scripts for generating Supplementary Figures S1–S6.

---

## Scripts

| Script | Figure | Description |
|--------|--------|-------------|
| `FigS1_PTGDS_inflection_validation.R` | Fig. S1 | Statistical validation of biphasic PTGDS trajectory (segmented regression, LOESS, Davies' test) |
| `FigS2_module_trajectories.R`         | Fig. S2 | Module-level trajectory analysis (Purinergic/Ca²⁺, NF-κB, IEG, Ferroptosis vs PTGDS/NDUFS1) |
| `FigS4_neuronal_subtypes.R`           | Fig. S4 | Neuronal subtype vulnerability analysis |
| `FigS6_crossplatform_validation.R`    | Fig. S6 | Cross-platform validation |

---

## Data note — FigS1

`SEAAD_processed_data.RData` contains cell-level expression data (~1.3M nuclei)
derived from the SEA-AD dataset and **cannot be redistributed** under the
[Allen Brain Atlas data use policy](https://portal.brain-map.org/explore/policies).

Instead, a pre-computed bin-level summary is provided:

```
data/SEA-AD/FigS1_bin_means.csv   ← 9 rows, one per pseudo-progression bin
```

This CSV contains all columns required by `FigS1_PTGDS_inflection_validation.R`:
`bin`, `PTGDS_mean`, `PTGDS_sd`, `n_cells`, `PTGDS_se`, `PTGDS_ma`

To **regenerate** this CSV from raw SEA-AD data:
1. Download the H5AD file from [Allen Brain Cell Atlas](https://portal.brain-map.org/atlases-and-data/rnaseq/human-mtg-10x_sea-ad)
2. Run `analysis/01_SEAAD_data_processing.R`
3. Uncomment the `write.csv()` line at the bottom of `FigS1_PTGDS_inflection_validation.R`

---

## Usage

```r
# Fig. S1
Rscript analysis/FigS/FigS1_PTGDS_inflection_validation.R

# Fig. S2
Rscript analysis/FigS/FigS2_module_trajectories.R
```

Output files are written to `output/FigS1/` and `output/FigS2/` by default.
