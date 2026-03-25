Generates Figure 1 panels A–E from the SEA-AD H5AD dataset.

## Usage
```r
Rscript analysis/Fig1_SEAAD/Fig1_panels.R
```

## Input
| File | Path | Description |
|------|------|-------------|
| SEA-AD H5AD | `data/SEA-AD/SEAAD_MTG_RNAseq_final-nuclei.h5ad` | 1.3M nuclei, 84 donors |
| Tissue overview | `data/Fig1/SEA-AD_tissue_overview.jpg` | Panel A image |
| BioRender schematic | `data/Fig1/Fig1A_schematic.png` | Panel B image |

## Output
- `output/Fig1/Fig1_final.png` — 300 dpi PNG
- `output/Fig1/Fig1_final.tiff` — 300 dpi TIFF, LZW (journal submission)
