# Fig. 2 — Zebrafish MCI Model Analysis

Cross-species validation of the PTGDS inflection using a reversible
zebrafish (*Danio rerio*) MCI model and timing-dependent pharmacological modulation.

---

## Scripts

| Script | Description |
|--------|-------------|
| `Fig2A_H_panels.R` | Generates individual panels A–H (PNG, 600 dpi) |
| `Fig2_combine.R`   | Assembles panels into the final figure (PNG + TIFF, 300 dpi) |

---

## Panel descriptions

| Panel | Content |
|-------|---------|
| A | Experimental timeline: triple-stressor MCI induction (LPS + D-galactose + HFD); early (6–14 dpf) vs late (15–21 dpf) BXP-101 intervention |
| B | Red ball avoidance at 14 dpf (6 groups, n = 10/group, mean ± SD) |
| C | BLBP (TRITC) fluorescence intensity — reactive gliosis (WIF, n = 10/group) |
| D | Nestin (GFP) fluorescence intensity — neural progenitor state (WIF, n = 10/group) |
| E | Representative whole-mount immunofluorescence images (3 × 3 grid) |
| F | Longitudinal *ptgdsb.1/2* qPCR at 14 and 21 dpf (n = 3–4 replicates/group) |
| G | ELISA: TNF-α and IL-6 protein levels at 14 dpf |
| H | Molecular recovery schematic: NF-κB suppression → preserved *ptgdsb* capacity → reduced TNF-α/IL-6 → NGFR re-expression |

---

## Usage

```r
# Step 1: Generate individual panels
Rscript analysis/Fig2_zebrafish/Fig2A_H_panels.R

# Step 2: Assemble final figure
Rscript analysis/Fig2_zebrafish/Fig2_combine.R
```

Output files are written to `output/Fig2/` by default.
Override with the `FIG2_OUT` environment variable:

```bash
export FIG2_OUT=/your/output/path
Rscript analysis/Fig2_zebrafish/Fig2A_H_panels.R
Rscript analysis/Fig2_zebrafish/Fig2_combine.R
```

---

## Input data

| File | Path | Description |
|------|------|-------------|
| `behavior_14dpf.csv` | `data/zebrafish/` | Red ball avoidance (6 groups × replicates) |
| `WIF_14dpf.csv`      | `data/zebrafish/` | BLBP/Nestin fluorescence intensity (n = 10/group) |
| `qPCR_14dpf.csv`     | `data/zebrafish/` | Multi-gene qPCR panel at 14 dpf |
| `qPCR_21dpf.csv`     | `data/zebrafish/` | *ptgdsb.1/2* qPCR at 21 dpf |
| `ELISA_14dpf.csv`    | `data/zebrafish/` | TNF-α and IL-6 protein (pg/ml) |
| WIF images (9 PNG)   | `data/zebrafish/WIF/` | Confocal images for Panel E |

---

## Requirements

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "cowplot", "magick"))
```

R >= 4.3.2
