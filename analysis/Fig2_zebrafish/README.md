# Fig. 2 — Cross-species Conservation of Biphasic PTGDS Dynamics

BXP-101 rescues cognitive and molecular phenotypes in a reversible zebrafish MCI model,
validating timing-dependent pharmacological intervention at the PTGDS phase boundary.

## Scripts

| File | Description |
|------|-------------|
| `Fig2A_I_panels.R` | Generate individual panels A–I (no labels) |
| `Fig2_combine.R` | Assemble panels with labels into final figure |

## Panels

| Panel | Description | Data file |
|-------|-------------|-----------|
| A | Experimental timeline (triple-stressor + BXP-101 early/late) | — (schematic) |
| B | Red ball avoidance at 14 dpf — MCI 22.7% reduction; BXP-101 rescue | `behavior_14dpf.csv` |
| C | BLBP fluorescence intensity — reactive gliosis (WIF) | `WIF_14dpf.csv` |
| D | Nestin fluorescence intensity — neural progenitor state (WIF) | `WIF_14dpf.csv` |
| E | Representative WIF images — 3×3 grid (MERGE / BLBP / Nestin) | `data/zebrafish/WIF/*.png` |
| F | *ptgdsb.1/2* longitudinal qPCR at 14 and 21 dpf | `qPCR_14dpf.csv`, `qPCR_21dpf.csv` |
| G | ELISA: TNF-α (−45%) and IL-6 (−75%) at 14 dpf | `ELISA_14dpf.csv` |
| H | *ptgdsb.1* vs *ngfr* inverse correlation | `qPCR_14dpf.csv` |
| I | Molecular recovery schematic — NF-κB → ptgdsb → LCN2 → NGFR | — (schematic) |

## Data files (`data/zebrafish/`)

```
behavior_14dpf.csv    Red ball avoidance (14 dpf), 6 groups, n=10/group
WIF_14dpf.csv         BLBP (TRITC) and Nestin (GFP) fluorescence intensity
qPCR_14dpf.csv        ptgdsb.1/2, ngfr, bdnf, lcn2 relative expression (14 dpf)
qPCR_21dpf.csv        ptgdsb.1/2 relative expression (21 dpf)
ELISA_14dpf.csv       TNF-alpha and IL-6 protein levels (pg/ml)
WIF/                  Representative fluorescence images (9 PNG files)
```

## WIF images (`data/zebrafish/WIF/`)

9 confocal images — 3 groups × 3 channels:

| File | Group | Channel |
|------|-------|---------|
| `Vehicle_MERGE.png` | Control | Merged |
| `Vehicle_RFP_BLBP.png` | Control | BLBP (TRITC) |
| `Vehicle_GFP_Nestin.png` | Control | Nestin (GFP) |
| `MCI_MERGE.png` | MCI | Merged |
| `MCI_RFP_BLBP.png` | MCI | BLBP (TRITC) |
| `MCI_GFP_Nestin.png` | MCI | Nestin (GFP) |
| `BXP_101_MERGE.png` | MCI+BXP-101 | Merged |
| `BXP_101_RFP_BLBP.png` | MCI+BXP-101 | BLBP (TRITC) |
| `BXP_101_GFP_Nestin.png` | MCI+BXP-101 | Nestin (GFP) |

Antibodies: anti-BLBP (Abcam ab32423, TRITC); anti-Nestin clone 10C2 (Merck MAB5326, AF488).  
n=10 larvae/group; 14 dpf; ImageJ fixed ROI quantification.

## Key results

- Red ball avoidance: MCI 54.5% vs Control 70.5% (p=0.020); BXP-101 0.6 μg/ml rescue to 67.6% (p=0.016)
- BLBP: MCI +51% vs Control (p<0.01); BXP-101 0.4 μg/ml −47% vs MCI (p<0.001)
- TNF-α: BXP-101 0.6 μg/ml −45% vs MCI; IL-6: −75% vs MCI (ELISA, 14 dpf)
- *ptgdsb.1/2* compensatory induction at 21 dpf mirrors SEA-AD CPS 0.1–0.46 phase

## Usage

```r
# Set paths (defaults are relative to repo root)
Sys.setenv(ZEBRAFISH_DATA = "data/zebrafish")
Sys.setenv(WIF_IMG_DIR    = "data/zebrafish/WIF")
Sys.setenv(FIG2_OUT       = "output/Fig2")

# Step 1 — individual panels
source("analysis/Fig2_zebrafish/Fig2A_I_panels.R")

# Step 2 — assemble
source("analysis/Fig2_zebrafish/Fig2_combine.R")
```

## Groups

| Label | Description |
|-------|-------------|
| Control | Vehicle (no stressor) |
| MCI | Triple stressor: LPS 5 μg/L + D-galactose 0.2 mg/L + HFD 10% w/w (6–14 dpf) |
| MCI+Donepezil | MCI + Donepezil 10 μM (6–14 dpf) |
| MCI+BXP-101 0.3/0.4/0.6 μg/ml | MCI + BXP-101 early intervention (6–14 dpf) |

> **Note:** Locomotor activity data (14 dpf) and late-intervention group (15–21 dpf)
> are pending and will be added upon receipt from Zefit Inc. (CRO).
