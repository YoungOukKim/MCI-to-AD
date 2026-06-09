# Fig. 4 — Network Pharmacology of BXP-101

| Script | Output | Description |
|--------|--------|-------------|
| `Fig4_network_pharmacology.R` | `output/Fig4/` | **Fig. 4A–D**: BXP-101 target network and in silico docking |

> PTGDS evolutionary conservation (**Supplementary Fig. S5**) lives in
> `analysis/FigS/FigS5_PTGDS_conservation.R`.

---

### Panels
| Panel | Description |
|-------|-------------|
| A | Venn diagram — target overlap among atractylodin, honokiol, wedelolactone |
| B | PPI network — top 12 hub targets (STRING v12 + Cytoscape v3.9) |
| C | Pathway enrichment — top GO/KEGG terms (Metascape, −log₁₀ FDR) |
| D | Molecular docking binding affinities (AutoDock Vina, exhaustiveness=32) |

### Key hub targets
NF-κB p65 (PDB: 1NFI), GSK3B (PDB: 1Q3D), PTGS2 (PDB: 5KIR)

### Required data files (`data/network_pharmacology/`)
```
targets_atractylodin.csv    # SwissTargetPrediction / PharmMapper output
targets_honokiol.csv
targets_wedelolactone.csv
ppi_network.tsv             # STRING export (node1, node2, combined_score)
enrichment_results.csv      # Metascape output (Term, FDR, Count, Category)
docking_scores.csv          # AutoDock Vina summary (Compound, Target, PDB, Affinity_kcal_mol)
```

---

## Usage

```r
Sys.setenv(NET_DATA = "data/network_pharmacology")
Sys.setenv(FIG4_OUT = "output/Fig4")
source("analysis/Fig4_network_pharmacology/Fig4_network_pharmacology.R")
```

## R package requirements

```r
install.packages(c("ggplot2","dplyr","readr","tidyr",
                   "ggVennDiagram","igraph","ggraph"))
```
