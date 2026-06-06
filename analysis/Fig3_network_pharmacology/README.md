# Fig. 3 & Fig. S3 — Network Pharmacology and PTGDS Conservation

## Contents

| Script | Output | Description |
|--------|--------|-------------|
| `Fig3_network_pharmacology.R` | `output/Fig3/` | Fig. 3A–D: BXP-101 target network and docking |
| `FigS3_PTGDS_conservation.R` | `output/FigS3/` | Fig. S3A–D: PTGDS evolutionary conservation |

---

## Fig. 3 — Network Pharmacology of BXP-101

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

## Fig. S3 — PTGDS Evolutionary Conservation

### Panels
| Panel | Description |
|-------|-------------|
| S3A | Phylogenetic tree — ptgdsb.1/2 (zebrafish), PTGDS (human), Ptgds (mouse) |
| S3B | Protein domain structure — Lipocalin domain (IPR022271), catalytic residues |
| S3C | NF-κB (RELA) binding motifs in promoters |
| S3D | Pairwise amino acid identity heatmap (BLASTP) |

### Key results (paper)
- ptgdsb.1 vs human PTGDS: **85.6% amino acid identity**
- Structural homology: **RMSD < 1.3 Å** (AlphaFold2 models)
- Lipocalin domain (IPR022271): conserved across all 4 species
- NF-κB motifs: PTGDS promoter n=7, LCN2 promoter n=5, ptgdsb.1 n=5

### Upstream tools
| Step | Tool | URL |
|------|------|-----|
| Sequences | NCBI Entrez | Gene IDs: 336492, 558122, 5730, 19215 |
| Alignment | MAFFT --auto | https://mafft.cbrc.jp |
| Phylogeny | IQ-TREE -m TEST -bb 1000 | https://iqtree.org |
| Domains | InterProScan IPR022271 | https://www.ebi.ac.uk/interpro |
| Motifs | JASPAR MA0107.1 (RELA), ≥80% | https://jaspar.genereg.net |

### Required data files (`data/conservation/`)
```
aligned.fasta.treefile      # IQ-TREE output (Newick format)
```
> S3B–D panels use values hard-coded from paper results and
> do not require additional input files.

---

## Usage

```r
# Fig. 3
Sys.setenv(NET_DATA = "data/network_pharmacology")
Sys.setenv(FIG3_OUT = "output/Fig3")
source("analysis/Fig3_network_pharmacology/Fig3_network_pharmacology.R")

# Fig. S3
Sys.setenv(CONSERVATION_DATA = "data/conservation")
Sys.setenv(FIGS3_OUT         = "output/FigS3")
source("analysis/Fig3_network_pharmacology/FigS3_PTGDS_conservation.R")
```

## R package requirements

```r
install.packages(c("ggplot2","dplyr","readr","tidyr",
                   "ggVennDiagram","igraph","ggraph","ape"))
# Bioconductor
BiocManager::install(c("ggtree","Biostrings","ggmsa"))
```
