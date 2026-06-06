# ==============================================================================
# Fig3_network_pharmacology.R
#
# Purpose : Reproduce Fig. 3 panels A-D
#           Network pharmacology and molecular docking of BXP-101 components
#           (atractylodin, honokiol, wedelolactone)
#
# Panels:
#   A = Venn diagram — compound target overlap
#   B = PPI network — top 12 hub targets (STRING + Cytoscape)
#   C = Pathway/GO enrichment bar chart (-log10 FDR)
#   D = Molecular docking binding affinity summary
#
# Upstream tools (run before this script):
#   1. SwissTargetPrediction (https://swisstargetprediction.ch)
#      → targets_atractylodin.csv / targets_honokiol.csv / targets_wedelolactone.csv
#   2. STRING v12  (https://string-db.org)  score >= 0.4, top 20 interactors
#      → ppi_network.tsv  (node1, node2, combined_score)
#   3. Metascape   (https://metascape.org)
#      → enrichment_results.csv  (Term, FDR, Count, Category)
#   4. AutoDock Vina  exhaustiveness=32
#      PDB: NF-kB p65 1NFI | GSK3B 1Q3D | PTGS2 5KIR
#      → docking_scores.csv  (Compound, Target, PDB, Affinity_kcal_mol)
#
# Paper : Kim Y et al. (2025). https://github.com/YoungOukKim/MCI-to-AD
# Usage : Rscript analysis/Fig3_network_pharmacology/Fig3_network_pharmacology.R
# Requires: ggplot2, ggVennDiagram, igraph, ggraph, dplyr, readr
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(tidyr)
  library(ggVennDiagram)   # install.packages("ggVennDiagram")
  library(igraph); library(ggraph)
})

DATA_DIR <- Sys.getenv("NET_DATA", unset = "data/network_pharmacology")
OUT_DIR  <- Sys.getenv("FIG3_OUT", unset = "output/Fig3")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
DPI <- 600; WHITE <- "white"

BASE <- theme_classic(base_size = 14) +
  theme(plot.margin = margin(10, 10, 10, 10))

placeholder <- function(msg)
  ggplot() + annotate("text", x=.5, y=.5, label=msg,
    color="gray55", size=4.5, hjust=0.5) + theme_void()

# ── Helper: safe CSV load ────────────────────────────────────────────────────
load_csv <- function(fname) {
  path <- file.path(DATA_DIR, fname)
  if (!file.exists(path)) { message("  [Missing] ", fname); return(NULL) }
  read_csv(path, show_col_types = FALSE)
}

# ==============================================================================
# PANEL A — Venn diagram
# ==============================================================================
message(">>> Panel A: Venn diagram ...")

get_targets <- function(fname)
  load_csv(fname) %>% { if (is.null(.)) character(0) else pull(., Gene_Symbol) %>% unique() %>% na.omit() }

target_list <- list(
  Atractylodin  = get_targets("targets_atractylodin.csv"),
  Honokiol      = get_targets("targets_honokiol.csv"),
  Wedelolactone = get_targets("targets_wedelolactone.csv")
)

pA <- if (all(lengths(target_list) > 0)) {
  ggVennDiagram(target_list, label_alpha = 0,
                set_color = c("#4472C4","#ED7D31","#70AD47")) +
    scale_fill_gradient(low="#F7FBFF", high="#2171B5") +
    theme(legend.position = "none")
} else placeholder("Venn diagram\n(supply targets_*.csv in data/network_pharmacology/)")

ggsave(file.path(OUT_DIR,"Fig3A_venn.png"), pA,
       width=6, height=5, dpi=DPI, bg=WHITE)
message("  Saved Fig3A_venn.png")

# ==============================================================================
# PANEL B — PPI network (top 12 hubs)
# ==============================================================================
message(">>> Panel B: PPI network ...")

ppi_df <- load_csv("ppi_network.tsv") %>%
  { if (!is.null(.)) rename_with(tolower, .cols = everything()) else NULL }

pB <- if (!is.null(ppi_df)) {
  g    <- graph_from_data_frame(ppi_df[,1:2], directed = FALSE)
  hubs <- names(sort(degree(g), decreasing=TRUE))[1:min(12, vcount(g))]
  g12  <- induced_subgraph(g, which(V(g)$name %in% hubs))
  V(g12)$type <- ifelse(V(g12)$name %in% c("GSK3B","NFKB1","RELA","PTGS2"),
                        "key","other")
  ggraph(g12, layout="fr") +
    geom_edge_link(alpha=0.35, color="gray60", width=0.8) +
    geom_node_point(aes(color=type, size=degree(g12)), alpha=0.9) +
    geom_node_text(aes(label=name), repel=TRUE,
                   size=3.5, fontface="bold") +
    scale_color_manual(values=c(key="#E41A1C", other="#4472C4"),
                       guide="none") +
    scale_size_continuous(range=c(4,12), guide="none") +
    theme_graph(base_size=13) +
    theme(plot.margin=margin(10,10,10,10))
} else placeholder("PPI network\n(supply ppi_network.tsv in data/network_pharmacology/)")

ggsave(file.path(OUT_DIR,"Fig3B_PPI.png"), pB,
       width=7, height=6, dpi=DPI, bg=WHITE)
message("  Saved Fig3B_PPI.png")

# ==============================================================================
# PANEL C — Pathway enrichment
# ==============================================================================
message(">>> Panel C: Pathway enrichment ...")

enrich <- load_csv("enrichment_results.csv")
# Expected columns: Term, FDR, Count, Category

pC <- if (!is.null(enrich)) {
  df <- enrich %>%
    arrange(FDR) %>% slice_head(n=15) %>%
    mutate(Term = factor(Term, levels=rev(Term)),
           neg_logP = -log10(FDR))
  ggplot(df, aes(x=neg_logP, y=Term, fill=Category)) +
    geom_bar(stat="identity", alpha=0.85,
             color="gray30", linewidth=0.3) +
    geom_vline(xintercept=-log10(0.05),
               linetype="dashed", color="red", linewidth=0.6) +
    scale_fill_manual(values=c(
      "GO BP"="4472C4","KEGG"="#ED7D31","Reactome"="#70AD47")) +
    labs(x=expression(-log[10](FDR)), y=NULL, fill=NULL) +
    BASE + theme(axis.text.y=element_text(size=11),
                 legend.position="bottom")
} else placeholder("Pathway enrichment\n(supply enrichment_results.csv)")

ggsave(file.path(OUT_DIR,"Fig3C_enrichment.png"), pC,
       width=9, height=6, dpi=DPI, bg=WHITE)
message("  Saved Fig3C_enrichment.png")

# ==============================================================================
# PANEL D — Docking binding affinity
# ==============================================================================
message(">>> Panel D: Docking scores ...")

dock <- load_csv("docking_scores.csv")
# Expected columns: Compound, Target, PDB, Affinity_kcal_mol

pD <- if (!is.null(dock)) {
  dock %>%
    mutate(Compound = factor(Compound,
             levels=c("Wedelolactone","Honokiol","Atractylodin")),
           Label = paste0(Target,"\n(",PDB,")")) %>%
  ggplot(aes(x=Label, y=abs(Affinity_kcal_mol), fill=Compound)) +
    geom_bar(stat="identity", position="dodge",
             alpha=0.85, color="gray30", linewidth=0.3) +
    scale_fill_manual(values=c(Wedelolactone="#375623",
                               Honokiol="#70AD47",
                               Atractylodin="#C5E0B4")) +
    labs(x=NULL, y="Binding affinity (|kcal/mol|)", fill=NULL) +
    BASE + theme(legend.position="bottom",
                 axis.text.x=element_text(angle=30, hjust=1))
} else placeholder("Docking affinity\n(supply docking_scores.csv)\nPose images from AutoDock Vina supplied separately")

ggsave(file.path(OUT_DIR,"Fig3D_docking.png"), pD,
       width=7, height=5, dpi=DPI, bg=WHITE)
message("  Saved Fig3D_docking.png")

message("\n", strrep("=",60))
message("  Fig. 3 panels saved to: ", OUT_DIR)
message(strrep("=",60))
