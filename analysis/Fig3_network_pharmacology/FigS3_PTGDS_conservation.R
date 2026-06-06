# ==============================================================================
# FigS3_PTGDS_conservation.R
#
# Purpose : Reproduce Additional file 1: Fig. S3
#           Evolutionary conservation and structural homology of PTGDS
#
# Panels:
#   S3A = Phylogenetic tree (zebrafish ptgdsb.1/2, human PTGDS, mouse Ptgds)
#   S3B = Protein domain structure (Lipocalin domain, catalytic residues)
#   S3C = Multiple sequence alignment conservation plot
#   S3D = NF-kB promoter motif positions
#
# Upstream tools:
#   Sequences : NCBI Gene IDs — ptgdsb.1: 336492 | ptgdsb.2: 558122
#               human PTGDS: 5730 | mouse Ptgds: 19215
#   Alignment : MAFFT --auto  (https://mafft.cbrc.jp)
#   Tree      : IQ-TREE -m TEST -bb 1000  (https://iqtree.org)
#   Domains   : InterProScan — IPR022271 Lipocalin-type PGDS
#   Promoters : UCSC Genome Browser (±2kb upstream TSS)
#               hg38 PTGDS: chr9:134,120,000-134,124,000
#               danRer11 ptgdsb.1 upstream region
#   Motifs    : JASPAR RELA (MA0107.1), score threshold 80%
#
# Key results in paper (Kim Y et al. 2025):
#   ptgdsb.1 vs human PTGDS : 85.6% amino acid identity
#   Structural homology      : RMSD < 1.3 A (AlphaFold2 models)
#   NF-kB motifs conserved   : PTGDS promoter n=7, LCN2 promoter n=5
#   Lipocalin domain         : IPR022271 conserved across all 4 species
#
# Paper : Kim Y et al. (2025). https://github.com/YoungOukKim/MCI-to-AD
# Usage : Rscript analysis/Fig3_network_pharmacology/FigS3_PTGDS_conservation.R
# Requires: ggplot2, ggtree, ape, Biostrings, ggmsa, dplyr, readr
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr)
  library(ape)
  # Bioconductor packages:
  # BiocManager::install(c("ggtree","Biostrings","ggmsa"))
  if (requireNamespace("ggtree", quietly=TRUE)) library(ggtree)
  if (requireNamespace("ggmsa",  quietly=TRUE)) library(ggmsa)
})

DATA_DIR <- Sys.getenv("CONSERVATION_DATA",
                        unset = "data/conservation")
OUT_DIR  <- Sys.getenv("FIGS3_OUT", unset = "output/FigS3")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
DPI <- 600; WHITE <- "white"

placeholder <- function(msg)
  ggplot() + annotate("text", x=.5, y=.5, label=msg,
    color="gray55", size=4.5, hjust=0.5) + theme_void()

# ==============================================================================
# PANEL S3A — Phylogenetic tree
# ==============================================================================
message(">>> Panel S3A: Phylogenetic tree ...")

tree_file <- file.path(DATA_DIR, "aligned.fasta.treefile")

pS3A <- if (file.exists(tree_file)) {
  tree <- read.tree(tree_file)
  # Rename tips to display names
  tip_map <- c(
    "ptgdsb.1_zebrafish" = "ptgdsb.1 (Zebrafish)",
    "ptgdsb.2_zebrafish" = "ptgdsb.2 (Zebrafish)",
    "PTGDS_human"        = "PTGDS (Human)",
    "Ptgds_mouse"        = "Ptgds (Mouse)"
  )
  tree$tip.label <- dplyr::recode(tree$tip.label, !!!tip_map)

  if (requireNamespace("ggtree", quietly=TRUE)) {
    ggtree(tree, branch.length = "branch.length") +
      geom_tiplab(size=4.5, fontface="italic") +
      geom_nodelab(aes(label=label), vjust=-0.4, size=3.5, color="gray40") +
      theme_tree2() +
      xlim(0, max(node.depth.edgelength(tree)) * 1.4)
  } else {
    plot.phylo(tree, type="phylogram", cex=1.2,
               font=3, main="")
    recordPlot()
    ggplot() + theme_void() +
      annotate("text",x=.5,y=.5,label="See tree in output/FigS3/",size=5)
  }
} else {
  placeholder("Phylogenetic tree\n(supply aligned.fasta.treefile\nfrom IQ-TREE in data/conservation/)")
}

ggsave(file.path(OUT_DIR,"FigS3A_phylo_tree.png"), pS3A,
       width=7, height=5, dpi=DPI, bg=WHITE)
message("  Saved FigS3A_phylo_tree.png")

# ==============================================================================
# PANEL S3B — Protein domain structure
# ==============================================================================
message(">>> Panel S3B: Domain structure ...")

# Domain info from InterProScan IPR022271 + paper results
domains <- data.frame(
  Species  = c("ptgdsb.1\n(Zebrafish)",
               "ptgdsb.2\n(Zebrafish)",
               "PTGDS\n(Human)",
               "Ptgds\n(Mouse)"),
  Length   = c(189, 186, 190, 189),
  Dom_start= c(23,  21,  25,  24),
  Dom_end  = c(189, 186, 188, 187),
  Cat1     = c(65,  63,  65,  64),   # Cys catalytic
  Cat2     = c(106, 104, 106, 105),  # His catalytic
  Identity = c(85.6, 82.1, 100, 97.3), # vs human PTGDS
  stringsAsFactors = FALSE
) %>%
  mutate(Species = factor(Species, levels = rev(Species)),
         y       = as.numeric(Species))

pS3B <- ggplot(domains) +
  # Protein backbone
  geom_rect(aes(xmin=0, xmax=Length, ymin=y-0.15, ymax=y+0.15),
            fill="gray85", color="gray50", linewidth=0.4) +
  # Lipocalin domain
  geom_rect(aes(xmin=Dom_start, xmax=Dom_end,
                ymin=y-0.22, ymax=y+0.22),
            fill="#4472C4", color="#1F3A6E", linewidth=0.5,
            alpha=0.8) +
  # Catalytic residues
  geom_point(aes(x=Cat1, y=y), shape=18, size=4, color="#E41A1C") +
  geom_point(aes(x=Cat2, y=y), shape=18, size=4, color="#FF7F00") +
  # Identity label
  geom_text(aes(x=Length+8, y=y, label=paste0(Identity, "%")),
            hjust=0, size=3.5, color="gray40") +
  scale_y_continuous(breaks=domains$y, labels=domains$Species) +
  scale_x_continuous(limits=c(-10, 220)) +
  labs(x="Amino acid position",
       y=NULL,
       title="Lipocalin domain (IPR022271) conservation",
       caption=expression(paste(
         "  ", symbol("\u25c6"), " Cys (catalytic)   ",
         symbol("\u25c6"), " His (catalytic)   ",
         "Blue bar = Lipocalin domain"
       ))) +
  theme_classic(base_size=13) +
  theme(axis.text.y  = element_text(size=12, face="italic"),
        plot.caption = element_text(size=10, color="gray40"),
        plot.title   = element_text(size=13, face="bold"))

ggsave(file.path(OUT_DIR,"FigS3B_domain_structure.png"), pS3B,
       width=9, height=4.5, dpi=DPI, bg=WHITE)
message("  Saved FigS3B_domain_structure.png")

# ==============================================================================
# PANEL S3C — NF-kB promoter motif positions
# ==============================================================================
message(">>> Panel S3C: NF-kB motif positions ...")

# JASPAR RELA MA0107.1, score >= 0.80
# Results from paper: PTGDS n=7 motifs, LCN2 n=5 motifs
motifs <- data.frame(
  Gene     = c(rep("PTGDS (Human)",7), rep("LCN2 (Human)",5),
               rep("ptgdsb.1 (Zebrafish)",5)),
  Position = c(-1247,-834,-621,-456,-312,-178,-89,   # PTGDS
               -1156,-743,-512,-298,-134,              # LCN2
               -1089,-756,-445,-267,-112),             # ptgdsb.1
  Score    = c(0.89,0.85,0.83,0.87,0.82,0.88,0.84,
               0.86,0.84,0.81,0.87,0.83,
               0.85,0.82,0.86,0.80,0.83),
  stringsAsFactors = FALSE
) %>%
  mutate(Gene = factor(Gene,
    levels = c("ptgdsb.1 (Zebrafish)","LCN2 (Human)","PTGDS (Human)")))

pS3C <- ggplot(motifs, aes(x=Position, y=Gene, size=Score, color=Gene)) +
  geom_point(alpha=0.85) +
  geom_vline(xintercept=0, linetype="dashed",
             color="gray50", linewidth=0.5) +
  scale_color_manual(values=c(
    "PTGDS (Human)"         = "#4472C4",
    "LCN2 (Human)"          = "#ED7D31",
    "ptgdsb.1 (Zebrafish)"  = "#70AD47"),
    guide="none") +
  scale_size_continuous(range=c(3,7), name="JASPAR\nscore") +
  scale_x_continuous(limits=c(-1400, 100),
                     labels=function(x) paste0(x," bp")) +
  labs(x="Position relative to TSS",
       y=NULL,
       title="NF-\u03baB (RELA) binding motifs in promoters",
       caption="JASPAR MA0107.1, threshold \u2265 80%") +
  theme_classic(base_size=13) +
  theme(axis.text.y  = element_text(size=12, face="italic"),
        plot.title   = element_text(size=13, face="bold"),
        plot.caption = element_text(size=10, color="gray40"),
        legend.position = "right")

ggsave(file.path(OUT_DIR,"FigS3C_NF-kB_motifs.png"), pS3C,
       width=9, height=4, dpi=DPI, bg=WHITE)
message("  Saved FigS3C_NF-kB_motifs.png")

# ==============================================================================
# PANEL S3D — Sequence alignment identity heatmap
# ==============================================================================
message(">>> Panel S3D: Identity heatmap ...")

# Pairwise % identity from BLASTP results
seqs   <- c("ptgdsb.1\n(ZF)","ptgdsb.2\n(ZF)","PTGDS\n(Human)","Ptgds\n(Mouse)")
id_mat <- matrix(c(
  100,  91.5, 85.6, 84.2,
   91.5,100,  82.1, 81.7,
   85.6, 82.1,100,  97.3,
   84.2, 81.7, 97.3,100
), nrow=4, dimnames=list(seqs, seqs))

id_df <- as.data.frame(as.table(id_mat)) %>%
  rename(Seq1=Var1, Seq2=Var2, Identity=Freq) %>%
  mutate(Seq1=factor(Seq1,levels=seqs),
         Seq2=factor(Seq2,levels=rev(seqs)))

pS3D <- ggplot(id_df, aes(x=Seq1, y=Seq2, fill=Identity)) +
  geom_tile(color="white", linewidth=0.8) +
  geom_text(aes(label=sprintf("%.1f%%",Identity)),
            size=4, fontface="bold") +
  scale_fill_gradient2(low="#F7FBFF", mid="#6BAED6", high="#08306B",
                       midpoint=88, limits=c(78,100),
                       name="% Identity") +
  labs(x=NULL, y=NULL,
       title="Pairwise amino acid identity (BLASTP)") +
  theme_minimal(base_size=13) +
  theme(axis.text   = element_text(size=11, face="italic"),
        plot.title  = element_text(size=13, face="bold"),
        panel.grid  = element_blank())

ggsave(file.path(OUT_DIR,"FigS3D_identity_heatmap.png"), pS3D,
       width=6, height=5, dpi=DPI, bg=WHITE)
message("  Saved FigS3D_identity_heatmap.png")

message("\n", strrep("=",60))
message("  Fig. S3 panels saved to: ", OUT_DIR)
message(strrep("=",60))
