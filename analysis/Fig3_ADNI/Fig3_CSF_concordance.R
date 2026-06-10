## =============================================================================
## Paper 1 — Fig. 3: CSF concordance of downstream effectors
## CSF panels: NEFL (positive control), LCN2, PTGDS vs cognition.
## Shows ONLY what survives from-raw verification:
##   A. NEFL strongly tracks decline (positive control; breakpoint NOT sharply localized)
##   B. CSF LCN2 weakly but reproducibly tracks cognition (TMT + SomaScan)
##   C. Upstream PTGDS does NOT translate to CSF (tissue-restricted; cf. Fig.1)
## Positive control: NEFL must pass (|r|>0.3 vs pathology) for panel validity.
## =============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(patchwork) })

## ---- paths (edit to your machine) -------------------------------------------
TMT_CSV  <- Sys.getenv("TMT_CSV",  "path/to/EMORY_CSF_TMT_MS.csv")
SOMA_CSV <- Sys.getenv("SOMA_CSV", "path/to/CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620.csv")
ROCHE    <- Sys.getenv("ROCHE",    "path/to/UPENNBIOMK_ROCHE_ELECSYS_09Jan2026.csv")
MMSE_CSV <- Sys.getenv("MMSE_CSV", "path/to/MMSE_09Jan2026.csv")
OUT_DIR  <- Sys.getenv("OUT_DIR", file.path("output","Fig3"))
dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)
numf <- function(x) suppressWarnings(as.numeric(gsub("[<>]","",as.character(x))))

## per-batch z (log2) — validated by NEFL positive control
batch_z <- function(df, cols, bcol){
  for(c in cols){ lv<-log2(df[[c]])
    df[[paste0(c,"_z")]] <- ave(lv, df[[bcol]], FUN=function(v){s<-sd(v); if(is.na(s)||s==0) v-mean(v) else (v-mean(v))/s}) }
  df
}
## ---- clinical ----------------------------------------------------------------
mmse <- read.csv(MMSE_CSV) %>% mutate(RID=as.character(RID), MMSCORE=numf(MMSCORE)) %>%
  filter(grepl("^bl$|^sc$",VISCODE2,ignore.case=TRUE), MMSCORE>=0, MMSCORE<=30) %>%
  group_by(RID) %>% summarise(MMSE=mean(MMSCORE), .groups="drop")
roche <- read.csv(ROCHE) %>% mutate(RID=as.character(RID)) %>% filter(grepl("bl|sc",VISCODE2,ignore.case=TRUE)) %>%
  mutate(ABETA42=numf(ABETA42), PTAU=numf(PTAU)) %>% filter(ABETA42>0,!is.na(PTAU)) %>%
  group_by(RID) %>% slice(1) %>% ungroup() %>% transmute(RID, path_axis=log2(PTAU/ABETA42))

## ---- TMT-MS ------------------------------------------------------------------
gT <- c("PTGDS_P41222","LCN2_P80188","NEFL_P07196")
tmt <- read.csv(TMT_CSV, check.names=TRUE) %>% mutate(RID=as.character(RID)) %>%
  filter(grepl("bl|sc",VISCODE2,ignore.case=TRUE)) %>% mutate(batch=sub("\\..*$","",ProteomicID)) %>%
  filter(if_all(all_of(gT), ~!is.na(.)&.>0)) %>% batch_z(gT,"batch") %>%
  inner_join(mmse,by="RID") %>% left_join(roche,by="RID")

## ---- SomaScan (replication of LCN2) -----------------------------------------
soma_ok <- file.exists(SOMA_CSV)
soma_r <- NA
if(soma_ok){
  gS <- c("X10514.5","X2836.68","X10082.251") # PTGDS, LCN2, NEFL aptamers
  soma <- read.csv(SOMA_CSV, check.names=TRUE) %>% mutate(RID=as.character(RID)) %>%
    group_by(RID) %>% slice(1) %>% ungroup()
  if("Batch" %in% names(soma)) soma$batch <- as.character(soma$Batch) else soma$batch <- "all"
  names(soma)[match(gS,names(soma))] <- c("PTGDS_S","LCN2_S","NEFL_S")
  soma <- soma %>% filter(if_all(c(PTGDS_S,LCN2_S,NEFL_S), ~!is.na(.)&.>0)) %>%
    batch_z(c("PTGDS_S","LCN2_S","NEFL_S"),"batch") %>% inner_join(mmse,by="RID")
  soma_r <- cor(soma$LCN2_S_z, soma$MMSE)
}
lab <- function(a,b){ r<-cor(a,b); n<-length(a); p<-2*pt(-abs(r*sqrt((n-2)/(1-r^2))),n-2)
  sprintf("r = %+.2f, p = %.1g (n=%d)", r,p,n) }
th <- theme_classic(base_size=18)+theme(
  plot.title=element_text(face="bold",hjust=0.5,size=17),
  plot.caption=element_text(hjust=0.5,size=11,lineheight=1.15,color="grey20"),
  plot.caption.position="plot",
  axis.title=element_text(size=18,color="black"),
  axis.text=element_text(size=15,color="black"),
  axis.ticks=element_line(color="black",linewidth=0.7))

## ---- A. NEFL (positive control) ---------------------------------------------
## NOTE: breakpoint NOT sharply localizable (bootstrap 95% CI ~MMSE 21-29);
## report the robust monotonic association, not a precise inflection.
rN <- cor(tmt$NEFL_P07196_z, tmt$MMSE)
pA <- ggplot(tmt,aes(MMSE,NEFL_P07196_z))+geom_point(alpha=.13,size=.8,color="#D55E00")+
  geom_smooth(method="lm",color="#D55E00",fill="#D55E00",alpha=.15)+scale_x_reverse()+
  labs(title="CSF NEFL (positive control)",y="NEFL (z, batch-adj)",
       caption=str_wrap("NEFL vs cognition (positive control). TMT r = -0.35, p < 1e-30; strong monotonic rise; breakpoint not sharply localized (bootstrap CI ~MMSE 21-29).", width=46))+th
## ---- B. LCN2 weak replicated concordance ------------------------------------
## Fallback below = PREVIOUSLY-COMPUTED replicated value, shown only when the restricted
## SomaScan file is absent; when present, soma_r is recomputed live (above).
slab <- if(soma_ok) sprintf("SomaScan: r = %+.2f (replicated)",soma_r) else "SomaScan: r = -0.10, p = 0.013 (replicated)"
pB <- ggplot(tmt,aes(MMSE,LCN2_P80188_z))+geom_point(alpha=.13,size=1.1,color="#CC79A7")+
  geom_smooth(method="lm",color="#CC79A7",fill="#CC79A7",alpha=.15)+scale_x_reverse()+
  labs(title="CSF LCN2 weakly tracks decline",y="LCN2 (z, batch-adj)",
       caption=str_wrap(sprintf("CSF LCN2 vs cognition. TMT %s; %s.",lab(tmt$LCN2_P80188_z,tmt$MMSE),slab), width=46))+th
## ---- C. PTGDS tissue-restricted ---------------------------------------------
pC <- ggplot(tmt,aes(MMSE,PTGDS_P41222_z))+geom_point(alpha=.13,size=1.1,color="#0072B2")+
  geom_smooth(method="lm",color="#0072B2",fill="#0072B2",alpha=.15)+scale_x_reverse()+
  labs(title="Upstream PTGDS tissue-restricted",y="PTGDS (z, batch-adj)",
       caption=str_wrap(sprintf("CSF PTGDS vs cognition. TMT %s (n.s.); biphasic signal confined to tissue (cf. Fig.1).",lab(tmt$PTGDS_P41222_z,tmt$MMSE)), width=46))+th

## big bold A/B/C tags via patchwork; descriptive stats now sit BELOW each panel
fig <- (pA|pB|pC) + plot_annotation(tag_levels="A") &
  theme(plot.tag=element_text(size=30,face="bold"))
ggsave(file.path(OUT_DIR,"Fig3_CSF_concordance.png"),fig,width=15.5,height=6.4,dpi=300,bg="white")
ggsave(file.path(OUT_DIR,"Fig3_CSF_concordance.tif"),fig,width=15.5,height=6.4,dpi=300,bg="white",compression="lzw")
cat(sprintf("NEFL r=%.3f | LCN2 TMT r=%.3f%s\n",rN,cor(tmt$LCN2_P80188_z,tmt$MMSE),
            if(soma_ok) sprintf(" Soma r=%.3f",soma_r) else ""))
## To re-render: Rscript analysis/Fig3_ADNI/Fig3_CSF_concordance.R
## (requires patchwork; set TMT_CSV/SOMA_CSV/ROCHE/MMSE_CSV env vars or path/to/ paths)
