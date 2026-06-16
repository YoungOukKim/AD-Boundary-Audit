# =============================================================================
# external_marker_panel_audit.R
# Adjudicates a PANEL of "promising" AD fluid markers with the P3 4-axis Class
# logic, using PUBLISHED independent-cohort directions (brain tissue vs CSF).
# NOT a fresh per-cell run; every direction is sourced (see refs[]).
# Output: output/tables/external_markers/panel_verdict.csv
# =============================================================================
# brain_dir / csf_dir in {up, down, flat, up_ns}; region_conflict for HMOX1.
ev <- read.csv(text=
"gene,csf_dir,brain_dir,brain_mods,breadth,region_conflict,source
PTGDS,na,up,RNA+protein,P1 ROSMAP/Banner ext-proteomics,no,P1+companion
HMOX1,na,up,RNA,MTG-up/EC-down,yes,P3
NPTX1,down,flat,RNA+protein,4 cohorts/6 regions,no,Dai2026/Lao2026/Xiao2017
NPTXR,down,flat,RNA+protein,4 cohorts/6 regions,no,Dai2026/Lao2026/Xiao2017
NPTX2,down,down,RNA+protein,4 cohorts/6 regions,no,Lao2026/Xiao2017
GFAP,up,up,RNA+protein,multi-cohort/region,no,reactive-astro/plasma
CHI3L1,up,up,RNA+protein,multi-cohort/region,no,Lananna2020/CSF-meta
TREM2,up,up,RNA(brain)+sTREM2(CSF),multi-cohort,no,DAM/preclinical
NRGN,up,down,brain/EV down; CSF up,multi-cohort,no,postmortem/Simoa
SNAP25,up,down,brain/EV down; CSF up,multi-cohort,no,EV-down/CSF-up
VAMP2,up,down,brain down; CSF up,multi-cohort,no,CSF-synapse/EV
VGF,down,down,RNA+protein,multi-cohort,no,Lao2026/granin-proteoform
SCG2,up_ns,down,brain down; CSF up/ns,multi-cohort,no,Lao2026/CSF-mixed
CNTN2,down,down,CSF+brain down,2 cohorts+postmortem,no,Amsterdam-DC", stringsAsFactors=FALSE)


# ---- honest per-marker evidence strength (NOT all equally certain) ----
EVID <- c(GFAP="robust", NPTX2="robust", TREM2="robust", NRGN="robust", SNAP25="robust",
          VGF="moderate-robust", VAMP2="moderate", CHI3L1="moderate(stage-inconsistent)",
          NPTX1="moderate(brain=inference)", NPTXR="moderate(brain=inference)",
          CNTN2="thin(2 cohorts,n=14 brain)", SCG2="thin(CSF up/ns contested)",
          PTGDS="own-program(P1 ext-validated)", HMOX1="own-program(P3 region-conflict)")
BASIS <- "literature-encoded endpoint directions (Control vs AD / module), NOT a trajectory engine run"

sign <- function(d) ifelse(d=="up",1, ifelse(d=="down",-1, ifelse(d=="up_ns",1, 0)))
classify <- function(r){
  if (r$region_conflict=="yes")
    return(c("Class I / region-specific artifact","REFUTED (direction flips across regions)"))
  bm <- r$brain_dir!="flat"
  if (!bm && r$csf_dir!="na")
    return(c("Class I / modality-specific","modality-limited: CSF signal, no brain transition"))
  if (bm && r$csf_dir=="na")                       # tissue-anchored (no CSF axis tested)
    return(c("Class II (tissue)","brain-replicating transition"))
  sb <- sign(r$brain_dir); sc <- sign(r$csf_dir)
  if (bm && sc!=0 && sb==sc)
    return(c("Class II / concordant","CONFIRMED: brain + CSF same direction, cross-axis"))
  if (bm && r$brain_dir=="down" && r$csf_dir=="up")
    return(c("Class II (tissue) / modality-INVERSE","brain loss; CSF rise = release readout (do NOT read as brain up)"))
  if (bm && r$brain_dir=="down" && r$csf_dir=="up_ns")
    return(c("Class II (tissue) / modality-divergent","brain loss; CSF up/ns -> modality axis unresolved"))
  c("indeterminate","-")
}
res <- t(sapply(seq_len(nrow(ev)), function(i){ c <- classify(ev[i,]); c }))
out <- cbind(ev[,c("gene","csf_dir","brain_dir","breadth","source")],
             class=res[,1], framework_read=res[,2], evidence=EVID[ev$gene], basis=BASIS)
dir.create("output/tables/external_markers", recursive=TRUE, showWarnings=FALSE)
write.csv(out, "output/tables/external_markers/panel_verdict.csv", row.names=FALSE)

cat("\n=== PANEL ADJUDICATION: promising AD markers through the 4-axis framework ===\n")
print(out[,c("gene","csf_dir","brain_dir","class","evidence")], row.names=FALSE)
cat("\n--- regime counts ---\n")
reg <- ifelse(grepl("concordant|tissue\\)$|^Class II \\(tissue\\)",out$class) & !grepl("INVERSE|divergent",out$class),"Class II concordant (CONFIRMED)",
       ifelse(grepl("INVERSE",out$class),"Class II modality-INVERSE (reinterpreted)",
       ifelse(grepl("divergent",out$class),"modality-divergent",
       ifelse(grepl("modality-specific",out$class),"Class I modality-specific",
       ifelse(grepl("artifact",out$class),"Class I region artifact (REFUTED)","other")))))
print(as.data.frame(table(regime=reg)), row.names=FALSE)
cat("\nCAVEAT: literature PRIOR (endpoint directions), not a fresh engine run; evidence varies per\n")
cat("marker (evidence column). Independent computational arbiter = SEA-AD donor-level run.\n")
cat("\nTakeaway: the panel does NOT uniformly confirm. The framework sorts 'promising'\n")
cat("markers into generalizable (concordant Class II), brain-loss-with-inverse-CSF,\n")
cat("modality-specific (CSF-only), and region-specific artifact regimes.\n")
