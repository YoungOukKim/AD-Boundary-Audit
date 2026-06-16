# =============================================================================
# external_marker_crosscohort.R
# INDEPENDENT-COHORT adjudication of the externally-proposed CSF prognostic
# markers NPTX1 / NPTXR / NPTX2 with the P3 4-axis Class logic.
#
# This is NOT a fresh per-cell engine run (the matched per-cell run is
# external_marker_audit.R on SEA-AD neurons). Gated AD snRNA atlases
# (ROSMAP/Mathys/Green/Leng/cellxgene) are not reachable from the sandbox, so
# here we apply the SAME Class I/II + modality-axis decision rules to PUBLISHED
# group-level evidence from cohorts INDEPENDENT of SEA-AD. Every value is sourced.
#
# Pre-registered prediction (external_marker_audit.R header):
#   NPTX1, NPTXR -> brain ~flat / CSF-only  => Class I / modality-specific
#   NPTX2        -> brain downward transition => Class II
#
# Sources (all independent of SEA-AD):
#  [Dai2026]  Dai et al., Nat Commun 2026; CANDI+DDI CSF, n=635. CSF NPTX1/NPTXR
#             drop with AD, predict MCI->dementia <3y (prognostic).
#  [Lao2026]  Lao et al., bioRxiv 2025.10.17.683150 (v2, 2026); 4 cohorts (JHU,
#             Banner, Irvine, Northwestern), MTG bulk RNA n=575 + PRM-MS n=135.
#             NPTX2 graded Control>MCI>AD (CvAD p<=1e-4; MCIvAD p<=1e-2; protein
#             CvAD p<=1e-3). NPTX1/NPTXR cluster in the STABLE synaptic-machinery
#             module (1/1/2), separate from the declining NPTX2 module (1/2/1).
#  [Xiao2017] Xiao et al., eLife 6:e23798; NPTX2 down in 6 brain regions (RNA+
#             protein); NPTX1/NPTXR "preserved in AD brain despite CSF reduction".
# =============================================================================
ev <- data.frame(
  gene          = c("NPTX1","NPTXR","NPTX2"),
  csf_dir       = c("down","down","down"),          # [Dai2026]/[Xiao2017]
  csf_prognostic= c("yes","yes","yes"),             # CSF predicts decline
  brain_rna_dir = c("flat","flat","down"),          # [Lao2026]/[Xiao2017]
  brain_prot_dir= c("flat","flat","down"),          # [Xiao2017]/[Lao2026]
  brain_module  = c("synaptic-machinery(stable)","synaptic-machinery(stable)",
                    "activity-dependent(declining)"),               # [Lao2026]
  brain_cohorts = c(4,4,4), brain_regions = c(6,6,6),               # [Lao2026]+[Xiao2017]
  cognition_linked_brain = c("n.s.","n.s.","yes"),  # NPTX2 loss absent in asymptomatic AD [Lao2026]
  stringsAsFactors = FALSE)

# ---- P3 decision rules (same as manuscript Class I/II + modality axis) -------
classify <- function(r){
  brain_moves <- r$brain_rna_dir!="flat" | r$brain_prot_dir!="flat"
  modality_concordant <- brain_moves && r$csf_dir!="flat"
  if (!brain_moves && r$csf_dir!="flat")
    return(c(class="Class I / modality-specific",
             axis="MODALITY: CSF signal without a brain-tissue transition",
             verdict="CSF-prognostic but NOT a brain transcriptional boundary"))
  if (modality_concordant)
    return(c(class="Class II",
             axis="cross-modality + cross-region + cross-cohort concordant",
             verdict="brain-replicating transition (true generalizable marker)"))
  c(class="indeterminate", axis="-", verdict="-")
}
out <- cbind(ev[,c("gene","csf_dir","brain_rna_dir","brain_prot_dir",
                   "brain_module","brain_cohorts","brain_regions")],
             t(sapply(seq_len(nrow(ev)), function(i) classify(ev[i,]))))
pred <- c(NPTX1="Class I / modality-specific", NPTXR="Class I / modality-specific",
          NPTX2="Class II")
out$prereg_prediction <- pred[out$gene]
out$match <- ifelse(out$class==out$prereg_prediction, "MATCH", "MISS")

dir.create("output/tables/external_markers", recursive=TRUE, showWarnings=FALSE)
write.csv(out, "output/tables/external_markers/crosscohort_verdict.csv", row.names=FALSE)
cat("\n=== INDEPENDENT CROSS-COHORT ADJUDICATION (NPTX family) ===\n")
print(out[,c("gene","csf_dir","brain_rna_dir","brain_module","class","prereg_prediction","match")],
      row.names=FALSE)
cat(sprintf("\nPrediction outcome: %s  (%d/%d markers)\n",
    ifelse(all(out$match=="MATCH"),"CONFIRMED","PARTIAL"), sum(out$match=="MATCH"), nrow(out)))
cat("Independent of SEA-AD: 4 brain cohorts (Lao2026) + 6 regions (Xiao2017) + CANDI/DDI CSF (Dai2026).\n")
cat("Matched per-cell datapoint = external_marker_audit.R on SEA-AD neurons (your run -> synaptic_audit.csv).\n")
