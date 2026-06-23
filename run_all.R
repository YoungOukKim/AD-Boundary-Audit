# =============================================================================
# run_all.R  -- reproduce Paper 3 (v39) figures & supplementary tables in R.
# Run from repo root:   Rscript run_all.R
#
# FULLY REPRODUCED from bundled data (verified against the manuscript):
#   Figure 2, Figure 5, Figure 6, Supplementary Fig S1; Supp Tables S1, S2, S3.
# ANALYSIS CHECK:
#   PTGDS_vertex_check.R  (reconciles PTGDS vertex 0.47 vs bin-level 0.56)
# REQUIRE ADDITIONAL DATA (scaffolds; see each script header for the CSV needed):
#   Figure 3, Figure 4, Supplementary Figs S2, S3.
# =============================================================================
src <- function(p) { message("\n=== ", p, " ==="); tryCatch(source(p),
  error = function(e) message("  [skipped] ", conditionMessage(e))) }

# compute Figure 5b adjusted-Cox HRs first (writes data/fig5b_prognosis_computed.csv; skipped if ADNI data absent)
src("R/analysis/Figure5b_prognosis_cox.R")
src("R/analysis/celltype_generalization_audit.R")
src("R/figures/Figure2_consensus.R")
src("R/figures/Figure5_brainCSF_prognosis.R")
src("R/figures/Figure6_driver_signatures.R")
src("R/figures/SuppFig_celltype_generalization.R")
src("R/figures/SuppFig_global_correction_control.R")
src("R/figures/SuppFigS1_pancellular.R")
src("R/figures/SuppFig_EC_nonreplication.R")
src("R/tables/make_supp_tables.R")
message("\n=== R/analysis/PTGDS_vertex_check.R ==="); source("R/analysis/PTGDS_vertex_check.R")
message("\n=== R/analysis/CSF_TMT_detection_audit.R ==="); source("R/analysis/CSF_TMT_detection_audit.R")
# scaffolds (run if their input CSVs are present)
src("R/figures/Figure3_generalization.R")
src("R/figures/Figure3b_CSF.R")
src("R/figures/Figure3_combined.R")   # combined portrait Figure 3 (a + b) embedded in the manuscript body
src("R/figures/Figure4_metabolic_artifact.R")
src("R/figures/SuppFigS2_GWAS.R")
src("R/figures/SuppFigS3_pseudobulk.R")
message("\nDone. Outputs in output/figures and output/tables.")
if (file.exists("R/analysis/external_marker_crosscohort.R")) try(source("R/analysis/external_marker_crosscohort.R"))
if (file.exists("R/analysis/external_marker_panel_audit.R")) try(source("R/analysis/external_marker_panel_audit.R"))
if (file.exists("data-external/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad")) try(source("R/analysis/external_marker_audit.R"))
if (file.exists("output/tables/external_markers/sea_ad_panel_verdict.csv")) try(source("R/analysis/panel_finalize.R"))
