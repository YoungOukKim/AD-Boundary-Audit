Cell-type generalization audit (manuscript Discussion).

binmeans_{Astrocyte,Microglia}.csv  9 CPS bins x 44 markers (per-marker z-score
  bin means). Astrocyte = reactive-astrocyte 44-panel; Microglia = independent
  microglia-specific 44-panel. Generated from the SEA-AD MTG h5ad by
  R/analysis/celltype_audit_from_raw.R (controlled-access input, data-external/).

Reproduce the boundary calls from these bundled files (no external data):
  Rscript R/analysis/celltype_generalization_audit.R
  -> output/tables/celltype_audit/multimethod_{tag}.csv + generalization_summary.txt
  -> identical 9-algorithm audit localizes a concordant early transition
     (Ward primary CPS ~0.20) in BOTH panels => framework is panel-/cell-type-agnostic.

Committed reference outputs are in output/tables/celltype_audit/.
