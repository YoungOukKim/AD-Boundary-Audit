Sensitivity of the boundary consensus to the SEA-AD specimen-metadata
correction published 2026-06-25 on the Sage Bionetworks AD Knowledge Portal.

  Run:  python3 analysis_py/correction_sensitivity_20260625.py
  Out:  output/tables/SuppTable_S7_correction_sensitivity.csv
        output/tables/SuppTable_S7_correction_sensitivity_readme.txt

WHY THIS IS PYTHON AND NOT R
----------------------------
The 19-component engine is re-implemented here and VERIFIED against the two
reference outputs already committed to this repository:
  * output/tables/celltype_audit/multimethod_Astrocyte.csv  -> 19/19 exact
  * output/tables/celltype_audit/multimethod_Microglia.csv   -> 19/19 exact
  * output/tables/SuppTable_S2_consensus.csv                 -> 16/19 exact
The three components that differ in the last comparison are Segmented_PC1..3
and the cause is documented below (.HAVE_SEG). An R re-implementation should be
added and cross-checked once an R installation is available.

KNOWN DEFECT FOUND WHILE WRITING THIS SCRIPT  (fix before submission)
--------------------------------------------------------------------
boundary_detection.R branches on `.HAVE_SEG <- requireNamespace("segmented")`.
With the package present, Segmented_PC1/2/3 = 0.665898 / 0.654355 / 0.416294;
without it, the 200-point fallback grid gives 0.835621 / 0.280384 / 0.243573.
That moves 3 of 19 components and changes EVERY region count reported in
Supp Table S2:
    Primary(CPS~0.21)      4 -> 5
    Substructure           4 -> 5
    Secondary(0.39-0.42)   2 -> 1
    Late split(0.54-0.67)  9 -> 7
    (and one component falls beyond the late split)
Supp Table S2 was generated WITH the package; the cell-type generalization
outputs were generated WITHOUT it. `segmented` must be declared a pinned hard
dependency (or the fallback removed) so the table is deterministic.
None of the three affected components lies in the Primary group under the
with-package branch, so the conclusions below are branch-independent; they are
reported separately for the 16 branch-independent components.

WHAT THE CORRECTION ACTUALLY SAYS
---------------------------------
Columns in the change log are "correct" then "org." (the original, wrong value).
Three MTG rows exist; only two can touch an MTG RNA-seq object:

  M2AT_220121_206_D01  donor        correct H21.33.047  org. H18.30.002  M/M  MTG/MTG
      -> M2AT = ATAC modality. Cannot affect an RNA-seq analysis. Its original
         (wrong) donor H18.30.002 is a Reference-series donor and is absent
         from the 84-donor cohort in any case.
  M2TX_210727_207_E01  donor;tissue correct H20.33.016/DFC  org. H21.33.037/MTG  F/F
  M1TX_210723_159_G01  donor;tissue correct H21.33.037/MTG  org. H20.33.016/DFC  F/F
      -> a clean two-way swap between two TRANSCRIPTOMIC specimens.

Consequence for a pre-correction MTG RNA object: the nuclei labelled
H21.33.037 are in fact the DFC tissue of donor H20.33.016, whose own CPS is
0.616531 (bin [0.6,0.7)); and H21.33.037's true MTG specimen was filed under
DFC and therefore excluded. Both specimens are female, so there is no sex
confound. CPS is a donor-level property, so the substituted specimen re-enters
bin 1 at the same axis position, CPS 0.180319.

STRUCTURE OF THE TEST
---------------------
Bins are FIXED CPS intervals, not donor quantiles, so removing a donor cannot
move a bin edge or move any other donor between bins; exactly one row of the
9 x 44 matrix changes. That row changes in two separable ways:

 (A) LOCATION -- exact, needs no expression data.
     centre' = (centre - w*CPS_d)/(1-w); every detector reports a bin centre or
     a midpoint of two adjacent centres, so the shift in the reported
     breakpoint is exact.
 (B) FRAGILITY -- the minimum perturbation radius r* over ALL directions in
     R^44 that moves the Ward primary call, inverted into the donor deviation
     ||mu_bin1 - mu_contaminant|| it would require.
 (C) MECHANISM -- the contaminant is a mid/late-trajectory sample sitting in
     the lowest bin, so removing it must move bin 1 AWAY from the bulk, whereas
     every breaking direction points TOWARD the bulk. Verified empirically by
     using each of the 8 other bins in turn as the contaminant's profile proxy.
