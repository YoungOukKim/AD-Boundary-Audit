# The earlier four-axis reproducibility audit

This directory holds the code, derived data and committed outputs of the earlier
study, "A reproducibility-audit framework for generalizable versus
dataset-specific molecular transition boundaries in Alzheimer's disease",
exactly as they were deposited. Nothing here has been edited; the files were
moved into this directory so that the two studies' figures, which share the
names `Figure1.png` through `Figure5.png`, cannot be confused.

The repository state as originally deposited, with these files at the top level,
is tagged **`v1.0-four-axis-audit`**:

```bash
git checkout v1.0-four-axis-audit
```

## Running it

The scripts use paths relative to their own root, so run from *this* directory
rather than from the repository root:

```bash
cd previous-study
Rscript run_all.R
```

Requires R ≥ 4.3. Figures 2, 5, 6 and 7, Supplementary Figure S1 and
Supplementary Tables S1–S3 reproduce from the bundled data; Figures 3 and 4 and
some other analyses need controlled-access inputs and run as scaffolds when
those inputs are absent. `FIGURE_MAPPING.txt` maps each figure to its script.

## What the present study found about this one

The study in the parent directory audits the boundary reported here and
**withdraws the specificity claim** made for it. Two findings bear directly on
this deposit and a reader should know them before reusing it:

- The published specificity value of p < 0.001 for the CPS 0.207 boundary
  compares a Ward clustering boundary against a null accumulated from a
  principal-component breakpoint. Scored on the same statistic, 753 of 1,000
  expression-matched random panels reproduce the boundary.
- The nine-algorithm consensus is substantially redundant. Across 300 matched
  random panels the 16 branch-independent components return a median of 10
  distinct values per panel, and Ward and k-means return the same first boundary
  in 83.3% of panels.

A third point concerns the engine itself: three components derived from the R
`segmented` package return different values depending on whether that package
is installed, so the deposited tables here are not self-consistent across
installations for those three. The present study excludes them.

## Files that look duplicated and are not

Three pairs of byte-identical files in this directory are deliberate and should
be left alone:

| Files | Why they match |
|---|---|
| `data/fig5b_prognosis_computed.csv`, `data/fig5b_prognosis_manuscript.csv` | the computed Cox hazard ratios and the manuscript's values; that they are identical is the check, not an accident, and `R/figures/Figure5_brainCSF_prognosis.R` falls back to the second when the first is absent |
| `data/fig5b_prognosis_manuscript.csv`, `output/tables/Figure5b_cox_HR.csv` | input and the output written from it by `R/analysis/Figure5b_prognosis_cox.R` |
| `data/consensus_breakpoints.csv`, `output/tables/SuppTable_S2_consensus.csv` | input and the supplementary table `R/tables/make_supp_tables.R` writes from it |
