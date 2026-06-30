# A reproducibility-audit framework for generalizable versus dataset-specific molecular transition boundaries in Alzheimer's disease

![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg) ![Status](https://img.shields.io/badge/status-preprint%20in%20deposition-orange) ![R](https://img.shields.io/badge/R-%E2%89%A5%204.3-276DC3)

Analysis code, derived data and figure/table outputs for a permutation-controlled, four-axis reproducibility audit of trajectory-based molecular staging in Alzheimer's disease. The framework asks whether a molecular transition boundary localized in one brain region behaves as an intrinsic, brain-wide event, and adjudicates candidate anchor markers as generalizable (Class II) or dataset-specific (Class I).

The manuscript is the source of truth; bundled CSVs reproduce the manuscript figures/tables in R (ggplot2). `Rscript run_all.R` from the repo root regenerates everything reproducible from bundled data; analyses requiring controlled-access inputs are run as scaffolds when those inputs are present (see Reproducibility).

## Authors

YoungOuk Kim¹† (corresponding), WooMyung Heo¹†, Se Jin Park²†, YoungChul Kim¹, Ye Eun Cho²

¹ BioXP Research Institute, Donghae, Gangwon-do, Republic of Korea
² Department of Food Biotechnology and Environmental Science, Kangwon National University, Chuncheon 24341, Republic of Korea

† These authors contributed equally to this work.
Corresponding author: YoungOuk Kim — yo.kim@bioxp.biz (BioXP Research Institute)

## Abstract

Molecular staging of Alzheimer's disease (AD) increasingly defines transition boundaries along single-cell pseudo-progression trajectories, yet whether such boundaries reproduce across brain regions, cohorts and molecular modalities is largely untested. We developed a permutation-controlled framework that combines nine boundary-detection algorithms with a fixed reactive-astrocyte marker panel and applied it across independent datasets. In the Seattle Alzheimer's Disease Brain Cell Atlas middle temporal gyrus, it localized a transition that was robust across algorithms and recovered in most cell types, but that did not generalize: its leading marker was attenuated or absent in prefrontal cortex, entorhinal cortex and cerebrospinal fluid, and an apparent cross-region conservation of glial metabolic genes reflected a global-expression offset rather than a shared program. The same framework nonetheless certified an externally validated marker (astrocytic PTGDS) as reproducible across regions and modalities, showing that it distinguishes generalizable anchors from dataset-specific ones rather than rejecting all signals. We provide this four-axis audit—algorithmic consensus, region, cohort and modality—as a transferable, code-available control for trajectory-based staging in AD and other progressive proteinopathies, and suggest such auditing before a boundary is read as a biological stage.

## Repository structure

```
AD-Boundary-Audit/
├── README.md                      # this file (GitHub landing page)
├── FIGURE_MAPPING.txt             # authoritative figure/Supp ↔ script map (manuscript numbering)
├── LICENSE                        # MIT (code); shared data CC BY 4.0
├── .gitignore
├── run_all.R                      # end-to-end orchestrator (Rscript run_all.R from repo root)
│
├── R/
│   ├── analysis/                  # boundary engine + audits
│   │   ├── boundary_detection.R           # standalone 9-algorithm / 19-component engine (set.seed 42)
│   │   ├── synthetic_benchmark.R          # Figure 7: ground-truth benchmark (self-contained synthetic data)
│   │   ├── celltype_audit_from_raw.R      # regenerate bin-means from SEA-AD h5ad (controlled input)
│   │   ├── celltype_generalization_audit.R# boundary calls from bundled bin-means
│   │   ├── verify_correction_consistency.R# HK-flat & metabolic conservation under genome-wide correction
│   │   ├── external_marker_audit.R / _generic.R / _crosscohort.R / _panel_audit.R  # marker adjudication
│   │   ├── panel_finalize.R               # merge per-cell + cross-cohort verdicts
│   │   ├── PTGDS_vertex_check.R           # reconciles PTGDS vertex 0.47 vs bin-level segmented 0.56
│   │   ├── CSF_TMT_detection_audit.R
│   │   └── Figure5b_prognosis_cox.R       # adjusted Cox HRs (needs ADNI; writes fig5b CSV)
│   ├── figures/                   # one script per figure / Supp figure (+ utils.R)
│   ├── tables/                    # make_supp_tables.R
│   └── assets/                    # Figure1_hero.svg (schematic; not R-generated)
│
├── data/                          # bundled derived CSVs (no raw / controlled data)
│   ├── celltype_audit/            # binmeans_{Astrocyte,Microglia}.csv (+ README)
│   └── metabolic_correction/      # panel_correction_z.csv
├── data-external/                 # controlled-access inputs — NOT redistributed (README only)
│
└── output/                        # committed reference outputs
    ├── figures/                   # Figure2–7 (Figure7_benchmark.png), Supp figures (.png/.tif); _manuscript_reference/
    └── tables/                    # Table1–3, Supp Tables S1–S3, external_markers/, celltype_audit/, benchmark/ (synthetic: partA/global/lomo.csv + summary.txt)
```

`FIGURE_MAPPING.txt` is the authoritative figure↔script map and follows the final manuscript numbering.

Filenames are intentionally not renamed to avoid breaking hard-coded paths. Note: `data/drivers_master_event_corrected.csv` is the primary CPS ≈ 0.21 transition ("master event" was the pre-revision label; the manuscript text uses "primary").

## Figure → script

| Manuscript | Content | Script → output |
|---|---|---|
| Fig. 1 | Reproducibility-audit framework (schematic) | `R/assets/Figure1_hero.svg` (not R-generated) |
| Fig. 2 | Multi-algorithm consensus transition in MTG | `R/figures/Figure2_consensus.R` → `Figure2.png` |
| Fig. 3 | Region/cohort/modality non-generalization | `R/figures/Figure3_combined.R` (components: `Figure3_generalization.R` + `Figure3b_CSF.R`) → `Figure3.png` |
| Fig. 4 | Metabolic conservation is a global-expression artefact | `R/figures/Figure4_metabolic_artifact.R` → `Figure4.png` |
| Fig. 5 | Brain–CSF discordance; no cognition-independent prognosis | `R/figures/Figure5_brainCSF_prognosis.R` (Cox: `R/analysis/Figure5b_prognosis_cox.R`) → `Figure5.png` |
| Fig. 6 | Driver signatures and within-MTG robustness | `R/figures/Figure6_driver_signatures.R` → `Figure6.png` |
| Fig. 7 | Synthetic ground-truth benchmark (operating characteristics) | `R/analysis/synthetic_benchmark.R` → `output/figures/Figure7_benchmark.png` |

## Supplementary → script

| Manuscript | Content | Script |
|---|---|---|
| Supp. Fig. S1 | CPS-0.21 transition recovered across MTG cell types | `R/figures/SuppFigS1_pancellular.R` |
| Supp. Fig. S2 | AD-GWAS convergence on the CPS-0.207 boundary (microglia) | `R/figures/SuppFigS2_GWAS.R` |
| Supp. Fig. S3 | Housekeeping negative control for the global-expression correction | `R/figures/SuppFigS3_global_correction_control.R` |
| Supp. Fig. S4 | Audit generalizes across cell types (microglia panel) | `R/figures/SuppFigS4_celltype_generalization.R` |
| Supp. Fig. S5 | Donor-level (pseudobulk) testing is required | `R/figures/SuppFigS5_pseudobulk.R` |
| Supp. Fig. (EC) | Entorhinal non-replication (supports Fig. 3) | `R/figures/SuppFig_EC_nonreplication.R` |
| Supp. Tables S1–S3 | Panel composition / consensus breakpoints / drivers | `R/tables/make_supp_tables.R` |

## Data availability

| Dataset | Access | Link |
|---|---|---|
| SEA-AD snRNA-seq (MTG, DLPFC; 1.3M nuclei) | Public | https://portal.brain-map.org |
| Entorhinal cortex snRNA-seq (Leng et al. 2021) | Public (GEO) | GSE147528 |
| ADNI CSF proteomics (TMT-MS / SomaScan) | Registration + DUA | https://adni.loni.usc.edu |
| Emory CSF TMT-MS columns | Header only (bundled) | `data/EMORY_CSF_TMT_MS_columns.txt` |

Controlled-access inputs live in `data-external/` and are not redistributed (see `data-external/README.txt`): the SEA-AD MTG `.h5ad`, ADNI `.rda` tables, and Emory CSF TMT-MS subject-level matrix. Only derived, de-identified summaries needed to reproduce figures are bundled under `data/` (e.g. `csf_marker_rho.csv`, `csf_mmse_binned.csv`, bin-means). Figure 2/3 CSF ρ were computed from ADNI SomaScan (HMOX1 SeqId 17398-55) + MMSE; raw ADNI matrices are not re-hosted. Users must obtain restricted data independently.

## Key findings

| Finding | Value | Source |
|---|---|---|
| Primary reactive-astrocyte transition (MTG astrocytes) | CPS 0.207; permutation 0/1000 random 44-gene panels | `data/consensus_breakpoints.csv` |
| Cross-cell-type recovery (MTG) | 7 of 8 cell types | Supp. Fig. S1 |
| Region test (HMOX1, Class I) | MTG Δz −0.69 → A9 ≈ −0.06 (1/8 cell types) | Table 1 |
| Cohort test (Leng EC, n = 10) | panel reversed; FTH1 +0.54 (opposite to MTG) | Table 1; `data/testA/B_edgeR_*_EC.csv` |
| Modality test (HMOX1 in CSF) | undetected (TMT-MS) / ρ = −0.11 (SomaScan) | Table 1; `data/csf_marker_rho.csv` |
| Reactive marker GFAP (re-anchor) | CSF ρ = −0.17, P = 4×10⁻⁶ | `data/csf_marker_rho.csv` |
| PTGDS positive control (Class II) | reproducible across regions/species/proteomics | external validation (companion study) |
| AD-GWAS convergence | microglia P = 0.029; astrocytes P = 0.22 | `data/P3_gwas_convergence_results.csv` |
| Metabolic "conservation" | artefact: only 8 of 21 genes same-sign after global-expression correction (13 divergent) | `data/fig4_metabolic.csv` |
| Prognosis | no value beyond cognition + amyloid–tau (adjusted Cox) | `output/tables/Figure5b_cox_HR.csv` |
| Audit operating characteristics (synthetic ground truth) | sensitivity 100%, specificity 94%; each artefact rejected primarily by a distinct axis | `R/analysis/synthetic_benchmark.R` → `output/tables/benchmark/summary.txt` |

The 44-marker panel deliberately excludes PTGDS, LCN2 and MAPT for independent boundary detection; PTGDS enters only as an externally validated positive control. The boundaries are cross-sectional staging constructs and carry no prognostic value beyond cognition and amyloid–tau status.

## Reproducibility

R ≥ 4.3. Figure packages: `ggplot2`, `dplyr`, `tidyr`, `scales`, `patchwork`, `ggrepel`. Analysis extras: `segmented`, `cluster`, `edgeR` (EC pseudobulk). See script headers.

Deterministic consensus engine (`set.seed(42)`); `segmented` uses `n.boot = 0` (record the installed version).

```bash
# from repo root
Rscript run_all.R
```

Reproduced from bundled data (verified against the manuscript): Figures 2, 5, 6; Supp Fig S1; Supp Tables S1–S3, and the EC non-replication that supports Fig. 3. Require additional/controlled inputs (run as scaffolds when present): Figures 3, 4 and Supp Figs S2, S3 input regeneration, plus the from-raw audits in `R/analysis/` that read `data-external/`. Outputs are written to `output/figures` and `output/tables`. Figure 7 (`R/analysis/synthetic_benchmark.R`) is fully self-contained—it generates its own synthetic data with `set.seed(42)` and therefore reproduces the benchmark (sensitivity 100%, specificity 94%; each artefact class rejected primarily by a distinct axis) on any machine without external inputs, writing `output/figures/Figure7_benchmark.png` and `output/tables/benchmark/` (~25 min single-core).

Prognosis values are final and computed. `output/tables/Figure5b_cox_HR.csv`, `data/fig5b_prognosis_computed.csv` and `data/fig5b_prognosis_manuscript.csv` all carry the confirmed baseline-deduped Cox output (n = 1,105 / 527 events; MMSE-threshold HR 3.79, molecular 0.99 unadjusted / 0.92 adjusted), which matches the manuscript Results, Methods and Fig. 5b caption. The manuscript fit was performed in Python (statsmodels `PHReg`); `R/analysis/Figure5b_prognosis_cox.R` is an equivalent R/`survival` re-implementation provided for full in-repo reproducibility. PTGDS vertex 0.47 (donor-level quadratic) vs bin-level segmented 0.56 is reconciled by `R/analysis/PTGDS_vertex_check.R`.

## Citation

Kim YO, Heo WM, Park SJ, Kim YC, Cho YE (2026). *A reproducibility-audit framework for generalizable versus dataset-specific molecular transition boundaries in Alzheimer's disease.* Preprint/DOI to be added upon deposition.

## Competing interests

W.M.H. is Chief Executive Officer of BioXP, Inc., and Y.O.K. is Director of its Research Institute; BioXP, Inc. provided research support for this study. W.M.H., Y.O.K. and Y.C.K. are affiliated with BioXP, Inc. and therefore declare a competing financial interest. BioXP, Inc. has filed a patent application (Republic of Korea, filed 17 April 2026), with W.M.H., Y.O.K. and Y.C.K. as inventors, relating to the use of PTGDS and LCN2 as biomarkers for the diagnosis of Alzheimer's disease; both markers are evaluated in this study. S.J.P. and Y.E.C. declare no competing interests.

## License

Code is released under the MIT License. Openly shared derived-data files are released under CC BY 4.0. Restricted source datasets (SEA-AD, ADNI) remain governed by their own data-use agreements.

## Contact

YoungOuk Kim (corresponding author) — BioXP Research Institute, Republic of Korea — yo.kim@bioxp.biz
