# Auditing molecular transition boundaries on binned disease axes

![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)
![Python](https://img.shields.io/badge/Python-3.13-3776AB)
![Status](https://img.shields.io/badge/status-manuscript%20in%20submission-orange)

Code, derived data and figures for a five-gate audit of molecular transition
boundaries detected on binned disease-progression axes, benchmarked against
ground truth and applied to a boundary in SEA-AD astrocytes that **we reported
ourselves**.

The study withdraws its own earlier specificity claim. A fixed 44-marker
astrocyte panel in SEA-AD middle temporal gyrus returns a transition at CPS
0.207. Audited here, that boundary survives dropping the 23 markers the assay
cannot adjudicate, but one donor of 84 moves it past the donor gate, and 753 of
1,000 expression-matched random panels land within 0.05 CPS of it. The published
value below 0.001 arose because the observation was a clustering boundary while
the null was accumulated from a principal-component breakpoint; scored like for
like the value is 0.753. The cause is structural: reading a boundary at the
first change of cluster label puts it on the first midpoint for two-thirds of
random nine-bin labellings, and for every draw on this cohort's grid.

The same audit run on neuropathological staging returns a different position, so
it is not a procedure that rejects everything.

## Authors

YoungOuk Kim<sup>1</sup>† (corresponding), WooMyung Heo<sup>1</sup>†,
Se Jin Park<sup>2</sup>†, YoungChul Kim<sup>1</sup>, Ye Eun Cho<sup>2</sup>

<sup>1</sup> BioXP Research Institute, Donghae, Gangwon-do, Republic of Korea
<sup>2</sup> Department of Food Biotechnology and Environmental Science, Kangwon
National University, Chuncheon 24341, Republic of Korea

† These authors contributed equally.
Corresponding author: YoungOuk Kim — yo.kim@bioxp.biz

## Quick start

```bash
python run_all.py              # every table, then every figure
python run_all.py --tables     # tables only
python run_all.py --figures    # figures only, from the tables already present
python run_all.py --korean     # also the Korean-language figures
```

Requires Python ≥ 3.13 with numpy, scipy, pandas, scikit-learn and matplotlib.
The versions this study was run under are pinned in `requirements.txt`. Nothing
has to be downloaded: every input the pipeline reads is in `analysis_py/inputs`.

## One table, one writer

Each table in `analysis_py/data` is written by exactly one script, and
`run_all.py` runs them in the only order in which they may be run. This is not
housekeeping. An earlier version of this pipeline had two scripts writing one
file name under different normalisation conventions, and running them in the
wrong order silently changed a reported result — the same class of defect this
study reports. The constraint is now enforced by the layout rather than by a
reader noticing it.

| Script | Writes |
|---|---|
| `detection.py` | `detection_astrocyte`, `hmox1_donors`, `bin1_donors` |
| `final.py` | `observed_components`, `reconstruction_fragility`, `nulls`, `restriction`, `independence`, `independence_distinct`, `leave_one_donor`, `bootstrap`, `bin_composition` |
| `conditions.py` | `boundary_conditions` |
| `staging.py` | `pathology_axes`, `pathology_peraxis`, `pathology_engine`, `pathology_null` |
| `benchmark.py` | `benchmark`, `benchmark_meta` |
| `region_grid.py` | `region_grid_bias` |

`engine.py` is the boundary-detection engine: a Python transcription of the R
engine used in the earlier study, validated against the committed outputs of
that study before it was used here. `paths.py` resolves every path from the
repository root, so the pipeline runs wherever it is cloned.

## Repository layout

```
AD-Boundary-Audit/
├── run_all.py              this study: entry point
├── analysis_py/            this study: code
│   ├── engine.py           19-component boundary-detection engine
│   ├── detection.py  final.py  conditions.py  staging.py
│   ├── benchmark.py  region_grid.py
│   ├── fig1.py … fig5.py   figures; style.py, make_ko_figs.py
│   ├── inputs/             donor-level inputs (read, never written)
│   ├── data/               tables the pipeline writes
│   └── figures/            figures the pipeline writes
│
├── run_all.R               earlier study: entry point
├── R/                      earlier study: code
├── data/  output/          earlier study: data and committed outputs
└── FIGURE_MAPPING.txt      earlier study: figure ↔ script map
```

## Inputs and what is not redistributed

`analysis_py/inputs` holds donor-level tables derived from the SEA-AD middle
temporal gyrus snRNA-seq release. The release itself is controlled-access and is
**not** redistributed here; obtain it from [sea-ad.org](https://sea-ad.org) or
the Sage Bionetworks AD Knowledge Portal
(`SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad`).

| File | What it is |
|---|---|
| `donor_mean_Astrocyte.csv` | per-donor mean log-normalised expression, astrocytes, 84 donors with a CPS |
| `donor_sumsq_Astrocyte.csv` | matching sums of squares |
| `detection_by_donor_Astrocyte.csv` | per-donor, per-gene detection rate |
| `detection_controls.csv` | the same detection audit across cell types, for the ambient band |
| `donor_meta.csv` | donor metadata: CPS, staging, covariates |
| `cps.xlsx` | the three published CPS fits compared in the axis-uncertainty analysis |
| `SEAAD_44panel_bin_means.csv` | the 44-marker bin matrix as deposited by the earlier study, used as the reference the reconstruction is checked against. Byte-identical to `data/SEAAD_44panel_bin_means.csv`; kept here under the same name so the pipeline runs standalone and so one file does not travel under two names |

The single extraction pass that produces the first five from the h5ad is
documented in the manuscript's Methods. It needs the controlled-access object
and is therefore not part of `run_all.py`.

## Reproducibility notes

- All random procedures use seed 42; the engine's k-means uses ten restarts.
- The `segmented` R package changes three components' values depending on
  whether it is installed. Those components are excluded from every analysis
  reported here, and the Python engine reproduces the remaining 16 against the
  earlier study's committed outputs.
- `run_all.py` regenerates every table in `analysis_py/data` and every figure in
  `analysis_py/figures` from `analysis_py/inputs` alone.

## The earlier study

The four-axis reproducibility audit that this study audits is preserved intact:
entry point `run_all.R`, code in `R/`, data in `data/`, committed outputs in
`output/`, and the state as deposited tagged `v1.0-four-axis-audit`.

Readers should note that the present study withdraws the specificity claim made
for the CPS 0.207 boundary in that work, and reports that its nine-algorithm
consensus is substantially redundant: across 300 matched random panels the 16
branch-independent components return a median of 10 distinct values per panel,
and Ward and k-means return the same first boundary in 83.3% of panels.

## License

MIT for code; shared derived data CC BY 4.0. See `LICENSE`.
