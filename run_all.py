#!/usr/bin/env python3
"""Regenerate every table and every figure of the boundary-audit study.

    python run_all.py              tables, then figures
    python run_all.py --tables     tables only
    python run_all.py --figures    figures only, from the tables already in
                                   analysis_py/data
    python run_all.py --korean     Korean-language figures as well

Each table in analysis_py/data is written by exactly one script, and the order
below is the order in which the scripts may be run. This is not housekeeping:
an earlier version of this pipeline had two scripts writing one file name under
different normalisation conventions, which is the same class of defect this
study reports, so the constraint is now enforced by the layout rather than by a
reader noticing it.

Inputs live in analysis_py/inputs and are donor-level tables derived from the
SEA-AD release; no controlled-access file is redistributed here. Outputs are
written to analysis_py/data and analysis_py/figures.

Table owners
    detection.py    detection_astrocyte, hmox1_donors, bin1_donors
    final.py        observed_components, reconstruction_fragility, nulls,
                    restriction, independence, independence_distinct,
                    leave_one_donor, bootstrap, bin_composition
    conditions.py   boundary_conditions
    staging.py      pathology_axes, pathology_peraxis, pathology_engine
    benchmark.py    benchmark, benchmark_meta
    region_grid.py  region_grid_bias
"""
import pathlib
import subprocess
import sys
import time

HERE = pathlib.Path(__file__).resolve().parent / "analysis_py"
TABLES = ["detection.py", "final.py", "conditions.py", "staging.py",
          "benchmark.py", "region_grid.py"]
FIGURES = ["fig1.py", "fig2.py", "fig3.py", "fig4.py", "fig5.py"]
KOREAN = ["make_ko_figs.py"]


def run(script):
    t0 = time.time()
    print(f"\n=== {script} " + "=" * max(0, 60 - len(script)), flush=True)
    r = subprocess.run([sys.executable, str(HERE / script)], cwd=HERE)
    if r.returncode:
        sys.exit(f"{script} failed with exit code {r.returncode}")
    print(f"--- {script} ok in {time.time() - t0:.1f}s", flush=True)


if __name__ == "__main__":
    args = sys.argv[1:]
    scripts = []
    if "--figures" not in args:
        scripts += TABLES
    if "--tables" not in args:
        scripts += FIGURES
    if "--korean" in args:
        scripts += KOREAN
    for s in scripts:
        run(s)
    print("\nall done")
