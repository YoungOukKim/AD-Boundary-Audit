Figure1.png in this directory is NOT the Figure 1 of the current manuscript.

The figures in output/ were produced by the R scripts in R/figures/ for the
earlier version of this study, whose claim was that a multi-algorithm consensus
localizes a reproducible astrocyte transition. The current manuscript withdraws
that claim, and its Figures 1-5 are produced by a separate Python pipeline
(fig1.py ... fig5.py, with analysis.py, benchmark.py, conditions.py, final.py
and style.py) that is not yet deposited here. The five PNGs in this directory
and the five figures in the manuscript therefore differ.

An earlier version of this file stated that Figure 1 was built from
R/assets/Figure1_hero.svg. That is incorrect: the SVG in R/assets/ is an older
schematic using the earlier Class I / Class II vocabulary and cannot produce the
five-gate figure the manuscript uses.

Before deposition, the Python pipeline that regenerates the manuscript's figures
and numbers must be added to this repository, and these stale outputs removed or
clearly marked, so that the Data and code availability statement is true.
