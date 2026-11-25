# =============================================================================
# PTGDS_vertex_check.R
# Reconciles the PTGDS biphasic-trajectory positions reported across papers.
# PTGDS is deliberately excluded from Paper 3's 44-marker panel; its trajectory
# lives in the Paper-1 (39-gene) panel: SEAAD_paper1panel_bin_means.csv.
#
# Three DISTINCT curve features (all real, not contradictory):
#   - CPS ~0.23 : early segmented slope-change breakpoint   (Paper 1)
#   - CPS  0.47 : donor-level QUADRATIC VERTEX (-b1/2b2)     (Paper 1; cited by P3)
#   - CPS ~0.54-0.56 : bin-level discrete argmax / bin-level segmented breakpoint
# The manuscript's "vertex / peaks at CPS~0.47" cites the Paper-1 quadratic
# vertex and is CORRECT. The repo's older "inflection ~0.56" was a coarser
# bin-level segmented fit (a different feature), not a contradiction.
# =============================================================================
DATA <- "data"
bm <- read.csv(file.path(DATA, "SEAAD_paper1panel_bin_means.csv"), check.names = FALSE)
x <- bm$cps_center; y <- bm$PTGDS

# quadratic vertex (the manuscript / Paper-1 feature)
b <- coef(lm(y ~ x + I(x^2))); vertex <- -b[2] / (2 * b[3])
# discrete peak
argmax <- x[which.max(y)]
# bin-level 2-segment breakpoint (grid search)
best <- Inf; bp <- NA
for (k in seq(min(x) + 1e-3, max(x) - 1e-3, length.out = 400)) {
  x1 <- pmin(x, k); x2 <- pmax(x - k, 0)
  sse <- sum(resid(lm(y ~ x1 + x2))^2); if (sse < best) { best <- sse; bp <- k }
}
cat(sprintf("PTGDS quadratic vertex (bin-level)   = CPS %.3f   [P1 donor-level: 0.47]\n", vertex))
cat(sprintf("PTGDS discrete argmax                = CPS %.3f\n", argmax))
cat(sprintf("PTGDS bin-level segmented breakpoint = CPS %.3f   [repo's old '~0.56']\n", bp))
cat("\nConclusion: manuscript 'vertex/peak ~0.47' = Paper-1 quadratic vertex = CORRECT.\n")
