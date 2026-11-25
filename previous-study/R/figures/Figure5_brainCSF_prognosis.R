# =============================================================================
# Figure5_brainCSF_prognosis.R
# Figure 5. Brain-CSF discordance and absence of cognition-independent prognosis.
#   (a) HMOX1 collapses at the boundary in brain (MTG astrocytes) across CPS bins,
#       and the brain delta z vs CSF rho discordance.
#   (b) Cognition-threshold survival model (circular) vs molecular panel score
#       (null after adjustment).
# Inputs : SEAAD_44panel_bin_means.csv      (HMOX1 trajectory, panel a-left)
#          fig5b_prognosis_manuscript.csv   (manuscript forest values, panel b)
# NOTE   : the bundled cox_HR_results_new.csv (HR 2.70) is the EARLIER prognostic
#          cohort; the manuscript Figure 5b reports the adjusted null analysis
#          (MMSE-threshold 3.79; molecular 0.99 unadj / 0.92 adj; baseline-deduped Cox, n=1,105/527), provided in
#          fig5b_prognosis_manuscript.csv. CSF rho (-0.11) is from ADNI SomaScan.
# Output : output/figures/Figure5.png
# =============================================================================
source("R/figures/utils.R")

# ---- (a-left) HMOX1 mean z across CPS bins ----
bm <- read.csv(file.path(DATA, "SEAAD_44panel_bin_means.csv"), check.names = FALSE)
hm <- data.frame(cps = bm$cps_center, z = bm$HMOX1)
f4a1 <- ggplot(hm, aes(cps, z)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = 0.21, linetype = "dashed", colour = UP, linewidth = 0.6) +
  geom_line(colour = DOWN, linewidth = 1) +
  geom_point(colour = DOWN, size = 2.2) +
  labs(title = "HMOX1 collapses at the boundary",
       x = "CPS", y = "HMOX1 mean z") + theme_paper +
  labs(tag = "a") + theme(plot.tag = element_text(face = "bold", size = 13.7), plot.tag.position = c(0.02, 0.97))

# ---- (a-right) brain delta z vs CSF rho discordance ----
disc <- data.frame(
  comp = factor(c("Brain\nsnRNA \u0394z","CSF\nSomaScan \u03c1"),
                levels = c("CSF\nSomaScan \u03c1","Brain\nsnRNA \u0394z")),
  val  = c(-0.60, -0.11),
  col  = c(DOWN, GREY))
f4a2 <- ggplot(disc, aes(val, comp, fill = col)) +
  geom_col(width = 0.6) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  annotate("text", x = -0.32, y = 1, label = "TMT-MS: undetected", size = 2.6, colour = "grey30") +
  scale_fill_identity(guide = "none") +
  labs(title = "Brain\u2013CSF discordance", x = "effect with progression", y = NULL) +
  theme_paper

# ---- (b) prognosis forest (manuscript adjusted analysis) ----
fb_file <- if (file.exists(file.path(DATA, "fig5b_prognosis_computed.csv")))
             "fig5b_prognosis_computed.csv" else "fig5b_prognosis_manuscript.csv"
message("Figure 5b source: ", fb_file, "  (computed = real Cox; manuscript = placeholder)")
pf <- read.csv(file.path(DATA, fb_file), check.names = FALSE)
pf$label <- gsub("\\\\n", "\n", pf$label)
pf$label <- factor(pf$label, levels = rev(pf$label))
pal <- c(circular = UP, molecular = GREY, molecular_adj = col_mci)
f4b <- ggplot(pf, aes(HR, label, colour = group)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.18, linewidth = 0.7) +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("%.2f [%.2f-%.2f]", HR, CI_lower, CI_upper)),
            vjust = -1.1, size = 3, colour = "black") +
  scale_x_log10(breaks = c(1,2,3,4,5)) +
  scale_colour_manual(values = pal, guide = "none") +
  labs(title = "Boundaries carry no prognosis independent of cognition",
       x = "hazard ratio (95% CI)  \u2014  endpoint: \u22653-pt MMSE decline or MMSE \u226423",
       y = NULL) +
  theme_paper + theme(axis.text.y = element_text(size = 8)) +
  labs(tag = "b") + theme(plot.tag = element_text(face = "bold", size = 13.7), plot.tag.position = c(0.02, 0.97))

fig4 <- (f4a1 | f4a2) / f4b + plot_layout(heights = c(1, 1))
save_fig(fig4, "Figure5.png", w = 11, h = 8)

# optional: export panels separately (for manuscript 2-image layout)
if (Sys.getenv("EXPORT_PANELS") == "1") {
  ggsave(file.path(OUTF, "Figure5_panelA.png"), (f4a1 | f4a2), width = 11, height = 5.324, dpi = 300, bg = "white")
  ggsave(file.path(OUTF, "Figure5_panelB.png"), f4b,          width = 11, height = 4.078, dpi = 300, bg = "white")
  message(">>> exported Figure4_panelA/B")
}
