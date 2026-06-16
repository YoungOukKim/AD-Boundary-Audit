# =============================================================================
# Figure3b_CSF.R
# Figure 3b. CSF cross-modality test (ADNI SomaScan 7K vs baseline MMSE).
#   (left)  HMOX1 and GFAP CSF abundance (z) across baseline-MMSE bins.
#   (right) Spearman rho of HMOX1/GFAP/FKBP5 vs MMSE; HMOX1 absent in TMT-MS.
# Derived from ADNI SomaScan (HMOX1 SeqId 17398-55) + MMSE; raw ADNI data are
# controlled-access and NOT redistributed, so only the derived summaries ship:
#   csf_marker_rho.csv, csf_mmse_binned.csv
# (recompute: log2 RFU, per-protein z, Spearman vs earliest MMSE; n~708)
# Output: output/figures/Figure3_panelB.png
# =============================================================================
source("R/figures/utils.R")
bn <- read.csv(file.path(DATA, "csf_mmse_binned.csv"), check.names = FALSE)
rho <- read.csv(file.path(DATA, "csf_marker_rho.csv"), check.names = FALSE)
names(rho)[1] <- "marker"

# (left) trajectories across MMSE bins (worse ->)
bn$bin <- factor(bn$bin, levels = c("29-30","27-28","24-26","21-23","<=20"))
bl <- pivot_longer(bn[, c("bin","HMOX1","GFAP")], -bin, names_to = "marker", values_to = "z")
pL <- ggplot(bl, aes(bin, z, colour = marker, group = marker)) +
  geom_hline(yintercept = 0, colour = "grey85") +
  geom_line(linewidth = 1) + geom_point(size = 2.2) +
  scale_colour_manual(values = c(HMOX1 = UP, GFAP = DOWN), name = NULL) +
  labs(title = "CSF cross-modality test", x = "baseline MMSE (worse \u2192)",
       y = "CSF abundance (z, SomaScan)") +
  theme_paper + theme(legend.position = c(0.18, 0.85))

# (right) Spearman rho bars
rho$marker <- factor(rho$marker, levels = c("HMOX1","GFAP","FKBP5"))
rho$star <- ifelse(rho$P < 1e-3, "***", ifelse(rho$P < 1e-2, "**", ifelse(rho$P < 0.05, "*", "n.s.")))
rho$col <- c(HMOX1 = UP, GFAP = DOWN, FKBP5 = GREY)[as.character(rho$marker)]
pR <- ggplot(rho, aes(marker, rho, fill = col)) +
  geom_col(width = 0.62) + geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_text(aes(label = star, vjust = ifelse(rho > 0, -0.4, 1.2)), size = 4) +
  scale_fill_identity() +
  labs(title = "\u03c1 (SomaScan); HMOX1 absent in TMT-MS",
       x = NULL, y = "Spearman \u03c1 vs MMSE") + theme_paper

pL <- pL + labs(tag = "b") + theme(plot.tag = element_text(face = "bold", size = 13.8), plot.tag.position = c(0.03, 0.96))
save_fig((pL | pR), "Figure3_panelB.png", w = 11, h = 4.56)
