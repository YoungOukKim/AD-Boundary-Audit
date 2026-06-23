# =============================================================================
# Figure3_combined.R
# Figure 3 (combined, portrait) for the manuscript body:
#   (a) Panel-wide boundary delta-z does not generalize MTG -> A9  (44 markers)
#   (b) CSF cross-modality test: (left) HMOX1/GFAP across MMSE bins,
#                                (right) Spearman rho vs MMSE
# Panel a is placed full-width on top with 7 pt marker labels so the 44-row
# bar stays legible (>=7 pt) at page size; panel b sits below. Panel logic is
# identical to Figure3_generalization.R (a) and Figure3b_CSF.R (b); only the
# layout, label sizes and canvas change.
# Output: output/figures/Figure3.png
# =============================================================================
source("R/figures/utils.R")

# ---- (a) MTG vs A9 dodged bar (44 markers) ----
d  <- read.csv(file.path(DATA, "fig3a_region_dz.csv"), check.names = FALSE)
dl <- pivot_longer(d, c("MTG","A9"), names_to = "region", values_to = "dz")
dl$region <- factor(dl$region, levels = c("MTG","A9"))
ord <- d$marker[order(d$MTG)]                      # HMOX1 bottom, FKBP5 top
dl$marker <- factor(dl$marker, levels = ord)
f3a <- ggplot(dl, aes(dz, marker, fill = region)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.72) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = c(MTG = "#5E4FA2", A9 = GREY), name = NULL) +
  scale_y_discrete(expand = c(0.012, 0.012)) +
  labs(title = "Boundary \u0394z does not generalize MTG \u2192 A9",
       x = "boundary \u0394z at CPS 0.21", y = NULL) +
  theme_paper +
  theme(axis.text.y = element_text(size = 7),
        legend.position = c(0.87, 0.06),
        legend.background = element_rect(fill = "white", colour = "grey80", linewidth = 0.3),
        legend.key.size = unit(0.42, "cm"))

# ---- (b) CSF cross-modality ----
bn <- read.csv(file.path(DATA, "csf_mmse_binned.csv"), check.names = FALSE)
rho <- read.csv(file.path(DATA, "csf_marker_rho.csv"), check.names = FALSE)
names(rho)[1] <- "marker"
bn$bin <- factor(bn$bin, levels = c("29-30","27-28","24-26","21-23","<=20"))
bl <- pivot_longer(bn[, c("bin","HMOX1","GFAP")], -bin, names_to = "marker", values_to = "z")
pL <- ggplot(bl, aes(bin, z, colour = marker, group = marker)) +
  geom_hline(yintercept = 0, colour = "grey85") +
  geom_line(linewidth = 1) + geom_point(size = 2.2) +
  scale_colour_manual(values = c(HMOX1 = UP, GFAP = DOWN), name = NULL) +
  labs(title = "CSF cross-modality test", x = "baseline MMSE (worse \u2192)",
       y = "CSF abundance (z, SomaScan)") +
  theme_paper + theme(legend.position = c(0.18, 0.85),
                      plot.title = element_text(size = 11))

rho$marker <- factor(rho$marker, levels = c("HMOX1","GFAP","FKBP5"))
rho$star <- ifelse(rho$P < 1e-3, "***", ifelse(rho$P < 1e-2, "**", ifelse(rho$P < 0.05, "*", "n.s.")))
rho$col <- c(HMOX1 = UP, GFAP = DOWN, FKBP5 = GREY)[as.character(rho$marker)]
pR <- ggplot(rho, aes(marker, rho, fill = col)) +
  geom_col(width = 0.62) + geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_text(aes(label = star, vjust = ifelse(rho > 0, -0.4, 1.2)), size = 4) +
  scale_fill_identity() +
  labs(title = "\u03c1 vs MMSE", subtitle = "HMOX1 absent in TMT-MS",
       x = NULL, y = "Spearman \u03c1 vs MMSE") +
  theme_paper + theme(plot.title = element_text(size = 11),
                      plot.subtitle = element_text(size = 8.5))

# ---- portrait assembly: a full-width on top, b (left|right) below ----
f3a <- f3a + labs(tag = "a") +
  theme(plot.tag = element_text(face = "bold", size = 13), plot.tag.position = c(0.015, 0.985))
pL  <- pL + labs(tag = "b") +
  theme(plot.tag = element_text(face = "bold", size = 13), plot.tag.position = c(0.02, 0.97))
fig3 <- f3a / (pL | pR) + plot_layout(heights = c(3.0, 1))
save_fig(fig3, "Figure3.png", w = 6.5, h = 8.8)
