# =============================================================================
# Figure3_generalization.R  (panel a; panel b = Figure3b_CSF.R)
# Figure 3a. Panel-wide boundary delta z does not generalize MTG -> A9.
# Input : fig3a_region_dz.csv  (marker, MTG, A9)
# Output: output/figures/Figure3_panelA.png
# =============================================================================
source("R/figures/utils.R")
d <- read.csv(file.path(DATA, "fig3a_region_dz.csv"), check.names = FALSE)
dl <- pivot_longer(d, c("MTG","A9"), names_to = "region", values_to = "dz")
dl$region <- factor(dl$region, levels = c("MTG","A9"))
ord <- d$marker[order(d$MTG)]                     # HMOX1 bottom, FKBP5 top
dl$marker <- factor(dl$marker, levels = ord)
f2a <- ggplot(dl, aes(dz, marker, fill = region)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.72) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = c(MTG = "#5E4FA2", A9 = GREY), name = NULL) +
  labs(title = "Boundary \u0394z does not generalize MTG \u2192 A9",
       x = "boundary \u0394z at CPS 0.21", y = NULL) +
  theme_paper + theme(axis.text.y = element_text(size = 5.5),
                      legend.position = c(0.85, 0.12))
f2a <- f2a + labs(tag = "a") + theme(plot.tag = element_text(face = "bold", size = 9.9), plot.tag.position = c(0.02, 0.98))
save_fig(f2a, "Figure3_panelA.png", w = 6.5, h = 7.78)
