# =============================================================================
# Figure4_metabolic_artifact.R
# Figure 4. Apparent metabolic conservation is a global-expression artifact.
#   (a) Uncorrected MTG vs A9 dz (apparent conservation, same-sign count)
#   (b) A9 near-uniform negative shift (sorted A9 dz)
#   (c) After global-expression correction (same-sign count rebalanced)
# Input : fig4_metabolic.csv (gene, MTG_dz, MTG_corrected, A9_dz, A9_corrected)
# Output: output/figures/Figure4.png
# =============================================================================
source("R/figures/utils.R")
m <- read.csv(file.path(DATA, "fig4_metabolic.csv"), check.names = FALSE)
ss <- function(x, y) sum(sign(x) == sign(y))
key <- c("SLC16A1","SLC16A3","SLC16A7","SLC2A1","SLC2A3","GAPDH","PKM","ENO2",
         "LDHA","LDHB","HK1","PGK1","MAPT","MBP","MOG","PLP1","ALDOA")
sc <- function(xx, yy, ttl, n) {
  d <- data.frame(x = m[[xx]], y = m[[yy]], g = m$gene)
  d$lab <- ifelse(d$g %in% key, d$g, NA)
  ggplot(d, aes(x, y)) +
    geom_abline(linetype = "dashed", colour = "grey60") +
    geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
    geom_point(colour = "#C0392B", size = 2) +
    geom_text_repel(aes(label = lab), size = 2.4, max.overlaps = 20, na.rm = TRUE) +
    labs(title = ttl, subtitle = sprintf("same-sign %d/21", n),
         x = "MTG \u0394z", y = "A9 \u0394z") + theme_paper
}
pa <- sc("MTG_dz","A9_dz", "Uncorrected (apparent conservation)", ss(m$MTG_dz, m$A9_dz))
pc <- sc("MTG_corrected","A9_corrected", "Global-expression corrected", ss(m$MTG_corrected, m$A9_corrected))
mb <- m[order(m$A9_dz), ]; mb$gene <- factor(mb$gene, levels = mb$gene)
pb <- ggplot(mb, aes(A9_dz, gene)) +
  geom_col(fill = "#2C6FBB", width = 0.7) + geom_vline(xintercept = 0, linewidth = 0.4) +
  labs(title = "A9: near-uniform negative shift", x = "A9 \u0394z", y = NULL) +
  theme_paper + theme(axis.text.y = element_text(size = 6))
save_fig((pa | pb | pc) + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold", size = 16.2)), "Figure4.png", w = 13, h = 5.2)
