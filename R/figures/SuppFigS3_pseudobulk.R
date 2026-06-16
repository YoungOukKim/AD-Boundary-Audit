# =============================================================================
# SuppFigS3_pseudobulk.R
# Supplementary Fig. S3. Donor-level (pseudobulk) testing is required.
#   (a) DEGs at FDR<0.05 vs nuclei: cell-level (red) vs pseudobulk (blue)
#   (b) inflation = cell-DEG / pseudobulk-DEG (sorted; median line)
#   (c) cell-level DEGs biased to high expression (mean expr DEG vs all)
# Input : P3_pseudobulk_vs_celllevel.csv
# Output: output/figures/SuppFigure_S3.png
# =============================================================================
source("R/figures/utils.R")
d <- read.csv(file.path(DATA, "P3_pseudobulk_vs_celllevel.csv"), check.names = FALSE)
d$infl <- d$deg_cell / pmax(d$deg_pb, 1)
la <- pivot_longer(d, c("deg_cell","deg_pb"), names_to = "method", values_to = "deg")
la$method <- factor(la$method, labels = c("cell-level (nucleus)","pseudobulk (donor)"))
pa <- ggplot(la, aes(n_used, deg, colour = method)) +
  geom_point(size = 2) + scale_y_log10() +
  scale_colour_manual(values = c("cell-level (nucleus)" = UP, "pseudobulk (donor)" = DOWN), name = NULL) +
  labs(title = "Cell-level inflates DEG counts", x = "nuclei per cell type", y = "DEGs at FDR<0.05") +
  theme_paper + theme(legend.position = c(0.7, 0.18))
ds <- d[order(d$infl), ]; ds$idx <- seq_len(nrow(ds))
pb <- ggplot(ds, aes(idx, infl)) +
  geom_col(fill = "#7E57C2", width = 0.8) +
  geom_hline(yintercept = median(d$infl), linetype = "dashed") +
  annotate("text", x = 4, y = median(d$infl) * 1.25, size = 2.8,
           label = sprintf("median %.1f\u00d7", median(d$infl))) +
  labs(title = "Inflation = cell-DEG / pseudobulk-DEG", x = "cell type (sorted)", y = "inflation") + theme_paper
above <- sum(d$meanexpr_cellDEG > d$meanexpr_all)
pc <- ggplot(d, aes(meanexpr_all, meanexpr_cellDEG)) +
  geom_abline(linetype = "dashed", colour = "grey60") +
  geom_point(colour = UP, size = 2) +
  labs(title = "DEGs biased to high expression",
       subtitle = sprintf("%d/%d types above diagonal", above, nrow(d)),
       x = "mean expr, all tested genes", y = "mean expr, cell-level DEGs") + theme_paper
save_fig((pa | pb | pc) + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold", size = 13)), "SuppFigure_S3.png", w = 13, h = 3.98)
