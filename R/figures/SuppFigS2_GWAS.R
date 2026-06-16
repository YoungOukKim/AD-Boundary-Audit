# =============================================================================
# SuppFigS2_GWAS.R   (matches manuscript caption)
# Supplementary Fig. S2. AD-GWAS risk genes converge on the CPS-0.207 boundary
# in microglia but not astrocytes.
#   (a) Per-gene boundary effect (dz): astrocytes (x) vs microglia (y); dashed
#       identity line; most genes above the line (larger in microglia).
#   (b) Observed mean |dz| of the risk-gene set (red diamond) vs matched null
#       (band: mean +/-1 and +/-1.96 s.d.): astro P=0.22, micro P=0.029.
# Inputs: P3_gwas_pergene_dz.csv, P3_gwas_convergence_results.csv
# Output: output/figures/SuppFigure_S2.png
# =============================================================================
source("R/figures/utils.R")
pg <- read.csv(file.path(DATA, "P3_gwas_pergene_dz.csv"), check.names = FALSE)
names(pg)[1] <- "gene"
lab <- c("MS4A6A","GRN","CASS4","MEF2C","SORT1","PLCG2","INPP5D","APOE","BIN1","HLA-DRB1","ABCA7","CR1")
pg$lab <- ifelse(pg$gene %in% lab, pg$gene, NA)
pa <- ggplot(pg, aes(Astro, `Micro-PVM`)) +
  geom_abline(linetype = "dashed", colour = "grey55") +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(colour = "#5E4FA2", size = 2, alpha = 0.85) +
  geom_text_repel(aes(label = lab), size = 2.6, max.overlaps = 20, na.rm = TRUE) +
  labs(title = "Per-gene AD-GWAS boundary effect",
       x = "astrocyte \u0394z", y = "microglia \u0394z") + theme_paper

cv <- read.csv(file.path(DATA, "P3_gwas_convergence_results.csv"), check.names = FALSE)
cv$cell <- factor(c("Astrocytes","Microglia (Micro-PVM)")[match(cv$cell_type, c("Astro","Micro-PVM"))],
                  levels = c("Microglia (Micro-PVM)","Astrocytes"))
cv$plab <- sprintf("P = %.3f", cv$emp_P)
pb <- ggplot(cv, aes(y = cell)) +
  geom_errorbarh(aes(xmin = null_mean - 1.96*null_sd, xmax = null_mean + 1.96*null_sd), height = 0.12, colour = "grey70", linewidth = 0.7) +
  geom_errorbarh(aes(xmin = null_mean - null_sd, xmax = null_mean + null_sd), height = 0.22, colour = "grey45", linewidth = 1.4) +
  geom_point(aes(x = null_mean), shape = 124, size = 4, colour = "grey30") +
  geom_point(aes(x = mean_abs_dz), shape = 18, size = 4.5, colour = "#C0392B") +
  geom_text(aes(x = mean_abs_dz, label = plab), vjust = -1.2, size = 3, colour = "#C0392B") +
  labs(title = "Observed vs matched null",
       x = "mean |\u0394z| of risk-gene set", y = NULL) + theme_paper

save_fig((pa | pb) + plot_layout(widths = c(1, 1.1)) + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold", size = 13)), "SuppFigure_S2.png", w = 13, h = 4.85)
