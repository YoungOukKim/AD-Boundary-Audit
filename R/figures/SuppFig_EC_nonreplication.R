# =============================================================================
# SuppFig_EC_nonreplication.R
# EC (Leng 2021) non-replication evidence underlying Figure 3 / Table S5.
#   (a) EC astrocyte panel logFC: FTH1 rises (opposite to MTG iron-storage
#       decline) while TFRC/GFAP are flat; HMOX1 not testable (detection floor).
#   (b) Positive control (GFAP-high vs GFAP-low astrocytes): canonical Leng
#       homeostatic genes are downregulated -> pipeline correctly recovers
#       reactive de-homeostasis.
# Inputs : testA_edgeR_panel_EC.csv, testB_edgeR_GFAPhigh_EC.csv
# Output : output/figures/SuppFig_EC_nonreplication.png
# =============================================================================
source("R/figures/utils.R")

# ---- (a) EC astrocyte panel ----
a <- read.csv(file.path(DATA, "testA_edgeR_panel_EC.csv"))
as <- a[a$celltype == "Astro" & a$testable, ]
as$dir <- ifelse(as$logFC > 0, "up", "down")
pa <- ggplot(as, aes(reorder(gene, logFC), logFC, fill = dir)) +
  geom_col(width = 0.7) + coord_flip() +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = c(up = UP, down = DOWN), guide = "none") +
  labs(title = "a  EC astrocytes: panel opposite to MTG (FTH1\u2191)",
       x = NULL, y = "EC \u0394logFC") + theme_paper

# ---- (b) positive control ----
b <- read.csv(file.path(DATA, "testB_edgeR_GFAPhigh_EC.csv"))
b$homeo <- ifelse(b$is_leng_homeo == TRUE | b$is_leng_homeo == "TRUE", "Leng homeostatic", "other")
pb <- ggplot(b, aes(reorder(gene, logFC), logFC, fill = homeo)) +
  geom_col(width = 0.7) + coord_flip() +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = c("Leng homeostatic" = DOWN, "other" = GREY), name = NULL) +
  labs(title = "b  Positive control (GFAP-high): homeostatic genes down",
       x = NULL, y = "GFAP-high vs -low \u0394logFC") +
  theme_paper + theme(legend.position = "bottom")

fig <- (pa | pb) + plot_annotation(
  title = "EC (Leng 2021) non-replication and pipeline positive control",
  theme = theme(plot.title = element_text(face = "bold", size = 13)))
save_fig(fig, "SuppFig_EC_nonreplication.png", w = 13, h = 5.5)
