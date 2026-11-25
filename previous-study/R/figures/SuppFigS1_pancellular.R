# =============================================================================
# SuppFigS1_pancellular.R
# Supplementary Fig. S1. Cross-cell-type concordance of the CPS-0.21 transition.
#   (a) Median boundary position (CPS) per cell type: master event (~0.21) vs
#       MCI-stage (~0.5), with dashed reference lines at 0.21 and 0.42.
#   (b) Number of supporting algorithms (of up to ~11) per cell type, per boundary.
# Input : SEAAD_multicell_summary.csv  (per-cell-type boundary summary)
# Output: output/figures/SuppFigure_S1.png
# =============================================================================
source("R/figures/utils.R")
ORANGE <- "#FF7F0E"; BLUE <- "#1F77B4"

s <- read.csv(file.path(DATA, "SEAAD_multicell_summary.csv"), check.names = FALSE)

# ---- (a) median boundary position per cell type ----
da <- s %>% transmute(cell_type,
                      `master event` = median_breakpoint_master_event,
                      `MCI-stage`     = median_breakpoint_MCI_stage) %>%
  pivot_longer(-cell_type, names_to = "boundary", values_to = "cps")
ord_a <- s$cell_type[order(s$median_breakpoint_MCI_stage)]     # Excitatory bottom -> Oligo top
da$cell_type <- factor(da$cell_type, levels = ord_a)
da$boundary  <- factor(da$boundary, levels = c("MCI-stage","master event"))
pa <- ggplot(da, aes(cps, cell_type, colour = boundary)) +
  geom_vline(xintercept = 0.21, linetype = "dashed", colour = ORANGE) +
  geom_vline(xintercept = 0.42, linetype = "dashed", colour = BLUE) +
  geom_point(size = 3) +
  scale_colour_manual(values = c("master event" = ORANGE, "MCI-stage" = BLUE), name = NULL) +
  labs(title = "Median boundary position per cell type", x = "CPS", y = NULL) +
  theme_paper + theme(legend.position = "bottom")

# ---- (b) algorithm support per cell type (alphabetical; Astrocyte at bottom) ----
db <- s %>% transmute(cell_type,
                      `master event` = n_algorithms_master_event,
                      `MCI-stage`     = n_algorithms_MCI_stage) %>%
  pivot_longer(-cell_type, names_to = "boundary", values_to = "n")
db$cell_type <- factor(db$cell_type, levels = rev(sort(unique(s$cell_type))))
db$boundary  <- factor(db$boundary, levels = c("MCI-stage","master event"))
pb <- ggplot(db, aes(n, cell_type, fill = boundary)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.7) +
  scale_fill_manual(values = c("master event" = ORANGE, "MCI-stage" = BLUE), name = NULL) +
  labs(title = "Algorithm support per cell type", x = "n algorithms", y = NULL) +
  theme_paper + theme(legend.position = "bottom")

fig <- (pa | pb) + plot_annotation(
  title = "Cross-cell-type concordance of the CPS-0.21 transition (SEA-AD)",
  tag_levels = "a", theme = theme(plot.title = element_text(face = "bold", size = 13))) &
  theme(plot.tag = element_text(face = "bold", size = 17.3))
save_fig(fig, "SuppFigure_S1.png", w = 14, h = 6)
