# =============================================================================
# Figure2_consensus.R
# Figure 2. A multi-algorithm consensus localizes a reproducible reactive-
#           astrocyte transition in SEA-AD MTG.
#   (a) Top drivers of the CPS~0.21 transition        [drivers_master_event_corrected.csv]
#   (b) Drivers of the secondary MCI-stage (CPS~0.42) [drivers_MCI_phase_corrected.csv]
#   (c) Marker x CPS-bin mean-z heatmap, both transitions [SEAAD_44panel_bin_means.csv]
# Output: output/figures/Figure2.png
# =============================================================================
source("R/figures/utils.R")

# ---- (a) & (b): top-12 drivers by |delta z| ----
driver_bar <- function(csv, title, xlab) {
  d <- read.csv(file.path(DATA, csv), check.names = FALSE)
  d <- d[order(-abs(d$Delta_z)), ][1:12, ]
  d$dir <- ifelse(d$Delta_z > 0, "up", "down")
  ggplot(d, aes(x = reorder(Gene, abs(Delta_z)), y = Delta_z, fill = dir)) +
    geom_col(width = 0.72) + coord_flip() +
    scale_fill_manual(values = c(up = UP, down = DOWN), guide = "none") +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    labs(title = title, x = NULL, y = xlab) + theme_paper
}
f1a <- driver_bar("drivers_master_event_corrected.csv",
                  "Top drivers at CPS\u22480.21", expression(Delta*z~"@0.21"))
f1b <- driver_bar("drivers_MCI_phase_corrected.csv",
                  "Top drivers at CPS\u22480.42 (secondary)", expression(Delta*z~"@0.42"))

# ---- (c): marker x CPS-bin heatmap ----
bm <- read.csv(file.path(DATA, "SEAAD_44panel_bin_means.csv"), check.names = FALSE)
markers <- setdiff(names(bm), c("cps_bin", "cps_center", "n_cells"))
# row order: exact manuscript order, top -> bottom (HMOX1 at the very bottom)
ms_top2bottom <- c("FKBP5","GFAP","BCL2","MERTK","SOX9","SERPINA3","TP53","STAT3",
  "OCLN","AQP4","NFE2L2","SOCS3","FTH1","FTL","TREM2","VCAM1","C3","IL1B","STAT1",
  "CLU","TNF","SELE","FBLN5","CXCL12","GDNF","CR1","CD33","CXCL16","VWF","PECAM1",
  "CASP3","IL6","MFGE8","NGF","NFKB1","ICAM1","NEFL","VIM","CCL2","APP","BDNF",
  "APOE","TFRC","HMOX1")
ms_top2bottom <- c(ms_top2bottom, setdiff(markers, ms_top2bottom))
ord <- rev(ms_top2bottom)   # ggplot puts levels[1] at bottom -> HMOX1 bottom, FKBP5 top
long <- bm %>%
  select(all_of(c("cps_bin", markers))) %>%
  pivot_longer(-cps_bin, names_to = "marker", values_to = "z") %>%
  mutate(marker = factor(marker, levels = ord),
         binlab = factor(cps_bin, levels = bm$cps_bin,
                         labels = c("0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5","0.5-0.6",
                                    "0.6-0.7","0.7-0.8","0.8-0.9","0.9-1.0")))
lim <- max(abs(long$z), na.rm = TRUE)
f1c <- ggplot(long, aes(binlab, marker, fill = z)) +
  geom_tile() +
  geom_vline(xintercept = 1.5, linetype = "dashed", linewidth = 0.6) +
  scale_fill_gradientn(colours = RdBu, limits = c(-lim, lim), name = "mean z") +
  labs(title = "Marker \u00d7 CPS-bin mean z (both transitions)", x = "CPS bin", y = NULL) +
  theme_paper +
  theme(axis.text.y = element_text(size = 6),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        panel.grid = element_blank())

fig1 <- (f1a / f1b | f1c) + plot_layout(widths = c(1, 1.15)) +
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold", size = 16.1))
save_fig(fig1, "Figure2.png", w = 13, h = 8)
