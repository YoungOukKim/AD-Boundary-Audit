# =============================================================================
# Figure6_driver_signatures.R
# Figure 6. Driver signatures and within-MTG robustness.
#   Distinct driver signatures of the two transitions:
#   per-marker boundary effect (delta z) at CPS~0.21 (x) vs CPS~0.42 (y).
# Inputs: drivers_master_event_corrected.csv, drivers_MCI_phase_corrected.csv
# Output: output/figures/Figure6.png
# =============================================================================
source("R/figures/utils.R")

m <- read.csv(file.path(DATA, "drivers_master_event_corrected.csv"))[, c("Gene","Delta_z")]
c <- read.csv(file.path(DATA, "drivers_MCI_phase_corrected.csv"))[, c("Gene","Delta_z")]
names(m)[2] <- "dz021"; names(c)[2] <- "dz042"
d <- merge(m, c, by = "Gene")

# label the notable markers (largest |effect| on either axis + named drivers)
key <- c("HMOX1","GFAP","FKBP5","BCL2","SOX9","MERTK","TFRC","APOE","FTH1","SERPINA3","NEFL")
d$lab <- ifelse(d$Gene %in% key, d$Gene, NA)

f5 <- ggplot(d, aes(dz021, dz042)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.4) +
  geom_point(size = 2.4, colour = "#5E4FA2", alpha = 0.85) +
  geom_text_repel(aes(label = lab), size = 3, max.overlaps = 30,
                  segment.colour = "grey60", segment.size = 0.3, na.rm = TRUE) +
  labs(title = "Distinct driver signatures of the two transitions",
       x = "driver \u0394z at CPS\u22480.21",
       y = "driver \u0394z at CPS\u22480.42") +
  theme_paper

save_fig(f5, "Figure6.png", w = 7.2, h = 6.2)
