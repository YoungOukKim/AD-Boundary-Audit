# =============================================================================
# utils.R  -- shared theme, colours, paths, helpers for Paper 3 figure scripts
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(scales)
  library(patchwork); library(ggrepel)
})

# ---- paths (run from repo root) ----
DATA <- "data"
OUTF <- "output/figures"
dir.create(OUTF, showWarnings = FALSE, recursive = TRUE)

# ---- colours (match manuscript: up=red, down=blue; diverging RdBu) ----
UP   <- "#B2182B"   # up / induction
DOWN <- "#2166AC"   # down / decrement
GREY <- "#9E9E9E"
col_master <- "#E64A19"   # CPS~0.21 boundary
col_mci    <- "#1E88E5"   # CPS~0.42 boundary
RdBu <- c("#2166AC","#4393C3","#92C5DE","#D1E5F0","#F7F7F7",
          "#FDDBC7","#F4A582","#D6604D","#B2182B")

# ---- clean theme ----
theme_paper <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        axis.line = element_line(colour = "black", linewidth = 0.4),
        plot.title = element_text(face = "bold", size = 12, hjust = 0),
        axis.title = element_text(size = 10),
        legend.title = element_text(size = 9),
        legend.text  = element_text(size = 8))

save_fig <- function(p, file, w, h, dpi = 300) {
  ggsave(file.path(OUTF, file), p, width = w, height = h, dpi = dpi, bg = "white")
  message(">>> wrote ", file.path(OUTF, file))
}
