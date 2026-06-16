# =============================================================================
# make_supp_tables.R
# Reproduces Supplementary Tables S1-S3 (CSV) from the bundled data, in the
# manuscript's citation order:
#   S1 = 44-marker panel composition   [panel_composition.csv]
#   S2 = per-algorithm consensus        [consensus_breakpoints.csv]
#   S3 = per-marker drivers, CPS~0.21   [drivers_master_event_corrected.csv]
# Output: output/tables/SuppTable_S{1,2,3}_*.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
DATA <- "data"; OUTT <- "output/tables"
dir.create(OUTT, showWarnings = FALSE, recursive = TRUE)

# ---- S1: panel composition, grouped by category ----
panel <- read.csv(file.path(DATA, "panel_composition.csv"), check.names = FALSE)
s1 <- panel %>% group_by(Category) %>%
  summarise(Markers = paste(Marker, collapse = ", "), .groups = "drop")
write.csv(s1, file.path(OUTT, "SuppTable_S1_panel.csv"), row.names = FALSE)

# ---- S2: per-algorithm consensus breakpoints ----
con <- read.csv(file.path(DATA, "consensus_breakpoints.csv"), check.names = FALSE)
write.csv(con, file.path(OUTT, "SuppTable_S2_consensus.csv"), row.names = FALSE)

# ---- S3: per-marker drivers of the primary CPS~0.21 transition ----
dr <- read.csv(file.path(DATA, "drivers_master_event_corrected.csv"), check.names = FALSE)
dr <- dr[order(-abs(dr$Delta_z)), ]
s3 <- data.frame(
  Gene = dr$Gene,
  `delta_z (boundary)`      = sprintf("%+.3f", dr$Delta_z),
  `mean z, pre-boundary bin` = sprintf("%+.3f", dr$Bin1_z),
  `mean z, post-boundary bin`= sprintf("%+.3f", dr$Bin2_z),
  check.names = FALSE)
write.csv(s3, file.path(OUTT, "SuppTable_S3_drivers.csv"), row.names = FALSE)

message(">>> wrote SuppTable_S1/S2/S3 (", nrow(s1), " categories, ",
        nrow(con), " consensus rows, ", nrow(s3), " drivers)")
