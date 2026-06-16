# =============================================================================
# panel_finalize.R  -- merge within-atlas SEA-AD per-cell directions with the
# independent cross-cohort literature verdict; flag agreement; final table.
# Run AFTER external_marker_audit.R (needs sea_ad_panel_verdict.csv).
# Output: output/tables/external_markers/final_panel_verdict.csv
# =============================================================================
D <- "output/tables/external_markers"
cc <- read.csv(file.path(D,"panel_verdict.csv"), stringsAsFactors=FALSE)       # cross-cohort (lit)
sa <- read.csv(file.path(D,"sea_ad_panel_verdict.csv"), stringsAsFactors=FALSE) # SEA-AD per-cell
m <- merge(cc, sa[,c("gene","cell_type","spearman_rho","sea_ad_dir","shape","boundary_cps","bin_stability")],
           by="gene", all.x=TRUE)
norm <- function(d) ifelse(grepl("up",d),"up", ifelse(grepl("down",d),"down",
                    ifelse(grepl("flat",d),"flat","na")))
m$agreement <- ifelse(is.na(m$sea_ad_dir), "pending(run SEA-AD)",
                ifelse(norm(m$brain_dir)==m$sea_ad_dir, "CONFIRMED in SEA-AD",
                       paste0("DIVERGENT (lit ", m$brain_dir, " vs SEA-AD ", m$sea_ad_dir, ")")))
m <- m[,c("gene","cell_type","csf_dir","brain_dir","sea_ad_dir","spearman_rho","shape",
          "bin_stability","class","framework_read","agreement","source")]
write.csv(m, file.path(D,"final_panel_verdict.csv"), row.names=FALSE)
cat("\n=== FINAL PANEL VERDICT (cross-cohort + SEA-AD per-cell) ===\n")
print(m[,c("gene","cell_type","brain_dir","sea_ad_dir","class","agreement")], row.names=FALSE)
ok <- sum(m$agreement=="CONFIRMED in SEA-AD"); tot <- sum(!grepl("pending",m$agreement))
if (tot>0) cat(sprintf("\nWithin-atlas agreement: %d/%d markers match the literature brain direction.\n", ok, tot))
