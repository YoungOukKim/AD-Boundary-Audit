# =============================================================================
# CSF_TMT_detection_audit.R
# Confirms which 44-panel markers are detected in ADNI Emory CSF TMT-MS, and
# that HMOX1 is ABSENT (supports the Figure 2b annotation "HMOX1 absent in
# TMT-MS"). NOTE: the rho-vs-MMSE curves in Figure 2b are from ADNI SomaScan
# joined to MMSE (HMOX1 rho -0.11, GFAP -0.17), which this TMT-MS file lacks.
# The full TMT-MS matrix is ADNI controlled-access and is NOT redistributed
# here; only the column header is shipped (EMORY_CSF_TMT_MS_columns.txt).
# =============================================================================
DATA <- "data"
full <- file.path(DATA, "EMORY_CSF_TMT_MS.csv")
if (file.exists(full)) {
  h <- names(read.csv(full, nrows = 1, check.names = FALSE))
} else {
  h <- readLines(file.path(DATA, "EMORY_CSF_TMT_MS_columns.txt"))
}
panel <- c("HMOX1","GFAP","FKBP5","MERTK","TFRC","APOE","NEFL","SERPINA3","CLU",
           "AQP4","VCAM1","C3","STAT3","BCL2","SOX9","FTH1","FTL","TREM2")
detected <- panel[sapply(panel, function(g) any(startsWith(toupper(h), paste0(g, "_"))))]
cat("Detected panel markers in TMT-MS:\n  ", paste(detected, collapse = ", "), "\n")
cat("HMOX1 detected? ", "HMOX1" %in% detected, "  (manuscript: absent in TMT-MS)\n")
