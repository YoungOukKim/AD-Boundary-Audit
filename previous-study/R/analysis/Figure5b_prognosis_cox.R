# =============================================================================
# Figure5b_prognosis_cox.R   [RUN ON D: — ADNIMERGE2 .rda + Emory TMT-MS]
#
# Computes the Figure 5b adjusted-Cox prognosis that the manuscript currently
# reports (3.79 / 0.99 / 0.92; baseline-deduped Cox, n=1,105/527 events) and that Figure 5b merely PLOTS. NOTE: the manuscript fit used Python statsmodels PHReg; this R/survival script is an equivalent re-implementation — confirm it reproduces these values on the controlled inputs.
# Writes a forest CSV that Figure5_brainCSF_prognosis.R reads -> real reproduction.
#
# DATA (your D:):
#   MMSE.rda                      -> longitudinal MMSE  (total col MMSCORE)
#   UPENNBIOMK_ROCHE_ELECSYS.rda  -> CSF Abeta42 / pTau (adjustment)
#   <Emory TMT-MS csv>            -> molecular PC1 of CSF panel proteins
#
# ENDPOINT: first visit with MMSE <= 23 OR MMSE <= (baseline MMSE - 3).
#   (1) baseline MMSE <= 26                       -> HR  [circular]
#   (2) molecular PC1 (TMT-MS panel)              -> HR/s.d. [unadjusted]
#   (3) PC1 + baseline MMSE + Abeta42 + pTau      -> HR/s.d. [adjusted]
#
# RUN:  Rscript R/analysis/Figure5b_prognosis_cox.R          (needs: survival)
# =============================================================================
suppressMessages(library(survival))

## ===========================================================================
## CONFIG  — edit paths/panel only   (RUN FROM REPO ROOT: setwd to AD-Boundary-Audit)
## ===========================================================================
ADNI_DIR <- "data-external/ADNIMERGE2"   # controlled-access ADNI .rda dir; set to your local path
PROT_CSV <- "data-external/EMORY_CSF_TMT_MS.csv"   # Emory TMT-MS (PTID/RID + protein cols); set to your local path
OUT_FOREST <- "data/fig5b_prognosis_computed.csv"
OUT_FULL   <- "output/tables/Figure5b_cox_HR.csv"

# CSF panel proteins (gene_UNIPROT) to build molecular PC1 from; missing ones dropped.
PANEL <- c("GFAP_P14136","SERPINA3_P01011","VIM_P08670","AQP4_P55087",
           "MERTK_Q12866","FKBP5_Q13451","C3_P01024","STAT3_P40763",
           "APOE_P02649","CLU_P10909","TREM2_Q9NZC2","FTH1_P02794",
           "FTL_P02792","TFRC_P02786","NEFL_P07196","BCL2_P10415")
## ===========================================================================

ld   <- function(f){ e <- new.env(); load(file.path(ADNI_DIR, f), envir = e); get(ls(e)[1], envir = e) }
pick <- function(df, cands, what){
  hit <- cands[cands %in% colnames(df)]
  if (!length(hit)) stop(sprintf("'%s' not found. tried {%s}\n  available: %s",
       what, paste(cands, collapse=", "), paste(colnames(df), collapse=", ")))
  hit[1]
}

## ---- load + show columns (so wrong guesses are obvious) -------------------
MMSE <- ld("MMSE.rda")
CSF  <- ld("UPENNBIOMK_ROCHE_ELECSYS.rda")
prot <- read.csv(PROT_CSV, check.names = FALSE, stringsAsFactors = FALSE)
cat("== MMSE cols ==\n"); print(colnames(MMSE))
cat("== CSF  cols ==\n"); print(colnames(CSF))

m_rid  <- pick(MMSE, c("RID"), "MMSE RID")
m_sc   <- pick(MMSE, c("MMSCORE","MMSE","MMSETOTAL","MMTOTAL"), "MMSE total")
m_date <- pick(MMSE, c("EXAMDATE","VISDATE","USERDATE"), "MMSE date")
c_rid  <- pick(CSF,  c("RID"), "CSF RID")
c_ab   <- pick(CSF,  c("ABETA42","ABETA","A_BETA_42"), "CSF Abeta42")
c_tau  <- pick(CSF,  c("PTAU","PTAU181","PTAU_"), "CSF pTau")
p_rid  <- pick(prot, c("RID"), "TMT-MS RID")

## ---- baseline + endpoint from longitudinal MMSE --------------------------
MMSE[[m_date]] <- as.Date(as.character(MMSE[[m_date]]))
MMSE <- MMSE[!is.na(MMSE[[m_rid]]) & !is.na(MMSE[[m_sc]]) & !is.na(MMSE[[m_date]]), ]
MMSE <- MMSE[order(MMSE[[m_rid]], MMSE[[m_date]]), ]
surv <- do.call(rbind, lapply(split(MMSE, MMSE[[m_rid]]), function(d){
  base_date <- d[[m_date]][1]; base_mmse <- d[[m_sc]][1]
  d$months  <- as.numeric(d[[m_date]] - base_date) / 30.44
  d$ev      <- (d[[m_sc]] <= 23) | (d[[m_sc]] <= base_mmse - 3)
  fu <- d[d$months > 0, ]; if (!nrow(fu)) return(NULL)
  k <- which(fu$ev)
  if (length(k)) { t <- fu$months[k[1]]; e <- 1L } else { t <- max(fu$months); e <- 0L }
  data.frame(RID = d[[m_rid]][1], time = t, event = e,
             mmse_bl = base_mmse, mmse_le26 = as.integer(base_mmse <= 26))
}))

## ---- baseline CSF Abeta42 / pTau -----------------------------------------
CSF[[c_ab]]  <- suppressWarnings(as.numeric(as.character(CSF[[c_ab]])))
CSF[[c_tau]] <- suppressWarnings(as.numeric(as.character(CSF[[c_tau]])))
if ("EXAMDATE" %in% colnames(CSF)) CSF <- CSF[order(CSF[[c_rid]], as.Date(as.character(CSF$EXAMDATE))), ]
csf_bl <- do.call(rbind, lapply(split(CSF, CSF[[c_rid]]), function(d){
  d <- d[!is.na(d[[c_ab]]) | !is.na(d[[c_tau]]), , drop = FALSE]; if (!nrow(d)) return(NULL)
  data.frame(RID = d[[c_rid]][1], amyloid = d[[c_ab]][1], tau = d[[c_tau]][1])
}))

## ---- molecular PC1 from TMT-MS panel -------------------------------------
## ---- dedup proteomics to ONE baseline row per subject (earliest visit) ----
if ("VISCODE2" %in% colnames(prot)) {                       # date-format-agnostic baseline pick
  .vc <- toupper(as.character(prot$VISCODE2))
  .mo <- suppressWarnings(as.integer(sub("^[^0-9]*", "", .vc)))  # M06->6, BL->NA
  .rk <- ifelse(.vc %in% c("BL","SC","SCMRI","M00"), 0L, ifelse(is.na(.mo), 999L, .mo))
  prot <- prot[order(prot[[p_rid]], .rk), ]
}
prot <- prot[!duplicated(prot[[p_rid]]), ]
cat(sprintf("TMT-MS deduped to baseline: %d unique subjects\n", nrow(prot)))

have <- intersect(PANEL, colnames(prot))
if (length(have) < 3) stop("only ", length(have), " panel proteins in TMT-MS; edit PANEL")
cat("panel PC1 from (", length(have), "): ", paste(have, collapse=", "), "\n")
Xp <- scale(as.matrix(sapply(prot[, have], function(z) suppressWarnings(as.numeric(z)))))
Xp[is.na(Xp)] <- 0
prot_pc <- data.frame(RID = prot[[p_rid]],
                      mol_pc1 = as.numeric(scale(prcomp(Xp, center = FALSE, scale. = FALSE)$x[, 1])))

## ---- merge (= cohort with MMSE + CSF + proteomics) & Cox -----------------
df <- merge(merge(surv, csf_bl, by = "RID"), prot_pc, by = "RID")
df <- df[is.finite(df$time) & df$time > 0, ]
cat(sprintf("\nCohort n = %d,  events = %d\n\n", nrow(df), sum(df$event)))

getHR <- function(fit, term){ ci <- summary(fit)$conf.int
  c(HR = unname(ci[term,"exp(coef)"]), lo = unname(ci[term,"lower .95"]), hi = unname(ci[term,"upper .95"])) }
m1 <- coxph(Surv(time, event) ~ mmse_le26, data = df)
m2 <- coxph(Surv(time, event) ~ mol_pc1, data = df)
m3 <- coxph(Surv(time, event) ~ mol_pc1 + mmse_bl + amyloid + tau, data = df)
a <- getHR(m1,"mmse_le26"); b <- getHR(m2,"mol_pc1"); cc <- getHR(m3,"mol_pc1")

forest <- data.frame(
  label = c("MMSE-threshold boundary (<=26)\n  predictor is cognition; endpoint MMSE-circular",
            "Molecular boundary (PC1), unadjusted",
            "Molecular boundary (PC1)\n  adjusted: baseline MMSE, amyloid, tau"),
  HR = round(c(a["HR"],b["HR"],cc["HR"]),2),
  CI_lower = round(c(a["lo"],b["lo"],cc["lo"]),2),
  CI_upper = round(c(a["hi"],b["hi"],cc["hi"]),2),
  group = c("circular","molecular","molecular_adj"), stringsAsFactors = FALSE)
cat("Figure 5b — computed:\n"); print(forest)                       # show BEFORE writing
dir.create(dirname(OUT_FOREST), showWarnings = FALSE, recursive = TRUE) # ensure data/
dir.create(dirname(OUT_FULL),   showWarnings = FALSE, recursive = TRUE) # ensure output/tables/
try(write.csv(forest, OUT_FOREST, row.names = FALSE), silent = TRUE)
try(write.csv(forest, OUT_FULL,   row.names = FALSE), silent = TRUE)
cat("\n[mol unadj p]", signif(summary(m2)$coefficients["mol_pc1","Pr(>|z|)"],3),
    "  [mol adj p]", signif(summary(m3)$coefficients["mol_pc1","Pr(>|z|)"],3), "\n")
cat(">>> wrote", OUT_FOREST, "(Figure5_brainCSF_prognosis.R reads this)\n")
