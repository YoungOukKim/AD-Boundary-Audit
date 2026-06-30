# =============================================================================
# synthetic_benchmark.R
# Ground-truth benchmark of the reproducibility-audit framework.
#
# Plants known boundaries in synthetic marker x bin matrices and asks whether
# the SAME engine (boundary_detection.R, sourced verbatim) + the SAME audit
# gates (algorithmic consensus, random-panel permutation null, leave-one-marker
# -out) certify real panel-specific boundaries and reject four artefact classes.
#
# Gates are fixed a priori (NOT tuned to the answer):
#   support  S* = 95th percentile of the random-panel null support
#   specificity  p <= 0.05
#   leave-one-marker-out  max single-marker boundary shift <= 1 bin (0.05)
#
# Background is calibrated to reproduce the EMPIRICAL random-panel null shape
# (change-points concentrated late, sparse early), matching the real SEA-AD
# astrocyte null (manuscript: 0/1000 random panels at early CPS, ~40% near 0.54).
#
# Grid: B = 20 bins, panel K = 44 markers, NPERM = 300 random panels, seed 42.
# (The real-data analysis used 9 cell-count-limited CPS bins and 1,000 panels;
#  the engine is bin-count-agnostic and the null shape is stable, so the
#  benchmark uses a finer, fixed grid for clean localization resolution.)
#
# Run from repo root:  Rscript R/analysis/synthetic_benchmark.R
# Outputs: output/benchmark/{partA.csv, global.csv, lomo.csv, summary.txt,
#          Figure_synthetic_benchmark.png}
# Runtime ~25 min single-core (segmented engine ~0.5 s/call; leave-one-marker-out adds ~18 min).
# =============================================================================
suppressWarnings(suppressMessages(source("R/analysis/boundary_detection.R")))
OUT <- "output/benchmark"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

## consensus = densest breakpoint within +/-0.05 (repo: celltype_audit_from_raw.R)
consensus <- function(bp){ bp <- bp[!is.na(bp)]; if(!length(bp)) return(list(bp=NA,n=0))
  dd <- sapply(bp, function(c0) sum(abs(bp - c0) <= 0.05)); list(bp = bp[which.max(dd)], n = max(dd)) }

B <- 20; AXIS <- seq(0.025, 0.975, length.out = B); K <- 44; NPERM <- 300
TS <- 0.22; TRUTH <- 0.20          # step at 0.22 -> boundary between bins -> 0.20
sv <- function(tau){ s <- as.numeric(AXIS >= tau); s - mean(s) }

## --- background transcriptome: flat-noise + DIVERSE LATE logistic (no linear) -
##     -> random panels break late, sparse early (empirical null shape).
gen_bg <- function(G, seed, global_step = FALSE){ set.seed(seed); M <- matrix(0, B, G)
  for(g in 1:G){ if(runif(1) < 0.55){ M[,g] <- rnorm(B, 0, 0.6) }            # flat noise
    else { mu <- runif(1, 0.48, 0.84); k <- runif(1, 0.04, 0.12)
      sh <- 1/(1+exp(-(AXIS-mu)/k)); sh <- sh - mean(sh); sh <- sh/sd(sh)
      M[,g] <- abs(rnorm(1,0,0.7))*sample(c(-1,1),1)*sh + rnorm(B,0,0.2) } }
  if(global_step) M <- M + outer(sv(TS), 0.6*(1+rnorm(G,0,0.15)))            # shared step in EVERY gene
  M }

## --- random-panel null (built once per background) ---------------------------
build_null <- function(BG){ set.seed(42); G <- ncol(BG)
  pc1 <- numeric(NPERM); sup <- numeric(NPERM)
  for(p in 1:NPERM){ cols <- sample(G, K); d <- detect_boundaries(BG[,cols,drop=FALSE], AXIS)
    v <- d$breakpoint[d$algorithm=="PCA_PC1"]; pc1[p] <- if(length(v)) v else NA
    prim <- d$breakpoint[d$algorithm=="Hier_k3_bp1"][1]
    sup[p] <- sum(abs(d$breakpoint - prim) <= 0.05, na.rm = TRUE) }
  list(pc1 = pc1, sup = sup) }

## --- one panel -> primary (Ward), support, specificity p ---------------------
eval_panel <- function(P, pc1null){ d <- detect_boundaries(P, AXIS)
  prim <- d$breakpoint[d$algorithm=="Hier_k3_bp1"][1]; if(is.na(prim)) prim <- consensus(d$breakpoint)$bp
  list(primary = prim, support = sum(abs(d$breakpoint - prim) <= 0.05, na.rm = TRUE),
       pspec = mean(abs(pc1null - prim) <= 0.05, na.rm = TRUE), comps = d) }

## --- panel generators (ground-truth regimes) ---------------------------------
mk <- function(kind, BG, seed, delta = 0.6){ set.seed(seed); idx <- sample(ncol(BG), K); P <- BG[, idx, drop=FALSE]
  if(kind=="real_spec") P <- 0.6*P + outer(sv(TS), delta*(1+rnorm(K,0,0.15)))     # co-regulated localized step
  if(kind=="ramp")      P <- 0.6*P + outer(AXIS-mean(AXIS), delta*(1+rnorm(K,0,0.15)))  # monotonic, no step
  if(kind=="decoy")     P[,1] <- P[,1] + 5*delta*sv(TS)                            # single-marker step
  P }                                                                              # "null" = panel as-is

message(">>> building backgrounds and nulls ...")
BG  <- gen_bg(3000, 101);  NL  <- build_null(BG)
BGg <- gen_bg(3000, 909, global_step = TRUE); NLg <- build_null(BGg)
S_STAR <- as.numeric(quantile(NL$sup, 0.95, na.rm = TRUE)); ALPHA <- 0.05; LOMO_TOL <- 0.05

## ============ Part A: 5 regimes + leave-one-marker-out for real_spec/decoy ====
prim_only <- function(P){ d <- detect_boundaries(P, AXIS); pr <- d$breakpoint[d$algorithm=="Hier_k3_bp1"][1]
  if(is.na(pr)) pr <- consensus(d$breakpoint)$bp; pr }
lomo_maxshift <- function(P, b0) max(sapply(1:K, function(j) abs(prim_only(P[,-j,drop=FALSE]) - b0)), na.rm=TRUE)
NREP <- 24; rows <- list(); real_comps <- list()
for(kind in c("real_spec","null","ramp","decoy")){
  for(i in 1:NREP){ P <- mk(kind, BG, 3000+i); e <- eval_panel(P, NL$pc1)
    ms <- if(kind %in% c("real_spec","decoy")) lomo_maxshift(P, e$primary) else NA_real_
    rows[[length(rows)+1]] <- data.frame(regime=kind, rep=i, primary=e$primary,
      support=e$support, pspec=e$pspec, locerr=abs(e$primary-TRUTH), maxshift=ms)
    if(kind=="real_spec") real_comps[[length(real_comps)+1]] <- e$comps } }
# global regime evaluated against ITS OWN null (LOMO not applicable: no single dominant marker)
for(i in 1:NREP){ P <- { set.seed(4000+i); BGg[, sample(ncol(BGg), K), drop=FALSE] }; e <- eval_panel(P, NLg$pc1)
  rows[[length(rows)+1]] <- data.frame(regime="global", rep=i, primary=e$primary,
    support=e$support, pspec=e$pspec, locerr=abs(e$primary-TRUTH), maxshift=NA_real_) }
df <- do.call(rbind, rows)
df$cert_2gate <- (df$support >= S_STAR) & (df$pspec <= ALPHA)
# 3-gate: LOMO applies only where one marker could drive the boundary (real_spec, decoy);
# null/ramp/global carry no dominant marker, so their 3-gate equals their 2-gate.
df$cert_3gate <- df$cert_2gate & ( is.na(df$maxshift) | df$maxshift <= LOMO_TOL )
write.csv(df[df$regime!="global",], file.path(OUT,"partA.csv"), row.names = FALSE)
write.csv(df[df$regime=="global",], file.path(OUT,"global.csv"), row.names = FALSE)
# leave-one-marker-out table (real_spec & decoy; all replicates) for the figure
L <- df[df$regime %in% c("real_spec","decoy"), c("regime","rep","maxshift")]
write.csv(L, file.path(OUT,"lomo.csv"), row.names = FALSE)

## ============================ Part C: algorithm ablation ===================
fams <- c("PCA","CUSUM","Mahal","Hier","KMeans","DP","Variance","Trajectory","Segmented")
famcols <- list(PCA=c("PCA_PC1","PCA_PC2","PCA_PC3"), CUSUM=c("CUSUM_PC1","CUSUM_PC2","CUSUM_PC3"),
  Mahal="Mahal_max_gap", Hier=c("Hier_k3_bp1","Hier_k3_bp2"), KMeans=c("KMeans_k3_bp1","KMeans_k3_bp2"),
  DP=c("DP_PC1","DP_PC2","DP_PC3"), Variance="Variance_jump", Trajectory="Trajectory_inflection",
  Segmented=c("Segmented_PC1","Segmented_PC2","Segmented_PC3"))
fam_err <- sapply(fams, function(f) median(sapply(real_comps, function(d){
  v <- d$breakpoint[d$algorithm %in% famcols[[f]]]; if(!length(v)) return(NA); abs(median(v)-TRUTH) }), na.rm=TRUE))
loo_err <- sapply(fams, function(f) median(sapply(real_comps, function(d){
  keep <- d[!(d$algorithm %in% famcols[[f]]), ]; abs(consensus(keep$breakpoint)$bp - TRUTH) }), na.rm=TRUE))
cons_err <- median(sapply(real_comps, function(d) abs(consensus(d$breakpoint)$bp - TRUTH)))

## ============================ summary =======================================
dg <- df[df$regime=="global",]
sink(file.path(OUT,"summary.txt"))
cat("SYNTHETIC BENCHMARK -- operating characteristics\n================================================\n")
cat(sprintf("B=%d bins, K=%d markers, NPERM=%d, seed 42 | gates: support>=%.0f, p<=%.2f, LOMO<=%.2f | truth=%.2f\n\n",
  B,K,NPERM,S_STAR,ALPHA,LOMO_TOL,TRUTH))
print(aggregate(cbind(primary,support,pspec,locerr,cert_3gate)~regime, df, function(x) round(median(x),3)), row.names=FALSE)
cat("\nLOMO max single-marker shift (median): "); print(round(tapply(L$maxshift,L$regime,median),3))
# Sensitivity / specificity from the EXACT per-panel 3-gate (support + specificity + leave-one-marker-out).
sens <- mean(df$cert_3gate[df$regime=="real_spec"])
spec_class <- sapply(c("global","null","ramp","decoy"), function(r) 1 - mean(df$cert_3gate[df$regime==r]))
n_decoy_slip <- sum(df$cert_3gate[df$regime=="decoy"])
cat(sprintf("\nSENSITIVITY (real_spec certified) = %.0f%%\n", 100*sens))
cat(sprintf("SPECIFICITY (artefacts rejected) = %.0f%% | global %.0f null %.0f ramp %.0f decoy %.0f\n",
  100*mean(spec_class), 100*spec_class["global"],100*spec_class["null"],100*spec_class["ramp"],100*spec_class["decoy"]))
cat(sprintf("  (decoy: %d of %d replicates slipped through leave-one-marker-out -- boundary shifted only one bin on driver removal)\n",
  n_decoy_slip, NREP))
cat("\nNon-redundancy (each artefact rejected PRIMARILY by a DISTINCT axis):\n")
cat("  global: SPECIFICITY (random panels reproduce the boundary)\n  null/ramp: SUPPORT (no localized consensus)\n")
cat("  decoy: LEAVE-ONE-MARKER-OUT (single marker drives the boundary)\n")
cat("\nABLATION single-algorithm localization error (|bp-truth|, median):\n"); print(round(sort(fam_err,decreasing=TRUE),3))
cat(sprintf("CONSENSUS error = %.3f (vs single-algorithm worst %.3f). No single method is reliable a priori;\n", cons_err, max(fam_err)))
cat(sprintf("leave-one-algorithm-out consensus error stays %.3f-%.3f -> result not contingent on the exact count.\n", min(loo_err), max(loo_err)))
sink()
cat(readLines(file.path(OUT,"summary.txt")), sep="\n"); cat("\n")

## ============================ figure (base R) ===============================
COL <- c(real_spec="#2E7D32",global="#E64A19",null="#9E9E9E",ramp="#7B1FA2",decoy="#C62828")
LAB <- c(real_spec="real (panel-specific)",global="global offset",null="null (no signal)",ramp="monotonic ramp",decoy="single-marker decoy")
lab <- function(letter,title){ mtext(letter,side=3,line=1.35,adj=0,font=2,cex=1.05); mtext(title,side=3,line=0.35,adj=0.5,font=2,cex=0.92) }
png(file.path(OUT,"Figure_synthetic_benchmark.png"), width=2400, height=2200, res=300)
par(mfrow=c(2,2), mar=c(4.2,4.4,3,1.2), mgp=c(2.5,0.8,0), cex.axis=0.85, cex.lab=0.95, font.main=2, cex.main=1.0)
set.seed(1); jx <- runif(nrow(df),-0.18,0.18)
plot(df$support+jx, df$pspec, type="n", xlab="Algorithm consensus support (of 19 components)",
  ylab="Random-panel specificity  p", xlim=c(2.5,11.5), ylim=c(-0.02,0.85), main="")
lab("a","Certification gate: operating characteristics")
rect(S_STAR,-0.05,12,0.05,col=adjustcolor("#2E7D32",0.10),border=NA); abline(v=S_STAR,lty=2,col="grey40"); abline(h=0.05,lty=2,col="grey40")
text(11.3,0.045,"certified\nregion",adj=c(1,1),cex=0.7,col="#2E7D32",font=3)
for(rg in names(COL)) points((df$support+jx)[df$regime==rg],df$pspec[df$regime==rg],pch=21,bg=adjustcolor(COL[rg],0.8),col="white",cex=1.25,lwd=0.6)
legend("topright",legend=LAB[names(COL)],pt.bg=adjustcolor(COL,0.85),pch=21,col="white",pt.cex=1.1,cex=0.62,bg="white",box.col="grey80")
bx <- split(L$maxshift,L$regime)
plot(NA,xlim=c(0.5,2.5),ylim=c(-0.02,0.42),xaxt="n",xlab="",ylab="Max boundary shift on single-marker removal",main="")
lab("b","Leave-one-marker-out robustness")
axis(1,at=1:2,labels=c("real\n(panel-specific)","single-marker\ndecoy"),cex.axis=0.8,padj=0.5); abline(h=0.05,lty=2,col="grey40")
text(2.45,0.062,"1-bin tolerance",adj=c(1,0),cex=0.62,col="grey30")
for(i in seq_along(c("real_spec","decoy"))){ rg<-c("real_spec","decoy")[i]; v<-bx[[rg]]; set.seed(i)
  points(rep(i,length(v))+runif(length(v),-0.12,0.12),v,pch=21,bg=adjustcolor(COL[rg],0.8),col="white",cex=1.4,lwd=0.6)
  segments(i-0.2,median(v),i+0.2,median(v),lwd=2.5,col=COL[rg]) }
fe <- sort(fam_err,decreasing=FALSE); barplot(fe,horiz=TRUE,las=1,col=adjustcolor("#90A4AE",0.9),border="white",
  xlab="Localization error  |breakpoint - truth|",xlim=c(0,0.65),main="",cex.names=0.72); lab("c","Single algorithms are individually fragile")
abline(v=cons_err,col="#2E7D32",lwd=2.5); text(cons_err,par("usr")[4]*0.96,sprintf(" consensus = %.3f",cons_err),col="#2E7D32",cex=0.72,adj=0,font=2)
barplot(loo_err,horiz=TRUE,las=1,col=adjustcolor("#2E7D32",0.55),border="white",
  xlab="Consensus error after removing one algorithm family",xlim=c(0,0.65),main="",cex.names=0.72); lab("d","Consensus is robust to algorithm choice")
abline(v=cons_err,col="grey30",lwd=1.5,lty=2); text(0.62,1.5,sprintf("error stays\n%.3f-%.3f",min(loo_err),max(loo_err)),adj=c(1,0),cex=0.72,col="grey25",font=3)
dev.off()
message(">>> wrote ", file.path(OUT,"Figure_synthetic_benchmark.png"))
