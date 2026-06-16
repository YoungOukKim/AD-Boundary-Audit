# =============================================================================
# celltype_audit_from_raw.R  (RUN FROM REPO ROOT)
# Generates data/celltype_audit/binmeans_{tag}.csv and output/tables/celltype_audit/*
# from the controlled-access SEA-AD MTG h5ad (place in data-external/).
# MANUSCRIPT CLAIM = the concordant LOCALIZATION (CPS ~0.20, Discussion).
# The PC1-piecewise permutation block is an optional diagnostic, NOT a manuscript claim.  — SELF-CONTAINED, clean A-vs-Microglia head-to-head
# -----------------------------------------------------------------------------
# Runs the IDENTICAL 9-algorithm / 19-component consensus + permutation audit on
# TWO cell-type-specific panels through the SAME background-pool null, so the
# panel-specificity numbers are directly comparable:
#   - Astrocyte  : reactive-astrocyte 44-panel  (positive control, expect ~0.207, specific)
#   - Microglia  : microglia-specific 44-panel
# Embedded engine = your boundary_detection.R (validated: reproduces astrocyte
# primary 0.207 via Hier/k-means/variance-jump, DP-PC2 0.388). segmented warnings
# are suppressed (they only touch the version-sensitive Segmented_PC rows).
#
# OUTPUT (paste the summary + 2 multimethod files back):
#   <OUT_DIR>/headtohead_summary.txt
#   <OUT_DIR>/multimethod_Astrocyte.csv , multimethod_Microglia.csv
#   <OUT_DIR>/binmeans_Astrocyte.csv    , binmeans_Microglia.csv
# =============================================================================
suppressPackageStartupMessages({ library(rhdf5); library(data.table) })

## ====================== CONFIG ======================
H5AD_PATH <- "data-external/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad"  # controlled-access SEA-AD h5ad
OUT_DIR   <- "output/tables/celltype_audit"   # relative; run from repo root
N_PERM <- 1000L; BG_POOL_SIZE <- 3000L; MIN_DETECT <- 0.02; X_IS_LOGNORM <- TRUE
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

ASTRO_PANEL <- c("FKBP5","GFAP","BCL2","MERTK","SOX9","SERPINA3","TP53","STAT3","OCLN","AQP4",
 "NFE2L2","SOCS3","FTH1","FTL","TREM2","VCAM1","C3","IL1B","STAT1","CLU","TNF","SELE","FBLN5",
 "CXCL12","GDNF","CR1","CD33","CXCL16","VWF","PECAM1","CASP3","IL6","MFGE8","NGF","NFKB1",
 "ICAM1","NEFL","VIM","CCL2","APP","BDNF","APOE","TFRC","HMOX1")
MICRO_PANEL <- c("P2RY12","P2RY13","CX3CR1","TMEM119","CSF1R","SELPLG","GPR34","OLFML3",
 "TREM2","TYROBP","APOE","CST7","LPL","CD9","ITGAX","CLEC7A","SPP1","GPNMB",
 "CD33","MS4A6A","MS4A4A","INPP5D","ABI3","PLCG2","BIN1","PICALM",
 "C1QA","C1QB","C1QC","C3","IL1B","TNF","NFKB1","CCL3","STAT1",
 "CTSB","CTSD","CTSS","GRN","MERTK","GAS6","FTL","FTH1","IRF8")

## =================== EMBEDDED ENGINE (your boundary_detection.R) ============
piecewise_fit <- function(x, y) {
  if (length(x) < 4 || sd(y, na.rm=TRUE) < 1e-10) return(list(bp=NA)); best_bp<-NA; best_sse<-Inf
  for (bp in seq(x[2]+0.02, x[length(x)-1]-0.02, length.out=50)) {
    left<-x<=bp; right<-x>bp; if (sum(left)<2||sum(right)<2) next
    tryCatch({ ll<-lm(y[left]~x[left]); rr<-lm(y[right]~x[right]); sse<-sum(residuals(ll)^2)+sum(residuals(rr)^2)
      if (sse<best_sse){best_sse<-sse;best_bp<-bp} }, error=function(e) NULL) }
  list(bp=best_bp) }
cusum_bp <- function(v,x) x[which.max(abs(cumsum(v-mean(v))))]
dp_changepoint <- function(x,y,min_seg=2){ n<-length(x); bi<-NA; bs<-Inf
  for (k in (min_seg+1):(n-min_seg)) tryCatch({ l<-lm(y[1:k]~x[1:k]); r<-lm(y[(k+1):n]~x[(k+1):n])
    s<-sum(residuals(l)^2)+sum(residuals(r)^2); if(s<bs){bs<-s;bi<-k} },error=function(e) NULL)
  if(is.na(bi)) return(NA); (x[bi]+x[bi+1])/2 }
.cluster_boundaries <- function(cl,a){ b<-c(); for(i in 2:length(cl)) if(cl[i]!=cl[i-1]) b<-c(b,(a[i-1]+a[i])/2); b }
.HAVE_SEG <- requireNamespace("segmented", quietly=TRUE)
detect_boundaries <- function(expr_bin, axis_vals) {
  expr_bin<-as.matrix(expr_bin); stopifnot(nrow(expr_bin)==length(axis_vals))
  res<-list(); pca<-prcomp(expr_bin,center=TRUE,scale.=FALSE); npc<-min(3,ncol(pca$x))
  for(k in 1:npc){f<-piecewise_fit(axis_vals,pca$x[,k]); if(!is.na(f$bp)) res[[paste0("PCA_PC",k)]]<-f$bp}
  for(k in 1:npc) res[[paste0("CUSUM_PC",k)]]<-cusum_bp(pca$x[,k],axis_vals)
  gaps<-sapply(2:nrow(expr_bin),function(i) sqrt(sum((expr_bin[i,]-expr_bin[i-1,])^2)))
  gp<-(axis_vals[-length(axis_vals)]+axis_vals[-1])/2; res[["Mahal_max_gap"]]<-gp[which.max(gaps)]
  hc<-hclust(dist(expr_bin,method="euclidean"),method="ward.D2"); b3<-.cluster_boundaries(cutree(hc,k=3),axis_vals)
  if(length(b3)>=1) res[["Hier_k3_bp1"]]<-b3[1]; if(length(b3)>=2) res[["Hier_k3_bp2"]]<-b3[2]
  set.seed(42); bk<-.cluster_boundaries(kmeans(expr_bin,centers=3,nstart=10)$cluster,axis_vals)
  if(length(bk)>=1) res[["KMeans_k3_bp1"]]<-bk[1]; if(length(bk)>=2) res[["KMeans_k3_bp2"]]<-bk[2]
  for(k in 1:npc){bp<-dp_changepoint(axis_vals,pca$x[,k]); if(!is.na(bp)) res[[paste0("DP_PC",k)]]<-bp}
  vd<-diff(apply(expr_bin,1,var)); if(length(vd)>0){j<-which.max(abs(vd)); res[["Variance_jump"]]<-(axis_vals[j]+axis_vals[j+1])/2}
  if(length(axis_vals)>=4){d0<-sapply(1:nrow(expr_bin),function(i) sqrt(sum((expr_bin[i,]-expr_bin[1,])^2)))
    res[["Trajectory_inflection"]]<-axis_vals[which.max(abs(diff(diff(d0))))+1]}
  for(k in 1:npc) tryCatch({
    if(.HAVE_SEG){ d<-data.frame(A=axis_vals,PC=pca$x[,k])
      sg<-suppressWarnings(segmented::segmented(lm(PC~A,data=d),seg.Z=~A,psi=stats::median(axis_vals)))
      res[[paste0("Segmented_PC",k)]]<-sg$psi[1,"Est."]
    } else { y<-pca$x[,k]; bs<-Inf;ps<-NA; for(kk in seq(axis_vals[2],axis_vals[length(axis_vals)-1],length.out=200)){
      x1<-pmin(axis_vals,kk);x2<-pmax(axis_vals-kk,0);s<-sum(resid(lm(y~x1+x2))^2);if(s<bs){bs<-s;ps<-kk}}; res[[paste0("Segmented_PC",k)]]<-ps}
  },error=function(e) NULL)
  d<-data.frame(algorithm=names(res),breakpoint=unlist(res),row.names=NULL); d[order(d$breakpoint),] }
consensus <- function(bp){ bp<-bp[!is.na(bp)]; dd<-sapply(bp,function(c0)sum(abs(bp-c0)<=0.05)); list(bp=bp[which.max(dd)],n=max(dd)) }

## =================== shared obs ===================
message(">>> reading obs ...")
sc_idx<-h5read(H5AD_PATH,"obs/Subclass"); sc_cat<-h5read(H5AD_PATH,"obs/__categories/Subclass")
cell_type<-sc_cat[sc_idx+1L]; cps<-as.numeric(h5read(H5AD_PATH,"obs/Continuous Pseudo-progression Score"))
va<-as.character(h5read(H5AD_PATH,"var/_index"))
indptr<-h5read(H5AD_PATH,"X/indptr",bit64conversion="double"); nC<-length(indptr)-1L
zf<-function(v){s<-sd(v); if(is.na(s)||s==0) return(rep(0,length(v))); (v-mean(v))/s}

read_cells <- function(keep, genes){            # fast block read, keep-filtered
  gim<-match(genes,va)-1L
  E<-matrix(0,sum(keep),length(genes),dimnames=list(NULL,genes)); ro<-integer(nC); ro[keep]<-seq_len(sum(keep))
  blk<-200000L
  for(s0 in seq(1L,nC,by=blk)){ e0<-min(s0+blk-1L,nC); if(!any(keep[s0:e0])) next
    sp<-indptr[s0]; cnt<-indptr[e0+1L]-sp; if(cnt<=0) next
    ci<-h5read(H5AD_PATH,"X/indices",start=sp+1L,count=cnt,bit64conversion="double")
    cd<-h5read(H5AD_PATH,"X/data",start=sp+1L,count=cnt)
    for(i in s0:e0){ if(!keep[i]) next; a<-indptr[i]-sp+1L; b<-indptr[i+1L]-sp; if(b<a) next
      h<-match(gim,ci[a:b]); ok<-!is.na(h); if(any(ok)) E[ro[i],ok]<-cd[a:b][h[ok]] }
    cat(sprintf("\r    read %5.1f%%",e0/nC*100)); flush.console() }
  cat("\n"); if(!X_IS_LOGNORM) E<-log1p(E); E }

run_celltype <- function(tag, keep, panel) {
  message(sprintf(">>> [%s] cells=%s", tag, format(sum(keep),big.mark=",")))
  set.seed(42)
  pin<-panel[panel%in%va]; bg<-sample(setdiff(va,pin), min(BG_POOL_SIZE,length(setdiff(va,pin))))
  E<-read_cells(keep, c(pin,bg)); mcps<-cps[keep]
  Ez<-apply(E,2,zf); cb<-cut(mcps,breaks=seq(0.1,1.0,0.1),include.lowest=TRUE,right=FALSE)
  binm<-function(cols){ ag<-aggregate(Ez[,cols,drop=FALSE],list(b=cb),mean); ag<-ag[order(ag$b),]
    list(M=as.matrix(ag[,-1,drop=FALSE]), center=as.numeric(tapply(mcps,cb,mean)[as.character(ag$b)]),
         n=as.integer(table(cb)[as.character(ag$b)]), bin=as.character(ag$b)) }
  pb<-binm(pin); axis<-pb$center
  det<-colMeans(E[,bg,drop=FALSE]>0); bgk<-bg[det>=MIN_DETECT]; bgB<-binm(bgk)$M
  write.csv(data.frame(cps_bin=pb$bin,cps_center=pb$center,n_cells=pb$n,pb$M,check.names=FALSE),
            file.path("data/celltype_audit",paste0("binmeans_",tag,".csv")), row.names=FALSE)
  obs<-detect_boundaries(pb$M, axis); write.csv(obs, file.path(OUT_DIR,paste0("multimethod_",tag,".csv")), row.names=FALSE)
  primary<-round(obs$breakpoint[obs$algorithm=="Hier_k3_bp1"],3)        # Ward early boundary (= published primary)
  supp<-sum(abs(obs$breakpoint-primary)<=0.05,na.rm=TRUE)               # methods agreeing within +/-0.05 of primary
  # PUBLISHED specificity statistic: PC1-piecewise breakpoint of each random panel
  k<-length(pin); nullpc1<-numeric(N_PERM)
  for(p in seq_len(N_PERM)){ cols<-sample(seq_len(ncol(bgB)),k); rp<-detect_boundaries(bgB[,cols,drop=FALSE],axis)
    v<-rp$breakpoint[rp$algorithm=="PCA_PC1"]; nullpc1[p]<-if(length(v)) v else NA
    if(p%%200==0){cat(sprintf("\r    [%s] perm %d/%d",tag,p,N_PERM));flush.console()} }
  cat("\n")
  spec_n<-sum(abs(nullpc1-primary)<=0.05,na.rm=TRUE)                    # near primary  (LOW = panel-specific)
  late_n<-sum(abs(nullpc1-0.54)<=0.05,na.rm=TRUE)                       # near late 0.54 (sanity; astro ~40%)
  rm(E); gc()
  list(tag=tag, n=sum(keep), npanel=k, bg=length(bgk), primary=primary, supp=supp,
       spec_n=spec_n, spec=spec_n/N_PERM, late_n=late_n, late=late_n/N_PERM,
       nmin=round(min(nullpc1,na.rm=TRUE),2), nmed=round(median(nullpc1,na.rm=TRUE),2),
       nmax=round(max(nullpc1,na.rm=TRUE),2), axis=axis, obs=obs)
}

A <- run_celltype("Astrocyte", cell_type=="Astrocyte", ASTRO_PANEL)
M <- run_celltype("Microglia", grepl("Micro",cell_type), MICRO_PANEL)

sink(file.path(OUT_DIR,"headtohead_summary.txt"))
cat("HEAD-TO-HEAD: 4-AXIS AUDIT, IDENTICAL PERMUTATION NULL\n")
cat("======================================================\n")
cat(sprintf("N_PERM=%d  BG_POOL=%d  MIN_DETECT=%.2f\n\n", N_PERM, BG_POOL_SIZE, MIN_DETECT))
fmt<-function(r) sprintf("%-11s n=%6d panel=%2d bg=%4d | primary(Ward) CPS=%.3f (%2d/19 methods) | null PC1-piecewise %.2f-%.2f(med %.2f) | specificity@primary=%d/%d (%.3f) | late@0.54=%d/%d (%.3f)\n",
  r$tag, r$n, r$npanel, r$bg, r$primary, r$supp, r$nmin, r$nmax, r$nmed, r$spec_n, N_PERM, r$spec, r$late_n, N_PERM, r$late)
cat(fmt(A)); cat(fmt(M))
cat("\nStatistic = PUBLISHED PC1-piecewise null (manuscript: astro 0/1000 @0.21, ~40% @0.54).\n")
cat("specificity@primary LOW => boundary NOT reproduced by random panels => panel-specific.\n\n")
cat("--- Astrocyte per-algorithm ---\n"); print(A$obs,row.names=FALSE)
cat("\n--- Microglia per-algorithm ---\n"); print(M$obs,row.names=FALSE)
sink()
cat(readLines(file.path(OUT_DIR,"headtohead_summary.txt")), sep="\n")
message("\n>>> DONE. Send me headtohead_summary.txt (+ the two multimethod CSVs).")
