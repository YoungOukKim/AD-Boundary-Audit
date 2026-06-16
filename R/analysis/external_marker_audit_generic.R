# =============================================================================
# external_marker_audit_generic.R   (RUN FROM REPO ROOT, on the data machine)
# Cohort-GENERIC, donor-level, trajectory-aware panel audit. Same corrected
# method as external_marker_audit.R, but parameterised so the SAME engine runs
# on an INDEPENDENT cohort (e.g., Leng 2021 entorhinal cortex, GSE147528, Braak
# axis) as well as SEA-AD (CPS axis). This is how the corrected method gets an
# independent-cohort run.
#
# DIRECTION is bin-FREE: Spearman(donor mean expr vs donor progression).
# SHAPE (monotonic / non-monotonic) + donor-binned boundary + bin-sensitivity
# (K in {6,8,9,10,12}) reported. Output tagged by COHORT.
#
# ---- EDIT THIS BLOCK PER COHORT ----
H5AD       <- "data-external/Leng2021_EC_GSE147528.h5ad"   # your file
COHORT     <- "LengEC"                                       # tag for outputs
PROG_OBS   <- "Braak"            # progression axis obs key (e.g. "Braak"; SEA-AD: "Continuous Pseudo-progression Score")
PROG_TYPE  <- "ordinal"          # "ordinal" (Braak labels -> int) or "continuous" (use as-is)
DONOR_OBS  <- "SampleID"         # donor/sample identity obs key
SUBCLASS   <- "clusterCellType"  # cell-type obs key (SEA-AD: "Subclass")
MIN_PROG   <- NA                 # e.g. 0.1 for CPS; NA to keep all
# ------------------------------------
suppressMessages(library(rhdf5))
OUTD<-"output/tables/external_markers"; dir.create(OUTD,recursive=TRUE,showWarnings=FALSE)
## ---- boundary-detection engine embedded inline (self-contained; no source needed) ----
# boundary_detection.R -- 9-algorithm / 19-component consensus engine (set.seed(42)).
# Sourced by audit/figure scripts. Standalone; no side effects.
.HAVE_SEG <- requireNamespace("segmented", quietly = TRUE)
## --- engine (boundary_detection.R, embedded verbatim; set.seed(42)) ---------
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
## ---- end engine ----
FLAT_RHO<-0.20; PCUT<-0.05; MIN_CELLS_DONOR<-10; KSET<-c(6,8,9,10,12)
MAP <- list(Astrocyte=c("GFAP","CHI3L1","PTGDS","HMOX1"),
            Microglia=c("TREM2"),
            Neuron=c("NPTX1","NPTX2","NPTXR","NRGN","SNAP25","VAMP2","VGF","SCG2","CNTN2"))
NONNEU<-"Astro|Micro|PVM|Oligo|OPC|Endo|VLMC|Peri|SMC|Fibro|Immune|Macro"

.cat<-function(key){p<-paste0("obs/",key);cp<-paste0("obs/__categories/",key)
  tryCatch({co<-as.integer(h5read(H5AD,p));ca<-as.character(h5read(H5AD,cp));ca[co+1L]},
           error=function(e) as.character(h5read(H5AD,p)))}
roman2int<-function(s){m<-c(I=1,II=2,III=3,IV=4,V=5,VI=6); ifelse(s %in% names(m), m[s], suppressWarnings(as.numeric(s)))}
to_num<-function(v){ if(PROG_TYPE=="continuous") return(as.numeric(v))
  x<-gsub("[^0-9IVX]","",toupper(as.character(v)))           # strip 'Braak ', etc.
  num<-suppressWarnings(as.numeric(x)); rom<-roman2int(x)
  ifelse(!is.na(num),num,rom) }

sub<-.cat(SUBCLASS); donor<-.cat(DONOR_OBS); prog<-to_num(.cat(PROG_OBS))
genes_all<-as.character(h5read(H5AD,"var/_index"))
indptr<-as.numeric(h5read(H5AD,"X/indptr")); ncell<-length(indptr)-1L
allg<-unique(unlist(MAP)); gidx<-match(allg,genes_all); names(gidx)<-allg
if(any(is.na(gidx))) cat("NOT in var:",paste(allg[is.na(gidx)],collapse=", "),"\n"); gidx<-gidx[!is.na(gidx)]
sel_pre <- !is.na(prog) & (if(is.na(MIN_PROG)) TRUE else prog>=MIN_PROG) &
       (grepl("Astro",sub,ignore.case=TRUE)|grepl("Micro|PVM",sub,ignore.case=TRUE)|!grepl(NONNEU,sub,ignore.case=TRUE))
pidx0<-gidx-1L; G<-names(gidx)                      # 0-based panel gene idx
indptr<-h5read(H5AD,"X/indptr",bit64conversion="double"); nC<-length(indptr)-1L
sel<-which(sel_pre); kset<-logical(nC);kset[sel]<-TRUE; rmap<-integer(nC);rmap[sel]<-seq_along(sel)
M<-matrix(0,length(sel),length(pidx0)); bs<-100000L; cat("reading expression ...\n")
for(s0 in seq(1L,nC,by=bs)){e0<-min(s0+bs-1L,nC);sp<-indptr[s0];cnt<-indptr[e0+1L]-sp;if(cnt<=0)next
  ci<-h5read(H5AD,"X/indices",start=sp+1L,count=cnt,bit64conversion="double")
  cd<-h5read(H5AD,"X/data",start=sp+1L,count=cnt)
  for(k in which(kset[s0:e0])){g<-s0+k-1L;a<-indptr[g]-sp+1L;b<-indptr[g+1L]-sp;if(b<a)next
    h<-match(pidx0,ci[a:b]);ok<-!is.na(h);M[rmap[g],ok]<-cd[a:b][h[ok]]}
  cat(sprintf("  ...%d/%d\r",e0,nC))}; cat("\n")
colnames(M)<-G; sub_s<-sub[sel]; donor_s<-donor[sel]; prog_s<-prog[sel]

shape_flag<-function(v,ctr){im<-which.max(v);it<-which.min(v);n<-length(v);intr<-function(i)i>1&&i<n
  if(intr(im)&&v[im]>max(v[1],v[n])+0.15) return(sprintf("non-monotonic(peak@%.2f)",ctr[im]))
  if(intr(it)&&v[it]<min(v[1],v[n])-0.15) return(sprintf("non-monotonic(trough@%.2f)",ctr[it]))
  "monotonic"}
binsign<-function(dx,dp,K){br<-seq(min(dp),max(dp),length.out=K+1);ctr<-(br[-1]+br[-(K+1)])/2
  bin<-cut(dp,br,include.lowest=TRUE,labels=FALSE);z<-scale(dx)[,1];m<-tapply(z,bin,mean,na.rm=TRUE)
  cc<-ctr[as.integer(names(m))];mid<-median(dp)
  sign(mean(m[cc>=mid],na.rm=TRUE)-mean(m[cc<mid],na.rm=TRUE))}

verd<-list();traj<-list()
for(ct in names(MAP)){
  ctm<-if(ct=="Neuron") !grepl(NONNEU,sub_s,ignore.case=TRUE) else if(ct=="Astrocyte") grepl("Astro",sub_s,ignore.case=TRUE) else grepl("Micro|PVM",sub_s,ignore.case=TRUE)
  genes<-intersect(MAP[[ct]],colnames(M)); if(sum(ctm)<30||!length(genes)) next
  d<-donor_s[ctm]; nC2<-table(d); ok_d<-names(nC2)[nC2>=MIN_CELLS_DONOR]
  if(length(ok_d)<5){cat(ct,": donors<5, skip\n"); next}
  Dm<-rowsum(M[ctm,genes,drop=FALSE],d)[ok_d,,drop=FALSE]/as.numeric(nC2[ok_d])
  Dp<-tapply(prog_s[ctm],d,function(x) x[1])[ok_d]
  K0<-min(9,length(unique(Dp))); br<-seq(min(Dp),max(Dp),length.out=K0+1); ctr<-(br[-1]+br[-(K0+1)])/2
  bin<-cut(Dp,br,include.lowest=TRUE,labels=FALSE); Z<-scale(Dm)
  M<-t(sapply(1:K0,function(bb) colMeans(Z[bin==bb,,drop=FALSE],na.rm=TRUE)))
  if(length(genes)==1) M<-matrix(M,ncol=1,dimnames=list(NULL,genes))
  ob<-tryCatch(detect_boundaries(M,ctr),error=function(e)NULL)
  prim<-if(!is.null(ob)) round(ob$breakpoint[ob$algorithm=="Hier_k3_bp1"][1],3) else NA
  for(g in genes){
    ct_t<-suppressWarnings(cor.test(Dm[,g],Dp,method="spearman")); rho<-unname(ct_t$estimate); pv<-ct_t$p.value
    dir<-if(!is.na(rho)&&abs(rho)>=FLAT_RHO&&pv<PCUT) ifelse(rho>0,"up","down") else "flat"
    st<-sapply(KSET,function(K) tryCatch(binsign(Dm[,g],Dp,K),error=function(e)NA))
    verd[[length(verd)+1]]<-data.frame(cohort=COHORT,gene=g,cell_type=ct,n_donors=nrow(Dm),n_cells=sum(ctm),
      spearman_rho=round(rho,3),spearman_p=signif(pv,3),dir=dir,
      shape=shape_flag(M[,g],ctr),boundary=prim,
      bin_stability=paste0(sum(st==sign(rho),na.rm=TRUE),"/",length(KSET)),stringsAsFactors=FALSE)
    traj[[length(traj)+1]]<-data.frame(cohort=COHORT,gene=g,cell_type=ct,prog_bin=round(ctr,3),mean_z=round(M[,g],3))
  }
}
V<-do.call(rbind,verd)
write.csv(V,file.path(OUTD,paste0(COHORT,"_panel_verdict.csv")),row.names=FALSE)
write.csv(do.call(rbind,traj),file.path(OUTD,paste0(COHORT,"_panel_trajectory.csv")),row.names=FALSE)
cat(sprintf("\n=== %s donor-level trajectory audit (progression=%s) ===\n",COHORT,PROG_OBS))
print(V[,c("gene","cell_type","spearman_rho","spearman_p","dir","shape","boundary","bin_stability")],row.names=FALSE)
cat("\nDirection bin-FREE (Spearman vs progression). Compare to panel_verdict.csv (literature) and SEA-AD run.\n")
