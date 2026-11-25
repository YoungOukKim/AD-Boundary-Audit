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

