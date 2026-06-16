# =============================================================================
# verify_correction_consistency.R   (RUN FROM REPO ROOT)
# Runs DIRECTLY on your canonical correction output -- no h5ad, no generation.
# The offset column in this file is already the genome-wide ("mean across all
# genes") correction your pipeline computed (Methods), so this settles whether
# BOTH hold under that ONE correction:
#   (a) housekeeping flat after correction, (b) metabolic conservation drops.
# Input: data/metabolic_correction/panel_correction_z.csv
#   columns region, cell_class, panel(HK|METAB), gene, bin, z_pre, offset, z_post
# =============================================================================
IN<-"data/metabolic_correction/panel_correction_z.csv"; CC<-"Astro"
d<-read.csv(IN,check.names=FALSE); d<-d[d$cell_class==CC,]
cat(sprintf("offset range (small & stable => genome-wide mean): %.3f .. %.3f\n",
            min(d$offset),max(d$offset)))
eff<-function(panel,col){ s<-d[d$panel==panel,]; g<-unique(s$gene)
  f<-function(reg) sapply(g,function(x){v<-s[s$gene==x & s$region==reg,]
      mean(v[[col]][v$bin>=0.45])-mean(v[[col]][v$bin<=0.35])})
  data.frame(gene=g,MTG=f("MTG"),A9=f("A9")) }
cons<-function(e) sum(abs(e$MTG)>0.10 & abs(e$A9)>0.10 & sign(e$MTG)==sign(e$A9))
absm<-function(e) c(MTG=mean(abs(e$MTG)),A9=mean(abs(e$A9)))
for(p in c("HK","METAB")){ pre<-eff(p,"z_pre"); post<-eff(p,"z_post")
  cat(sprintf("\n%s (n=%d): mean|dz| pre[MTG %.3f A9 %.3f] post[MTG %.3f A9 %.3f] | conserved %d -> %d\n",
      p,nrow(pre),absm(pre)["MTG"],absm(pre)["A9"],absm(post)["MTG"],absm(post)["A9"],cons(pre),cons(post))) }
hk<-eff("HK","z_post"); mp<-eff("METAB","z_pre"); mq<-eff("METAB","z_post")
hk_ok<-all(absm(hk)<0.10); met_ok<-cons(mq)<cons(mp)
cat(sprintf("\nVERDICT: (a) HK flat after = %s | (b) METAB conservation drops = %s\n",
            ifelse(hk_ok,"YES","NO"),ifelse(met_ok,"YES","NO")))
cat(if(hk_ok&&met_ok) ">>> PASS: both hold under one genome-wide correction. Section sound.\n"
    else ">>> CHECK: see which condition failed above.\n")
