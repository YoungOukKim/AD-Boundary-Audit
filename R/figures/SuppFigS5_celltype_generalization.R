# =============================================================================
# SuppFigS4_celltype_generalization.R   (ggplot, matches manuscript theme_paper)
# The identical multi-algorithm consensus audit on an independent microglia
# 44-panel localizes a transition concordant with the reactive-astrocyte panel.
#   (a) 19-component boundary localization, both panels (primary cluster ~0.20)
#   (b) panel PC1 trajectory with localized early transition
# Inputs: data/celltype_audit/binmeans_{tag}.csv, output/.../multimethod_{tag}.csv
# Output: output/figures/SuppFigure_S5.{png,tif}
# =============================================================================
source("R/figures/utils.R")
CTD<-"data/celltype_audit"; MM<-"output/tables/celltype_audit"
prim<-c(Astrocyte=0.207, Microglia=0.201)
bpf<-function(t){x<-read.csv(file.path(MM,paste0("multimethod_",t,".csv"))); x$panel<-t
  x$prim<-prim[t]; x$incl<-abs(x$breakpoint-prim[t])<=0.05; x}
binf<-function(t){b<-read.csv(file.path(CTD,paste0("binmeans_",t,".csv")),check.names=FALSE)
  list(cps=b$cps_center,M=as.matrix(b[,!(names(b)%in%c("cps_bin","cps_center","n_cells"))]))}
pc1<-function(M){p<-prcomp(M,center=TRUE,scale.=FALSE)$x[,1]; sign(cor(seq_along(p),p))*p}

bp<-rbind(bpf("Astrocyte"),bpf("Microglia"))
bp$panel<-factor(bp$panel,levels=c("Microglia","Astrocyte"))
lab<-data.frame(panel=factor(c("Astrocyte","Microglia"),levels=c("Microglia","Astrocyte")),
  prim=prim, txt=c(sprintf("CPS 0.207  (%d/19)",sum(bpf("Astrocyte")$incl)),
                   sprintf("CPS 0.201  (%d/19)",sum(bpf("Microglia")$incl))))
pa<-ggplot(bp,aes(breakpoint,panel))+
  geom_vline(data=lab,aes(xintercept=prim),colour=col_master,linetype="dashed",linewidth=0.6)+
  geom_point(aes(colour=incl),shape=124,size=5,stroke=1.1)+
  geom_text(data=lab,aes(x=prim,y=panel,label=txt),colour=col_master,hjust=-0.06,vjust=-1.2,size=3.1)+
  scale_colour_manual(values=c(`TRUE`=DOWN,`FALSE`=GREY),
    labels=c(`TRUE`="primary cluster (within 0.05 CPS)",`FALSE`="other components"),name=NULL)+
  scale_x_continuous(limits=c(0,1))+
  labs(title="Consensus boundary localization (19 algorithm components)",
       x="Continuous Pseudo-progression Score (CPS)",y=NULL)+
  theme_paper+theme(legend.position="top",legend.text=element_text(size=8))

mg<-binf("Microglia"); as<-binf("Astrocyte")
df<-rbind(data.frame(cps=mg$cps,pc1=pc1(mg$M),panel="Microglia panel"),
          data.frame(cps=as$cps,pc1=pc1(as$M),panel="Astrocyte panel"))
df$panel<-factor(df$panel,levels=c("Microglia panel","Astrocyte panel"))
pb<-ggplot(df,aes(cps,pc1,colour=panel,linetype=panel))+
  geom_vline(xintercept=prim["Microglia"],colour=col_master,linetype="dashed",linewidth=0.6)+
  geom_line(linewidth=1)+geom_point(size=1.7)+
  scale_colour_manual(values=c(`Microglia panel`=DOWN,`Astrocyte panel`=GREY),name=NULL)+
  scale_linetype_manual(values=c(`Microglia panel`="solid",`Astrocyte panel`="32"),name=NULL)+
  scale_x_continuous(limits=c(0,1))+
  labs(title="Panel trajectory and localized early transition",
       x="Continuous Pseudo-progression Score (CPS)",y="Panel PC1 (a.u.)")+
  theme_paper+theme(legend.position=c(0.015,0.99),legend.justification=c(0,1),
                    legend.background=element_rect(fill=alpha("white",0.7),colour=NA),legend.text=element_text(size=8))

fig<-(pa/pb)+plot_layout(heights=c(0.85,1.1))+
  plot_annotation(tag_levels="a")& theme(plot.tag=element_text(face="bold",size=16.1))
save_fig(fig,"SuppFigure_S5.png",13,11.595)
tiff(file.path(OUTF,"SuppFigure_S5.tif"),width=13,height=11.595,units="in",res=600,compression="lzw"); print(fig); dev.off()
