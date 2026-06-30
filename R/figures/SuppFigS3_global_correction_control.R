# =============================================================================
# SuppFigS3_global_correction_control.R   (ggplot, matches theme_paper)
# Housekeeping negative control: 14 constitutive HK genes stay FLAT before/after
# the genome-wide correction -> it removes only the shared offset, not signal.
#   (a) HK trajectories pre/post   (b) effect magnitude vs signal-gene scale
# Input: data/metabolic_correction/panel_correction_z.csv
# Output: output/figures/SuppFigure_S3.{png,tif}
# =============================================================================
source("R/figures/utils.R")
IN<-"data/metabolic_correction/panel_correction_z.csv"; CC<-"Astro"; REG<-"MTG"
d<-read.csv(IN,check.names=FALSE); d<-d[d$cell_class==CC & d$panel=="HK" & d$region==REG,]
long<-rbind(data.frame(gene=d$gene,cps=d$bin,z=d$z_pre, state="before correction"),
            data.frame(gene=d$gene,cps=d$bin,z=d$z_post,state="after correction"))
long$state<-factor(long$state,levels=c("before correction","after correction"))
pa<-ggplot(long,aes(cps,z,group=interaction(gene,state),colour=state))+
  annotate("rect",xmin=0.1,xmax=1,ymin=0.40,ymax=0.62,fill=UP,alpha=0.07)+
  annotate("rect",xmin=0.1,xmax=1,ymin=-0.62,ymax=-0.40,fill=UP,alpha=0.07)+
  annotate("text",x=0.55,y=0.53,label="signal-gene range (~0.48)",colour=UP,size=3)+
  geom_hline(yintercept=0,colour="grey80",linewidth=0.3)+
  geom_line(linewidth=0.5,alpha=0.65)+
  scale_colour_manual(values=c(`before correction`=GREY,`after correction`=DOWN),name=NULL)+
  coord_cartesian(ylim=c(-0.62,0.62))+
  labs(title="Housekeeping trajectories",x="CPS",y="z-score (per CPS bin)")+
  theme_paper+theme(legend.position=c(0.02,0.04),legend.justification=c(0,0),
                    legend.background=element_rect(fill=alpha("white",0.7),colour=NA),legend.text=element_text(size=8))
dz<-function(col){e<-sapply(unique(d$gene),function(g){v<-d[d$gene==g,];v<-v[order(v$bin),]
  mean(v[[col]][v$bin>=0.45])-mean(v[[col]][v$bin<=0.35])});mean(abs(e))}
bars<-data.frame(grp=factor(c("HK\nbefore","HK\nafter","signal\ngenes"),
                            levels=c("HK\nbefore","HK\nafter","signal\ngenes")),
                 val=c(dz("z_pre"),dz("z_post"),0.48),fill=c("before","after","signal"))
pb<-ggplot(bars,aes(grp,val,fill=fill))+geom_col(width=0.62)+
  geom_hline(yintercept=0.10,linetype="dotted",colour="grey45")+
  annotate("text",x=0.55,y=0.125,label="conserved-call threshold 0.10",hjust=0,size=2.7,colour="grey40")+
  geom_text(aes(label=ifelse(val>0.4,"~0.48",sprintf("%.3f",val))),vjust=-0.5,size=3.1)+
  scale_fill_manual(values=c(before=GREY,after=DOWN,signal=UP),guide="none")+
  coord_cartesian(ylim=c(0,0.55))+
  labs(title="Effect magnitude",x=NULL,y=expression("mean |"*Delta*"z| across transition"))+
  theme_paper
fig<-(pa|pb)+plot_layout(widths=c(1.35,1))+
  plot_annotation(tag_levels="a")& theme(plot.tag=element_text(face="bold",size=16.1))
save_fig(fig,"SuppFigure_S3.png",13,5.417)
tiff(file.path(OUTF,"SuppFigure_S3.tif"),width=13,height=5.417,units="in",res=600,compression="lzw"); print(fig); dev.off()
