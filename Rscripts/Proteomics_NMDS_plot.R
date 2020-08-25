library(vegan)
library(fields)
library(MASS)
library(ggplot2)
library(ggforce)

cols<-c("#ED8636","#0433FF")

#####  Read the data for RTE protoemics  
rte_bio<-read.csv('Data/RTE_BioRep12.csv', header=T)
rte_bio1=t(rte_bio)
rte_bio2<-rte_bio1[-c(1:2),]
rte_bio3<-apply(rte_bio2,2,as.numeric)

#run NMDS
rtebio1<-metaMDS(rte_bio3)
ordiplot(rtebio1,type= "text", display="sites" )


## Plot the results
rte1 = as.data.frame(scores(rtebio1))
rte1$site<-c(rep("N", times=3),rep("O", times=3),rep("N", times=3),rep("O", times=3))
rte1$Origin<-c(rep("Nearshore", times=6),rep("Offshore", times=6))

ggplot(rte1, aes(x=NMDS1,y=NMDS2,color=Origin,shape=site)) + 
    geom_point(size=3.5) + theme_bw() +
    theme(legend.text = element_text(size=12), legend.title = element_text(size=13))+
    scale_shape_manual(values=c(19,17), name="Transplant site",labels=c("Nearshore", "Offshore"))+
    xlab("Axis 1") + ylab("Axis 2")+
    theme(axis.text=element_text(size=10))+
    scale_color_manual(values=paste0(cols,"CC"))+
    guides(shape = guide_legend(order = 2),color=guide_legend(order=1))+
    geom_ellipse(aes(x0=-0.06, y0=-0.11, a=.12, b=.12, angle=10), color=cols[1])+
    geom_ellipse(aes(x0=0, y0=0.095, a=.33, b=0.06, angle=207), color=cols[2])+
    theme(axis.text = element_text(size=12), axis.title = element_text(size=13))

ggsave("./Output/RTE_nmds.ellipse2.pdf", height = 5, width = 6.8)

### legends inside the plot for PNAS ####
ggplot(rte1, aes(x=NMDS1,y=NMDS2,color=Origin,shape=site)) + 
    geom_point(size=3.5) + theme_bw() +
    theme(legend.text = element_text(size=12), legend.title = element_text(size=13))+
    scale_shape_manual(values=c(19,17), name="Transplant site",labels=c("Nearshore", "Offshore"))+
    xlab("Axis 1") + ylab("Axis 2")+
    theme(axis.text=element_text(size=10))+
    scale_color_manual(values=paste0(cols,"CC"))+
    theme(legend.position = "none")+ 
    geom_ellipse(aes(x0=-0.06, y0=-0.11, a=.12, b=.12, angle=10), color=cols[1])+
    geom_ellipse(aes(x0=0, y0=0.095, a=.33, b=0.06, angle=207), color=cols[2])+
    theme(axis.text = element_text(size=12), axis.title = element_text(size=13))+
    annotate(geom="text", x=.18, y=.24, label="Origin",size=3.5, hjust=0)+
    annotate("point", x = .19, y =0.22 , colour = cols[1], size=3.5)+
    annotate(geom="text", x=.21, y=.22, label="Nearshore",color =cols[1], size=3.5, hjust=0)+
    annotate("point",  x = .19, y =0.2 ,colour = cols[2], size=3.5)+
    annotate(geom="text", x=.21, y=.2, label="Offshore",color =cols[2], size=3.5, hjust=0)+
    annotate(geom="text", x=.18, y=-0.18, label="Transplant site",size=3.5, hjust=0)+
    geom_text(aes(x=.19, y=-0.22), label = "▲", size = 3.5, family = "HiraKakuPro-W3", color=1)+
    annotate(geom="text", x=.21, y=-0.22, label="Offshore",color =1, size=3.5, hjust=0)+
    annotate("point", x = .19, y =-0.2 , colour = 1, size=3.5)+
    annotate(geom="text", x=.21, y=-.2, label="Nearshore",color =1, size=3.5, hjust=0)
ggsave("Output/RTE_nmds.legendInside.pdf", height = 5, width = 5)
