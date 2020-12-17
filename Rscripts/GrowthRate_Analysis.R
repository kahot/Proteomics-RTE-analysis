#Growth Rate analysis
library(Rcmdr)
library(ggplot2)
library(reshape2)
library(plotrix)

cols<-c("#ED8636","#0433FF")

######
#read the data
GR<-read.csv('Data/Rel_growth_gram.csv')
colnames(GR)[3:14]<-gsub("\\.", " ", colnames(GR)[3:14])

grm<-melt(GR, id_vars=c("Origin","Sample"))

#Test the assumptions
leveneTest((value)^2 ~ Sample, data=grm) #2.2e-16 ***
shapiro.test(grm$value)   #2.2e-16
# can't use parametric tests

GR2<-GR[,c(1:2,14 )]
grm2<-melt(GR2, id_vars=c("Origin","Sample"))
leveneTest(log(value) ~ Sample, data=grm2) #0.006212 **

shapiro.test(log(grm2$value)) p-value = 0.4455

#Run Wilcoxon Rank Sum Test
wilcox.test(GR$`Day 77`[GR$Origin=="N"],GR$`Day 77`[GR$Origin=="O"], "greater")
#W = 797, p-value = 1.532e-08

#######
se<-aggregate(GR[3:14], GR['Origin'], std.error)
me<-aggregate(GR[3:14], GR['Origin'], mean)

SE<-melt(se)
ME<-melt(me)
colnames(ME)[2:3]<-c("Day", "Growth")
ME$SE<-SE$value

ggplot(ME, aes(x=Day, y=Growth, color=Origin))+
    geom_point(position=position_dodge(width=0.2),size =2)+
    scale_color_manual(values=cols,  labels=c("N-coral","O-coral"))+
    geom_errorbar(aes(ymin=Growth-SE, ymax=Growth+SE), width=.2, size=.2,position=position_dodge(width=0.2))+
    geom_line(aes(x=Day, y=Growth, group=Origin),position=position_dodge(width=0.2), stat='identity')+
    theme_bw()+
    labs(x="")+
    ylab("Relative growth (g)")+
    theme(axis.text.x = element_text(size=12, color=1,angle=45, hjust=1),
          axis.text.y = element_text(color=1),axis.title.y = element_text(size=13))+
    theme(panel.grid.major.x = element_blank(),legend.title = element_blank())+
    annotate("text", x=11.3, y=.03,label= expression(~italic(P)~" = 1.5 x"~10^-8), size=3.5)+
    annotate("point", x = 1, y =.46, colour = cols[1], size=3)+
    annotate(geom="text", x=2, y=.46, label="N-coral",size=3)+
    annotate("point", x = 1, y =.43 , color=cols[2], size=3)+
    annotate(geom="text", x=2, y=.43, label="O-coral",size=3)+
    theme(legend.position="none")
ggsave("Output/Growth_plot.pdf", width =7 ,height = 3.7,useDingbats=FALSE)
