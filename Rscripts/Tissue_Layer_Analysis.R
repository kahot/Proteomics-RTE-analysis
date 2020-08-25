library(Rcmdr)
library(plotrix)
library(ggplot2)
library(ggthemes)

#Read the data
after<-read.csv("Data/tissue_thickness_A.csv", header=T) 

#Check the assumptions for parametric testing
leveneTest(Thickness ~ Treatment, data=after)
#        Df F value Pr(>F)
#group   3  0.8642 0.4606
#      196   

#"Shapiro-Wilk normality test"
shapiro.test(log(after$Thickness)) 
#W = 0.98698, p-value = 0.06337

##### Comparing the treatments  #####
####  nested 2-way anova   
mod1<-aov(log(Thickness)~Site*Origin/Indiv,data=after)
summary(mod1)
#                   Df Sum Sq Mean Sq F value   Pr(>F)    
#Site                1  0.306  0.3063   33.27 3.44e-08 ***
#Origin              1  0.503  0.5031   54.64 5.25e-12 ***
#Site:Origin         1  0.269  0.2686   29.17 2.08e-07 ***
#Site:Origin:Indiv  16  8.342  0.5214   56.63  < 2e-16 ***
#Residuals         180  1.657  0.0092 

TukeyHSD(x=mod1,'Site:Origin', conf.level=0.95)
#$`Site:Origin`
#                diff         lwr         upr     p adj
#O:N-N:N  0.004979642 -0.04478639  0.05474567 0.9938582
#N:O-N:N -0.173599999 -0.22336603 -0.12383397 0.0000000
#O:O-N:N -0.022043741 -0.07180977  0.02772229 0.6599050
#N:O-O:N -0.178579641 -0.22834567 -0.12881361 0.0000000
#O:O-O:N -0.027023383 -0.07678941  0.02274265 0.4959223
#O:O-N:O  0.151556258  0.10179023  0.20132229 0.0000000

#Rewriting to Genotyp:Site style
#NO-NN  0.9938582
#ON-NN  0.0000000
#NN-OO  0.6599050
#NO-ON  0.0000000
#OO-NO  0.4959223
#OO-ON  0.0000000






#######
### Plot the results ####

cols<-c("#ED8636","#0433FF")

nn<-subset(after, Treatment=="NN")
no<-subset(after, Treatment=="NO")
on<-subset(after, Treatment=="ON")
oo<-subset(after, Treatment=="OO")

Origin<-c("Nearshore","Nearshore","Offshore","Offshore")
stats2<-data.frame(Site=c("Nearshore","Offshore","Nearshore","Offshore"),
                   Origin=c("Nearshore","Nearshore","Offshore","Offshore"),
                   Type=c("NN","ON","NO","OO"), 
                   Tissue=c(mean(nn$Thickness), mean(no$Thickness),mean(on$Thickness),mean(oo$Thickness)), 
                   SE=c(std.error(nn$Thickness), std.error(no$Thickness), std.error(on$Thickness),std.error(oo$Thickness)))


ggplot(stats2, aes(x=Site,y=Tissue,color=factor(Origin),fill=factor(Origin))) +
    geom_errorbar(position=position_dodge(width=0.1),mapping=aes(ymin=Tissue-SE, ymax=Tissue+SE,color=factor(Origin)), width=.18)+
    scale_color_manual(values=paste0(cols))+
    geom_segment(data=stats2, mapping=aes(x=1, y=stats2[1,4], xend=2, yend=stats2[2,4]), size=.8, color=cols[1])+
    geom_segment(data=stats2, mapping=aes(x=2.05, y=stats2[4,4], xend=1.05, yend=stats2[3,4]), size=.8,  color='#0433FF')+
    geom_point(position=position_dodge(width=0.1),size=9, shape=21, color=c("#ED8636","#ED8636","#0433FF","#0433FF"))+
    scale_fill_manual(values=c("#E8AE80","#4E74FF"))+
    theme(plot.margin = unit(c(0,0,0,0), "lines"))+
    labs(x="Transplant site",y="Tissue layer thickness (mm)") +theme_bw()+theme(legend.position="none")+scale_y_continuous(limits = c(2.4, 3.2))+
    theme(axis.text.x =element_text(size=12, face="bold",color=c("#ED8636","#0433FF")), panel.grid.major.x = element_blank(), axis.title = element_text(size=14))+
    annotate(geom="text", x=1.05, y=3.07, label="a",color ='black', size=4)+
    annotate(geom="text", x=1.1,  y=2.55, label="b",color ='black', size=4, parse=TRUE, hjust=0)+
    annotate(geom="text", x=1.18,  y=2.55, label=paste("(","italic(P)", " < 0.001)"),color ='black', size=3.5, parse=TRUE, hjust=0)+
    annotate(geom="text", x=2.08, y=3.1, label="a",color ='black', size=4)+
    annotate(geom="text", x=2.11,  y=2.92, label="a",color ='black', size=4)+
    annotate(geom="text", x=0.975, y=stats2[1,4], label="N",color ='black', size=6)+
    annotate(geom="text", x=1.025, y=stats2[3,4], label="O",color ='white', size=6)+
    annotate(geom="text", x=1.975, y=stats2[2,4], label="N",color ='black', size=6)+
    annotate(geom="text", x=2.025, y=stats2[4,4], label="O",color ='white', size=6)
ggsave(filename="Output/Tissue.pdf",width=4, height=3.8 )



                        
                        

