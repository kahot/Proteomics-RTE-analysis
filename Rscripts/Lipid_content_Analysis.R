library(Rcmdr)
library(plotrix)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(grid)

cols<-c("#ED8636","#0433FF")
#Read the data
lipids<-read.csv("Data/lipids.csv", header=T, stringsAsFactors = F) 

#Check the assumptions for parametric testing
leveneTest(Lipids ~ Treatment, data=lipids)
#      Df F value Pr(>F)
#group  1  0.1585 0.6928
#      38    
shapiro.test(lipids$Lipids)
#W = 0.94888, p-value = 0.06942

#after RTE analysis
after<-lipids[lipids$Treatment=="After",]
after$Treat2<-paste0(after$Origin,after$Site)
### ANOVA
a1 <- aov(Lipids~Treat2,data=after)
summary(a1)
#            Df  Sum Sq  Mean Sq F value Pr(>F)  
#Treatment    3 0.02502 0.008339   3.151 0.0539 .
#Residuals   16 0.04234 0.002646                 

TukeyHSD(x=a1,conf.level=0.95)
#              diff         lwr         upr     p adj
#NO-NN  0.016651917 -0.07643267 0.109736503 0.9550729
#ON-NN  0.081475239 -0.01160935 0.174559824 0.0973067 
#OO-NN -0.008889362 -0.10197395 0.084195224 0.9925862
#ON-NO  0.064823321 -0.02826126 0.157907906 0.2315865
#OO-NO -0.025541279 -0.11862586 0.067543306 0.8601419
#OO-ON -0.090364600 -0.18344919 0.002719985 0.0586391  .



######### Before and After Comparison #################

### N->N  comparing before and after###
lipN1<-aov(Lipids~ Treatment, data=lipids[lipids$Origin=="N"&lipids$Site=="N",] )
summary(lipN1)
#            Df  Sum Sq  Mean Sq F value Pr(>F)
#Treatment    1 0.00469 0.004691   0.568  0.473
#Residuals    8 0.06606 0.008257   

### N->O before and after ###
lipN2<-aov(Lipids~ Treatment, data=lipids[lipids$Origin=="N"&lipids$Site=="O",] )
summary(lipN2)
#            Df  Sum Sq  Mean Sq F value Pr(>F)
#Treatment    1 0.00899 0.008991   0.996  0.348
#Residuals    8 0.07222 0.009028


### O->N  before and after ###
lipO1<-aov(Lipids~ Treatment, data=lipids[lipids$Origin=="O"&lipids$Site=="N",] )
summary(lipO1)
#            Df  Sum Sq Mean Sq F value  Pr(>F)   
#Treatment    1 0.06676 0.06676   33.67 0.000404 ***
#Residuals    8 0.01586 0.00198

### O->O  before and after ###
lipO2<-aov(Lipids~ Treatment, data=lipids[lipids$Origin=="O"&lipids$Site=="O",] )
summary(lipO2)
#            Df   Sum Sq Mean Sq F value   Pr(>F)    
#Treatment    1 0.013339 0.013339   64.72 0.0000419 ***
#Residuals    8 0.001649 0.000206


# Extracting before exp values
before<-lipids[c(1:5,21:25),]

lip.b<-aov(Lipids~ Origin, data=before )
summary(lip.b)
#             Df  Sum Sq Mean Sq F value Pr(>F)
#Origin       1 0.00373 0.003728   0.526  0.489
#Residuals    8 0.05672 0.007091 


#######################################################################
#Calculate variance of var(x/y) 
#var(x/y)= var(x)/y^2 + var(y)x^2/y^4
lipid_raw<-read.csv("Data/Lipids_raw_scaled.csv")
post<-lipid_raw[lipid_raw$Time=="Post",]

cal_lipid_variance<-function(x){ #x = a data.frame
    sum<- aggregate(x["total"], by=list(x$Treatment), mean)
    lip<- aggregate(x["lipid"], by=list(x$Treatment), mean)
    varT<-aggregate(x["total"], by=list(x$Treatment), var)  
    varL<-aggregate(x["lipid"], by=list(x$Treatment), var)  
    
    sum$Var.t<-varT$total
    sum$lipid<-lip$lipid
    sum$Var.l<-varL$lipid
    
    #x=lipid y=tissue
    sum$Var<-sum$Var.l/(sum$total^2)+sum$Var.t* sum$lipid^2 / (sum$total^4)
    sum$SEM<-sqrt(sum$Var/5)
    sum$Ratio<-sum$lipid/sum$total
    return(sum)
}


sum<-cal_lipid_variance(post)


##### Plot the reslts

sum$Site<-c("N","O","N","O")
sum$Origin<-c("N","N","O","O")
sum$Ratio<-sum$Ratio*100
sum$SEM<-sum$SEM*100
ggplot(sum, aes(x=Site,y=Ratio,color=factor(Origin),fill=factor(Origin)))+
    geom_errorbar(mapping=aes(ymin=Ratio-SEM, ymax=Ratio+SEM,color=factor(Origin)), position=position_dodge(width=0.1),width=.18)+
    geom_segment(data=sum, mapping=aes(x=1, y=sum[1,8], xend=2, yend=sum[2,8]), size=.7, color=cols[1])+
    geom_segment(data=sum, mapping=aes(x=2.05, y=sum[4,8], xend=1.05, yend=sum[3,8]), size=.7,  color=cols[2])+
    scale_color_manual(values=cols)+
    scale_fill_manual(values=c("#E8AE80","#4E74FF"))+
    theme(axis.text.x =element_text(size=12, color=1), panel.grid.major.x = element_blank(), axis.title = element_text(size=13))+
    theme(plot.margin = unit(c(0,0,0,0), "lines"))+
    geom_point(position=position_dodge(width=0.1), aes(fill=factor(Origin)),size=9, shape=21, color=c("#ED8636","#ED8636","#0433FF","#0433FF"))+
    labs(x="Transplant site",y="Tissue lipid content (w/w %)") +theme_bw()+theme(legend.position="none")+
    theme(axis.text.x =element_text(size=12, face="bold",color=c("#ED8636","#0433FF")), panel.grid.major.x = element_blank(), axis.title = element_text(size=14))+
    scale_x_discrete(labels=c("Nearshore", "Offshore"))+
    annotate(geom="text", x=1.1, y=29.5, label=paste0("a^{.}~(italic(P) == 0.054)"),color ='black', size=4,parse=TRUE,hjust=0)+
    annotate(geom="text", x=1.05, y=20.7, label="a",color ='black', size=4)+
    annotate(geom="text", x=2.05, y=22.8, label="a",color ='black', size=4)+
    annotate(geom="text", x=2.12, y=19.7, label="a",color ='black', size=4)+
    annotate(geom="text", x=0.975, y=sum[1,8], label="N",color ='black', size=6)+
    annotate(geom="text", x=1.025, y=sum[3,8], label="O",color ='white', size=6)+
    annotate(geom="text", x=1.975, y=sum[2,8], label="N",color ='black', size=6)+
    annotate(geom="text", x=2.025, y=sum[4,8], label="O",color ='white', size=6)
ggsave(filename="Output/Lipid.pdf",width=4, height=3.8 )


######################
### Before and after plot ###

#lipids<-read.csv("Data/lipids.csv", header=T, stringsAsFactors = F) 
before<-lipid_raw[lipid_raw$Time=="Initial",]
pre_sum<-cal_lipid_variance(before)

pre_sum$Site<-c("N","O")
pre_sum$Origin<-c("N","O")
pre_sum[,7:8]<-pre_sum[,7:8]*100
sum<-rbind(sum,pre_sum)
sum$Time<-c(rep("Post-exp", time=4), rep("Initial", times=2))

write.csv(sum, "Data/Lipids_summary.csv")


text_before <- textGrob("Initial", gp=gpar(fontsize=11))
text_after <- textGrob("Post-exp", gp=gpar(fontsize=11))


N_plot<-
    ggplot(sum[sum$Origin=="N",], aes(x=Time,y=Ratio,color=Site))+
    scale_x_discrete(labels=c("Initial", "Post-exp"))+
    geom_errorbar(mapping=aes(ymin=Ratio-SEM, ymax=Ratio+SEM,color=factor(Site)), position=position_dodge(width=0.2),width=.18)+
    geom_segment(data=sum, mapping=aes(x=1,y=sum[5,8], xend=1.97, yend=sum[2,8]), size=.5,  color=cols[2],arrow=arrow(length = unit(0.3, "cm")))+
    geom_segment(data=sum, mapping=aes(x=1,y=sum[5,8], xend=1.88, yend=sum[1,8]), size=.6, color=cols[1],arrow=arrow(length = unit(0.3, "cm")))+
    labs(x="",y="Tissue lipid content (w/w %)") +
    theme_bw()+theme(legend.position="none")+
    scale_color_manual(values=c(cols, cols[1]))+
    scale_fill_manual(values=c("#E8AE80"))+
    ylim(6.4,33.9)+
    theme(plot.margin = unit(c(1,0.5,1,1), "lines"))+
    geom_point(aes(x=1.95, y=sum[1,8]),col=cols[1],bg=cols[1], pch=21, size=5)+
    geom_point(aes(x=2.05   , y=sum[2,8]),col="#0433FF",bg=cols[1], pch=21, size=5)+
    geom_point(aes(x=1,    y=sum[5,8]),col=cols[1], bg=cols[1], pch=21, size=5)+
    theme(axis.ticks.x = element_blank(), axis.text.x=element_text(color=1, size=12),
          panel.grid.major.x = element_blank())+
    annotate("point", x = 2.1, y =6.7 , colour = cols[1], size=3)+
    annotate(geom="text", x=2.4, y=6.7, label=expression("N" %->% "N"),size=3)+
    annotate("point", x = 2.1, y =8.5 , col="#0433FF",bg=cols[1], pch=21, size=3)+
    annotate(geom="text", x=2.4, y=8.5, label=expression("N" %->% "O"),size=3)
N_plot   

##       
O_plot<-
    ggplot(sum[sum$Origin=="O",], aes(x=Time,y=Ratio,color=Site))+
    geom_errorbar(mapping=aes(ymin=Ratio-SEM, ymax=Ratio+SEM,color=factor(Site)), position=position_dodge(width=0.1),width=.18)+
    geom_segment(data=sum, mapping=aes(x=1,y=sum[6,8], xend=1.97, yend=sum[4,8]), size=.5,  color=cols[2],arrow=arrow(length = unit(0.3, "cm")))+
    geom_segment(data=sum, mapping=aes(x=1,y=sum[6,8], xend=1.93, yend=sum[3,8]), size=.6, color=cols[1],arrow=arrow(length = unit(0.3, "cm")))+
    labs(x="",y="") +
    theme_bw()+theme(legend.position="none")+
    scale_color_manual(values=c(cols, cols[2]))+
    scale_fill_manual(values=c("#4E74FF"))+
    #geom_point(position=position_dodge(width=0.1), aes(fill=factor(Origin)),size=5, shape=21, color=c("#ED8636",cols[2],"#ED8636"))+
    ylim(6.4,33.9)+
    theme(plot.margin = unit(c(1,1,1,.5), "lines"))+
    geom_point(aes(x=2.02, y=sum[4,8]),col=cols[2],bg=cols[2], pch=21, size=5)+
    geom_point(aes(x=1.98   , y=sum[3,8]),col=cols[1],bg=cols[2], pch=21, size=5)+
    geom_point(aes(x=1,    y=sum[6,8]),col=cols[2], bg=cols[2], pch=21, size=5)+
    theme(axis.text.x = element_text(color=1, size=12),axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank())+
    annotation_custom(text_before, xmin=1, xmax=1, ymin=-5, ymax=-5)+
    annotation_custom(text_after, xmin=2, xmax=2, ymin=-5, ymax=-5)+
    annotation_custom(text_before, xmin=1, xmax=1, ymin=-5, ymax=-5)+
    annotation_custom(text_after, xmin=2, xmax=2, ymin=-5, ymax=-5)+
    annotate("point", x = 2.1, y =6.7 , colour = "#0433FF", size=3)+
    annotate(geom="text", x=2.4, y=6.7, label=expression("O" %->% "O"),size=3)+
    annotate("point", x = 2.1, y =8.5 , bg=cols[2],col=cols[1], pch=21, size=3.5)+
    annotate(geom="text", x=2.4, y=8.5, label=expression("O" %->% "N"),size=3)

O_plot 

pdf("Output/LLipids_before.after.pdf",width=7, height=3.7)
grid.arrange(N_plot,O_plot,nrow=1)
dev.off()

