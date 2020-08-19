library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(plyr)
#library(gplots)
cols<-c("#ED8636","#0433FF")

## Read RTE adjustedNSAF values for 12 samples
adjnsaf<-read.csv('Data/RTE_adjNSAF.csv',stringsAsFactors = F, row.names = 1)


##########################
#Look at detoxification proteins (ave. loffold change is larger in NO but abundance is higher in NN)
### look at adjNSAF values for detox proteins
a<-read.table("RTE_BP/17.detox.txt")
a<-as.vector(a$V1)
de_prot<-filter(adjnsaf, PROTID %in% a)
write.csv(de_prot,'Output/detox_nsaf.csv')


### look at adjNSAF values for detox proteins
#Read the file with shorter names for GO terms
deto<-read.csv('Output/detox_nsaf2.csv', stringsAsFactors = F, row.names = 1)
colnames(deto)[4:15]<-c(rep("NN", times=3),rep("NO", times=3),rep("ON", times=3),rep("OO", times=3))
samples<-c("NN","NO","ON","OO")
de.mean<-deto[,2:3]
for (i in 1:4){
        id<-samples[i]
        nums<-which(colnames(deto)==id)
        de.mean[,paste(id)]<-apply(deto[nums],1,mean)
}

detox<-melt(de.mean, id.vars=c("Names","Protein"))
detox$Origin<-''
detox$Site<-''

for (i in 1:nrow(detox)){
        if (detox$variable[i]=="NN") {detox$Origin[i]<-"N"; detox$Site[i]<-"N"}
        else if (detox$variable[i]=="NO") {detox$Origin[i]<-"N"; detox$Site[i]<-"O"}
        else if (detox$variable[i]=="ON") {detox$Origin[i]<-"O"; detox$Site[i]<-"N"}
        else   { detox$Origin[i]<-"O"; detox$Site[i]<-"O"}
}
colnames(detox)[4]<-"NSAF"
write.csv(detox,"Output/detox_NSAF_summary.csv")
proteins<-deto$Names

#Create RTE-plots for each proteins
for (i in 1:length(proteins)){
        df<-detox[detox$Name==proteins[i],]    
        ggplot(df, aes(x=Site,y=NSAF,color=Origin, fill=Origin))+
                geom_point(size=5)+
                scale_color_manual(values=cols)+
                theme_bw()+
                ggtitle(df$Name[1])+
                theme(plot.title= element_text(size = 11))+
                labs(x="",y="NSAF") +
                theme(legend.position="none")+
                geom_segment(data=df, mapping=aes(x=1, y=df[df$Site=="N"&df$Origin=="N","NSAF"], xend=2, yend=df[df$Site=="O"&df$Origin=="N","NSAF"]), size=.8, color=cols[1])+
                geom_segment(data=df, mapping=aes(x=1, y=df[df$Site=="N"&df$Origin=="O","NSAF"], xend=2, yend=df[df$Site=="O"&df$Origin=="O","NSAF"]), size=.8,  color='#0433FF')+
                theme(axis.text.x =element_text(size=12, color="black"), panel.grid.major.x = element_blank(), axis.title= element_text(size=11))
        ggsave(paste0("Output/Detox_plot/Detox_nonStd_",i,".pdf"), height=2.5, width=2.5)
        
}


#######   Create HEATMAPs ##########
#convert to matrix for heatmap
#standardize the de.mean matrix

deDF<-de.mean[,c(2:6)]
rownames(deDF)<-deDF$Names
deDF<-deDF[,-1]
dt<-scale(t(deDF))
deDF<-t(dt)
write.csv(deDF,"Output/detox_DF.csv")


detox.m<-melt(deDF)
colnames(detox.m)<-c("Protein","Treatment","NSAF")

#Order the treatments:
detox.m$Treatment<-factor(detox.m$Treatment, levels=c("NO","NN","ON","OO"))

#order the proteins -> skip
#deDF<-data.frame(deDF)
#deDF2<-deDF[order(deDF$ON, decreasing=F),]

#de_1<-deDF2[16:27,]
#de_2<-deDF2[1:15,]
#de_2<-de_2[order(de_2$OO),]

#de_comb<-rbind(de_2, de_1)
#protorder<-rownames(de_comb)


### Read the protein order 
deDF<-read.csv("Output/detox_DF_ordered.csv", row.names = 1)
protorder<-rownames(deDF)

detox.m$Protein<-factor(detox.m$Protein, levels=protorder)
detox.m$origin<-c(rep("N coral", times=54), rep("O coral", times=54))


detox.m$Treatment.label<-mapvalues(detox.m$Treatment,c("NO","NN","ON","OO"), 
                                   c("N\u2192O","N\u2192N","O\u2192N","O\u2192O") )

#write.csv(detox.m,"Output/RTE/detox.m.csv")

ggplot(detox.m, aes(x=Treatment.label, y=Protein, fill=NSAF)) + geom_tile()+
    scale_fill_distiller(palette = "RdBu", name="Standardized\nNSAF")+
    theme(axis.text.x = element_text(angle=45, hjust=1))+
    theme(legend.title = element_text(size=10), axis.title=element_blank())+
    geom_vline(xintercept=2.5, color="white")+
    facet_grid(~origin, scales="free_x",space="free_x")+
    theme(panel.spacing.x=unit(0, "lines"), strip.background.x = element_rect(color="gray50"))
ggsave("Output/Heatmap_detox.pdf", width = 5,height = 4.5)


#save deto.melt files
write.csv(detox.m,"Output/Detox.melt.file.csv")

