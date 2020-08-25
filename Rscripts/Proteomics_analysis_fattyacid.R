library(dplyr)
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(plyr)
cols<-c("#ED8636","#0433FF")
#library(Rcmdr)


#Create a NSAF dataframe for GO term

adjnsaf<-read.csv('Data/RTE_adjNSAF.csv',stringsAsFactors = F, row.names = 1)
a<-read.table("RTE_BP/2.beta_oxidation2.txt")
a<-as.vector(a$V1)
a_prot<-filter(adjnsaf, PROTID %in% a)
###
protIDs<-read.csv("Data/ProteinID_2.11.17_sorted.csv", stringsAsFactors = F)
protIDs$Name<-gsub(" OS.+",'', protIDs$Name )
a_prot$Protein<-apply(a_prot["PROTID"], 1, function(x) if(length(protIDs$Name[protIDs$PROTID==x])>0) paste0(protIDs$Name[protIDs$PROTID==x]," (",x,")") else x)
#remove the row m.3482 since they are all zero values
a_prot<-a_prot[-which(a_prot$PROTID=="m.3482"),]

write.csv(a_prot,'Output/betaoxidation.csv')


####### read the saved file #######
prot<-read.csv('Output/betaoxidation.csv', stringsAsFactors = F, row.names = 1)

#rename the long names
prot$Protein<-gsub(' mitochondrial','',prot$Protein)
prot$Protein<-gsub('coenzyme','CoA', prot$Protein)
prot$Protein<-gsub('dehydrogenase','DH', prot$Protein)
prot$Protein<-gsub('specific','sp',prot$Protein)


#
createHeatmapDF<-function(nsaf){
        row.names(nsaf)<-nsaf$Protein
        #get the means for each sample
        nsaf<-nsaf[,-c(1,ncol(nsaf))]
        colnames(nsaf)<-c(rep("NN", times=9),rep("NO", times=9),rep("ON", times=8),rep("OO", times=9))
        
        samples<-c("NN","NO","ON","OO")
        nsaf2<-nsaf[,c(1,10)]
        for (i in 1:4){
                id<-samples[i]
                nums<-which(colnames(nsaf)==id)
                nsaf2[,paste(id)]<-apply(nsaf[nums],1,mean)
        }
        
        #standardize
        dt<-scale(t(nsaf2))
        DF<-t(dt)
        
        return(DF)
}
        
meltDF<-function(df){
        DF2<-data.frame(df)
        DF2$Protein<-rownames(df)
        DFm<-reshape2::melt(DF2)
        colnames(DFm)[2:3]<-c("Treatment","NSAF")
        DFm$Treatment<-factor(DFm$Treatment, levels=c("NO","NN","ON","OO"))
        return(DFm)
        
}        


#create the matrix/data frame for heatmap plotting
df1<-createHeatmapDF(prot) 

dfm<-meltDF(df1)

#add origin
dfm$origin<-c(rep("N coral", times=nrow(dfm)/2), rep("O coral", times=nrow(dfm)/2))

#Order the treatments:
dfm$Treatment<-factor(dfm$Treatment, levels=c("NO","NN","ON","OO"))

#Order the proteins
df2<-data.frame(df1)
df2<-df2[order(df2$ON),]
ordered<-rownames(df2)     

dfm$Protein<-factor(dfm$Protein, levels=ordered)

dfm$Treatment.label<-mapvalues(dfm$Treatment,c("NO","NN","ON","OO"), 
                                   c("N\u2192O","N\u2192N","O\u2192N","O\u2192O") )

ggplot(dfm, aes(x=Treatment.label, y=Protein, fill=NSAF)) + geom_tile()+
        scale_fill_distiller(palette = "RdBu", name="Standardized\nNSAF")+
                theme(axis.text.x = element_text(angle=45, hjust=1))+
        theme(legend.title = element_text(size=10), axis.title=element_blank())+
        facet_grid(~origin, scales="free",space="free")+
        theme(panel.spacing.x=unit(0, "lines"),strip.background.x = element_rect(color="gray50"))
ggsave("Output/betaoxidation_Heatmap_stdNSAF.pdf", width = 6,height = 2.5)



### Combine with detox proteins
detox.m<-read.csv("Output/Detox.melt.file.csv", row.names = 1, stringsAsFactors = F)
deDF<-read.csv("Output/detox_DF_ordered.csv", row.names = 1)
protorder<-rownames(deDF)
detox.m$Protein<-factor(detox.m$Protein, levels=protorder)

prot1<-detox.m
prot1$group<-"Detoxification"

dfm$Protein<-factor(dfm$Protein, levels=ordered)
prot2<-dfm
prot2$group<-"Lipid oxidation"

protC<-rbind(prot1,prot2)
write.csv(protC, "Output/Detox_Lipid_nsafMelt.csv")


protC$Treatment.label<-factor(protC$Treatment.label, levels= c("N\u2192O","N\u2192N","O\u2192N","O\u2192O"))

unique(protC$Treatment.label)
ggplot(protC, aes(x=Treatment.label, y=Protein, fill=NSAF)) + geom_tile()+
        scale_fill_distiller(palette = "RdBu", name="Standardized\nNSAF")+
        theme(axis.text.x = element_text(angle=45, hjust=1))+
        theme(axis.text.y = element_text(size=8))+
        theme(legend.title = element_text(size=10), axis.title=element_blank())+
        facet_grid(group~origin, scales="free",space="free")+
        theme(panel.spacing.x=unit(0, "lines"),strip.background.x = element_rect(color="gray50"))
ggsave("Output/Dotox_Lipid_heatmap.pdf", width=6, height=6)
#save as Height 580, width=480 using Export->eps (this will not show the arrow)
# you can copy image and open with Preveiw -> save as






######################################################
# direction plots for individual proteins

#don't standardize:
createHeatmapDF2<-function(nsaf){
        row.names(nsaf)<-nsaf$Protein
        #get the means for each sample
        nsaf<-nsaf[,-c(1,ncol(nsaf))]
        colnames(nsaf)<-c(rep("NN", times=9),rep("NO", times=9),rep("ON", times=8),rep("OO", times=9))
        
        samples<-c("NN","NO","ON","OO")
        nsaf2<-nsaf[,c(1,10)]
        for (i in 1:4){
                id<-samples[i]
                nums<-which(colnames(nsaf)==id)
                nsaf2[,paste(id)]<-apply(nsaf[nums],1,mean)
        }
        
        return(nsaf2)
}

lipi<-createHeatmapDF2(prot) 
lipi.m<-meltDF(lipi)


lipi.m$Site<-  c(rep("N", times=nrow(lipi.m)/4), rep("O", times=nrow(lipi.m)/4),rep("N", times=nrow(lipi.m)/4), rep("O", times=nrow(lipi.m)/4))
lipi.m$Origin<-c(rep("N", times=nrow(lipi.m)/2), rep("O", times=nrow(lipi.m)/2))

proteins<-prot$Protein

#Create RTE-plots for each proteins
for (i in 1:length(proteins)){
        df<-lipi.m[lipi.m$Protein==proteins[i],]    
        ggplot(df, aes(x=Site,y=NSAF,color=Origin, fill=Origin))+
                geom_point(size=5)+
                scale_color_manual(values=cols)+
                theme_bw()+
                ggtitle(df$Protein[1])+
                theme(plot.title= element_text(size = 7))+
                labs(x="Transplant site",y="NSAF") +
                theme(legend.position="none")+
                geom_segment(data=df, mapping=aes(x=1, y=df[df$Site=="N"&df$Origin=="N","NSAF"], xend=2, yend=df[df$Site=="O"&df$Origin=="N","NSAF"]), size=.8, color=cols[1])+
                geom_segment(data=df, mapping=aes(x=1, y=df[df$Site=="N"&df$Origin=="O","NSAF"], xend=2, yend=df[df$Site=="O"&df$Origin=="O","NSAF"]), size=.8,  color='#0433FF')+
                theme(axis.text.x =element_text(size=12, color="black"), panel.grid.major.x = element_blank(), axis.title= element_text(size=11))
        ggsave(paste0("Output/Lipid_plot/Lipid_nonStd",i,".pdf"), height=2.5, width=2.6)
        
}



### #2.lipid abundance comparison from Qspec
Ncoral<-read.csv("Data/Qspec_Ncorals.csv", stringsAsFactors = F) 
Ocoral<-read.csv("Data/Qspec_Ocorals.csv", stringsAsFactors = F) 
cross<-read.csv("Data/Qspec_NOON.csv", stringsAsFactors = F)
backt<-read.csv("Data/Qspec_back.csv", stringsAsFactors = F)
Nsite<-read.csv("Data/Qspec_Nsite.csv", stringsAsFactors = F)
Osite<-read.csv("Data/Qspec_Osite.csv", stringsAsFactors = F)


lip<-list()
lipNsite<-Nsite[Nsite$Protein %in% prot$PROTID,]
lipOsite<-Osite[Osite$Protein %in% prot$PROTID,]

lipNcoral<-Ncoral[Ncoral$Protein %in% prot$PROTID,]
lipOcoral<-Ocoral[Ocoral$Protein %in% prot$PROTID,]

lipCross<-cross[cross$Protein %in% prot$PROTID,]
lipBack<-backt[backt$Protein %in% prot$PROTID,]

lips<-c("lipNsite",
        "lipOsite",
        "lipNcoral",
        "lipOcoral",
        "lipCross" ,
        "lipBack")
for (i in 1:6){
        lip[[i]]<-get(lips[i])
}

id="m.17984"
id="m.22274"
id="m.27714"
ids<-c("m.17984", id="m.22274", id="m.27714")

for (k in 1:3){
    id<-ids[k]
    
    stats<-data.frame(matrix(ncol=15, nrow=6))
    colnames(stats)[c(1:3,10:15)]<-c("treatment",colnames(lipNsite)[c(1:2,9:14)])
    for (i in 1:6){
        df<-lip[[i]]
        df<-df[df$Protein==id,]
        
        stats$treatment[i]<-paste0(lips[i]," ",colnames(df)[3],"-",colnames(df)[6])
        stats[i,2:16]<-df[1,1:14]
    }
    stats2<-stats[,c(1:2,10:15)]
    write.csv(stats2,paste0("Output/Lipid_abundance_comparison_",ids[k],".csv"))
}
