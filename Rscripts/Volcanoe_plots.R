#Create volcano plots of proteomic results
library(ggplot2)
library(calibrate)

cols<-c("#ED8636","#0433FF")

#1. Nearshore-site
n<-read.csv("Data/Qspec_Nsite.csv", stringsAsFactors = F) #Qspec result
n$DEFLINE<-gsub(" OS.+",'', n$DEFLINE)

#Find significant proteins 
nSig<-n[abs(n$LogFoldChange)>=0.5 & abs(n$Zstatistic)>=2,]
nSig<-nSig[!grepl("symbB1", nSig$Protein),]

#Add the two semi significant proteins m.6724 and m.17984
add<-n[n$Protein=="m.17984"|n$Protein=="m.6724",]
nSig<-rbind(nSig,add)

#Protein IDs for detox and lipid oxidataion (created from CompGO previously)
fat<-read.table("RTE_BP/2.beta_oxidation2.txt") 
fa<-as.vector(fat$V1)
deto<-read.table("RTE_BP/17.detox.txt")
de<-as.vector(deto$V1)

# Add other significant GO categories
# CC file from CompGO result
cc<-read.table("Data/NNNOupNO_CC.txt", fill = T, sep="\t", stringsAsFactors = F)
#fidn the proteins associated with the following GO terms
ccv<-unlist(strsplit(cc$V6[cc$V2=="clathrin coat"],","))
vc<-unlist(strsplit(cc$V6[cc$V2=="vesicle coat"],","))


#BP
bp<-read.table("Data/NNNOupNO_BP.txt", fill = T, sep="\t", stringsAsFactors = F)

#tRNA
trna<-unlist(strsplit(bp$V6[bp$V2=="tRNA aminoacylation"],","))

#MF (all overlapping with BP or CC. Do not need to run)

#create a volcano plot 
ggplot()+
        geom_point(data=n, aes(x=LogFoldChange, y=abs(Zstatistic)), color="gray", size=.6, shape=21)+
        theme_bw()+xlab("Log2 fold change")+ ylab('Significance (|Zstat|)')+
        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
        geom_hline(yintercept=2, color="gray70", linetype=2)+
        geom_point(data=nSig[nSig$Zstatistic>=2 & nSig$LogFoldChange >=0.5,], aes(x=LogFoldChange, y=Zstatistic), color=cols[2],shape=21, size=0.7)+
        geom_point(data=nSig[nSig$Zstatistic<=(-2) & nSig$LogFoldChange <=(-0.5),], aes(x=LogFoldChange, y=abs(Zstatistic)), color=cols[1],shape=21, size=0.7)+
        geom_point(data=nSig[nSig$Protein%in%trna,], aes(x=LogFoldChange, y=Zstatistic), color=cols[2],bg="red",shape=21, size=1.3)+
        geom_point(data=nSig[nSig$Protein%in%vc,], aes(x=LogFoldChange, y=Zstatistic), color=cols[2],bg="lawngreen",shape=21, size=1.3)+
        geom_point(data=nSig[nSig$Protein%in%ccv,], aes(x=LogFoldChange, y=Zstatistic), color=cols[2],bg="lawngreen",shape=21, size=1.3)+
        geom_point(data=nSig[nSig$Protein%in%fa,], aes(x=LogFoldChange, y=abs(Zstatistic)), color=cols[2],bg="lightpink",shape=21, size=1.3)+
        geom_point(data=nSig[nSig$Protein%in%de,], aes(x=LogFoldChange, y=abs(Zstatistic)), color=cols[2],bg="yellow",shape=21, size=1.3)+
        #geom_point(data=nSig[nSig$Protein=="m.6024",], aes(x=LogFoldChange, y=abs(Zstatistic)), shape=21, col=cols[1], bg="dodgerblue", cex=1.3)+
        annotate("point", x = 2.5, y = 21, color=cols[2],bg="red",shape=21, size=1.5)+
        annotate(geom="text", x=2.7, y=21, label="tRNA aminoacylation ",size=2, hjust=0)+
        annotate("point", x = 2.5, y = 20, color=cols[2],bg="lawngreen",shape=21, size=1.5)+
        annotate(geom="text", x=2.7, y=20, label="vesicle/clathrin coat",size=2, hjust=0)+
        annotate("point", x = 2.5, y = 19, color=cols[2],bg="lightpink",shape=21, size=1.5)+
        annotate(geom="text", x=2.7, y=19, label="lipid oxidation",size=2, hjust=0)+
        #annotate("point", x = 2.5, y = 18, color=cols[2],bg="paleturquoise1",shape=21, size=1.5)+
        #annotate(geom="text", x=2.7, y=18, label="lipid oxidation",size=2, hjust=0)+
        annotate("point", x = 2.5, y = 18, color=cols[2],bg="yellow",shape=21, size=1.5)+
        annotate(geom="text", x=2.7, y=18, label="detoxification",size=2, hjust=0)+
        #annotate("point", x = 2.5, y = 16, color=cols[2],bg="dodgerblue",shape=21, size=1.5)+
        #annotate(geom="text", x=2.7, y=16, label="Acid ceramidase",size=2, hjust=0)+
        annotate(geom="text", x=-3.5, y= 21, label="Significantly more abundant in: ", size=2, hjust=0)+
        annotate("point", x = -3.4, y= 20, color=cols[2],shape=21, size=1.5)+
        annotate(geom="text", x=-3.2, y= 20, label="Offshore coral ",size=2, hjust=0)+
        annotate("point", x = -3.4, y= 19, color=cols[1],shape=21, size=1.5)+
        annotate(geom="text", x=-3.2, y= 19, label="Nearshore coral",size=2, hjust=0)
        #geom_text(data=nSig[nSig$Zstatistic<=(-2) & nSig$LogFoldChange <=(-2),], aes(x=LogFoldChange, y=abs(Zstatistic),label=DEFLINE),hjust=0, vjust=0, size=2)
ggsave("Output/RTE/Nsite_volcano_plot.pdf", width =4 ,height = 3)        
        

######
#2. Offshore site

o<-read.csv("Data/Qspec_Osite.csv", stringsAsFactors = F)
o$DEFLINE<-gsub(" OS.+",'', o$DEFLINE)

#
oSig<-o[abs(o$LogFoldChange)>=0.5 & abs(o$Zstatistic)>=2,]
oSig<-oSig[!grepl("symbB1", oSig$Protein),]

#find significant proteins belong to portein translation:

t<-read.table("RTE_BP/tRNA.txt")
a1<-as.vector(t$V1)

oSig$tRNA<-apply(oSig["Protein"],1, function(x) if(x%in%a1) x=1 else x=0)

#find beta oxidataion protines
fat<-read.table("RTE_BP/2.beta_oxidation2.txt") #use the file that include m.6724
a2<-as.vector(fat$V1)

fatO<-oSig[oSig$Protein %in% a2,]
#write.csv(fatO, "Output/LipidOxidation_Osite.csv")

ggplot()+
    geom_point(data=o, aes(x=LogFoldChange, y=abs(Zstatistic)), color="gray", size=.6, shape=21)+
    theme_bw()+xlab("Log2 fold change")+ ylab('Significance (|Zstat|)')+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    geom_hline(yintercept=2, color="gray70", linetype=2)+
    geom_point(data=oSig[oSig$Zstatistic>=2 & oSig$LogFoldChange >=0.5,], aes(x=LogFoldChange, y=Zstatistic), color=cols[2],shape=21, size=0.7)+
    geom_point(data=oSig[oSig$Zstatistic<=(-2) & oSig$LogFoldChange <=(-0.5),], aes(x=LogFoldChange, y=abs(Zstatistic)), color=cols[1],shape=21, size=0.7)+
    geom_point(data=oSig[oSig$tRNA==1,], aes(x=LogFoldChange, y=Zstatistic), color=cols[2],bg="olivedrab1",shape=21, size=1.3)+
    annotate("point", x = 1.7, y = 46.4, color=cols[2],bg="olivedrab1",shape=21, size=1.5)+
    annotate(geom="text", x=1.9, y=47, label="tRNA aminoacylation \nfor protein translation ",size=2, hjust=0, vjust=1)+
    annotate(geom="text", x=-4.8, y= 47, label="Significantly more abundant in: ", size=2, hjust=0)+
    annotate("point", x = -4.7, y= 45, color=cols[2],shape=21, size=1.5)+
    annotate(geom="text", x=-4.5, y= 45, label="Offshore coral ",size=2, hjust=0)+
    annotate("point", x = -4.7, y= 43, color=cols[1],shape=21, size=1.5)+
    annotate(geom="text", x=-4.5, y= 43, label="Nearshore coral",size=2, hjust=0)

#ggsave("Output/Osite_volcano_plot.pdf", width =4 ,height = 3)        


#remove the outliers
ggplot()+
    geom_point(data=o, aes(x=LogFoldChange, y=abs(Zstatistic)), color="gray", size=.6, shape=21)+
    theme_bw()+xlab("Log2 fold change")+ ylab('Significance (|Zstat|)')+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    geom_hline(yintercept=2, color="gray70", linetype=2)+
    geom_point(data=oSig[oSig$Zstatistic>=2 & oSig$LogFoldChange >=0.5,], aes(x=LogFoldChange, y=Zstatistic), color=cols[2],shape=21, size=0.7)+
    geom_point(data=oSig[oSig$Zstatistic<=(-2) & oSig$LogFoldChange <=(-0.5),], aes(x=LogFoldChange, y=abs(Zstatistic)), color=cols[1],shape=21, size=0.7)+
    geom_point(data=oSig[oSig$tRNA==1,], aes(x=LogFoldChange, y=Zstatistic), color=cols[2],bg="olivedrab1",shape=21, size=1.3)+
    annotate("point", x = 1.7,    y= 21.8, color=cols[2],bg="olivedrab1",shape=21, size=1.5)+
    annotate(geom="text", x=1.9,  y= 22, label="tRNA aminoacylation \nfor protein translation ",size=2, hjust=0, vjust=1)+
    annotate(geom="text", x=-4.8, y= 22, label="Significantly more abundant in: ", size=2, hjust=0)+
    annotate("point", x = -4.7,   y= 21, color=cols[2],shape=21, size=1.5)+
    annotate(geom="text", x=-4.5, y= 21, label="Offshore coral ",size=2, hjust=0)+
    annotate("point", x = -4.7,   y= 20, color=cols[1],shape=21, size=1.5)+
    annotate(geom="text", x=-4.5, y= 20, label="Nearshore coral",size=2, hjust=0)+ylim(0,22)

ggsave("Output/Osite_volcano_plot.pdf", width =4 ,height = 3)        

