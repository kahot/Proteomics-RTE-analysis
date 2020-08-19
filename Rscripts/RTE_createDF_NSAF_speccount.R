
library(dplyr)
library(gplots)
library(DataCombine)


#read the abacus out put file
abacus<-read.csv("Data/ABACUS_All_output.csv", header=T)

colnames(abacus)[1:100]

#select Normalized NSAF value columns
nsafadj<-select(abacus,contains("_ADJNSAF"))
nsafadj<-cbind(abacus$PROTID, nsafadj)
colnames(nsafadj)[1]<-"PROTID"
colnames(nsafadj)<-gsub('X2016_AUG_01_KAHO_UWPR_QE','',colnames(nsafadj))
colnames(nsafadj)<-gsub("_ADJNSAF",'',colnames(nsafadj))
write.csv(nsafadj, file='Data/Abacus_nsafadj.csv')

##select total NSAF values
nsaftot<-select(abacus,contains("_TOTNSAF"))
nsaftot<-cbind(abacus$PROTID, nsaftot)
colnames(nsaftot)[1]<-"PROTID"
colnames(nsaftot)<-gsub('X2016_AUG_01_KAHO_UWPR_QE','',colnames(nsaftot))
colnames(nsaftot)<-gsub("_TOTNSAF",'',colnames(nsaftot))
write.csv(nsaftot, file='Data/Abacus_nsaftot.csv')




#select adjusted spectral count values
specadj<-select(abacus,contains("_NUMSPECSADJ"))
specadj<-cbind(abacus$PROTID, specadj)
colnames(specadj)[1]<-"PROTEID"
colnames(specadj)[3:ncol(specadj)]<-gsub('X2016_AUG_01_KAHO_UWPR_QE','',colnames(specadj)[3:ncol(specadj)])
write.csv(specadj, file='Data/Abacus_specadj.csv')


#select Normalized NSAF value columns
spectot<-select(abacus,contains("_NUMSPECSTOT"))
spectot<-cbind(abacus$PROTID, spectot)
colnames(spectot)[1]<-"PROTID"
colnames(spectot)[3:ncol(spectot)]<-gsub('X2016_AUG_01_KAHO_UWPR_QE','',colnames(spectot)[3:ncol(spectot)])
colnames(spectot)<-gsub("_NUMSPECSTOT",'',colnames(spectot))
write.csv(spectot, file='Data/Abacus_spectot.csv')

####
#Check the sample IDs
ids<-read.table("Data/sampleIDs.txt")
spcs<-read.csv("Data/Abacus_sample1-39_2016.csv")

for (i in 1:35){
        specT<-spcs[,c("PROTID", paste0("K",ids$V2[i]))]
        df<-spectot[,c("PROTID", paste0("_",ids$V1[i]))]
        DF<-merge(specT, df,by="PROTID")
        DF$Diff<-DF[,2]-DF[,3]
        ids$Check[i]<-length(DF$Diff[DF$Diff!=0])
        
}


#select RTE samples
ids$aba_id<-paste0("_", ids$V1)
rte_nsafadj<-nsafadj[,c("PROTID",paste(ids$aba_id))]
colnames(rte_nsafadj)<-c("PROTID",paste(ids$V3))
write.csv(rte_nsafadj, "Data/RTE_adjNSAF.csv")

rte_nsaftot<-nsaftot[,c("PROTID",paste(ids$aba_id))]
colnames(rte_nsaftot)<-c("PROTID",paste(ids$V3))
rte_nsaftot[rte_nsaftot$PROTID=="m.3482",]
#Stil zero????
write.csv(rte_nsaftot, "Data/RTE_totNSAF.csv")



#create an average for each sample
rte_nsaf<-rte_nsafadj[,1:2]
idlist<-unique(ids$V3)
for (i in 1:12){
        id<-idlist[i]
        nums<-which(colnames(rte_nsafadj)==id)
        rte_nsaf[,paste(id)]<-apply(rte_nsafadj[nums],1,mean)
}

write.csv(rte_nsaf,"Data/RTE_adjNSAF12.csv")




## Read RTE adjustedNSAF values for 12 samples
adjnsaf<-read.csv('Data/RTE_adjNSAF.csv',stringsAsFactors = F, row.names = 1)



#### 
setwd("~/Dropbox/R/LCMS/RTE_BP")
NSAFtotal<-list()

for (i in 1:14){
  go<-file.choose()
  a<-read.table(go)
  a<-as.vector(a$V1)
  a_prot<-filter(adjnsaf, PROTID %in% a)
  a_prot1<-a_prot[,-1]
  a_sum<-colSums(a_prot1)
  NSAFtotal[[i]]<-a_sum
 }

NSAFsummary2<-do.call(rbind,NSAFtotal)
write.csv(NSAFsummary2,'NSAFsummary2.csv')

