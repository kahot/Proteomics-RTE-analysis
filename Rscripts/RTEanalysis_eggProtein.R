#Look into the proteins related to reproduction (oocyte development)

### Egg protein
protIDs<-read.csv("Data/ProteinID_2.11.17_sorted.csv", stringsAsFactors = F)
protIDs$Name<-gsub(" OS.+",'', protIDs$Name )

eggProteins<-protIDs$PROTID[grep("egg", protIDs$Name,ignore.case = T)]
write.csv(eggProteins,"Output/Eggproteins.csv")


#Qsepc data
Ncoral<-read.csv("Data/Qspec_Ncorals.csv", stringsAsFactors = F) 
Ocoral<-read.csv("Data/Qspec_Ocorals.csv", stringsAsFactors = F) 
Nsite<-read.csv("Data/Qspec_Nsite.csv", stringsAsFactors = F)
Osite<-read.csv("Data/Qspec_Osite.csv", stringsAsFactors = F)

#Look for proteins with egg or vitellin in names

findProteins<-function(protlist){
    pList<-list()
    pList[[1]]<- Nsite[Nsite$Protein %in% protlist,]
    pList[[2]]<- Osite[Osite$Protein %in% protlist,]
    pList[[3]]<- Ncoral[Ncoral$Protein %in% protlist,]
    pList[[4]]<- Ocoral[Ocoral$Protein %in% protlist,]
    
    return(pList)
}

eggs<-protIDs$PROTID[grep("egg", protIDs$Name,ignore.case = T)]
vit<-protIDs$PROTID[grep("vitellogenin", protIDs$Name,ignore.case = T)]
yolk<-protIDs$PROTID[grep("yolk", protIDs$Name,ignore.case = T)]

oo<-c(eggs,vit, yolk)
oogenesis<-findProteins(oo)
treats<-c("Nsite","Osite","N-coral","O-coral")
for (i in 1:4){
    write.csv(oogenesis[[i]],paste0("Output/Oogenesis_proteins_",treats[i],".csv"))
}

