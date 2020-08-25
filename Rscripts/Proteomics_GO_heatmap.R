#Create a heatmap of enriched GO terms unique to different treatments

library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(plyr)
library(dplyr)
library(stringr)

#Read the file with GO terms with P-values and # of protens
go<-read.csv('Data/EnrichedGO_All_numbered.csv')

# Use log(p-values) for the heatmap
pval<-go[,c(2,4:8)]
colnames(pval)[2:6]<-c("treatment","NN","NO","ON","OO")
pval$GO.Name<-trimws(pval$GO.Name, "l")

#use melt function to create the data frame for a heatmap
pval.m<-melt(pval, id.vars=c("GO.Name", "treatment"))
colnames(pval.m)[3:4]<-c("Coral","Protein")

pval.m$Coral<-factor(pval.m$Coral, levels=c("NO","NN","ON","OO"))

#label for Facet_grid
labs<-list(
        "ON"= expression(O%->%N),
        'NN'= expression(N%->%N),
        'NO'= expression(N%->%O),
        "OO"= expression(O%->%O),
        "N-coral"= "N coral",
        "O-coral"= "O coral",
        "N-Site" = "N site",
        "O-Site"="O site",
        "CrossT" ="CrossT",
        "N"="N coral",
        "O"="O coral"
       )


plot_labeller <- function(variable,value){
        return(labs[value])
}


#Specify the order of different treatments
a<-c("ON","NN","NO","OO","N-coral","O-coral","N-Site","O-Site" ,"CrossT") 
pval.m$treatment<-factor(pval.m$treatment, levels=a)

#include arrows for labels
pval.m$Coral.label<-mapvalues(pval.m$Coral,c("NO","NN","ON","OO"), 
                              c("N\u2192O","N\u2192N","O\u2192N","O\u2192O") )

#Add origin and transplant site informaiton
pval.m$origin<-c(rep("N",  times=nrow(pval.m)/2), rep("O", times=nrow(pval.m)/2))
pval.m$site<-c(rep("Nsite", time=nrow(pval.m)/4), rep("Osite", time=nrow(pval.m)/4),rep("Nsite", time=nrow(pval.m)/4),rep("Osite", time=nrow(pval.m)/4))

#Format the y-axis (reverse the orders in each section) for the heatmap
max_width = max(nchar(as.character(pval$GO.Name)))
pval.m$label_padded = sprintf(paste0("%-", max_width, "s"),pval.m$GO.Name)
pval$padded<-sprintf(paste0("%-", max_width, "s"),pval$GO.Name)

goOrders<-c(paste(rev(pval$padded[1:4])),  paste(rev(pval$padded[5:19])),
            paste(rev(pval$padded[20:25])),paste(rev(pval$padded[26:51])),
            paste(rev(pval$padded[52:54])),paste(rev(pval$padded[55:60])),
            paste(rev(pval$padded[61:63])),paste(rev(pval$padded[64:65])),
            paste(rev(pval$padded[66:73])))

pval.m$label_padded<-factor(pval.m$label_padded, levels=goOrders)

ggplot(pval.m, aes(x=Coral.label, y=label_padded, fill=log10(Protein))) + geom_tile(color="gray90", size=0.05)+
    scale_fill_gradient2(low = "#E33244",high="white", midpoint=-1,name="P-value",breaks=c(-4,-3,-2,-1,0),labels=c(expression(10^-4),expression(10^-3),expression(10^-2),expression(10^-1),1))+
    theme_bw()+ theme(axis.text.x = element_text(angle=45, hjust=1))+
    theme(legend.title = element_text(size=10), axis.title=element_blank())+
    facet_grid(treatment~origin, scales="free",space="free", labeller=plot_labeller)+
    theme(strip.text.y = element_text(angle = 0), strip.background.x = element_rect(color="gray20"))+
    theme(axis.text.y=element_text(size=7, hjust=0))+
    theme(panel.spacing.x=unit(0, "lines"))+
    theme(legend.text.align = 0, legend.text = element_text(size=8))

ggsave("Output/LCMS/GOAll.pdf",width = 5.5,height =7.2)



