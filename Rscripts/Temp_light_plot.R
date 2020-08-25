# Plot the environmental parameter data from the HOBO logger
library(ggplot2)
library(gtable)
library(grid)

env<-read.csv("Data/RTE_templight.csv",stringsAsFactors = F)

n<-env[,c(1,2,4)]
o<-env[,c(1,3,5)]
n$Site<-"N"
o$Site<-"O"
colnames(n)[2:3]<-c("Temp","LUX")
colnames(o)[2:3]<-c("Temp","LUX")
ev<-rbind(n,o)

library(ggplot2)
library(ggthemes)
cols<-c("#ED8636","#0433FF")

dates<-ev$Date
for (i in 2:length(dates)){
    if (i%%30!=0) dates[i]<-'' 
    
}
breaks1<-seq(1, length(dates)/2, by=30)
temp<-ggplot(ev,aes(x=Date, y=Temp, color=Site))+
    theme(legend.position = "none", panel.grid.major.y = element_line(color="gray60", linetype=2,size=0.2))+ 
    geom_line(aes(x=Date, y=Temp, group=Site), linetype = 1, size=0.2)+
    scale_color_manual(values=cols) +
    ylab(expression("Temperature ("*~degree*C*")"))+theme_few()+
    xlab("")+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    theme(panel.background = element_rect(color="white"))+
    theme(legend.position = "none", panel.grid.major.y = element_line(color="gray60", linetype=2,size=0.2))+
    annotate("segment", x = 1000, xend = 1050, y = 28.5, yend = 28.5, colour = cols[1], size=.4)+
    annotate(geom="text", x=1140, y=28.5, label="  Nearshore",size=3)+
    annotate("segment", x = 1000, xend = 1050, y = 28.2, yend = 28.2, color=cols[2], size=.4)+
    annotate(geom="text", x=1140, y=28.2, label="Offshore",size=3)
 temp
light<-ggplot(ev,aes(x=Date, y=LUX, color=Site))+
    theme(legend.position = "none", panel.grid.major.y = element_line(color="gray60", linetype=2,size=0.2))+ 
    geom_line(aes(x=Date, y=LUX, group=Site), linetype = 1, size=0.2)+
    scale_color_manual(values=cols) +
    ylab("Light intensity (Lux)")+theme_few()+
    scale_x_discrete(labels=dates)+xlab("")+
    theme(axis.text.x = element_text(size =6, angle=90, hjust=1, color="gray20"), axis.ticks.x=element_blank())+
    theme(panel.background = element_rect(color="white"))+
    theme(legend.position = "none", panel.grid.major.y = element_line(color="gray60", linetype=2,size=0.2)) 


g1 <- ggplotGrob(temp)
g2 <- ggplotGrob(light)
g <- rbind(g1, g2, size = "first")
g$widths <- unit.pmax(g1$widths, g2$widths)
grid.newpage()
grid.draw(g)

ggsave(filename="Output/Temp_light.pdf",width =8, height =6)
