# Plot the environmental parameter data from the HOBO logger
library(ggplot2)
library(gtable)
library(grid)

env<-read.csv("Data/CommonGardenExp_Temp.csv",stringsAsFactors = F)
class(env$Date)
env$Date<-as.POSIXct(env$Date,format="%m/%d/%y %H:%M", tz="US/Hawaii")

env$day<-substr(env$Date, 1, 10)
dates<-unique(env$day)

temp<-data.frame(Date=dates)
for (i in 1:length(dates)){
    df<-env[env$day==dates[i],]
    temp$Max[i]<-max(df$Temp)
    temp$Min[i]<-min(df$Temp)
}

temp$DailyRange<-temp$Max-temp$Min
max(temp$DailyRange)
temp$Max30<-NA
temp$Min30<-NA
temp$MonthlyRange<-NA
for (i in 1:(nrow(temp)-30)){
    temp$Max30[i]<-max(temp$Max[i:(i+29)])
    temp$Min30[i]<-min(temp$Min[i:(i+29)])
    temp$MonthlyRange[i]<-temp$Max30[i]-temp$Min30[i]
    
}
   
    
