library(data.table)
library(gplots)
library(plotrix)

relabund <-read.table('chow.relabund.metadata', header=T,sep='\t',strip.white=TRUE)
relabund <- relabund[,c(2:4,6:ncol(relabund))]


#Cage - the cage that is used for the plot
# meta is the relabund plot that includes metadata
#lim is the limmit that an otu must reach at least once during the experiment to be included in the plot

relabund_plot<-function(cage,meta,lim){
  std <- function(x){sd(x)/sqrt(length(x))}
  focus <- meta[meta$Cage==cage,]
  days<- sort(unique(focus$Day))


  for(i in 1:length(days)){
    day = days[i]
    focus_day <- focus[focus$Day==day,]
    x<-apply(focus_day[,3:ncol(focus_day)],2,mean)
    y<-apply(focus_day[,3:ncol(focus_day)],2,std)
    if(i==1){
      means<- rbind(x)
      stderror <- rbind(y)
    }
    else {
      means<-rbind(means,x)  
      stderror<-rbind(stderror,y)
    }
  }
# Remove the OTU that never had an abundace over 0.01
  k<-1
  bad <- c()
  for(i in 2:ncol(means)){
    test<-which(means[,i]>lim)
    if(length(test)==0){
      bad[k]=i
      k=k+1
    }
  }

#Use Data.table to remove 'bad' columns there may be a better way to do this so that we can just used data frames
#or maybe just data tables to be consistant

  means <-data.table(means)
  stderror<-data.table(stderror)

  badnames = quote(names(means)[bad])

  means<-set(means,j=eval(badnames),value=NULL)

  badnames = quote(names(stderror)[bad])   # since means is changed above and the identity of badnames is quoted so it also changes
  stderror <-set(stderror, j=eval(badnames), value=NULL)

  no_days_means <- means[,Day:=NULL]
  no_days_ste <- stderror[,Day:=NULL]

#And now that everything is removed we can jump back in to the world of data frames since it is easier to handel plotting by
#just using the row numbers
  no_days_means<-data.frame(no_days_means)
  no_days_ste<-data.frame(no_days_ste)


#Plotting###
plot(1, type='n', xlim=c(0,days[day]), ylim=c(0,0.5), ylab="Relative Abundance", xlab="Day", xaxt="n", main = paste("Relative abundance in cage:",cage,sep=" "))

  for(i in 1:ncol(no_days_means)){
    points(days,no_days_means[,i],pch=i,col=i,type='l')
    plotCI(days,no_days_means[,i], uiw=no_days_ste[,i], add=T, col=i)
  }

}
