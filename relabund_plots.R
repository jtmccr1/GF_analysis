library(gplots)
library(plotrix)

#relabund <-read.table('simple.relabund.metadata', header=T,sep='\t',strip.white=TRUE)
#relabund <- relabund[,c(2:4,6:ncol(relabund))]


#Cage - the cage that is used for the plot
# meta is the relabund plot that includes metadata
#lim is the limmit that an otu must reach at least once during the experiment to be included in the plot

relabund_plot<-function(cage,meta,type,lim){
  #type='count' Debugging
  #cage='A'
  #meta=relabund
  #lim=2
  std <- function(x){sd(x)/sqrt(length(x))}
  focus <- meta[meta$Cage==cage,]
  days<- sort(unique(focus$Day))


  for(i in 1:length(days)){
    day = days[i]
    focus_day <- focus[focus$Day==day,]
    x<-apply(focus_day[,4:ncol(focus_day)],2,mean)
    y<-apply(focus_day[,4:ncol(focus_day)],2,std)
    if(i==1){
      means<- rbind(x)
      stderror <- rbind(y)
    } else {
      means<-rbind(means,x)  
      stderror<-rbind(stderror,y)
    }
  }
  means<-data.frame(means)  # Yeah the row names are duplicated but it is of no
  stderror<-data.frame(stderror) #concern for us. Right?
  good <- c()    # to be used in the following if statements
  if(type=='percent'){
# Remove the OTU that never had an abundace over lim
    k<-1
   
    for(i in 2:ncol(means)){
      test<-which(means[,i]>=lim)
      if(length(test)!=0){
        good[k]=i
        k=k+1
      }
    }
  
  }
    if(type=='count'){
      totalmeans<-sort(apply(means,2,mean),decreasing=TRUE)
      filtered_means<-totalmeans[1:lim]
      good<- which(names(means) %in% names(filtered_means))
    }

means <-means[,c(good)]

stderror <-stderror[,c(good)]
#why not take the good columns and not the bad ones then I don't need to use data tables 
#Use Data.table to remove 'bad' columns there may be a better way to do this so that we can just used data frames
#or maybe just data tables to be consistant






#Plotting###
plot(1, type='n', xlim=c(days[1],days[length(days)]), ylim=c(0,0.8), ylab="Relative Abundance", xlab="Day", main = paste("Relative abundance in cage:",cage,sep=" "))
  
if(is.null(ncol(means))==TRUE){
    size = 1
  } else {
    size = ncol(means)
  }
  for(i in 1:size){
    points(days,means[,i],pch=i,col=i,type='l')
    plotCI(days,means[,i], uiw=stderror[,i], add=T, col=i)
  }
#return(means)
#return(stderror)
#return(days)
list(means,stderror,days)
}
