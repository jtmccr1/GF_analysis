library(gplots)
library(plotrix)

trends<-function(Phylolevel,Inocs,lim){
file<- paste('~/Desktop/GF_data/colonus.relabund.tx.',Phylolevel,'.metadata', sep='')
meta <- read.table(file,header=T,sep='\t',strip.white=TRUE)
std <- function(x){sd(x)/sqrt(length(x))}
str2num<-function(x){as.numeric(as.character(x))}

#Inocs<-unique(meta$Inoc) # for comparing all the innocula
#Inocs<-c('H1','H2','H3','C1','C2','C3')

#=========Trim metadata to just metadata and relative abundance data==============================================================================
meta$X=NULL  #An extra colum of NA was read in removes this
meta<-meta[,c(2:6,8:ncol(meta))]  #just metadata and OTU

#========Trim everything down to just the inocula of interest=====================================================================================
meta<-meta[meta$Inoc %in% Inocs,]

#=========Remove OTUs that are never above lim in at least 2 samples===============================================================================
q<-1  
good <- c()    # to be used in the following if statements
for(h in 6:ncol(meta)){  # start at first OTU column
  test<-which(meta[,h]>=lim)
  if(length(test)>=2){
    good[q]=h  # Save the id of the good OTU
    q=q+1
  }
}
meta <-meta[,c(1:5,good)] # Trim meta down to the metadata and those OTUs that past the limit test

#==========Get info for each OTU for each inoc=================================================================================================== 

for(c in 6:ncol(meta)){  # highlight 1 OTU at a time there are 5 columns of meta data


  #+++++++++++++++set up Plotting+++++++++++++++++++++++++++++++++++++++++++
  plot(1, type='n', xlim=c(-1,30), ylim=c(0,1), ylab="Relative Abundance", xlab="Day", main = paste("Plot of",names(meta[c]),"at Phylolevel",Phylolevel))
  key<-c() # initiate key for legends
  llave<-1 # initiate counter for legend
  key_colors<-c()
  focus_otu<-meta[,c(1:5,c)]


  #=============Work with one Inocula ==================
    for(i in 1:length(Inocs)){  # Highlight one inocula at a time
     
      focus_inoc<-focus_otu[focus_otu$Inoc==Inocs[i],]
      Cages<-unique(focus_inoc$Cage)  # get the cages for this inocula
      
      #==============Take a cage from that Inocula========================    
      for(j in 1:length(Cages)){ # take each cage separately 
      
        focus_cage<-focus_inoc[focus_inoc$Cage==Cages[j],]
        if(unique(focus_cage$Cage)=='INOC'){ # There should only be one cage so unique should be length 1
          focus_cage$Day=0
        }
        Days<- unique(focus_cage$Day) # Get the days for this cage
      
        Days<- sort(str2num(Days)) # convert to number vector by function above and then sort
      
        
        mediocrist<-c() # Set up vector to collect mean values across days
        error<- c() # Set up vector to collect error values across days

        #==============Average over mice/day in cage and add to vector=========================================        
        
        for(k in 1:length(Days)){
          focus_day<-focus_cage[focus_cage$Day==Days[k],] # Get the mice we have in this cage for this day
          mediocrist[k]<-mean(focus_day[,6])  # should take the last column which is one housing the otu data
          error[k]<- std(focus_day[,6])
          
                  
      }#Close day

      #===============Add Data from cage to Plot=================================================
      
      color<-hsv(i/10+(j*2)/100) # i gives the tenth spot and j gives either + 0.00 or 0.05 since the most we have of cages is 2 
      
      points(Days,mediocrist,pch=16,col=color,type='b')  # color goes with Inocula(i) and cage(j)
      plotCI(Days,mediocrist, uiw=error, add=T, col=color)
      
      key[llave]<-paste(Inocs[i],Cages[j],sep='_')  #llave is spanish for key and serves as the counter here
      key_colors[llave]<-color
      llave<-llave+1   
    }#close Cage
    
  } #close Inoc
  
legend(25,1,key,pch=16,col=key_colors) # takes color from distriubution of j which corresponds to cages

} #Close OTU
}
  