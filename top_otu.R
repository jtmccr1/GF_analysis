library(gplots)
library(plotrix)

top_otu<-function(Phylolevel,Inoc,lim,Cages){
  file<- paste('~/Desktop/GF_data/colonus.relabund.tx.',Phylolevel,'.metadata', sep='')
  ref_file<-paste('~/Desktop/GF_data/colonus.reference.',Phylolevel,'.taxonomy',sep='')
  meta <- read.table(file,header=T,sep='\t',strip.white=TRUE)
  taxa<- read.table(ref_file,header=T,sep='\t',strip.white=TRUE)
  std <- function(x){sd(x)/sqrt(length(x))}
  str2num<-function(x){as.numeric(as.character(x))}
  get_name<-function(x){taxa$Taxonomy[taxa$OTU==x]}
  
  #Inocs<-unique(meta$Inoc) # for comparing all the innocula
  
  #=========Trim metadata to just metadata and relative abundance data==============================================================================
  meta$X=NULL  #An extra colum of NA was read in removes this
  meta<-meta[,c(2:6,8:ncol(meta))]  #just metadata and OTU
  
  #========Trim everything down to just the inocula of interest=====================================================================================
  focus_Inoc<-meta[meta$Inoc %in% Inoc,]
  
  #=========Take the top lim OTU from the Inocs=====================================================================================================
  ave_pop<-sort(apply(focus_Inoc[,6:ncol(focus_Inoc)],2,mean),decreasing=TRUE)  # average the relabunances of each OTU and sort in decreasing order
  passed_cutoff<-ave_pop[1:lim]   # Identify those OTU which pass the cutoff of being in the top lim OTU
  good<- which(names(focus_Inoc) %in% names(passed_cutoff))  # take the good colums and save
  
  
  focus_Inoc <-focus_Inoc[,c(1:5,good)] # Trim meta down to the metadata and those OTUs that past the limit test 
  
  if(missing(Cages)){
    Cages<-unique(focus_Inoc$Cage)  # get the cages for this inocula and otu
  }
  
  #==================Set up Plotting================================================================================================================
  #par(mar=c(5.1, 4.1, 4.1, 15.1), xpd=TRUE) # Set margins for legend outside of plot
  plot(1, type='n', xlim=c(-5,30), ylim=c(0,1.3), ylab="Relative Abundance", xlab="Day", main = paste("Plot of",Inoc,"Cage :",Cages,"Phylolevel",Phylolevel))
  key<-c() # initiate key for legends
  llave<-1 # initiate counter for legend
  key_colors<-c()
  colors<-c('powderblue','purple','red','rosybrown','royalblue','salmon','seagreen','seashell','siena','skyblue','slateblue','slategray')
    
  #==========Get info for each OTU ================================================================================================================
  
  for(c in 6:ncol(focus_Inoc)){  # highlight 1 OTU at a time there are 5 columns of meta data
      
      focus_otu<-focus_Inoc[,c(1:5,c)] # +1 because there is a column   /[1]Scientist/[2]Inoc/[3]Cage/[4]Mouse/[5]Day/[6]OTU
      
     #==============Take a cage from that Inocula========================    
      for(j in 1:length(Cages)){ # take each cage separately 
          
        focus_cage<-focus_otu[focus_otu$Cage==Cages[j],]
        if(unique(focus_cage$Cage)=='INOC'){ # There should only be one cage so unique should be length 1
          focus_cage$Day=-5
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
        
        otu_index<-which(names(meta) %in% names(focus_day[6])) #Bases color off of index in meta so colors should be maintained across Inocs
        color<-colors[otu_index-5] # c starts at 6 
        
        points(Days,mediocrist,pch=16,col=color,type='b')  
        plotCI(Days,mediocrist, uiw=error, add=T, col=color)
        
        key[llave]<-paste(get_name(names(focus_day[6])))  #llave is spanish for key and serves as the counter here
        key_colors[llave]<-color
        if(j==length(Cages)){ # move to another legend spot if all the cages have been plotted
        llave<-llave+1   
        }
      }#close Cage
      
    
    
    
  } #Close OTU
  legend("topright",key,pch=16,col=key_colors,xpd=TRUE) # takes color from distriubution of j which corresponds to cages
} #Close Function
