# This script finds ponds associated to cliffs

#Author: Marin Kneib
#Work address: Swiss Federal Research Institute WSL
#Email: marin.kneib@wsl.ch


###### clear entire workspace (excl. packages)
rm(list = ls())
gc()

###### Upload libraries
library(raster)
library(gsubfn) # For functions with several outputs
library(CircStats) # To get aspect mean & std
library(circular)
library(readxl) # Read tracking data
library(tibble)
library(rgdal)
library(ggplot2)

###### Define paths
root<-'path/to/root'
results_folder<-'path/to/results/folder'

####### Parameters
CliffPond_dist_thr<-10  # maximum distance between cliff & associated pond (in m)
max_ass<-10 # Max number of associated daughter cliffs / ponds
max_distCenter<-2000

##### Load cliff & pond info
load(paste(results_folder,'/Cliff_df_2009.RData',sep=''))  # Load cliff pixel information
load(paste(results_folder,'/Pond_df_2009.RData',sep=''))  # Load pond pixel information
load(paste(results_folder,'/Cliff_df_info_2009.RData',sep=''))  # Load cliff general information
load(paste(results_folder,'/Pond_df_info_2009.RData',sep=''))  # Load cliff general information
Cliff_px1<-Cliff_df
Cliff_gl1<-Cliff_df_info
Pond_px1<-Pond_df
Pond_gl1<-Pond_df_info

####### Find associated ponds
Ass_pond<-matrix(NA, nrow = n_cliff1, ncol = max_ass+1)

for (i in 1:n_cliff1){
  Ass_pond[i,1]<-as.numeric(paste(Cliff_gl1[[i]]$Rank))
  for (j in 1:n_pond1){
    # Only consider ponds less than 2km away
    dist_centr<-((as.numeric(paste(Cliff_gl1[[i]]$CenterX))-as.numeric(paste(Pond_gl1[[j]]$CenterX)))^2
                 +(as.numeric(paste(Cliff_gl1[[i]]$CenterY))-as.numeric(paste(Pond_gl1[[j]]$CenterY)))^2)^(1/2)
    if (dist_centr < max_distCenter){
      n_ass_pond<-0
      # Minimum distance between the 2 features
      minDist<-lapply(1:length(Cliff_px1[[i]]$x), function(x){
        +     spDistsN1(pts = cbind(Pond_px1[[j]]$x,Pond_px1[[j]]$y), pt = cbind(Cliff_px1[[i]]$x,Cliff_px1[[i]]$y)[x, ])})
      minDist<-lapply(1:length(Cliff_px1[[i]]$x), function(x){min(minDist[[x]])})
      minDist<-min(as.numeric(minDist))
      if (minDist < CliffPond_dist_thr & n_ass_pond < max_ass){
        n_ass_pond<-n_ass_pond+1
        Ass_pond[i,n_ass_pond+1]<-as.numeric(paste(Pond_gl1[[j]]$Rank))
      }
    }
  }
}

###### Output results
write.csv(Ass_pond, file = paste(results_folder,"/Associated_ponds.csv",sep=""))
  