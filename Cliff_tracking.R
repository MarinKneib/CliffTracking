# This script tracks cliffs given general & detailed information on the cliff at 2 specific dates

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

max_ass<-10 # Max number of associated daughter cliffs 
max_distCenter<-2000 # Max distance (to constrain the search)
max_backwasting<-20/365 # Max distancing rate (m.d-1)
max_aspect_std<-45 # Max aspect std for which the cliff is considered in its entirety
max_aspect_mean_diff<-30/365 # Max rotation rate (deg.d-1)
DATE<-c("2009-11-03","2010-10-26") # Dates of datasets

##### Load cliff info 
Date1<-as.Date(DATE[n],format='%Y-%m-%d')
Date2<-as.Date(DATE[n+1],format='%Y-%m-%d')
t_step<-as.double(Date2)-as.double(Date1)

# initial date
load(paste(results_folder,'/Cliff_df_2009.RData',sep=''))  # Load cliff pixel information for the 1st date
load(paste(results_folder,'/Cliff_df_info_2009.RData',sep=''))  # Load general cliff information for the 1st date
Cliff_px1<-Cliff_df_2009
Cliff_gl1<-Cliff_df_info_2009

# final date
load(paste(results_folder,'/Cliff_df_2010.RData',sep=''))  # Load cliff pixel information for the 2nd date
load(paste(results_folder,'/Cliff_df_info_2010.RData',sep=''))  # Load general cliff information for the 2nd date
Cliff_px2<-Cliff_df_2010
Cliff_gl2<-Cliff_df_info_2010

# number of cliffs
n_cliff1<-length(Cliff_px1)
n_cliff2<-length(Cliff_px2)

####### Load velocity maps (in m.yr-1)

vx<-raster(paste('path/to/velocity/data/vxx.tif',sep=''))
vy<-raster(paste('path/to/velocity/data/vy.tif',sep=''))

###### Displacement maps in between the 2 dates (in m)
dx<-vx*t_step/365
dy<-vy*t_step/365

####### Find daughter cliffs - initialize output matrix
Daughter_cliffs<-matrix(NA, nrow = n_cliff1, ncol = max_ass+1)

####### Algorithm description
# For all cliffs of date 1
# Find cliffs of date 2 which are not too far (within a few kms) - this is to constrain the search
# Of all these cliffs, only select the ones which are at a distance less than the max distanciation rate (*time interval) allowed 
# Take surface displacement into account at that point by adding surface motion of mother cliffs
# Add condition on the aspect:
#   If they have low aspect spread and if they have similar aspect mean (+/- 180 deg!!!), they are related
#   If one has a large spread (and the other one low), 
#     in the cliff with a large spread, take the pixels with similar aspect values to the cliff with a low spread
#     if some pixels with similar aspect values are close to the cliff with a low spread (+/- 180 deg!!!), they are related
#   If both have a large spread
#   In each feature take the pixels with aspect close to 0 degrees, then 90 degrees, then 180 degrees, then 270 degrees
#     if some pixels of the 2 features with similar aspect are close to eachother (+/- 180 deg!!!), the cliffs are related
#   Due to lack of DEM, it is impossible to know if the cliff has backwasted or not


for (i in 1:n_cliff1){
  Daughter_cliffs[i,1]<-as.numeric(paste(Cliff_gl1[[i]]$Rank))
  n_Daughter_cliff<-0
  # Get displacement of the mother cliff
  Centr_mother<-c(as.numeric(paste(Cliff_gl1[[i]]$CenterX)),as.numeric(paste(Cliff_gl1[[i]]$CenterY)))
  cell<-cellFromXY(dx,as.matrix(data.frame(t(Centr_mother))))
  mother_dx<-extract(dx,cell)
  mother_dy<-extract(dy,cell)
  
  for (j in 1:n_cliff2){
    # Only consider cliffs whose center is less than 2km away
    dist_centr<-((as.numeric(paste(Cliff_gl1[[i]]$CenterX))+mother_dx-as.numeric(paste(Cliff_gl2[[j]]$CenterX)))^2
                 +(as.numeric(paste(Cliff_gl1[[i]]$CenterY))+mother_dy-as.numeric(paste(Cliff_gl2[[j]]$CenterY)))^2)^(1/2)
    if (dist_centr < max_distCenter){
      # Minimum distance between the 2 features
      minDist<-lapply(1:length(Cliff_px1[[i]]$x), function(x){
        +     spDistsN1(pts = cbind(Cliff_px2[[j]]$x,Cliff_px2[[j]]$y), pt = cbind(Cliff_px1[[i]]$x+mother_dx,Cliff_px1[[i]]$y+mother_dy)[x, ])})
      minDist<-lapply(1:length(Cliff_px1[[i]]$x), function(x){min(minDist[[x]])})
      minDist<-min(as.numeric(minDist))
      # Consider only cliffs within a buffer corresponding to the max distancing rate
      if (minDist < max_backwasting*t_step){
        # Add conditions on aspect
        aspSTD1<-Cliff_gl1[[i]]$stdAspect # Aspect STD in degrees
        aspSTD2<-Cliff_gl2[[j]]$stdAspect
        aspM1<-Cliff_gl1[[i]]$meanAspect # Aspect mean in degrees
        aspM2<-Cliff_gl2[[j]]$meanAspect
        aspMean<-circ.mean(c(aspM2*pi/180,aspM1*pi/180))*180/pi
        aspM_diff<-180 - abs(180 - abs(aspM1 - aspM2) %% 360)
        n_valid<-0
        n_angle_good<-0
        n_angle_bad<-0
        if (aspSTD1 < max_aspect_std & aspSTD2 < max_aspect_std & (aspM_diff < max_aspect_mean_diff*t_step | aspM_diff > (180-max_aspect_mean_diff*t_step)) & n_Daughter_cliff < max_ass){
          n_Daughter_cliff<-n_Daughter_cliff+1
          Daughter_cliffs[i,n_Daughter_cliff+1]<-as.numeric(paste(Cliff_gl2[[j]]$Rank))
        }
        if (aspSTD1 < max_aspect_std & aspSTD2 > max_aspect_std & n_Daughter_cliff < max_ass){
          # Look at pixels in feature with similar aspect than other feature's main aspect
          aspMean<-aspM1
          for (l in 1:length(Cliff_px2[[j]]$Aspect)){
            diff_asp_aspMean<-180 - abs(180 - abs(aspMean - Cliff_px2[[j]]$Aspect[l]) %% 360)
            # Calculate distance if pixel aspect similar to the one of the first feature
            if (diff_asp_aspMean < max_aspect_mean_diff*t_step | (diff_asp_aspMean > (180-max_aspect_mean_diff*t_step))){
              minDist<-spDistsN1(pts = cbind(Cliff_px1[[i]]$x,Cliff_px1[[i]]$y), pt = cbind(Cliff_px2[[j]]$x[l],Cliff_px2[[j]]$y[l]))
              minDist<-min(minDist)
              # If at least one pixel with a not too high aspect rotation is close enough to previous cliff, then the cliffs are related
              if (minDist < max_backwasting*t_step){
                n_valid<-n_valid+1
              }
            }
          }
        }
        if (aspSTD1 > max_aspect_std & aspSTD2 < max_aspect_std & n_Daughter_cliff < max_ass){
          # Look at pixels in feature with similar aspect than other feature's main aspect
          aspMean<-aspM2
          for (k in 1:length(Cliff_px1[[i]]$Aspect)){
            diff_asp_aspMean<-180 - abs(180 - abs(aspMean - Cliff_px1[[i]]$Aspect[k]) %% 360)
            # Calculate distance if pixel aspect similar to the one of the first feature (+/- 180 deg)
            if (diff_asp_aspMean < max_aspect_mean_diff*t_step | (diff_asp_aspMean > (180-max_aspect_mean_diff*t_step))){
              minDist<-spDistsN1(pts = cbind(Cliff_px2[[j]]$x,Cliff_px2[[j]]$y), pt = cbind(Cliff_px1[[i]]$x[k],Cliff_px1[[i]]$y[k]))
              minDist<-min(minDist)
              # If at least one pixel with a not too high aspect rotation is close enough to previous cliff, then the cliffs are related
              if (minDist < max_backwasting*t_step){
                n_valid<-n_valid+1
              }
            }
          }
        }
        if (aspSTD1 > max_aspect_std & aspSTD2 > max_aspect_std & n_Daughter_cliff < max_ass){
          # Look at pixels in feature with similar aspect than other feature's aspect
          for (aspMean in c(0,90,180,270)){
            for (k in 1:length(Cliff_px1[[i]]$Aspect)){
              diff_asp_aspMean1<-180 - abs(180 - abs(aspMean - Cliff_px1[[i]]$Aspect[k]) %% 360)
              if (diff_asp_aspMean1 < max_aspect_std){
                for (l in 1:length(Cliff_px2[[j]]$Aspect)){
                  diff_asp_aspMean2<-180 - abs(180 - abs(aspMean - Cliff_px2[[j]]$Aspect[l]) %% 360)
                  # Calculate distance if pixel aspect similar to the one of the first feature
                  if (diff_asp_aspMean2 < max_aspect_std){
                    minDist<-spDistsN1(pts = cbind(Cliff_px2[[j]]$x[l],Cliff_px2[[j]]$y[l]), pt = cbind(Cliff_px1[[i]]$x[k],Cliff_px1[[i]]$y[k]))
                    if (minDist < max_backwasting*t_step){
                      n_valid<-n_valid+1
                    }
                  }
                }
              }
            }
          }
        }
        if (n_valid > 0){
          n_Daughter_cliff<-n_Daughter_cliff+1
          Daughter_cliffs[i,n_Daughter_cliff+1]<-as.numeric(paste(Cliff_gl2[[j]]$Rank))
        }
      }
    }
  }
}


###### Output results
write.csv(Daughter_cliffs, file = paste(results_folder,"/Daughter_cliffs.csv",sep=""))
