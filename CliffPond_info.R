# This script extracts necessary information for tracking of cliffs & ponds from shapefiles

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
path_func<-'path/2/function/called/in/script/Extract_cliff_info_noasp_func.r'
source(path_func)
root<-'path/2/root'
results_folder<-'path/2/results/folder'
dir.create(results_folder)
shapefiles<-'path/2/shapefiles/folder'
images_folder<-'path/2/coregistered/image/used/for/mapping/of/cliffs/and/ponds'

###### Upload cliff & pond shapefiles + image
img<-stack(paste(images_folder,'/Image.tif',sep=''))
cliff_shp<-shapefile(paste(shapefiles,'/Cliff.shp',sep=''))
pond_shp<-shapefile(paste(shapefiles,'/Pond.shp',sep=''))

###### Use function ExtractCliffInfo_noasp to derive data frames containing all necessary information to track cliffs & ponds
list[Cliff_df,Cliff_df_info,Pond_df,Pond_df_info]<-ExtractCliffInfo_noasp(img,cliff_shp,pond_shp)

###### Save files
save(Cliff_df, file = paste(results_folder,'/Cliff_df_2009.RData',sep=''))
save(Cliff_df_info, file = paste(results_folder,'/Cliff_df_info_2009.RData',sep=''))
save(Pond_df, file = paste(results_folder,'/Pond_df_2009.RData',sep=''))
save(Pond_df_info, file = paste(results_folder,'/Pond_df_info_2009.RData',sep=''))


