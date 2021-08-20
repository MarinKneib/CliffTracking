# Script to run cliff population dynamics model

#Author: Marin Kneib
#Work address: Swiss Federal Research Institute WSL
#Email: marin.kneib@wsl.ch


###### clear entire workspace (excl. packages)
rm(list = ls())
gc()
graphics.off()

###### Upload libraries
library(raster)
library(gsubfn) # For functions with several outputs
library(CircStats) # To get aspect mean & std
library(circular)
library(readxl) # Read tracking data
library(tibble)
library(rgdal)
library(ggplot2)
library(data.tree)
library(rgeos)
library(sp)
library(rgdal)
library(matrixStats)

###### Define paths
root<-'path/to/root'
results_folder<-'path/to/results/folder'
setwd(root)
source("PopDynMod_MainLoop.R") # Model function

###### Input parameters (example for Langtang - parameters calculated over full time series)

## LANGTANG
data_folder_L<-'N:/gebhyd/8_Him/Personal_folders/Marin/Remote_sensing/Cliff_tracking/RapidEye/Langtang/Tracking/Tracking_SM_2020_09_25'

# Death events (Ndeath = a*tot_pop+b + residuals + lognormal distr parameters + uncertainties)
a_death<-0.41
b_death<--23
std_res_death<-0.35
avg_res_death<--0.07
mu_death<-6.45
sigma_death<-0.69

# Size boundaries
max_size<-11350
min_size<-100

# Merge events (mean & std merge rate + distr: choose randomly from existing cliffs + mean & std area fraction indpdt of initial size)
Mean_nb_M_evts<-6.1
STD_nb_M_evts<-2.8
Mean_mother_ratio_M<-2.03
STD_mother_ratio_M<-0.06
a_mean_afrac_M<--0.13
b_mean_afrac_M<-1.17
a_sigma_afrac_M<--0.28
b_sigma_afrac_M<-2.78

# Split events (mean & std split rate + lognormal distr parameters + area fraction distribution: frac = b*size^a for mean & std)
Mean_nb_S_cliffs<-7.1
STD_nb_S_cliffs<-3.8
Mean_daughter_ratio_S<-2.06
STD_daughter_ratio_S<-0.12
a_Mean_afrac_S<--0.31
b_Mean_afrac_S<-2.67
a_STD_afrac_S<-0.01
b_STD_afrac_S<-0.46
Mean_size_init_S<-8.0
STD_size_init_S<-0.88

# Mixed events
Mean_nb_Mi_evts<-6.2
STD_nb_Mi_evts<-3.0
Mean_daughter_ratio_Mi<-2.46
STD_daughter_ratio_Mi<-0.34
Mean_mother_ratio_Mi<-2.48
STD_mother_ratio_Mi<-0.39
a_mean_afrac_Mi<-0.10
b_mean_afrac_Mi<--0.87
a_sigma_afrac_Mi<--0.03
b_sigma_afrac_Mi<-0.72
Mean_size_init_Mi<-7.4
STD_size_init_Mi<-1.01

# Size evolution (applied to all cliffs that do not die/merge/split)
a_Mean_afrac_R<--0.23
b_Mean_afrac_R<-1.7
a_STD_afrac_R<-0.051
b_STD_afrac_R<-0.19

# Birth events (mean & std birth rate + lognormal distr parameters for size)
Mean_birth_rate<-62.4
STD_birth_rate<-20.4

mu_birth<-6.42
sigma_birth<-0.71

# Glacier characterisitcs
Area_L<-11249486  # in m2

start_year_L<-2009

# Initial cliff population (year 2009) - cliff characteristics
load(paste(data_folder_L,"/Cliff_df_info_2009.RData",sep=""))
Init_pop_L<-Cliff_df_info

# remove what's not needed (just keep information on cliff area)
for (i in 1:length(Init_pop_L)){
  Init_pop_L[[i]]<-Init_pop_L[[i]]$Area
}

###### Parameter list
param_list_L<-c(a_death,b_death,std_res_death,avg_res_death,mu_death,sigma_death,max_size,min_size,Mean_nb_M_evts,
                STD_nb_M_evts,Mean_mother_ratio_M,STD_mother_ratio_M,a_mean_afrac_M,b_mean_afrac_M,a_sigma_afrac_M,b_sigma_afrac_M,
                Mean_nb_S_cliffs,STD_nb_S_cliffs,Mean_daughter_ratio_S,STD_daughter_ratio_S,a_Mean_afrac_S,b_Mean_afrac_S,a_STD_afrac_S,
                b_STD_afrac_S,Mean_size_init_S,STD_size_init_S,Mean_nb_Mi_evts,STD_nb_Mi_evts,Mean_mother_ratio_Mi,
                STD_mother_ratio_Mi,Mean_daughter_ratio_Mi,STD_daughter_ratio_Mi,a_mean_afrac_Mi,b_mean_afrac_Mi,
                a_sigma_afrac_Mi,b_sigma_afrac_Mi,Mean_size_init_Mi,STD_size_init_Mi,a_Mean_afrac_R,b_Mean_afrac_R,a_STD_afrac_R,
                b_STD_afrac_R,Mean_birth_rate,STD_birth_rate,mu_birth,sigma_birth)


## general model paramters (not site-specific)
n_sim<-200
# Length at which to run the model
t_step<-365         # in d
N_yrs<-10

######## Run model

output_list_L<-PopDynMod_MainLoop(param_list_L,n_sim,Init_pop_L,N_yrs)

# Model outputs (for plotting/saving)
tot_sim_nb_L<-output_list_L[1]
tot_sim_sz_L<-output_list_L[2]

