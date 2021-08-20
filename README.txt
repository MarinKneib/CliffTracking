These are the scripts used for the tracking of cliffs & associated ponds and modeling of the cliff population dynamics

The scripts are used in the following order:

CliffPond_info.r: 	Script to derive cliff & pond information from shapefile & raster data
					This script calls the function Extract_cliff_info_noasp_func.r that extracts the necessary information

Cliff_tracking.r: 	To track cliffs from one date to the other based on the cliff information

Associated_ponds.r:	To find ponds associated to cliffs based on the cliff & pond information

CliffModel_script.r:To run the cliff population dynamics model, initializing all the parameters and running
					the model contained in the PopDynMod_MainLoop.r function
					
					
					
					
Author: Marin Kneib
Work address: Swiss Federal Research Institute WSL
Email: marin.kneib@wsl.ch