# DATA MERGING BETWEEN MODALITIES ---------------------------------------------
# Description:  This script merges the NAPLS3 data across data modalities.
#               
# Author:       Roell, Lukas      
# Created:      2025/08/28
# License:      Creative Common: CC-BY

##### Load packages #####
# load necessary packages
library(pacman)
p_load(tidyr,dplyr,stringr,readxl,ggplot2,gridExtra,data.table)


# BASIC PREPARATIONS ---------------------------------------------------------------------------
# show working directory and ensure that the folder results including the data is located there
getwd()
# define input directory
in_dir = "results/napls/"
# define output directory
out_dir = "results/napls/03_NAPLS_ModalityMerging/"
# create output directory for this script
dir.create(out_dir, recursive = TRUE,showWarnings=F)
# read behavioral data
df_behav_wide = read.csv(paste0(in_dir,"01_NAPLS_BehavData_PreparationCleaning/napls_behavioral_wide_cleaned.csv"))
# read fMRI data
df_fc_long = read.csv(paste0(in_dir,"02_NAPLS_rsfMRIData_PreparationCleaning/napls_rsfMRI_FunctionalConnectivity_cleaned.csv"))
# convert FC data to wide format
df_fc_wide = df_fc_long %>% pivot_wider(id_cols = subject_id,names_from = session,values_from = -c(subject_id, session),names_sep = "__")


# MERGING DATA FRAMES --------------------------------------------------------------------------
# merge data frames by subject ID while keeping all subjects with behavioral data
df_napls_wide = merge(df_behav_wide,df_fc_wide,by="subject_id",all.x = T)


# DATA EXPORT --------------------------------------------------------------------------
# export data frame
write.csv(df_napls_wide,paste0(out_dir,"napls_MergedData_VarsOfInterest.csv"),row.names = F)

