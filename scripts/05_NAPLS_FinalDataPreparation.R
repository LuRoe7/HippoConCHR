# FINAL DATA PREPARATIONS ---------------------------------------------
# Description:  This script makes some final adjustments to the data prior to statistics
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
# define input directory for behavioral and functional connectivity data
in_dir = "results/napls/04_NAPLS_BehavData_PrepareMedication/"
# define output directory
out_dir = "results/napls/05_NAPLS_FinalDataPreparation/"
# create output directory for this script
dir.create(out_dir, recursive = TRUE,showWarnings=F)
# read behavioral data
df_napls_wide_raw = read.csv(paste0(in_dir,"napls_MergedData_IncludingMed.csv"))
# define numeric longitudinal variables
vars_long_numeric = colnames(df_napls_wide_raw)[grepl("__",colnames(df_napls_wide_raw))]
# bring to long format
df_napls_verylong = df_napls_wide_raw %>% pivot_longer(cols = all_of(vars_long_numeric),
                                                   names_to = c("outcome","session"),
                                                   names_sep = "__",
                                                   values_to = "value")
# define baseline variables
vars_base = c("subject_id","subjectkey","session","site","subjecttype","interview_age","sex","educat","primary_language","ethnic_group","wasi_IQ")
# bring to shorter long format
df_napls_long = df_napls_verylong %>% pivot_wider(id_cols = all_of(vars_base),names_from = outcome)


# RENAME COLUMNS ---------------------------------------------------------------------------
# shorten column names
colnames(df_napls_long) = colnames(df_napls_long) %>% str_replace("subject_id","subj") %>% str_replace("subjecttype","group") %>%
  str_replace("interview_age","age") %>% str_replace("educat","ey") %>% str_replace("primary_language","language") %>%
  str_replace("ethnic_group","ethnicity") %>% str_replace("wasi_IQ","IQ") %>% str_replace("fc_","fc") %>% str_replace("hip","hi") %>%
  str_replace("amyg","amy") %>% str_replace("sops_negative","NegS") %>% str_replace("sops_positive","PosS") %>%
  str_replace("sops_general","GenS") %>% str_replace("sops_disorganization","DisS") %>%
  str_replace("sops_p1_UnusualThought_DelusionalIdeas","P1") %>% str_replace("sops_p2_Suspciousness_PersecutoryIdeas","P2") %>%
  str_replace("sops_p3_GrandioseIdeas","P3") %>% str_replace("sops_p4_PerceptualAbnorm_Hallucinations","P4") %>%
  str_replace("sops_p5_DisorganizedCommunication","P5") %>% str_replace("sops_n1_SocialAnhedonia","N1") %>%
  str_replace("sops_n2_Avolition","N2") %>% str_replace("sops_n3_DecreasedEmotionExpression","N3") %>%
  str_replace("sops_n4_DecreasedExperienceOfEmotionsSelf","N4") %>% str_replace("sops_n5_DecreasedIdeationalRichness","N5") %>%
  str_replace("sops_n6_OccupationalFunctioning","N6") %>% str_replace("sops_d1_OddBehaviorAppearance","D1") %>%
  str_replace("sops_d2_BizarreThinking","D2") %>% str_replace("sops_d3_TroubleFocusAttention","D3") %>%
  str_replace("sops_d4_ImpairedHygiene","D4") %>% str_replace("sops_g1_SleepDisturbance","G1") %>%
  str_replace("sops_g2_DysphoricMood","G2") %>% str_replace("sops_g3_MotorDisturbance","G3") %>%
  str_replace("sops_g4_ImpairedStressTolerance","G4") %>% str_replace("calg_ts","DepS") %>% str_replace("gafrat","gaf") %>%
  str_replace("bacs_sc_raw","SyCod") %>% str_replace("hvltr_raw","VeMem")
# drop a few columns that are not required
df_napls_long = df_napls_long %>% select(-c(subjectkey,gafyear))
# convert to wide format again
df_napls_wide = df_napls_long %>% pivot_wider(id_cols = all_of(c("subj","site","group","sex","age","ey","language","ethnicity","IQ")),
                                               names_from = session,
                                               values_from = -c(subj,site,session,group,sex,age,ey,language,ethnicity,IQ),
                                               names_sep = "_")


# DATA EXPORT ---------------------------------------------------------------------------
# export data frame
write.csv(df_napls_wide,paste0(out_dir,"napls_MergedData_Unstandardized.csv"),row.names = F)

