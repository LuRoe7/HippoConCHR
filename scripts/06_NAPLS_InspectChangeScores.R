# INSPECTION OF CHANGE SCORES ---------------------------------------------
# Description:  This script creates plots to inspect the fluctuation of different variables over time.
#               
# Author:       Roell, Lukas      
# Created:      2025/08/28
# License:      Creative Common: CC-BY

##### Load packages #####
# load necessary packages
library(pacman)
p_load(readr,tidyr,dplyr,stringr,ggplot2,ggrepel,gridExtra,readxl,lme4,broom)


# BASIC PREPARATIONS ---------------------------------------------------------------------------
# show working directory and ensure that the folder results including the data is located there
getwd()
# define input directory for behavioral and functional connectivity data
in_dir = "results/napls/05_NAPLS_FinalDataPreparation/"
# define output directory
out_dir = "results/napls/06_NAPLS_InspectChangeScores/"
# create output directory for this script
dir.create(out_dir, recursive = TRUE,showWarnings=F)
# read data
df_napls_wide = read.csv(paste0(in_dir,"napls_MergedData_Unstandardized.csv"))
# make copy
df_inspect = df_napls_wide
# define numeric variables that change over time
vars_numeric_wide = colnames(df_napls_wide)[grepl("_",colnames(df_napls_wide))]
# convert data frame to long format
df_inspect_long = df_inspect %>% pivot_longer(cols = all_of(vars_numeric_wide),
                                              names_to = c("outcome","session"),
                                              names_sep = "_",
                                              values_to = "value")
# insert numeric time variable
df_inspect_long = df_inspect_long %>% mutate(month = case_when(session == "BL" ~ 0,
                                                               session == "M2" ~ 2,
                                                               session == "M4" ~ 4,
                                                               session == "M6" ~ 6,
                                                               session == "M8" ~ 8,
                                                               TRUE ~ NA_real_))
# select only enhanced subjects and healthy controls for plotting
df_inspect_long = subset(df_inspect_long,group=="E"|group=="C")
# rename group column entries
df_inspect_long$group = df_inspect_long$group %>% str_replace("C","HC")  %>% str_replace("E","CHR-P")
# mapping for nice titles
label_map = c(
  "PosS" = "Att. Positive Symptoms", "NegS" = "Negative Symptoms","DepS" = "Depressive Symptoms",
  "gaf" = "Levels of Functioning", "SyCod" = "Symbol Coding", "VeMem" = "Verbal Memory",
  "fchil" = "FC in left Hippocampus", "fchir" = "FC in right Hippocampus","fchih" = "FC between Hippocampi",
  "fchi" = "FC in Hippocampus")
# define vector with behavioral outcome variables
vars_behav = c("PosS","NegS","DepS","gaf","SyCod","VeMem")
# define vector with fMRI outcome variables
vars_fmri = c("fchi","fchil","fchir","fchih")


# INDIVIDUAL AND GROUP AVERAGE COURSES OF RAW VALUES OF VARIABLES ---------------------------------------------

##### Behavioral variables #####
# select only behavioral variables
df_inspect_behav = df_inspect_long %>% filter(outcome %in% vars_behav)
# rename those for plotting
df_inspect_behav = df_inspect_behav %>% mutate(outcome = str_replace_all(outcome, label_map))
# adjust factor levels
df_inspect_behav$outcome = factor(df_inspect_behav$outcome,
                                  levels=c("Att. Positive Symptoms","Negative Symptoms","Depressive Symptoms",
                                           "Levels of Functioning","Symbol Coding","Verbal Memory"))
# create a data frame with means
df_mean_behav = df_inspect_behav %>%
  group_by(month,outcome,group) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")
# create data frame for line drawings (without NAs)
df_lines_behav = df_inspect_behav %>%
  filter(!is.na(value))
# plot
p = ggplot(df_inspect_behav, aes(x=month, y=value)) +
  geom_line(data = df_lines_behav,aes(group = interaction(subj, outcome)),color = "grey", alpha = 0.5, linetype = "dashed") +
  geom_point(size=1,color="grey20") +
  geom_point(data = df_mean_behav, aes(x = month, y = mean_value,color=group), 
             shape = 15, size = 3) +
  geom_line(data = df_mean_behav, aes(x = month, y = mean_value,color=group)) +
  scale_color_manual(values = c("darkred","darkblue")) +
  theme_classic() +
  scale_x_continuous(labels = c("BL", "M2", "M4", "M6", "M8"),breaks = c(0,2,4,6,8), limits = c(-0.5,8.2)) +
  labs(x = "", y = "Variable Score") +
  facet_wrap(~ group + outcome, scales = "free_y",nrow=2) +
  theme(axis.text.x = element_text(size=14, color="black", face="bold"),
        axis.title.y = element_text(size=20, color="black", face="bold"),
        axis.text.y = element_text(size=14, color="black"),
        strip.text = element_text(size = 16, face = "bold"),
        legend.position = "none")
# export plot
ggsave(paste0(out_dir,"Napls_Behav_Course.jpeg"),plot=p,height=350,width=550,unit="mm",dpi=500)

##### Functional connectivity variables #####
# select only fmri variables
df_inspect_fmri = df_inspect_long %>% filter(outcome %in% vars_fmri)
# rename those for plotting
df_inspect_fmri = df_inspect_fmri %>% mutate(outcome = str_replace_all(outcome, label_map))
# adjust factor levels
df_inspect_fmri$outcome = factor(df_inspect_fmri$outcome,
                                 levels=c("FC in Hippocampus","FC in left Hippocampus","FC in right Hippocampus",
                                          "FC between Hippocampi"))
# create a data frame with means
df_mean_fmri = df_inspect_fmri %>%
  group_by(month,outcome,group) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")
# create data frame for line drawings (without NAs)
df_lines_fmri = df_inspect_fmri %>%
  filter(!is.na(value))
# plot
p = ggplot(df_inspect_fmri, aes(x=month, y=value)) +
  geom_line(data = df_lines_fmri,aes(group = interaction(subj, outcome)),color = "grey", alpha = 0.5, linetype = "dashed") +
  geom_point(size=1,color="grey20") +
  geom_point(data = df_mean_fmri, aes(x = month, y = mean_value,color=group), 
             shape = 15, size = 3) +
  geom_line(data = df_mean_fmri, aes(x = month, y = mean_value,color=group)) +
  scale_color_manual(values = c("darkred","darkblue")) +
  theme_classic() +
  scale_x_continuous(labels = c("BL", "M2", "M4", "M6", "M8"),breaks = c(0,2,4,6,8), limits = c(-0.5,8.2)) +
  labs(x = "", y = "FC Value") +
  facet_wrap(~ group + outcome, scales = "free_y",nrow=2) +
  theme(axis.text.x = element_text(size=14, color="black", face="bold"),
        axis.title.y = element_text(size=20, color="black", face="bold"),
        axis.text.y = element_text(size=14, color="black"),
        strip.text = element_text(size = 16, face = "bold"),
        legend.position = "none")
# export plot
ggsave(paste0(out_dir,"Napls_fMRI_Course.jpeg"),plot=p,height=350,width=550,unit="mm",dpi=500)


