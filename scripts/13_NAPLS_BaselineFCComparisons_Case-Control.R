# BASELINE COMPARISONS IN HIPPOCAMPAL FC  ---------------------------------------------
# Description:  This script compares baseline FC in the hippocampus between patients and controls and tests for baseline associations
# between symptoms and hippocampal FC
#               
# Author:       Roell, Lukas      
# Created:      2025/09/05
# License:      Creative Common: CC-BY

##### Load packages #####
# load necessary packages
library(pacman)
p_load(readr,tidyr,dplyr,stringr,ggplot2,ggrepel,gridExtra,ggExtra,readxl,writexl,lme4,broom,lavaan,viridis,glue,purrr,reshape2,
       MplusAutomation,tibble)


# BASIC PREPARATIONS -------------------------------------------------------------------------------------------------------------
# show working directory and ensure that the folder results including the data is located there
getwd()
# define input directory for behavioral and FC data
in_dir = "PATH/TO/DIRECTORY/"
# define output directory
out_dir = "PATH/TO/DIRECTORY/"
# create output directory for this script
dir.create(out_dir, recursive = TRUE,showWarnings=F)
# read data including behavioral and FC data and all slopes with imputed SEs
df_napls_mlr_imp = read.csv(paste0(in_dir,"Napls_BehavFCandSlopes_impSE.csv"))
# remove non-enhanced subjects
df_napls_mlr_imp = subset(df_napls_mlr_imp,group !=2)
# create vectors with behavioral variables
behavvars = c("NegS","DepS","gaf")
# create vector for hippocampal FC
neurovars = c("fchi")


# CASE-CONTROL COMPARISON IN FC  ---------------------
# Note: This code runs a multiple linear regression comparing baseline hippocampal FC between patients and controls.
# make copy
df_mlr = df_napls_mlr_imp
# convert group column
df_mlr$group = df_mlr$group %>% as.character() %>% str_replace("1","CHR-P") %>% str_replace("0","HC")
# convert sex column to factor
df_mlr$female = df_mlr$female %>% as.factor()
# convert site column to factor
df_mlr$site = df_mlr$site %>% as.factor()
# run regression
fit = lm(fchi_BL ~ group + age + female + site, data = df_mlr)
# extract summary
s = summary(fit)
# extract coefficients as a data frame
df_stats = as.data.frame(s$coefficients)
colnames(df_stats) = c("Estimate", "StdError", "tValue", "pValue")
# add confidence intervals
ci = confint(fit)
df_stats$CI_low = ci[, 1]
df_stats$CI_high = ci[, 2]
# add sample size
df_stats$N = nobs(fit)
# export test statistics
write.csv(df_stats,file.path(out_dir, "Napls_Case-Control_FCHippocampus_ResultsSummary.csv"), row.names = TRUE)
# visualize results
p = ggplot(df_mlr,aes(x=group,y=fchi_BL)) +
  geom_violin(aes(fill=group),alpha=0.6) +
  geom_boxplot(aes(fill=group),outlier.shape=NA,width=0.4) +
  geom_jitter(width=0.1,alpha=0.2) +
  scale_fill_manual(values=c("darkseagreen4", "grey60")) +
  theme_classic() +
  labs(x="",y = "FC within Hippocampus at Baseline") +
  theme(axis.text.x  = element_text(size = 16, color = "black",face = "bold"),
        axis.title.y = element_text(size = 16, color = "black", face = "bold"),
        axis.text.y  = element_text(size = 14, color = "black"),
        legend.position = "none")
# export plot
ggsave(paste0(out_dir,"Napls_Case-Control_FCHippocampus_ResultsSummary.jpeg"),
       plot=p,height=150,width=100,unit="mm",dpi=500)


# FC-SYMPTOMS ASSOCIATIONS AT BASELINE ---------------------
# Note: This code runs a multiple linear regression comparing baseline hippocampal FC between patients and controls.
# make copy
df_mlr = df_napls_mlr_imp
# convert group column
df_mlr$group = df_mlr$group %>% as.character() %>% str_replace("1","CHR-P") %>% str_replace("0","HC")
# convert sex column to factor
df_mlr$female = df_mlr$female %>% as.factor()
# convert site column to factor
df_mlr$site = df_mlr$site %>% as.factor()
# run regression loop in patients and controls
for (bvar in behavvars) {
  # run regression
  fit = lm(formula(paste0(paste0(bvar,"_BL")," ~ fchi_BL + group + fchi_BL*group + age + female + site")), data = df_mlr)
  # extract summary
  s = summary(fit)
  # extract coefficients as a data frame
  df_stats = as.data.frame(s$coefficients)
  colnames(df_stats) = c("Estimate", "StdError", "tValue", "pValue")
  # add confidence intervals
  ci = confint(fit)
  df_stats$CI_low = ci[, 1]
  df_stats$CI_high = ci[, 2]
  # add sample size
  df_stats$N = nobs(fit)
  # export test statistics
  write.csv(df_stats,file.path(out_dir, paste0("Napls_FCHippocampus-",bvar,"_ResultsSummary_Case-Control.csv")), row.names = TRUE)
}
# select only patients
df_mlr = subset(df_mlr,group == "CHR-P")
# run regression loop only in patients
for (bvar in behavvars) {
  # run regression
  fit = lm(formula(paste0(paste0(bvar,"_BL")," ~ fchi_BL + age + female + site")), data = df_mlr)
  # extract summary
  s = summary(fit)
  # extract coefficients as a data frame
  df_stats = as.data.frame(s$coefficients)
  colnames(df_stats) = c("Estimate", "StdError", "tValue", "pValue")
  # add confidence intervals
  ci = confint(fit)
  df_stats$CI_low = ci[, 1]
  df_stats$CI_high = ci[, 2]
  # add sample size
  df_stats$N = nobs(fit)
  # export test statistics
  write.csv(df_stats,file.path(out_dir, paste0("Napls_FCHippocampus-",bvar,"_ResultsSummary_OnlyPatients.csv")), row.names = TRUE)
}
# make copy
df_plot = df_mlr
# select variables of interest
df_plot = df_plot[c("subj","fchi_BL","NegS_BL","DepS_BL","gaf_BL")]
# visualize
df_plot = df_plot %>%
  pivot_longer(cols = c(NegS_BL, DepS_BL, gaf_BL),names_to = "scale",values_to = "value") %>%
  mutate(scale = recode(scale,
      NegS_BL = "Negative Symptoms",
      DepS_BL = "Depressive Symptoms",
      gaf_BL  = "Psychosocial Functioning"))
  

# visualize results
p = ggplot(df_plot,aes(x=fchi_BL,y=value)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_classic() +
  labs(x = "FC within Hippocampus at Baseline",y="Clinical Score at Baseline") +
  facet_wrap(~scale,scales="free") +
  theme(axis.title.x = element_text(size = 16, color = "black"),
        axis.text.x  = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 16, color = "black"),
        axis.text.y  = element_text(size = 14, color = "black"),
        strip.text = element_text(size=16,color="black",face="bold"))
p
# export plot
ggsave(paste0(out_dir,"Napls_FCHippocampus-Symptoms_ResultsSummary.jpeg"),
       plot=p,height=150,width=250,unit="mm",dpi=500)





















