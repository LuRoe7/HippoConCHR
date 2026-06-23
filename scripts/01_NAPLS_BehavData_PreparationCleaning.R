# DATA PREPARATION AND CLEANING OF BEHAVIORAL DATA ---------------------------------------------
# Description:  This script prepares and cleans the behavioral data from the NAPLS3 study.
#               
# Author:       Roell, Lukas      
# Created:      2025/08/28
# License:      Creative Common: CC-BY

##### Load packages #####
# load necessary packages
library(pacman)
p_load(readr,tidyr,dplyr,stringr,ggplot2,gridExtra,tibble)


# BASIC PREPARATIONS ---------------------------------------------------------------------------
# show working directory and ensure that the folder behavioral including the data is located there
getwd()
# define input directory
in_dir = "PATH/TO/DIRECTORY/"
# define output directory
out_dir = "PATH/TO/DIRECTORY/"
# create output directory for this script
dir.create(out_dir, recursive = TRUE,showWarnings=F)
# read behavioral data
df_wide_raw = read.csv(paste0(in_dir,"napls_behavioral_wide_old.csv"))
# add s to subject ID
df_wide_raw$src_subject_id = paste0("s",df_wide_raw$src_subject_id)
# insert site column based on first two numbers of subject ID
df_wide_raw = df_wide_raw %>% mutate(site = paste0("site", sprintf("%02d", as.integer(substr(src_subject_id, 2, 3))))) %>% 
  relocate(site, .before = 3)
# adjust entries in site column
df_wide_raw$site = df_wide_raw$site %>% str_replace("0","") %>% str_replace("site","site0")
# define behavioral variables of interest
vars = c("subjectkey","src_subject_id","site","subjecttype",
         "interview_age","sex","educat","primary_language","ethnic_group","wasi_IQ", # demographics + IQ
         colnames(df_wide_raw)[grep("weight_",colnames(df_wide_raw))], # weight over time
         colnames(df_wide_raw)[grep("height_",colnames(df_wide_raw))], # height over time
         colnames(df_wide_raw)[grep("sops_",colnames(df_wide_raw))], # scale of psychotic risk symptoms
         colnames(df_wide_raw)[grep("calg_",colnames(df_wide_raw))], # Calgary Depression = CDSS
         colnames(df_wide_raw)[grep("gaf",colnames(df_wide_raw))], # global functioning,
         colnames(df_wide_raw)[grep("bacs_sc_raw",colnames(df_wide_raw))], # cognition: symbol coding
         colnames(df_wide_raw)[grep("hvltr_raw",colnames(df_wide_raw))] # cognition: verbal memory
         )
# select those variables from data frame
df_wide = df_wide_raw[,vars]
# rename subject ID column
colnames(df_wide) = colnames(df_wide) %>% str_replace("src_subject_id","subject_id")
# convert interview age to years
df_wide$interview_age = df_wide$interview_age/12
# adjust entries in subjecttype columns with shorter names
df_wide$subjecttype = df_wide$subjecttype %>% str_replace("Control","C") %>% str_replace("Non-Enhanced","NE") %>% str_replace("Enhanced","E")
# create vector with categorical baseline variables
vars_categorical_base = c("sex","primary_language","ethnic_group")
# convert categorical baseline variables to strings
df_wide[vars_categorical_base] = lapply(df_wide[vars_categorical_base], as.character)
# define time-independent variables that do not change over measurement time points
vars_base = colnames(df_wide)[!grepl("__",colnames(df_wide))]
# define time-dependent variables that change over measurement time points
vars_long = colnames(df_wide)[grepl("__",colnames(df_wide))]
# select only numeric time-dependent variables
vars_long_numeric = vars_long[sapply(df_wide[vars_long], is.numeric)]
# convert data to very long format with subject, session and time-dependent-variables as columns
df_long = df_wide %>% pivot_longer(cols = all_of(vars_long_numeric),
                                   names_to = c("outcome","session"),
                                   names_sep = "__",values_to = "value")


# VISUALIZATION OF NUMERIC (LONGITUDINAL + BASELINE) AND CATEGORICAL DATA (BASELINE) -----------------------------------

##### Numeric longitudinal variables #####
# create vector with numeric longitudinal variables
vars_long_numeric_plot = unique(df_long$outcome)
# visualize numeric longitudinal variables to detect implausible values or erroneous data entries
p = ggplot(df_long,aes(x=session,y=value)) +
  geom_violin(fill = "lightblue", alpha = 0.5, trim = FALSE) +
  geom_boxplot(fill="dimgrey",outlier.shape = NA,width=0.4) +
  geom_jitter(size=1,color="grey",width = 0.1,alpha=0.3) +
  facet_wrap(~outcome,nrow=7,scales="free") +
  labs(x = "", y = "") +
  theme_classic() +
  theme(axis.title.x = element_text(size=10,color="black",face="bold"),
        axis.text.x = element_text(size=10,color="black"),
        axis.title.y = element_text(size=10,color="black",face="bold"),
        axis.text.y = element_text(size=10,color="black"),
        strip.text = element_text(face='bold'))
# export plot
ggsave(paste0(out_dir,"napls_behavioral_numeric_long.jpeg"),plot=p,height=450,width=500,unit="mm",dpi=500)

##### Numeric baseline variables #####
# subset data frame and keep only baseline session and one random longitudinal numeric variable to avoid artificial increase of sample size
df_long_bl = subset(df_long,session=="BL"&outcome=="weight") # random longitudinal variable = weight
# create vector with numeric baseline variables
vars_numeric_base_plot = c("interview_age","educat","wasi_IQ")
# visualize numeric baseline variables to detect implausible values or erroneous data entries
plot_list_numeric_base = list()
for (variab in vars_numeric_base_plot) {
  plot_list_numeric_base[[variab]] = local({
    variab = variab
    p = ggplot(df_long_bl, aes(x = unlist(df_long_bl[, "subjecttype"]), y = unlist(df_long_bl[, variab]))) +
      geom_violin(fill = "lightblue", alpha = 0.5, trim = FALSE) +
      geom_boxplot(fill = "dimgrey", outlier.shape = NA, width = 0.4) +
      geom_jitter(size = 1, color = "grey", width = 0.1, alpha = 0.3) +
      labs(x = "", y = variab) +
      theme_classic() +
      theme(
        axis.title.x = element_text(size = 10, color = "black", face = "bold"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title.y = element_text(size = 8, color = "black", face = "bold"),
        axis.text.y = element_text(size = 10, color = "black")
      )
  })
}
# arrange plots
plots_numeric_base = grid.arrange(grobs=plot_list_numeric_base ,nrow=1)
# export plot
ggsave(paste0(out_dir,"napls_behavioral_numeric_baseline.jpeg"),plot=plots_numeric_base,height=175,width=300,unit="mm",dpi=500)

##### Categorical baseline variables #####
# insert sample column for plotting
df_long_bl$sample = "napls"
# create vector with categorical baseline variables for plotting
vars_categorical_base_plot = c("site","sex","ethnic_group")
# visualize categorical columns to detect implausible values or erroneous data entries
plot_list_categorical = list()
for (variab in vars_categorical_base_plot) {
  plot_list_categorical[[variab]] = local({
    variab = variab
    p = ggplot(df_long_bl,aes(x=unlist(df_long_bl[,"sample"]),fill=unlist(df_long_bl[,variab]))) +
      geom_bar(position = "dodge") +
      labs(x = "", y = "", fill = variab) +
      theme_classic() +
      theme(axis.title.x = element_text(size=12,color="black",face="bold"),
            axis.text.x = element_text(size=10,color="black"),
            axis.title.y = element_text(size=12,color="black",face="bold"),
            axis.text.y = element_text(size=10,color="black"))
  })
}
# arrange plots with clinical scores
plots_categorical_base = grid.arrange(grobs=plot_list_categorical,nrow=1)
# export plot
ggsave(paste0(out_dir,"napls_behavioral_categorical_base.jpeg"),plot=plots_categorical_base,height=100,width=300,unit="mm",dpi=500)


# CLEANING DATA BASED ON VISUAL INSPECTION OF PLOTS ----------------------------------------------------------------------------
# Note: The created plots were inspected and implausible or lacking data was removed or converted properly, respectively.
# make copy of data frame
df_wide_clean = df_wide
# convert empty entries in ethnicity column to NA
df_wide_clean$ethnic_group[df_wide_clean$ethnic_group == ""] = NA
# remove BACS symbol coding scores below 10 which we assume to be rather a motivation problem than a real cognitive impairment (only 2 entries!)
df_wide_clean$bacs_sc_raw__BL[df_wide_clean$bacs_sc_raw__BL < 10] = NA
df_wide_clean$bacs_sc_raw__M8[df_wide_clean$bacs_sc_raw__M8 < 10] = NA
# remove very short height which we assume to be a typo or to be the weight (impossible to clarify for us)
df_wide_clean$height__BL[df_wide_clean$height__BL < 125] = NA


# EXPORT CLEANED DATA -----------------------------------------------------------------------------------------------------------
# export data frame
write_csv(df_wide_clean,paste0(out_dir,"napls_behavioral_wide_cleaned.csv"))

