# DEMOGRAPHICS ---------------------------------------------
# Description:  This script summarizes the demographics and compares them between CHR-P and HCs
#               
# Author:       Roell, Lukas      
# Created:      2025/08/28
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
df_demo = read.csv(paste0(in_dir,"Napls_BehavFCandSlopes_impSE.csv"))
# remove subjects that have non-enhanced psychosis risk, as only 1 of these subjects has valid data for FC at two time points
df_demo = subset(df_demo,group !=2)
# rename group column to CHR
df_demo = df_demo %>% rename(CHR = group)
# adjust subject column
df_demo$subj = paste0("s",as.character(df_demo$subj))
# convert site to factor
df_demo$site = as.factor(df_demo$site)
# convert sex to factor
df_demo$female = df_demo$female %>% as.character() %>% str_replace("1","F") %>% str_replace("0","M")
# adjust group column
df_demo$CHR = df_demo$CHR %>% as.character() %>% str_replace("1","CHR-P") %>% str_replace("0","HC")


# DEMOGRAPHICS: SUMMARY OF DEMOGRAPHIC MEASURES OF INTEREST -----------------------------------------------------------------------
# Note: Demographics are computed based on all subjects that had at least two measurement time points in negative symptoms and FC of
# hippocampus
# select subject ID, negative symptom change and FC hippocampus change
df_demo_tpl = df_demo[c("subj","NegS_CH","fchi_CH")]
# keep only complete cases that were considered in final analysis for negative symptoms
df_demo_tpl = df_demo_tpl[complete.cases(df_demo_tpl), ]
# save subjects as vector
subjects_demo = df_demo_tpl$subj
# select the subjects considered for final analysis
df_demo = df_demo[df_demo$subj %in% subjects_demo,]
# compute BMI
df_demo$BMI = df_demo$weight_BL/((df_demo$height_BL/100)^2)
# rename entries in language
df_demo = df_demo %>% mutate(language = case_when(
    str_detect(tolower(language), "english") ~ "english",
    str_detect(tolower(language), "spanish") ~ "spanish",
    TRUE ~ "other"
  ))
# select variables of interest and shorten variable names
df_demo = df_demo %>%
  select(subj, site, CHR, age, female, ey, language, ethnicity, IQ,
         BMI,cpz_BL,PosS_BL,NegS_BL,gaf_BL,VeMem_BL,SyCod_BL) %>%
  rename(ID = subj,Site = site,Group = CHR,Age = age,Sex = female,Education = ey,
         Language = language,Ethnicity = ethnicity,IQ = IQ,CPZ = cpz_BL,
         SOPSPos = PosS_BL,SOPSNeg = NegS_BL,GAF = gaf_BL,VerbalMemory = VeMem_BL,SymbolCoding = SyCod_BL)
# define which columns are numeric vs categorical
numeric_vars = c("Age", "Education", "IQ", "BMI","CPZ","SOPSPos", "SOPSNeg", "GAF","VerbalMemory","SymbolCoding")
categorical_vars = c("Site","Sex", "Language", "Ethnicity")

##### NUMERIC DEMOGRAPHICS #####
# summarise numeric demographics at baseline
df_numeric_summary = df_demo %>%
  group_by(Group) %>% # group by group
  summarise(across(all_of(numeric_vars),
                   list(mean = ~mean(.x, na.rm = TRUE),
                        sd = ~sd(.x, na.rm = TRUE),
                        n = ~sum(!is.na(.x))),
                   .names = "{.col}_{.fn}")) %>% # compute mean, SD, and n
  pivot_longer(-Group, names_to = c("Variable", ".value"), names_sep = "_") %>% # bring to long format
  mutate(summary = paste0(round(mean, 2), " (", round(sd, 2), ")")) %>% # add column with merged mean and SD
  select(Group, Variable, summary, n) # select columns of interest
# compute Wilcoxon p-values for each numeric variable
numeric_pvals = map_dfr(numeric_vars, function(var) {
  data_subset = df_demo %>% select(Group, all_of(var)) %>% filter(!is.na(.data[[var]]))
  p = wilcox.test(as.formula(paste(var, "~ Group")), data = data_subset)$p.value # compute Wilcox test
  tibble(Variable = var, p_value = p) # create data frame
})
# merge p-values into numeric summary
df_numeric_summary = df_numeric_summary %>% left_join(numeric_pvals, by = "Variable")


##### CATEGORICAL DEMOGRAPHICS #####
# summarize counts and proportions for categorical variables
df_categorical_summary = df_demo %>%
  select(Group, all_of(categorical_vars)) %>% # select variables
  pivot_longer(-Group, names_to = "Variable", values_to = "Value") %>% # convert to long format
  group_by(Group, Variable, Value) %>% # grou variables
  summarise(n = n(), .groups = "drop_last") %>% # compute N
  mutate(prop = n / sum(n)) %>%  # compute proportion as decimal
  group_by(Group, Variable) %>% # group again
  mutate(n_nonmissing = sum(n)) %>% # compute non missing points
  ungroup() %>% # ungroup 
  mutate(count_prop = paste0(n, " (", round(prop * 100, 1), "%)")) %>% # summarise number and proportion in one column
  select(Group, Variable, Value,count_prop, n_nonmissing) # select variables of interest
# compute Fisher tests
categorical_pvals = map_dfr(categorical_vars, function(var) {
  tab = table(df_demo$Group, df_demo[[var]])
  test = fisher.test(tab)
  tibble(Variable = var, p_value = test$p.value)
})
# merge p-values into categorical summary
df_categorical_summary = df_categorical_summary %>%
  left_join(categorical_pvals, by = "Variable")

##### GET NUMBER OF SUBJECTS WITH ANTIPSYCHOTICS ####
# count CHR-P subjects with a non-zero CPZ value
cpz_nonzero_count = df_demo %>%
  filter(Group == "CHR-P") %>% # keep only CHR-P subjects
  filter(!is.na(CPZ) & CPZ != 0) %>% # exclude missing and zero CPZ values
  nrow() # count remaining rows
# also get the total number of CHR-P subjects with any valid (non-NA) CPZ value, for context
cpz_total_chrp = df_demo %>%
  filter(Group == "CHR-P") %>% # keep only CHR-P subjects
  filter(!is.na(CPZ)) %>% # exclude missing CPZ values only
  nrow() # count remaining rows
# print result as a readable summary
cat("CHR-P subjects with non-zero CPZ:", cpz_nonzero_count,
    "out of", cpz_total_chrp, "with valid CPZ data\n")  # print both counts to console


##### FORMAT P VALUES #####
# define function that formats the p-values
format_p = function(p) {
  case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "p < 0.001***",
    p < 0.01  ~ "p < 0.01**",
    p < 0.05  ~ "p < 0.05*",
    TRUE      ~ paste0("p = ", signif(p, 3))
  )
}
# apply to numeric summary
df_numeric_summary = df_numeric_summary %>%
  mutate(p_value = format_p(p_value))
# apply to categorical summary
df_categorical_summary = df_categorical_summary %>%
  mutate(p_value = format_p(p_value))
# export data frames
write.csv(df_numeric_summary,paste0(out_dir,"Napls_Demographics_Numeric.csv"), row.names = FALSE)
write.csv(df_categorical_summary,paste0(out_dir,"Napls_Demographics_Categorical.csv"), row.names = FALSE)

