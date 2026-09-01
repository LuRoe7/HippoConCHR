# SINGLE SUBJECT REGRESSIONS ---------------------------------------------
# Description:  This script runs the final statistical analysis.
#               
# Author:       Roell, Lukas      
# Created:      2025/08/28
# License:      Creative Common: CC-BY

##### Load packages #####
# load necessary packages
library(pacman)
p_load(readr,tidyr,dplyr,stringr,ggplot2,ggrepel,gridExtra,ggExtra,readxl,writexl,tibble)


# BASIC PREPARATIONS -------------------------------------------------------------------------------------------------------------
# show working directory and ensure that the folder results including the data is located there
getwd()
# define input directory for behavioral and FC data
in_dir = "results/napls/05_NAPLS_FinalDataPreparation/"
# define output directory
out_dir = "results/napls/07_NAPLS_SingleSubjectRegressions/"
# create output directory for this script
dir.create(out_dir, recursive = TRUE,showWarnings=F)
# read merged data including clinical and behavioral variables
df_napls_wide = read.csv(paste0(in_dir,"napls_MergedData_Unstandardized.csv"))
# convert data frame to long very format with session and variable as columns
df_napls_long = df_napls_wide %>%
  pivot_longer(cols = matches(".*_(BL|M2|M4|M6|M8)$"),
               names_to = c("outcome", "session"),
               names_pattern = "(.*)_(BL|M2|M4|M6|M8)",
               values_to = "value")
# insert numeric session variable to data frame in long format
df_napls_long = df_napls_long %>% mutate(month = case_when(session == "BL" ~ 0,
                                                           session == "M2" ~ 2,
                                                           session == "M4" ~ 4,
                                                           session == "M6" ~ 6,
                                                           session == "M8" ~ 8))
# count available (=non-missing) time points per variable
df_timecounts = df_napls_long %>%
  group_by(subj, outcome) %>%
  summarise(nt = sum(!is.na(value)), .groups = "drop") %>%
  pivot_wider(names_from = outcome, values_from = nt, names_prefix = "nt_")
# add time point information to data frame in wide format
df_napls_wide = df_napls_wide %>% left_join(df_timecounts, by = "subj")
# create vectors with behavioral variables
behavvars = c("PosS","NegS","DepS","gaf","SyCod", "VeMem")
# create vector with medication variable
medvars = c("cpz")
# create vectors with FC variables to study
neurovars = c("fchi","fchil","fchir","fchih","fcsal","fcdmn","fcfpn","fclim","fcsmn","fcvis",
              "fctha","fcamy","fcbg","fcput","fcpal","fccau","fcacu")
# select variables of interest from long format data frame
df_napls_long_sub = df_napls_long %>% filter(outcome %in% c(behavvars,medvars,neurovars))


# SINGLE-SUBJECT REGRESSIONS TO COMPUTE INDIVIDUAL LINEAR CHANGE OVER TIME PER VARIABLE ----------------------------------------------------------------------------------------------------
# Note: All available data points per variable are used to compute individual slopes. This can lead to the situation that
# one subject has a slope for a clinical variable that is computed based on 4 time points, whereas the slope for the MRI variable is
# only computed based on 3 time points.

# run single subject regressions and extract statistics:
# beta, standard errors, r², mean squared error, variance in timepoints (x variable), variance across time in outcome values (y-variable)
df_mlrparams = df_napls_long_sub %>%
  group_by(subj, outcome) %>% # group by subjects and outcome
  filter(sum(!is.na(value)) >= 2) %>% # include only subjects with two or more measurement time points
  group_modify(~ {
    .x_non_na = .x %>% filter(!is.na(value))
    model = lm(value ~ month, data = .x_non_na) # run regression
    coefs = coef(summary(model)) # get coefficients
    est_int  = coefs["(Intercept)", "Estimate"] # extract intercept
    se_int   = coefs["(Intercept)", "Std. Error"] # extract standard error of intercept
    est_slope = coefs["month", "Estimate"] # extract slope
    se_slope  = coefs["month", "Std. Error"] # extract standard error of slope
    r2 = summary(model)$r.squared # compute r² of whole model
    sigma2 = summary(model)$sigma^2 # compute mean squared error of whole model
    var_month = var(.x_non_na$month)  # compute variance of timepoints
    var_value = var(.x_non_na$value)  # compute variance of value
    # create one-row result
    tibble(
      subj = .x$subj[1], # save subject ID
      outcome = .x$outcome[1], # save outcome
      int = est_int, # save intercept
      int_se = se_int, # save standard error of intercept
      slope = est_slope, # save slope
      slope_se = se_slope, # save standard error of slope
      r2 = r2, # save r² of model
      sigma2 = sigma2, # save MSE of model
      var_month = var_month, # save variance in available time points
      var_value = var_value # save variance in available outcome values per time points
    )
  }) %>%
  ungroup()
# convert data frame with model parameters to wide format
df_mlrparams_wide = df_mlrparams %>% pivot_wider(id_cols = subj,names_from = outcome,
                                                 values_from = c(int, int_se,slope,slope_se,r2,sigma2,var_month,var_value),
                                                 names_glue = "{outcome}_{.value}")
# merge with original wide data frame, while keeping also subjects with only one measurement time point in all variables
df_napls_mlr = merge(df_napls_wide, df_mlrparams_wide, by = "subj", all.x = TRUE)
# shorten column names
names(df_napls_mlr) = names(df_napls_mlr) %>% str_replace("_slope_se","_SE") %>% str_replace("_slope","_CH") %>% str_replace("_r2","_R2") %>%
  str_replace("_sigma2","_MSE") %>% str_replace("_var_month","_Vx") %>% str_replace("_var_value","_Vy")
# dummy code sex and group, but keep them as factors
df_napls_mlr$sex = df_napls_mlr$sex %>% str_replace("F","1") %>% str_replace("M","0") %>% as.factor()
df_napls_mlr$group = df_napls_mlr$group %>% str_replace("NE","2") %>% str_replace("E","1") %>% str_replace("C","0") %>% as.factor()
# rename both columns to female and CHR
df_napls_mlr = df_napls_mlr %>% rename(female = sex)
# remove the s from subject column and convert it to factor
df_napls_mlr$subj = df_napls_mlr$subj %>% str_replace("s","") %>% as.factor()
# create a dummy coded site column
df_napls_mlr = cbind(df_napls_mlr, model.matrix(~ site - 1, data = df_napls_mlr))
# select the new dummy coded site columns
site_columns = names(df_napls_mlr)[grepl("^sitesite", names(df_napls_mlr))] 
# adjust column names
colnames(df_napls_mlr) = colnames(df_napls_mlr) %>% str_replace("sitesite","site")
# rename entries in site and convert as factor
df_napls_mlr$site = df_napls_mlr$site %>% str_replace("site0","") %>% as.factor()
# convert missings to NA
df_napls_mlr[df_napls_mlr == "NaN"] = NA


# IMPUTING MISSING STANDARD ERRORS FOR SUBJECTS WITH ONLY TWO MEASUREMENT TIMEPOINTS ----------------------------------------------------------------------------------------------------
# Note: Because the accuracy of the estimation of linear change over time per subject in each outcome differs
# depending on the available time points per subject in each outcome, we aim to account for this. To this end, we need the
# standard errors of the linear change estimate from the single subject regressions.
# As subjects with only two measurement time points do not have a standard error, we impute them using the number of available time
# points, the variance across time and the variance in the outcome values as predictors.

# create output directory for imputation
dir.create(paste0(out_dir,"SE_Imputation/"), showWarnings = FALSE, recursive = TRUE)
# copy data framae
df_napls_mlr_imp = df_napls_mlr
# train regression model to predict standard erros (based on subjects with 3, 4, and 5 timepoints) and output relevant test statistics
for (outcome in c(behavvars, neurovars,medvars)) {
  # create formula
  formula_str = paste(paste0(outcome, "_SE"),"~",
                      paste0(outcome, "_CH"),"*",
                      paste0("nt_", outcome),"*",
                      paste0(outcome, "_Vx"),"*",
                      paste0(outcome, "_Vy"))
  # run regression model
  res = lm(as.formula(formula_str), data = df_napls_mlr_imp)
  # extract predictor-level statistics
  coef_table = summary(res)$coefficients
  conf_int = confint(res)
  pred_summary = data.frame(
    Predictor = rownames(coef_table),
    Beta = coef_table[, "Estimate"],
    CI_lower = conf_int[, 1],
    CI_upper = conf_int[, 2],
    p_value = coef_table[, "Pr(>|t|)"])
  # export predictor stats
  write.csv(pred_summary, file = paste0(out_dir,"SE_Imputation/ImputationSE_predictors_", outcome, ".csv"), row.names = FALSE)
  # extract model-level statistics
  model_summary = summary(res)
  r2 = model_summary$r.squared # r2
  adj_r2 = model_summary$adj.r.squared # adjusted r2
  fstat = model_summary$fstatistic # F stats
  model_p = pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE) # model p-value
  # count missing SEs and compute imputation stats
  outcome_se_var = paste0(outcome, "_SE")
  outcome_nt_var = paste0("nt_", outcome)
  rows_to_impute = which(is.na(df_napls_mlr_imp[[outcome_se_var]]) & df_napls_mlr_imp[[outcome_nt_var]] == 2)
  n_imputed = length(rows_to_impute) # number of imputed SEs
  n_total = sum(df_napls_mlr_imp[[outcome_nt_var]] > 1, na.rm = TRUE) # total number of subjects with more than 2 time points
  prop_imputed = ifelse(n_total > 0, n_imputed / n_total, NA) # proportion
  # extract rows to impute
  df_test = df_napls_mlr_imp[rows_to_impute, ]
  # predict SEs
  se_pred = predict(res, newdata = df_test)
  # count negative predicted SEs
  n_negative = sum(se_pred < 0)
  prop_negative = ifelse(n_imputed > 0, n_negative / n_imputed, NA)
  # impute values
  df_napls_mlr_imp[rows_to_impute, outcome_se_var] = se_pred
  # save model stats as data frame
  model_stats_df = data.frame(Outcome = outcome,
                              R2 = r2,
                              Adjusted_R2 = adj_r2,
                              Model_p_value = model_p,
                              N_total_with_nt_larger2 = n_total,
                              N_imputed = n_imputed,
                              Prop_imputed = prop_imputed,
                              N_negative_predicted = n_negative,
                              Prop_negative_predicted = prop_negative)
  # export model stats
  write.csv(model_stats_df, file = paste0(out_dir,"SE_Imputation/ImputationSE_model_", outcome, ".csv"), row.names = FALSE)
}
# define path to CSV files containing imputation results
imp_path = paste0(out_dir,"SE_Imputation/")
# list all csv files matching the pattern
imp_files = list.files(path = imp_path,pattern = "^ImputationSE_model_.*\\.csv$",full.names = TRUE)
# read them all and bind them row-wise
df_imp_all = imp_files %>% lapply(read.csv) %>% bind_rows()
# export data
write.csv(df_imp_all,file = paste0(out_dir,"ImputationSE_Summary.csv"), row.names = FALSE)


# EXPORT DATA FRAMES ----------------------------------------------------------------------------------------------------
# export data frame without imputed standard errors
write.csv(df_napls_mlr,file = paste0(out_dir,"Napls_BehavFCandSlopes_rawSE.csv"), row.names = FALSE)
# export data frame with imputed standard errors
write.csv(df_napls_mlr_imp,file = paste0(out_dir,"Napls_BehavFCandSlopes_impSE.csv"), row.names = FALSE)

