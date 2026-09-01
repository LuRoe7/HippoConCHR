# LINEARITY CHECKS ---------------------------------------------
# Description:  This script checks linearity of change scores.
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
in_dir = "results/napls/P1_HippoFC/07_NAPLS_SingleSubjectRegressions/"
# define output directory
out_dir = "results/napls/P1_HippoFC/17_NAPLS_SingleSubjectRegressions_LinearityCheck/"
# create output directory for this script
dir.create(out_dir, recursive = TRUE,showWarnings=F)
# read merged data including clinical and behavioral variables
df_napls_mlr_imp = read.csv(paste0(in_dir,"Napls_BehavFCandSlopes_impSE.csv"))
# remove subjects that have non-enhanced psychosis risk, as only 1 of these subjects has valid data for FC at two time points
df_napls_mlr_imp = subset(df_napls_mlr_imp,group !=2)
# rename group column to CHR
df_napls_mlr_imp = df_napls_mlr_imp %>% rename(CHR = group)
# select subjects with existent change score in negative symptoms and hippocampaus
df = df_napls_mlr_imp[complete.cases(df_napls_mlr_imp[c("NegS_CH","fchi_CH")]), ]
# select only CHR subjects
df = subset(df,CHR==1)
# create vectors with behavioral variables
behavvars = c("PosS","NegS","DepS","gaf","SyCod", "VeMem")
# create vector for hippocampal FC
neurovars = c("fchi")
# combine both variable vectors into one list to loop over
allvars = c(behavvars, neurovars)  # all variables to process: behavioral + neuro


# LINEARITY CHECKS -------------------------------------------------------------------------------------------------------------

##### COMPUTE AVERAGE R2 PER VARIABLE #####
# for each behavioral and neuro variable, compute the mean and SD of R2 across subjects, but only include subjects who have more than 2 time points (nt_ > 2) for that variable,
# since R2 from a regression with only 2 points is trivially 1 and uninformative
# pre-allocate an empty data frame to collect results across all variables
df_r2_summary = data.frame(
  variable  = character(),  # variable name (e.g. "NegS", "fchi")
  n         = numeric(),    # number of subjects included after nt_ filter
  mean_r2   = numeric(),    # average R2 across included subjects
  sd_r2     = numeric(),     # SD of R2 across included subjects
  min_r2    = numeric(),
  max_r2    = numeric()
)
# loop over each variable
for (var in allvars) {
  r2_col = paste0(var, "_R2")   # name of the R2 column for this variable (e.g. "NegS_R2")
  nt_col = paste0("nt_", var)   # name of the time point count column for this variable (e.g. "nt_NegS")
  # check that both columns exist in the data frame before proceeding
  if (!r2_col %in% names(df) | !nt_col %in% names(df)) {
    cat("Skipping", var, "— column(s) not found in df\n")  # warn if a column is missing
    next  # skip to the next variable
  }
  # keep only subjects with more than 2 time points for this variable
  # using > 2 excludes subjects where R2 = 1 by definition (only 2 points perfectly fit any line)
  df_sub = df[!is.na(df[[nt_col]]) & df[[nt_col]] > 2, ]  # filter by nt_ threshold, also drop NA nt_ rows
  # compute mean and SD of R2 across the filtered subjects, ignoring any remaining NAs in R2
  mean_r2 = mean(df_sub[[r2_col]], na.rm = TRUE)  # average R2 for this variable
  sd_r2   = sd(df_sub[[r2_col]],   na.rm = TRUE)  # SD of R2 for this variable
  min_r2   = min(df_sub[[r2_col]],   na.rm = TRUE)  # Minimum of R2 for this variable
  max_r2   = max(df_sub[[r2_col]],   na.rm = TRUE)  # Maximum of R2 for this variable
  n       = sum(!is.na(df_sub[[r2_col]]))          # number of subjects with non-missing R2 after filter
  # append results for this variable as a new row in the summary data frame
  df_r2_summary = rbind(df_r2_summary, data.frame(
    variable = var,      # variable name
    n        = n,        # sample size after nt_ > 2 filter
    mean_r2  = mean_r2,  # mean R2
    sd_r2    = sd_r2,     # SD of R2
    min_r2   = min_r2,
    max_r2   = max_r2
  ))
}
write.csv(df_r2_summary,
          paste0(out_dir, "Napls_AverageR2_PerVariable.csv"),
          row.names = FALSE)  # export without row index column

##### RESET TEST FOR LINEARITY PER VARIABLE AND SUBJECT #####
# the Ramsey RESET test checks whether adding polynomial terms of the fitted values
# (fitted^2, fitted^3) significantly improves the linear model.
# a non-significant result = no evidence against linearity.

library(lmtest)  # provides the resettest() function

# time points as numeric values matching the column suffixes
timepoints = c(0, 2, 4, 6, 8)           # BL=0, M2=2, M4=4, M6=6, M8=8
tp_labels  = c("BL","M2","M4","M6","M8") # column suffixes corresponding to each time point


##### LOOP OVER ALL VARIABLES #####
reset_results_list = list()  # empty list; one entry per variable

for (var in allvars) {
  
  nt_col  = paste0("nt_", var)   # column tracking how many non-missing time points each subject has
  tp_cols = paste0(var, "_", tp_labels)  # e.g. NegS_BL, NegS_M2, ... NegS_M8
  
  # check required columns exist
  if (!all(tp_cols %in% names(df)) | !nt_col %in% names(df)) {
    cat("Skipping", var, "— columns not found\n")
    next
  }
  
  # pre-allocate a results data frame for this variable: one row per subject
  df_var = data.frame(
    subj      = df$subj,
    variable  = var,
    nt        = df[[nt_col]],
    reset_F   = NA_real_,
    reset_df1 = NA_real_,
    reset_df2 = NA_real_,
    reset_p   = NA_real_
  )
  
  for (i in seq_len(nrow(df))) {
    
    # extract observed values and drop missing time points
    y_vals = as.numeric(df[i, tp_cols])
    df_i   = data.frame(time = timepoints, score = y_vals)
    df_i   = df_i[!is.na(df_i$score), ]  # drop time points with missing scores
    
    # use the actual row count in df_i as the authoritative time point count
    # nt_ column may differ if it was computed differently; this ensures consistency
    n_obs = nrow(df_i)  # true number of observed time points for this subject and variable
    
    # skip if fewer than 5 observed time points:
    # with power = 2 the RESET test adds 1 term to a 2-parameter model (intercept + time),
    # leaving 5 - 3 = 2 residual df — the minimum for a stable F-test.
    # with only 4 points residual df = 1, which produces unstable or NaN F-statistics
    # when the added polynomial term explains all remaining variance.
    if (n_obs < 4) next
    
    # fit the linear model
    lm_i = lm(score ~ time, data = df_i)
    
    # run RESET test with power = 2 only (adds fitted^2 as single additional term).
    # using power = 2:3 with few time points often yields NaN because the augmented model
    # consumes too many df, making the F-statistic undefined or negative.
    reset_i = tryCatch(
      resettest(lm_i, power = 2, type = "fitted"),  # single polynomial term — more stable with 5 points
      error   = function(e) NULL,   # return NULL on hard errors
      warning = function(w) NULL    # also catch warnings (e.g. NaN in pf) and skip rather than store bad values
    )
    
    # only store result if test ran cleanly and p-value is a valid number
    if (!is.null(reset_i) && !is.nan(reset_i$p.value)) {
      df_var$reset_F[i]   = reset_i$statistic
      df_var$reset_df1[i] = reset_i$parameter[1]
      df_var$reset_df2[i] = reset_i$parameter[2]
      df_var$reset_p[i]   = reset_i$p.value
    }
  }
  
  reset_results_list[[var]] = df_var
}


##### COMBINE ALL RESULTS INTO ONE DATA FRAME #####
df_reset_all = do.call(rbind, reset_results_list)  # stack results for all variables into one data frame
rownames(df_reset_all) = NULL                       # remove row names carried over from rbind


##### SUMMARIZE RESET TEST RESULTS PER VARIABLE #####
# for each variable, compute: how many subjects were tested, mean F, and proportion with p < .05
df_reset_summary = df_reset_all %>%
  filter(!is.na(reset_p)) %>%                         # keep only subjects where test was computable
  group_by(variable) %>%                              # summarize separately for each variable
  summarise(
    n_tested      = n(),                              # number of subjects for whom RESET test ran
    mean_F        = round(mean(reset_F), 3),          # average F-statistic across subjects
    sd_F          = round(sd(reset_F), 3),            # SD of F-statistic across subjects
    mean_p        = round(mean(reset_p), 3),          # average p-value across subjects
    pct_sig       = round(mean(reset_p < .05) * 100, 1)  # % of subjects showing significant non-linearity
  )


##### PRINT AND EXPORT #####
write.csv(df_reset_all,                                               # full per-subject results
          paste0(out_dir, "Napls_RESET_PerSubject.csv"),
          row.names = FALSE)

write.csv(df_reset_summary,                                           # summary across subjects per variable
          paste0(out_dir, "Napls_RESET_Summary.csv"),
          row.names = FALSE)





