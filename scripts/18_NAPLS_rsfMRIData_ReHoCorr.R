# DATA PREPARATION AND CLEANING OF RS-fMRI DATA ---------------------------------------------
# Description:  This script prepares the rs-fMRI data from the NAPLS3 study.
#               
# Author:       Roell, Lukas      
# Created:      2025/04/03
# License:      Creative Common: CC-BY

##### Load packages #####
# load necessary packages
library(pacman)
p_load(readr,tidyr,dplyr,stringr,readxl,ggplot2,gridExtra,scales,viridis)


# BASIC PREPARATIONS ---------------------------------------------------------------------------
# show working directory and ensure that the folder fmri including the data is located there
getwd()
# define input directory
in_dir = "results/napls/P1_HippoFC/02_NAPLS_rsfMRIData_PreparationCleaning/"
# define output directory
out_dir = "results/napls/P1_HippoFC/18_NAPLS_rsfMRIData_ReHoCorr/"
# create output directory for this script
dir.create(out_dir, recursive = TRUE,showWarnings=F)
# read fMRI data
df_fc = read_csv(paste0(in_dir,"napls_rsfMRI_FunctionalConnectivity_cleaned.csv"))
# rename subject column
df_fc = df_fc %>% rename(mri_id = subject_id)
# read ReHo data
load("data/napls/fmri/napls_reho_long_Brainnetome.RData")
# compute mean ReHo
df_reho$reho_hip = rowMeans(df_reho[c("rHippL","rHippR","cHippL","cHippR")],na.rm=T)
# merge by subject and session
df_fmri = merge(df_fc,df_reho,by=c("mri_id","session"),all.x=T)
# read behavioral data
df_behav = read_csv("results/napls/P1_HippoFC/01_NAPLS_BehavData_PreparationCleaning/napls_behavioral_wide_cleaned.csv")
# rename subject column
df_behav = df_behav %>% rename(mri_id = subject_id)
# select columns of interest
df_behav = df_behav[c("mri_id","subjecttype")]
# merge
df = merge(df_fmri, df_behav, by = "mri_id", all.x = TRUE)


# CORRELATION BETWEEN FC AND REHO OF HIPPOCAMPUS ---------------------------------------------------------------------------
# select only baseline
df_cor = subset(df,session=="BL"&subjecttype == "C")
cor_result = cor.test(df_cor$fc_hip, df_cor$reho_hip)  # run correlation test and store full output
# extract test statistics into a single-row data frame
df_cor_sum = data.frame(
  x       = "ReHo in Hippocampus",
  y       = "FC in Hippocampus",
  r       = cor_result$estimate,        # Pearson correlation coefficient
  t       = cor_result$statistic,       # t-statistic
  df      = cor_result$parameter,       # degrees of freedom (n - 2)
  p_value = cor_result$p.value,         # two-tailed p-value
  ci_low  = cor_result$conf.int[1],     # lower bound of 95% CI for r
  ci_high = cor_result$conf.int[2],     # upper bound of 95% CI for r
  row.names = NULL                      # suppress automatic row name
)
write.csv(df_cor_sum,
          paste0(out_dir, "Napls_FC_ReHo_Correlation.csv"),
          row.names = FALSE)

# TEST-RETEST RELIABILITY OF FC ---------------------------------------------------------------------------

##### SELECT BL AND NEXT AVAILABLE TIMEPOINT PER CONTROL SUBJECT #####
# for subjects with more than 2 sessions, retain only BL and the earliest follow-up available,
# prioritizing M2 > M4 > M6 > M8 as the second timepoint

# define session priority order for the follow-up selection
session_order = c("BL", "M2", "M4", "M6", "M8")  # BL is anchor; follow-up selected left to right

df_ctrl = df %>%
  filter(subjecttype == "C") %>%                   # keep only control subjects
  filter(session %in% session_order) %>%           # drop any unexpected session labels
  mutate(session = factor(session, levels = session_order))  # encode session as ordered factor for sorting

# for each subject, keep BL and the single earliest available follow-up
df_ctrl = df_ctrl %>%
  group_by(mri_id) %>%
  arrange(session, .by_group = TRUE) %>%           # sort sessions within subject by priority order
  filter(
    session == "BL" |                              # always keep BL
      session == first(session[session != "BL"])     # keep only the first non-BL session available
  ) %>%
  # exclude subjects who only have BL and no follow-up session at all
  filter(n() > 1) %>% 
  ungroup()


library(psych)  # for ICC(); install if needed: install.packages("psych")

##### COMPUTE TEST-RETEST RELIABILITY VIA ICC ACROSS TWO TIMEPOINTS #####

df_ctrl_wide <- df_ctrl %>%
  select(mri_id, session, fc_hip, reho_hip) %>%
  mutate(session = as.character(session)) %>%
  pivot_wider(
    id_cols     = mri_id,
    names_from  = session,
    values_from = c(fc_hip, reho_hip)
  )

icc_vars <- c("fc_hip")

df_retest <- data.frame(
  variable   = character(),
  n          = numeric(),
  icc        = numeric(),   # ICC(2,1): two-way random, absolute agreement
  ci_low     = numeric(),   # lower 95% CI
  ci_high    = numeric(),   # upper 95% CI
  f          = numeric(),   # F-statistic
  df1        = numeric(),
  df2        = numeric(),
  p_value    = numeric(),
  followup_sessions_used = character(),  # which follow-up sessions contributed data
  stringsAsFactors = FALSE
)

for (var in icc_vars) {
  
  bl_col        <- paste0(var, "_BL")
  followup_cols <- paste0(var, "_", c("M2", "M4", "M6", "M8"))
  followup_cols <- followup_cols[followup_cols %in% names(df_ctrl_wide)]
  
  if (!bl_col %in% names(df_ctrl_wide)) {
    cat("Skipping", var, "— BL column not found\n")
    next
  }
  
  # Collapse follow-up sessions: take earliest available per subject
  # NOTE: Variable retest interval is a limitation — consider restricting to one session
  df_ctrl_wide$followup_val <- apply(
    df_ctrl_wide[, followup_cols, drop = FALSE], 1,
    function(x) x[!is.na(x)][1]
  )
  
  # Track which sessions actually contributed follow-up values
  session_source <- apply(
    df_ctrl_wide[, followup_cols, drop = FALSE], 1,
    function(x) { idx <- which(!is.na(x))[1]; if (!is.na(idx)) names(x)[idx] else NA }
  )
  sessions_used <- paste(sort(unique(na.omit(session_source))), collapse = ", ")
  
  # Restrict to complete pairs
  df_pair <- df_ctrl_wide[
    !is.na(df_ctrl_wide[[bl_col]]) & !is.na(df_ctrl_wide$followup_val), 
    c(bl_col, "followup_val")
  ]
  
  # ICC requires a matrix with one column per timepoint (n_subjects x 2)
  icc_mat <- as.matrix(df_pair)
  
  # psych::ICC() returns all 6 forms; extract row 2 = ICC2 (two-way random, absolute agreement)
  icc_res  <- ICC(icc_mat, missing = FALSE, alpha = 0.05)
  icc2     <- icc_res$results[icc_res$results$type == "ICC2", ]
  
  df_retest <- rbind(df_retest, data.frame(
    variable               = var,
    n                      = nrow(df_pair),
    icc                    = icc2$ICC,
    ci_low                 = icc2$`lower bound`,
    ci_high                = icc2$`upper bound`,
    f                      = icc2$F,
    df1                    = icc2$df1,
    df2                    = icc2$df2,
    p_value                = icc2$p,
    followup_sessions_used = sessions_used,
    stringsAsFactors       = FALSE,
    row.names              = NULL
  ))
}

##### PRINT AND EXPORT #####
print(df_retest)

write.csv(df_retest,
          paste0(out_dir, "Napls_TestRetest_ICC_Controls.csv"),
          row.names = FALSE)



##### COMPUTE TEST-RETEST RELIABILITY VIA CORRELATION ACROSS TWO TIMEPOINTS #####
# for each variable, correlate the value at BL with the value at the follow-up session,
# across all control subjects who have both timepoints available

# pivot to wide format so BL and follow-up are in separate columns per subject
df_ctrl_wide = df_ctrl %>%
  select(mri_id, session, fc_hip, reho_hip) %>%   # keep only relevant columns
  mutate(session = as.character(session)) %>%       # convert factor back to character for pivot
  pivot_wider(
    id_cols     = mri_id,                          # one row per subject
    names_from  = session,                         # column names from session labels
    values_from = c(fc_hip, reho_hip)              # create BL and follow-up columns for each variable
  )

# define variables to test
icc_vars = c("fc_hip", "reho_hip")  # variables to correlate across timepoints

# pre-allocate results
df_retest = data.frame(
  variable = character(),
  n        = numeric(),    # number of subjects with both timepoints
  r        = numeric(),    # Pearson correlation coefficient
  t        = numeric(),    # t-statistic
  df       = numeric(),    # degrees of freedom
  p_value  = numeric(),    # two-tailed p-value
  ci_low   = numeric(),    # lower 95% CI for r
  ci_high  = numeric(),    # upper 95% CI for r
  stringsAsFactors = FALSE
)

for (var in icc_vars) {
  
  # identify BL column and all possible follow-up columns for this variable
  bl_col       = paste0(var, "_BL")                                   # e.g. fc_hip_BL
  followup_cols = paste0(var, "_", c("M2","M4","M6","M8"))            # all possible follow-up columns
  
  # keep only follow-up columns that actually exist after pivoting
  followup_cols = followup_cols[followup_cols %in% names(df_ctrl_wide)]
  
  # skip if BL column is missing entirely
  if (!bl_col %in% names(df_ctrl_wide)) {
    cat("Skipping", var, "— BL column not found\n")
    next
  }
  
  # for each subject, collapse all follow-up columns into one by taking the first non-NA value
  # this picks the earliest available follow-up (M2 > M4 > M6 > M8) since pivot_wider preserves order
  df_ctrl_wide$followup_val = apply(
    df_ctrl_wide[, followup_cols, drop = FALSE], 1,         # apply row-wise across follow-up columns
    function(x) x[!is.na(x)][1]                            # take first non-NA value = earliest follow-up
  )
  
  # keep only subjects with both BL and a follow-up value
  df_pair = df_ctrl_wide[!is.na(df_ctrl_wide[[bl_col]]) & !is.na(df_ctrl_wide$followup_val), ]
  
  # run correlation between BL and follow-up values
  cor_i = cor.test(df_pair[[bl_col]], df_pair$followup_val)  # Pearson correlation
  
  # store results
  df_retest = rbind(df_retest, data.frame(
    variable = var,
    n        = nrow(df_pair),           # number of subjects with complete pairs
    r        = cor_i$estimate,          # correlation coefficient
    t        = cor_i$statistic,         # t-statistic
    df       = cor_i$parameter,         # degrees of freedom (n - 2)
    p_value  = cor_i$p.value,           # two-tailed p-value
    ci_low   = cor_i$conf.int[1],       # lower 95% CI
    ci_high  = cor_i$conf.int[2],       # upper 95% CI
    stringsAsFactors = FALSE,
    row.names = NULL
  ))
}

##### PRINT AND EXPORT #####
write.csv(df_retest,
          paste0(out_dir, "Napls_TestRetest_Correlation_Controls.csv"),
          row.names = FALSE)
