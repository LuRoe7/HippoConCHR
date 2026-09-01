# FC AND SYMPTOM CHANGE PREDICTIONS  ---------------------------------------------
# Description:  This script runs a prediction of symptom course based on FC Change and vice versa using regression models.
# It considers up to four measurement time points to compute the early change of FC and correct both slopes for measurement error.
#               
# Author:       Roell, Lukas      
# Created:      2025/09/12
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
in_dir = "results/napls/P1_HippoFC/07_NAPLS_SingleSubjectRegressions/"
# define output directory
out_dir = "results/napls/P1_HippoFC/12b_NAPLS_FCandSymptomPredictions_LateChange_OnlyPatients/"
# create output directory for this script
dir.create(out_dir, recursive = TRUE,showWarnings=F)
# read data including behavioral and FC data and all slopes with imputed SEs
df_napls_mlr_imp = read.csv(paste0(in_dir,"Napls_BehavFCandSlopes_impSE.csv"))
# keep CHR individuals 
df_napls_mlr_imp = subset(df_napls_mlr_imp,group == 1)
# rename group column
df_napls_mlr_imp = df_napls_mlr_imp %>% rename(CHR = group)
# create vectors with behavioral variables
behavvars = c("NegS","DepS","gaf")
# create vector for hippocampal FC
neurovars = c("fchi")


# LATENT CHANGE PREDICTION: PREDICTING SYMPTOM CHANGE BASED ON PRIOR FC CHANGE ONLY IN PATIENTS  ---------------------
# run loop
for (bvar in behavvars) {
  for (nvar in neurovars) {
    
    ##### Preparations #####
    # select variables of interest
    df = df_napls_mlr_imp[c("subj","CHR","age","female","site02","site03","site04","site05","site06","site07","site08","site09",
                            paste0("nt_",bvar),paste0("nt_",nvar),
                            paste0(bvar,"_BL"),paste0(bvar,"_M2"),paste0(bvar,"_M4"),paste0(bvar,"_M6"),paste0(bvar,"_M8"),
                            paste0(nvar,"_BL"),paste0(nvar,"_M2"),paste0(nvar,"_M4"),paste0(nvar,"_M6"),paste0(nvar,"_M8"))]
    # select only subjects with at least 2 measurement timepoints in FC and clinical variable and baseline value in FC
    df = df[!is.na(df[[paste0(nvar, "_BL")]]) & df[[paste0("nt_",bvar)]] > 1 & df[[paste0("nt_",nvar)]] > 1, ]
    # reset rownames
    rownames(df) = NULL
    # map the month labels to their values
    months_map = c(BL = 0, M2 = 2, M4 = 4, M6 = 6, M8 = 8)
    # timepoints to consider per variable
    tps_n = c("BL","M2","M4","M6") # for early FC change
    tps_b = c("M2","M4","M6","M8") # for late clinical change
    # define column names based on these timepoints
    cols_orig_n = paste0(nvar, "_", tps_n)
    cols_orig_b = paste0(bvar, "_", tps_b)
    # define number of rows (= subjects)
    N = nrow(df)
    # insert columns to store some important data on early slope
    slope_n       = rep(NA_real_, N)  # early FC slope
    slope_n_se    = rep(NA_real_, N) # standard error of the early slope
    slope_n_nt    = rep(NA_integer_, N) # number of time points used to compute the early slope
    slope_n_var_m = rep(NA_real_, N) # variance across months
    slope_n_var_v = rep(NA_real_, N) # variance of FC values over time
    earliest_n_month= rep(NA_real_, N)  # earliest month used to compute the early slope
    earliest_n_tp   = rep(NA_character_, N) # label of earliest month to compute the early slope
    latest_n_month= rep(NA_real_, N)  # latest month used to compute the early slope
    latest_n_tp   = rep(NA_character_, N) # label of latest month to compute the early slope
    # insert columns to store some important data on late slope
    slope_b       = rep(NA_real_, N)  # late clinical slope
    slope_b_se    = rep(NA_real_, N) # standard error of the late slope
    slope_b_nt    = rep(NA_integer_, N) # number of time points used to compute the early slope
    slope_b_var_m = rep(NA_real_, N) # variance across months
    slope_b_var_v = rep(NA_real_, N) # variance of clinical values over time
    earliest_b_month= rep(NA_real_, N)  # earliest month used to compute the late slope
    earliest_b_tp   = rep(NA_character_, N) # label of earliest month to compute the late slope
    latest_b_month= rep(NA_real_, N)  # latest month used to compute the late slope
    latest_b_tp   = rep(NA_character_, N) # label of latest month to compute the late slope
    
    ##### Estimating early FC slopes and late clinical slopes for each subject #####
    # define maximal time points to consider for each variable to compute the early and late slopes
    max_pts = 4
    # loop over subjects to identify early FC slope and late clinical slope for each variable
    for (i in seq_len(N)) {
      # extract and preprocess data
      vals_n = df[i, cols_orig_n] # get subject-specific FC values per time point 
      colnames(vals_n) = colnames(vals_n) %>% str_replace(paste0(nvar,"_"),"") # rename columns to keep only session label
      vals_b = df[i, cols_orig_b]  # get subject-specific clinical values per time point
      colnames(vals_b) = colnames(vals_b) %>% str_replace(paste0(bvar,"_"),"") # rename columns to keep only session label
      # convert column names to months
      cols_n = colnames(vals_n) # save column names as character vector
      cols_b = colnames(vals_b) # save column names as character vector
      months_n = gsub("^BL$", "0", cols_n) # relabel BL to 0 for FC data
      months_n = gsub("^M", "", months_n) # extract month as integer for FC data
      months_n = as.numeric(months_n) # convert to numeric
      months_b = gsub("^M", "", cols_b) # extract month as integer for clinical data
      months_b = as.numeric(months_b) # convert to numeric
      # keep only non-missing values for both FC and clinical data
      vals_n_num = as.numeric(vals_n[1, ]) # save values as numeric vector
      vals_b_num = as.numeric(vals_b[1, ])  # save values as numeric vector
      keep_n = !is.na(vals_n_num) # remove missings
      keep_b = !is.na(vals_b_num) # remove missings
      # redefine FC timepoints only based on non-missing ones
      cols_n = cols_n[keep_n] # session names
      months_n = months_n[keep_n] # session months
      vals_n_num = vals_n_num[keep_n] # values of FC per session
      # redefine clinical timepoints only based on non-missing ones
      cols_b = cols_b[keep_b] # session names
      months_b = months_b[keep_b] # session months
      vals_b_num = vals_b_num[keep_b] # values of clinical variable per session
      # define lengths available
      available_n = length(vals_n_num)
      available_b = length(vals_b_num)
      # extract the maximums
      max_n = min(available_n, max_pts)
      max_b = min(available_b, max_pts)
      # initiziale empty objects that are used for the subsequent decision tree
      chosen_fc = NULL
      chosen_clin = NULL
      # initialize boolean object that indicates if a subject fulfilled any of the subsequent cases
      success = FALSE
      
      # CASE 1: Runs, if a subject has 3 available time points for FC and clinical variable -> will use BL, M2 and M4 for FC slope and M4, M6 and M8 for
      # clinical slope if available
      if (available_n >= 3 & available_b >= 3) {
        # create combinations of sessions
        combs_fc = combn(seq_along(cols_n), 3, simplify = FALSE) # create possible combinations between FC time points
        combs_clin = combn(seq_along(cols_b), 3, simplify = FALSE) # create combinations between clinical time points
        # ensure that baseline is considered for FC
        if (any(toupper(cols_n) == "BL")) {
          has_bl = sapply(combs_fc, function(i) any(toupper(cols_n[i]) == "BL"))
          combs_fc = c(combs_fc[has_bl], combs_fc[!has_bl])
        }
        # sort session combinations properly
        combs_fc = combs_fc[order(sapply(combs_fc, function(i) max(months_n[i])), decreasing = TRUE)]
        combs_clin = combs_clin[order(sapply(combs_clin, function(i) min(months_b[i])))]
        # loop over all clinical combinations
        for (ci in combs_clin) {
          first_clin = min(months_b[ci]) # get earliest month in this clinical combo
          for (fi in combs_fc) { # loop over all FC combinations
            last_fc = max(months_n[fi]) # get latest month in this FC combo
            inter = intersect(cols_n[fi], cols_b[ci]) # find overlapping time points between FC and clinical
            if (last_fc <= first_clin) { # only valid if FC ends before or at clinical start
              if (length(inter) == 0 || (length(inter) == 1 && last_fc == first_clin)) { # allow no overlap, or exactly one overlap if it's the boundary point
                chosen_fc = fi # store FC combo
                chosen_clin = ci # store clinical combo
                success = TRUE # mark as successful pairing
                break # exit FC loop
              }
            }
          }
          if (success) break # exit clinical loop if success found
        }
      }
      
      # CASE 2: Runs, if there is two or more time points for each variable and no match was found in case 1 -> tries to maximize the clinical time points and to
      # reduce the FC time point in order to find a fit between early FC and late clinical slope
      if (!success & available_n >= 2 & available_b >= 2) { # only run if case 1 failed and at least 2 FC and clinical points exist
        for (clin_k in seq(min(max_b, available_b), 2, by = -1)) { # loop over decreasing numbers of clinical points (from max down to 2)
          combs_clin = combn(seq_along(cols_b), clin_k, simplify = FALSE) # generate all clinical combinations of size clin_k
          combs_clin = combs_clin[order(sapply(combs_clin, function(i) min(months_b[i])))] # sort clinical combinations by earliest month
          for (ci in combs_clin) { # loop over clinical combinations
            first_clin = min(months_b[ci]) # earliest month in this clinical combo
            fc_candidates = which(months_n <= first_clin) # select FC points that occur before or at first clinical month
            bl_idx = which(toupper(cols_n) == "BL") # index of BL if present
            if (length(bl_idx) > 0 & !(bl_idx %in% fc_candidates)) { # if BL exists but is not yet included
              fc_candidates = sort(c(bl_idx, fc_candidates)) # force include BL
            }
            fc_candidates = fc_candidates[order(months_n[fc_candidates])] # sort FC candidates by month ascending
            if (length(fc_candidates) > max_n) { # keep at most max_n FC points
              fc_candidates = fc_candidates[1:max_n]
            }
            if (length(fc_candidates) >= 2) { # require at least 2 FC points
              chosen_fc = fc_candidates # store chosen FC indices
              chosen_clin = ci # store chosen clinical indices
              success = TRUE # mark success
              break # exit clinical loop
            }
          }
          if (success) break # exit loop if success achieved
        }
      }
      
      # CASE 3: Runs, if there is two or more time points for each variable and no match was found in case 2 -> tries to maximize the FC time points and to
      # reduce the clinical time point in order to find a fit between early FC and late clinical slope
      if (!success & available_n >= 2 & available_b >= 2) { # only run if case 1 and case 2 failed, and both FC and clinical have at least 2 points
        for (fc_k in seq(min(max_n, available_n), 2, by = -1)) { # loop over decreasing numbers of FC points (from max down to 2)
          combs_fc = combn(seq_along(cols_n), fc_k, simplify = FALSE) # generate all FC combinations of size fc_k
          if (any(toupper(cols_n) == "BL")) { # if BL exists among FC timepoints
            has_bl = sapply(combs_fc, function(i) any(toupper(cols_n[i]) == "BL")) # check which combinations include BL
            combs_fc = c(combs_fc[has_bl], combs_fc[!has_bl]) # prioritize combinations that contain BL
          }
          combs_fc = combs_fc[order(sapply(combs_fc, function(i) max(months_n[i])), decreasing = TRUE)] # sort FC combinations by latest month (descending)
          
          for (fi in combs_fc) { # loop over FC combinations
            last_fc = max(months_n[fi]) # get latest month in this FC combo
            
            for (clin_k in seq(min(max_b, available_b), 2, by = -1)) { # loop over decreasing numbers of clinical points
              combs_clin = combn(seq_along(cols_b), clin_k, simplify = FALSE) # generate all clinical combinations of size clin_k
              combs_clin = combs_clin[order(sapply(combs_clin, function(i) min(months_b[i])))] # sort clinical combinations by earliest month (ascending)
              
              for (ci in combs_clin) { # loop over clinical combinations
                first_clin = min(months_b[ci]) # earliest clinical month in this combo
                inter = intersect(cols_n[fi], cols_b[ci]) # check overlapping timepoints
                if (last_fc <= first_clin) { # ensure FC ends before or at clinical start
                  if (length(inter) == 0 || (length(inter) == 1 && last_fc == first_clin)) { # allow no overlap, or one overlap if boundary
                    chosen_fc = fi # store FC combo
                    chosen_clin = ci # store clinical combo
                    success = TRUE # mark success
                    break # exit clinical loop
                  }
                }
              }
              if (success) break # exit clinical loop if success achieved
            }
            if (success) break # exit FC loop if success achieved
          }
          if (success) break # exit loop if success achieved
        }
      }
      
      # CASE 4: Only runs if cases 1 to 3 were not successful -> takes two time points for each variable and checks temporal order (FC prior to clinical slope)
      if (!success & available_n >= 2 & available_b >= 2) { # only run if all previous steps failed and both FC and clinical have at least 2 points
        chosen_fc = seq_len(min(2, available_n)) # pick the first 2 FC points (or fewer if less available)
        chosen_clin = seq_len(min(2, available_b)) # pick the first 2 clinical points (or fewer if less available)
        if (max(months_n[chosen_fc]) > min(months_b[chosen_clin])) { # check temporal order: FC must not extend beyond first clinical
          chosen_fc = NULL # reset if invalid
          chosen_clin = NULL # reset if invalid
        }
      }
      
      # extract values of interest and fit slopes among other important parameters
      if (!is.null(chosen_fc) & !is.null(chosen_clin)) { # proceed only if valid FC and clinical selections exist
        # extract data
        fc_cols = cols_n[chosen_fc] # selected FC column names
        fc_months = months_n[chosen_fc] # selected FC months
        fc_vals = vals_n_num[chosen_fc] # selected FC numeric values
        clin_cols = cols_b[chosen_clin] # selected clinical column names
        clin_months = months_b[chosen_clin] # selected clinical months
        clin_vals = vals_b_num[chosen_clin] # selected clinical numeric values
        ord_fc = order(fc_months) # sort FC by month
        fc_cols = fc_cols[ord_fc]; fc_months = fc_months[ord_fc]; fc_vals = fc_vals[ord_fc] # reorder FC variables consistently
        ord_clin = order(clin_months) # sort clinical by month
        clin_cols = clin_cols[ord_clin]; clin_months = clin_months[ord_clin]; clin_vals = clin_vals[ord_clin] # reorder clinical variables consistently
        # compute the single subject regressions for early FC period
        if (length(fc_vals) >= 2) { # require at least 2 FC points
          fit_fc = lm(fc_vals ~ fc_months) # fit linear model of FC values over months
          slope_n[i] = coef(summary(fit_fc))["fc_months","Estimate"] # store slope estimate
          slope_n_se[i] = coef(summary(fit_fc))["fc_months","Std. Error"] # store slope SE
          slope_n_nt[i] = length(fc_vals) # number of FC points
          slope_n_var_m[i] = var(fc_months) # variance of months (predictor)
          slope_n_var_v[i] = var(fc_vals) # variance of values (outcome)
          earliest_n_month[i] = min(fc_months) # earliest FC month
          earliest_n_tp[i] = fc_cols[which.min(fc_months)] # column corresponding to earliest FC time point
          latest_n_month[i] = max(fc_months) # latest FC month
          latest_n_tp[i] = fc_cols[which.max(fc_months)] # column corresponding to latest FC time point
        }
        # compute the single subject regressions for late clinical period
        if (length(clin_vals) >= 2) { # require at least 2 clinical points
          fit_clin = lm(clin_vals ~ clin_months) # fit linear model of clinical values over months
          slope_b[i] = coef(summary(fit_clin))["clin_months","Estimate"] # store slope estimate
          slope_b_se[i] = coef(summary(fit_clin))["clin_months","Std. Error"] # store slope SE
          slope_b_nt[i] = length(clin_vals) # number of clinical points
          slope_b_var_m[i] = var(clin_months) # variance of months (predictor)
          slope_b_var_v[i] = var(clin_vals) # variance of values (outcome)
          earliest_b_month[i] = min(clin_months) # earliest clinical month
          earliest_b_tp[i] = clin_cols[which.min(clin_months)] # column corresponding to earliest clinical timepoint
          latest_b_month[i] = max(clin_months) # latest clinical month
          latest_b_tp[i] = clin_cols[which.max(clin_months)] # column corresponding to latest clinical timepoint
        }
      }
    }
    # add FC results to data frame
    df[[paste0(nvar, "_slope_early")]] = slope_n # store FC slope estimates
    df[[paste0(nvar, "_slope_early_se")]] = slope_n_se # store standard errors of FC slopes
    df[[paste0(nvar, "_slope_early_nt")]] = slope_n_nt # store number of FC points used
    df[[paste0(nvar, "_slope_early_var_month")]] = slope_n_var_m # store variance of FC months
    df[[paste0(nvar, "_slope_early_var_value")]] = slope_n_var_v # store variance of FC values
    df[[paste0(nvar, "_earliest_month")]] = earliest_n_month # store earliest FC month used
    df[[paste0(nvar, "_earliest_tp")]] = earliest_n_tp # store column name of earliest FC timepoint
    df[[paste0(nvar, "_latest_month")]] = latest_n_month # store latest FC month used
    df[[paste0(nvar, "_latest_tp")]] = latest_n_tp # store column name of latest FC timepoint
    # add clinical results to data frame
    df[[paste0(bvar, "_slope_late")]] = slope_b # store clinical slope estimates
    df[[paste0(bvar, "_slope_late_se")]] = slope_b_se # store standard errors of clinical slopes
    df[[paste0(bvar, "_slope_late_nt")]] = slope_b_nt # store number of clinical points used
    df[[paste0(bvar, "_slope_late_var_month")]] = slope_b_var_m # store variance of clinical months
    df[[paste0(bvar, "_slope_late_var_value")]] = slope_b_var_v # store variance of clinical values
    df[[paste0(bvar, "_earliest_month")]] = earliest_b_month # store earliest clinical month used
    df[[paste0(bvar, "_earliest_tp")]] = earliest_b_tp # store column name of earliest clinical timepoint
    df[[paste0(bvar, "_latest_month")]] = latest_b_month # store latest clinical month used
    df[[paste0(bvar, "_latest_tp")]] = latest_b_tp # store column name of latest clinical timepoint
    # replace NaN with NA
    df[df == "NaN"] = NA
    # remove subjects without slope in FC and clinical variable
    df = df[!is.na(df[paste0(nvar, "_slope_early")]) & !is.na(df[paste0(bvar, "_slope_late")]),]
    # remove subjects for which the latest FC time point and the earliest clinical timepoint do not match
    df = df[df[[paste0(nvar, "_latest_tp")]] == df[[paste0(bvar, "_earliest_tp")]], ]
    # remove subjects that do not have a baseline as earliest timepoint in FC
    df = subset(df,fchi_earliest_tp == "BL")
    # export this data
    write.csv(df, file = paste0(out_dir,"Napls_EarlySlopeFC_LateSlopeClin_", bvar, ".csv"), row.names = FALSE)
    
    
    ##### Imputation of standard errors #####
    # select variables of interest
    df = df[c("subj","CHR","age","female","site02","site03","site04","site05","site06","site07","site08","site09",
              paste0(nvar,"_slope_early"),paste0(nvar, "_slope_early_nt"),paste0(nvar, "_slope_early_se"),
              paste0(nvar, "_slope_early_var_month"),paste0(nvar, "_slope_early_var_value"),
              paste0(bvar,"_slope_late"),paste0(bvar, "_slope_late_nt"),paste0(bvar, "_slope_late_se"),
              paste0(bvar, "_slope_late_var_month"),paste0(bvar, "_slope_late_var_value"))]
    # create formula for clinical slope
    formula_bvar = paste(paste0(bvar, "_slope_late_se"),"~",
                         paste0(bvar, "_slope_late"),"+",
                         paste0(bvar, "_slope_late_nt"),"+",
                         paste0(bvar, "_slope_late_var_month"),"+",
                         paste0(bvar, "_slope_late_var_value"))
    # run regression model
    res = lm(as.formula(formula_bvar), data = df)
    # extract model-level statistics
    model_summary = summary(res)
    r2 = model_summary$r.squared # r2
    adj_r2 = model_summary$adj.r.squared # adjusted r2
    fstat = model_summary$fstatistic # F stats
    model_p = pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE) # model p-value
    # count missing SEs and compute imputation stats
    outcome_se_var = paste0(bvar, "_slope_late_se")
    outcome_nt_var = paste0(bvar, "_slope_late_nt")
    rows_to_impute = which(is.na(df[[outcome_se_var]]) & df[[outcome_nt_var]] == 2)
    n_imputed = length(rows_to_impute) # number of imputed SEs
    n_total = sum(df[[outcome_nt_var]] > 1, na.rm = TRUE) # total number of subjects with more than 2 time points
    prop_imputed = ifelse(n_total > 0, n_imputed / n_total, NA) # proportion
    # extract rows to impute
    df_test = df[rows_to_impute, ]
    # predict SEs
    se_pred = predict(res, newdata = df_test)
    # count negative predicted SEs
    n_negative = sum(se_pred < 0)
    prop_negative = ifelse(n_imputed > 0, n_negative / n_imputed, NA)
    # impute values
    df[rows_to_impute, outcome_se_var] = se_pred
    # save model stats as data frame
    model_stats_df = data.frame(Outcome = bvar,
                                R2 = r2,
                                Adjusted_R2 = adj_r2,
                                Model_p_value = model_p,
                                N_total_with_nt_larger2 = n_total,
                                N_imputed = n_imputed,
                                Prop_imputed = prop_imputed,
                                N_negative_predicted = n_negative,
                                Prop_negative_predicted = prop_negative)
    # export model stats
    write.csv(model_stats_df, file = paste0(out_dir,"ImputationSE_SyChangePrediction_ClinicalVariable", bvar, ".csv"), row.names = FALSE)
    
    # create formula for FC slope
    formula_nvar = paste(paste0(nvar, "_slope_early_se"),"~",
                         paste0(nvar, "_slope_early"),"+",
                         paste0(nvar, "_slope_early_nt"),"+",
                         paste0(nvar, "_slope_early_var_month"),"+",
                         paste0(nvar, "_slope_early_var_value"))
    # run regression model
    res = lm(as.formula(formula_nvar), data = df)
    # extract model-level statistics
    model_summary = summary(res)
    r2 = model_summary$r.squared # r2
    adj_r2 = model_summary$adj.r.squared # adjusted r2
    fstat = model_summary$fstatistic # F stats
    model_p = pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE) # model p-value
    # count missing SEs and compute imputation stats
    outcome_se_var = paste0(nvar, "_slope_early_se")
    outcome_nt_var = paste0(nvar, "_slope_early_nt")
    rows_to_impute = which(is.na(df[[outcome_se_var]]) & df[[outcome_nt_var]] == 2)
    n_imputed = length(rows_to_impute) # number of imputed SEs
    n_total = sum(df[[outcome_nt_var]] > 1, na.rm = TRUE) # total number of subjects with more than 2 time points
    prop_imputed = ifelse(n_total > 0, n_imputed / n_total, NA) # proportion
    # extract rows to impute
    df_test = df[rows_to_impute, ]
    # predict SEs
    se_pred = predict(res, newdata = df_test)
    # count negative predicted SEs
    n_negative = sum(se_pred < 0)
    prop_negative = ifelse(n_imputed > 0, n_negative / n_imputed, NA)
    # impute values
    df[rows_to_impute, outcome_se_var] = se_pred
    # save model stats as data frame
    model_stats_df = data.frame(Outcome = nvar,
                                R2 = r2,
                                Adjusted_R2 = adj_r2,
                                Model_p_value = model_p,
                                N_total_with_nt_larger2 = n_total,
                                N_imputed = n_imputed,
                                Prop_imputed = prop_imputed,
                                N_negative_predicted = n_negative,
                                Prop_negative_predicted = prop_negative)
    # export model stats
    write.csv(model_stats_df, file = paste0(out_dir,"ImputationSE_SyChangePrediction_FCVariable",nvar,"_",bvar, ".csv"), row.names = FALSE)
    
    
    ##### Latent variable regression #####
    # re-specify working directory which must not be too long for MPLUS
    wd = "/home/lukas/mplus"
    setwd(wd)
    # print process
    print(paste0("Running Single Indicator Approach for ",bvar, " and ",nvar,"..."))
    # shorten variable names of MRI variables for MPLUS
    nvar_mplus = nvar %>% str_replace("fc","")
    # create data frame with variables of interest
    df = df[c(paste0(bvar,"_slope_late"),paste0(nvar,"_slope_early"),paste0(nvar,"_slope_early_se"),paste0(bvar,"_slope_late_se"),
              "age","female","site02","site03","site04","site05","site06","site07","site08","site09")]
    # get average error variances across group from standard errors of both variables
    se_nvar = mean(df[[paste0(nvar,"_slope_early_se")]],na.rm=T)^2
    se_bvar = mean(df[[paste0(bvar,"_slope_late_se")]],na.rm=T)^2
    # drop standard error columns again
    df = df[ , !grepl("_slope_late_se$", names(df)) ]
    df = df[ , !grepl("_slope_early_se$", names(df)) ]
    # replace NA with -99
    df[is.na(df)] = -99
    
    ###### Two latent change variables #####
    # define the input and output files
    datafile = file.path(wd, paste0("NaplsData_", nvar, "_", bvar, "_SyChangePrediction_2LV.dat"))
    inpfile = file.path(wd, paste0("NaplsData_", nvar, "_", bvar, "_SyChangePrediction_2LV.inp"))
    # write output file
    write.table(df,datafile,quote = F,row.names = F, col.names = F,append = F, sep = "\t")
    # define model in MPLUS
    inp_text = paste0("
    Title:      LATENT VARIABLE REGRESSION

    Data:       File is ", datafile, ";
    
    Variable:
    names are
    ", bvar, "_CH
    ", nvar_mplus, "_CH
    age
    female
    site02
    site03
    site04
    site05
    site06
    site07
    site08
    site09;
    
    usevariables are
    ", bvar, "_CH
    ", nvar_mplus, "_CH
    age
    female
    site02
    site03
    site04
    site05
    site06
    site07
    site08
    site09;
    
    missing are all (-99);
  
    Analysis:   type = random;
                ESTIMATOR=MLR;
                ALGORITHM=INTEGRATION;
                STARTS=20;
                STITERATION=50000;
    
    Model:
    
    LV_", nvar_mplus, " by ", nvar_mplus, "_CH;
    ", nvar_mplus, "_CH@", format(se_nvar, scientific = FALSE), ";
    
    LV_", bvar, " by ", bvar, "_CH;
    ", bvar, "_CH@", format(se_bvar, scientific = FALSE), ";
    
    
    LV_", bvar, " on LV_", nvar_mplus, " age 
              female site02 site03 site04
               site05 site06 site07 site08
               site09;
    
    LV_", nvar_mplus, " WITH age female site02 site03 site04
                       site05 site06 site07 site08
                       site09;
    
    age WITH female site02 site03 site04 site05
                    site06 site07 site08
                    site09;
    
    female WITH site02 site03 site04 site05
                    site06 site07 site08
                    site09;
    
    site02 WITH site03 site04 site05
                    site06 site07 site08
                    site09;
    
    site03 WITH site04 site05 site06 site07 site08
                    site09;
    
    site04 WITH site05 site06 site07 site08
                    site09;
    
    site05 WITH site06 site07 site08
                    site09;
    
    site06 WITH site07 site08 site09;
    
    site07 WITH site08 site09;
    
    site08 WITH site09;
    
    Output:     tech1 tech4 stdyx sampstat;
    "
    )
    # write as input file for MPLUS
    writeLines(inp_text, con = inpfile)
    # run MPLUS
    system(paste("wine Mplus-8.6.exe", basename(inpfile)))
    
    ##### Summarize results from version with two latent change variables (others are kept in MPLUS format) #####
    # specify the output file
    outfile = file.path(wd, paste0("naplsdata_", tolower(nvar), "_", tolower(bvar), "_sychangeprediction_2lv.out"))
    # read the output
    model_out = readModels(outfile)
    # extract estimates
    df_res = as.data.frame(model_out$parameters$unstandardized)
    # compute z value and exact p values where applicable
    df_res$zval = NA
    df_res$pval_exact = NA
    valid = df_res$se > 0 & df_res$pval != 999
    df_res$zval[valid] = with(df_res[valid, ], est / se)
    df_res$pval_exact[valid] = with(df_res[valid, ], 2 * (1 - pnorm(abs(zval))))
    # mark significance
    df_res$p_sig = df_res$pval_exact < 0.05
    # export model stats
    write.csv(df_res, file = paste0("SingleIndicatorSyChangePredictions_PathsUnstandardized_",bvar,"_",nvar, "_2LV.csv"), row.names = FALSE)
    # extract standardized estimates
    df_res = as.data.frame(model_out$parameters$stdyx.standardized)
    # compute z value and exact p values where applicable
    df_res$zval = NA
    df_res$pval_exact = NA
    valid = df_res$se > 0 & df_res$pval != 999
    df_res$zval[valid] = with(df_res[valid, ], est / se)
    df_res$pval_exact[valid] = with(df_res[valid, ], 2 * (1 - pnorm(abs(zval))))
    # mark significance
    df_res$p_sig = df_res$pval_exact < 0.05
    # export model stats
    write.csv(df_res, file = paste0("SingleIndicatorSyChangePredictions_PathsStandardized_",bvar,"_",nvar, "_2LV.csv"), row.names = FALSE)
    # save model summary
    fit = model_out$summaries
    # select relevant meta and fit indices 
    df_model = data.frame(
      estimator       = if (is.null(fit$Estimator)) NA else fit$Estimator,
      n_obs           = if (is.null(fit$Observations)) NA else fit$Observations,
      n_parameters    = if (is.null(fit$Parameters)) NA else fit$Parameters,
      aic             = if (is.null(fit$AIC)) NA else fit$AIC,
      bic             = if (is.null(fit$BIC)) NA else fit$BIC,
      chi_sq          = if (is.null(fit$ChiSqM_Value)) NA else fit$ChiSqM_Value,
      chi_df          = if (is.null(fit$ChiSqM_DF)) NA else fit$ChiSqM_DF,
      chi_pval        = if (is.null(fit$ChiSqM_PValue)) NA else fit$ChiSqM_PValue,
      rmsea           = if (is.null(fit$RMSEA_Estimate)) NA else fit$RMSEA_Estimate,
      rmsea_ci_lower  = if (is.null(fit$RMSEA_90CI_LB)) NA else fit$RMSEA_90CI_LB,
      rmsea_ci_upper  = if (is.null(fit$RMSEA_90CI_UB)) NA else fit$RMSEA_90CI_UB,
      cfi             = if (is.null(fit$CFI)) NA else fit$CFI,
      tli             = if (is.null(fit$TLI)) NA else fit$TLI,
      srmr            = if (is.null(fit$SRMR)) NA else fit$SRMR,
      n_miss_patterns = if (!is.null(model_out$missing$univariateMissingDescriptives$nPatterns)) {
        model_out$missing$univariateMissingDescriptives$nPatterns
      } else { NA })
    # export stats
    write.csv(df_model, file = paste0("SingleIndicatorSyChangePredictions_Model_", bvar, "_", nvar, "_2LV.csv"), row.names = FALSE)
    # define and set old working directory
    old_wd = "/home/lukas/Desktop/LukasLinux/Projects/LongitMechanisms/"
    setwd(old_wd)
    # create sub-directory in output folder
    dir.create(paste0(out_dir,"MplusOutputs"), recursive = TRUE,showWarnings=F)
    # specify paths to output files
    res_files = list.files(path = wd,pattern = "\\.(csv|dat|inp|out)$",recursive = TRUE, full.names = TRUE)
    # copy to output sub-directory in output directory
    copied = file.copy(res_files,paste0(out_dir,"MplusOutputs"), overwrite = TRUE)
    # remove copied files from intermediate working directory
    file.remove(res_files[copied])
  }
}
  

##### Unstandardized results: Correction for multiple comparisons #####
# loop to extract the results for each MRI variable, perform a p-value correction across the six tests, and merge to one data frame 
df_si_list = list()
for (nvar in neurovars) {
  # shorten variable names of MRI variables for MPLUS and convert to upper case
  nvar_mplus = nvar %>% str_replace('fc','') %>% toupper()
  # specify paths to csv files with path statistics
  tmp_paths = list.files(path = paste0(out_dir,"MplusOutputs"),
                         pattern = sprintf("^SingleIndicatorSyChangePredictions_PathsUnstandardized_.*_%s\\_2LV.csv$", nvar),
                         full.names = TRUE)
  # create empty list
  df_list = list()
  # import CSV files in loop
  for (path in tmp_paths) {
    df = read.csv(path)
    df_filtered = subset(df,param==paste0("LV_",nvar_mplus) & paramHeader != 'Variances' & paramHeader != 'Means')
    df_list[[length(df_list) + 1]] = df_filtered
  }
  # create one data frame
  df_si_stats = do.call(rbind, df_list)
  # insert adjusted p-values and mark significance using different methods
  df_si_stats$p_fdr = p.adjust(df_si_stats$pval_exact,method="fdr")
  df_si_stats$p_fdr_sig = df_si_stats$p_fdr < 0.05
  df_si_stats$p_bonf = p.adjust(df_si_stats$pval_exact,method="bonferroni")
  df_si_stats$p_bonf_sig = df_si_stats$p_bonf < 0.05
  # add data frame to list
  df_si_list[[length(df_si_list) + 1]] = df_si_stats
}

##### Unstandardized results: Summary of main effect results #####
# create one data frame
df_si_stats_full = as.data.frame(do.call(rbind, df_si_list))
# export data
write.csv(df_si_stats_full,paste0(out_dir,"Napls_LatentVariableSyChangePrediction_UnstandardizedResultsSummary_MainEffBonf_FCHip_2LV.csv"), row.names = FALSE)
write_xlsx(df_si_stats_full,path = paste0(out_dir,"Napls_LatentVariableSyChangePrediction_UnstandardizedResultsSummary_MainEffBonf_FCHip_2LV.xlsx"))

##### Unstandardized results: Summary of full results #####
# list all csv containing full test stats of hippocampal FC
path_files_hip = list.files(path = paste0(out_dir,"MplusOutputs"),
                          pattern = "^SingleIndicatorSyChangePredictions_PathsUnstandardized_.*fchi\\_2LV.csv$",
                          full.names = TRUE)
# read them all and bind rowwise
df_path_hip = path_files_hip %>% lapply(read.csv) %>% bind_rows() %>% as.data.frame()
# replace 999 with NA
df_path_hip[df_path_hip == 999] = NA
# export data
write.csv(df_path_hip,file = paste0(out_dir,"Napls_LatentVariableSyChangePrediction_UnstandardizedResultsSummary_FullStats_FCHip_2LV.csv"), row.names = FALSE)
write_xlsx(df_path_hip,path = paste0(out_dir,"Napls_LatentVariableSyChangePrediction_UnstandardizedResultsSummary_FullStats_FCHip_2LV.xlsx"))

##### Standardized results: Correction for multiple comparisons #####
# loop to extract the results for each MRI variable, perform a p-value correction across the six tests, and merge to one data frame 
df_si_list = list()
for (nvar in neurovars) {
  # shorten variable names of MRI variables for MPLUS and convert to upper case
  nvar_mplus = nvar %>% str_replace('fc','') %>% toupper()
  # specify paths to csv files with path statistics
  tmp_paths = list.files(path = paste0(out_dir,"MplusOutputs"),
                         pattern = sprintf("^SingleIndicatorSyChangePredictions_PathsStandardized_.*_%s\\_2LV.csv$", nvar),
                         full.names = TRUE)
  # create empty list
  df_list = list()
  # import CSV files in loop
  for (path in tmp_paths) {
    df = read.csv(path)
    df_filtered = subset(df,param==paste0("LV_",nvar_mplus) & paramHeader != 'Variances' & paramHeader != 'Means')
    df_list[[length(df_list) + 1]] = df_filtered
  }
  # create one data frame
  df_si_stats = do.call(rbind, df_list)
  # insert adjusted p-values and mark significance using different methods
  df_si_stats$p_fdr = p.adjust(df_si_stats$pval_exact,method="fdr")
  df_si_stats$p_fdr_sig = df_si_stats$p_fdr < 0.05
  df_si_stats$p_bonf = p.adjust(df_si_stats$pval_exact,method="bonferroni")
  df_si_stats$p_bonf_sig = df_si_stats$p_bonf < 0.05
  # add data frame to list
  df_si_list[[length(df_si_list) + 1]] = df_si_stats
}

##### Standardized results: Summary of main effect results #####
# create one data frame
df_si_stats_full = data.frame(do.call(rbind, df_si_list))
# export data
write.csv(df_si_stats_full,paste0(out_dir,"Napls_LatentVariableSyChangePrediction_StandardizedResultsSummary_MainEffBonf_FCHip_2LV.csv"), row.names = FALSE)
write_xlsx(df_si_stats_full,path = paste0(out_dir,"Napls_LatentVariableSyChangePrediction_StandardizedResultsSummary_MainEffBonf_FCHip_2LV.xlsx"))

##### Standardized results: Summary of full results #####
# list all csv containing full test stats of hippocampal FC
path_files_hip = list.files(path = paste0(out_dir,"MplusOutputs"),
                            pattern = "^SingleIndicatorSyChangePredictions_PathsStandardized_.*fchi\\_2LV.csv$",
                            full.names = TRUE)
# read them all and bind rowwise
df_path_hip = path_files_hip %>% lapply(read.csv) %>% bind_rows() %>% as.data.frame()
# replace 999 with NA
df_path_hip[df_path_hip == 999] = NA
# export data
write.csv(df_path_hip,file = paste0(out_dir,"Napls_LatentVariableSyChangePrediction_StandardizedResultsSummary_FullStats_FCHip_2LV.csv"), row.names = FALSE)
write_xlsx(df_path_hip,path = paste0(out_dir,"Napls_LatentVariableSyChangePrediction_StandardizedResultsSummary_FullStats_FCHip_2LV.xlsx"))


##### Summary of model fit results #####
# run loop to extract their full test statistics on model fit results
df_list = list()
for (bvar in behavvars) {
  for (nvar in neurovars) {
    p = list.files(path = paste0(out_dir,"MplusOutputs"),
                   pattern = sprintf("^SingleIndicatorSyChangePredictions_Model_%s_%s\\_2LV.csv$",bvar,nvar),
                   full.names = TRUE)
    df = read.csv(p)
    df = df %>% add_column(behav=bvar,neuro=nvar,.before=1)
    df_list[[paste0(bvar,"+",nvar)]] = df
  }
}
# create one data frame
df_modelfit = as.data.frame(do.call(rbind, df_list))
# reset rownames
rownames(df_modelfit) = NULL
# export data
write.csv(df_modelfit,file = paste0(out_dir,"Napls_LatentVariableSyChangePrediction_ResultsSummary_ModelFit_FCHip_2LV.csv"), row.names = FALSE)
write_xlsx(df_modelfit, path = paste0(out_dir,"Napls_LatentVariableSyChangePrediction_ResultsSummary_ModelFit_FCHip_2LV.xlsx"))

