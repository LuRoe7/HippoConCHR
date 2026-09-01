# FC AND SYMPTOM CHANGE PREDICTIONS  ---------------------------------------------
# Description:  This script runs a prediction of symptom course based on FC Change and vice versa using regression models.
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
old_wd = "/home/lukas/Desktop/LukasLinux/Projects/LongitMechanisms/"
setwd(old_wd)
getwd()
# define input directory for behavioral and FC data
in_dir = "results/napls/P1_HippoFC/07_NAPLS_SingleSubjectRegressions/"
# define output directory
out_dir = "results/napls/P1_HippoFC/12a_NAPLS_FCandSymptomPredictions_EarlyChange_OnlyPatients/"
# create output directory for this script
dir.create(out_dir, recursive = TRUE,showWarnings=F)
# read data including behavioral and FC data and all slopes with imputed SEs (note that those are not used in this analysis)
df_napls_mlr_imp = read.csv(paste0(in_dir,"Napls_BehavFCandSlopes_impSE.csv"))
# keep only CHR individuals
df_napls_mlr_imp = subset(df_napls_mlr_imp,group == 1)
# rename group column
df_napls_mlr_imp = df_napls_mlr_imp %>% rename(CHR = group)
# create vectors with behavioral variables that were significant in main analysis
behavvars = c("NegS","DepS","gaf")
# create vector with neuro variables
neurovars = c("fchi")


# LATENT CHANGE PREDICTION: PREDICTING SYMPTOM CHANGE BASED ON PRIOR FC CHANGE  ---------------------
# Note: This code identifies the first two available measurement time points for FC and computes the slope between those.
# Then, it computes the subsequent slope of the clinical variables, performs the
# imputation of standard errors for subjects with only two measurement time points for the clinical variable. Lastly, it runs
# the latent variable regression with FC slope as manifest change variable and clinical slopes as latent change variables.
# Additionally, as a control analysis, it runs the regression also with a manifest change variable for clinical variables.
# run slope estimations, imputations and latent variable regression in one loop
for (bvar in behavvars) {
  for (nvar in neurovars) {
    ##### Slope estimations based on single subject regressions #####
    # select variables of interest
    df = df_napls_mlr_imp[c("subj","CHR","age","female","site02","site03","site04","site05","site06","site07","site08","site09",
                            paste0("nt_",bvar),paste0("nt_",nvar),
                            paste0(bvar,"_BL"),paste0(bvar,"_M2"),paste0(bvar,"_M4"),paste0(bvar,"_M6"),paste0(bvar,"_M8"),
                            paste0(nvar,"_BL"),paste0(nvar,"_M2"),paste0(nvar,"_M4"),paste0(nvar,"_M6"),paste0(nvar,"_M8"))]
    # select only subjects with at least 2 measurement timepoints in FC and clinical variable and baseline value in FC
    df = df[!is.na(df[[paste0(nvar, "_BL")]]) & df[[paste0("nt_",bvar)]] > 1 & df[[paste0("nt_",nvar)]] > 1, ]
    # reset rownames
    rownames(df) = NULL
    # define follow-up timepoints and the corresponding months
    timepoints = c("M2", "M4", "M6", "M8")
    months = c(2, 4, 6, 8)
    # compute slope per month from baseline for each follow-up timepoint
    baseline_col = paste0(nvar, "_BL") # save baseline column
    for (i in seq_along(timepoints)) {
      tp = timepoints[i] # define timepoint
      month = months[i] # define month
      followup_col = paste0(nvar, "_", tp) # define follow-up column
      slope_col = paste0(nvar, "_slope_", tp)  # new column for slope
      df[[slope_col]] = (df[[followup_col]] - df[[baseline_col]]) / month # compute slope
    }
    # initialize new columns
    df[[paste0(nvar, "_slope_first")]]      = NA
    df[[paste0(nvar, "_slope_first_end")]]  = NA
    df[[paste0(bvar, "_slope_post")]]       = NA
    df[[paste0(bvar, "_slope_post_se")]]    = NA
    df[[paste0(bvar, "_slope_post_var_month")]] = NA
    df[[paste0(bvar, "_slope_post_var_value")]]        = NA
    df[[paste0(bvar, "_slope_post_nt")]]          = NA
    df[[paste0(bvar, "_slope_post_start")]] = NA
    df[[paste0(bvar, "_slope_post_end")]]   = NA
    
    for (i in 1:nrow(df)) {
      # extract slopes for FC variable
      slopes = sapply(timepoints, function(tp) df[[paste0(nvar, "_slope_", tp)]][i])
      # get earliest non-missing slope for FC variable
      first_slope_idx = which(!is.na(slopes))[1]
      if (!is.na(first_slope_idx)) {
        df[[paste0(nvar, "_slope_first")]][i]     = slopes[first_slope_idx]
        slope_month                               = months[first_slope_idx]
        df[[paste0(nvar, "_slope_first_end")]][i] = slope_month
        
        # get values of clinical variable starting at or after slope_month
        bvar_values = sapply(timepoints, function(tp) df[[paste0(bvar, "_", tp)]][i])
        valid_idx   = which(!is.na(bvar_values) & months >= slope_month)
        
        if (length(valid_idx) >= 2 && months[valid_idx[1]] == slope_month) {
          # extract timepoints and values of clinical variable
          t_vals = months[valid_idx]
          v_vals = bvar_values[valid_idx]
          
          # run slope calculation for clinical variable
          fit = lm(v_vals ~ t_vals)
          coefs     = coef(summary(fit))
          est_slope = coefs["t_vals", "Estimate"]
          se_slope  = coefs["t_vals", "Std. Error"]
          
          var_month = var(t_vals)
          var_value = var(v_vals)
          
          df[[paste0(bvar, "_slope_post")]][i] = est_slope
          df[[paste0(bvar, "_slope_post_se")]][i] = se_slope
          df[[paste0(bvar, "_slope_post_var_month")]][i] = var_month
          df[[paste0(bvar, "_slope_post_var_value")]][i] = var_value
          df[[paste0(bvar, "_slope_post_nt")]][i] = length(valid_idx)
          df[[paste0(bvar, "_slope_post_start")]][i] = min(t_vals)
          df[[paste0(bvar, "_slope_post_end")]][i]   = max(t_vals)
        }
      }
    }
    # replace NaN with NA
    df[df == "NaN"] = NA
    # remove subjects without slope in clinical variable
    df = df[!is.na(df[paste0(bvar, "_slope_post")]),]
    # export data
    write.csv(df, file = paste0(out_dir,"Napls_EarlySlopeFC_LateSlopeClin_",nvar,"_",bvar, ".csv"), row.names = FALSE)
    
    
    ##### Imputation of standard errors #####
    # select variables of interest
    df = df[c("subj","CHR","age","female","site02","site03","site04","site05","site06","site07","site08","site09",
              paste0(nvar,"_slope_first"),paste0(bvar,"_slope_post"),
              paste0(bvar, "_slope_post_se"),paste0(bvar, "_slope_post_var_month"),paste0(bvar, "_slope_post_var_value"),
              paste0(bvar, "_slope_post_nt"))]
    # create formula
    formula_str = paste(paste0(bvar, "_slope_post_se"),"~",
                        paste0(bvar, "_slope_post"),"+",
                        paste0(bvar, "_slope_post_nt"),"+",
                        paste0(bvar, "_slope_post_var_month"),"+",
                        paste0(bvar, "_slope_post_var_value"))
    # run regression model
    res = lm(as.formula(formula_str), data = df)
    # extract model-level statistics
    model_summary = summary(res)
    r2 = model_summary$r.squared # r2
    adj_r2 = model_summary$adj.r.squared # adjusted r2
    fstat = model_summary$fstatistic # F stats
    model_p = pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE) # model p-value
    # count missing SEs and compute imputation stats
    outcome_se_var = paste0(bvar, "_slope_post_se")
    outcome_nt_var = paste0(bvar, "_slope_post_nt")
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
    write.csv(model_stats_df, file = paste0(out_dir,"ImputationSE_model_SyChangePred_",nvar,"_",bvar, ".csv"), row.names = FALSE)
    
    ##### Latent variable regression - manifest FC slope and latent clinical slope #####
    # re-specify working directory which must not be too long for MPLUS
    wd = "/home/lukas/mplus"
    setwd(wd)
    # print process
    print(paste0("Running Single Indicator Approach for ",bvar, " and ",nvar,"..."))
    # shorten variable names of MRI variables for MPLUS
    nvar_mplus = nvar %>% str_replace("fc","")
    # create data frame with variables of interest
    df = df[c(paste0(bvar,"_slope_post"),paste0(nvar,"_slope_first"),paste0(bvar,"_slope_post_se"),
              "age","female","site02","site03","site04","site05","site06","site07","site08","site09")]
    # replace NA with -99
    df[is.na(df)] = -99
    # get average error variances across group from standard errors of both variables
    se_bvar = mean(df[[paste0(bvar,"_slope_post_se")]],na.rm=T)^2
    # drop standard error columns again
    df = df[ , !grepl("_slope_post_se$", names(df)) ]
    # define the input and output files
    datafile = file.path(wd, paste0("NaplsData_", nvar, "_", bvar, "_SyChangePred.dat"))
    inpfile = file.path(wd, paste0("NaplsData_", nvar, "_", bvar, "_SyChangePred.inp"))
    # write output file
    write.table(df,datafile,quote = F,row.names = F, col.names = F,append = F, sep = "\t")
    # define the MPLUS input file containing the model
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
                ESTIMATOR = MLR;
                ALGORITHM=INTEGRATION;
                STARTS=20;
                STITERATION=50000;
    
    Model:
    
    LV_", bvar, " by ", bvar, "_CH;
    ", bvar, "_CH@", format(se_bvar, scientific = FALSE), ";
    
    LV_", bvar, " on ", nvar_mplus, "_CH
               age female site02 site03 site04
               site05 site06 site07 site08
               site09;
    
    ",
            
    nvar_mplus, "_CH WITH age female site02 site03 site04
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
    
    SAVEDATA: ESTIMATES = power_",nvar,"_",bvar,"_sychangepred.dat;
    Output:     tech1 tech4 stdyx sampstat;
    "
    )
    # write as input file for MPLUS
    writeLines(inp_text, con = inpfile)
    # run MPLUS
    system(paste("wine Mplus-8.6.exe", basename(inpfile)))
    # specify the output file
    outfile = file.path(wd, paste0("naplsdata_", tolower(nvar), "_", tolower(bvar), "_sychangepred.out"))
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
    write.csv(df_res, file = paste0("SingleIndicatorSyChangePred_PathsUnstandardized_",bvar,"_",nvar, ".csv"), row.names = FALSE)
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
    write.csv(df_res, file = paste0("SingleIndicatorSyChangePred_PathsStandardized_",bvar,"_",nvar, ".csv"), row.names = FALSE)
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
    write.csv(df_model, file = paste0("SingleIndicatorSyChangePred_Model_", bvar, "_", nvar, ".csv"), row.names = FALSE)
    
    ##### Latent variable regression (control analysis) - manifest FC slope and manifest clinical slope #####
    # define the input and output files
    datafile = file.path(wd, paste0("NaplsData_", nvar, "_", bvar, "_SyChangePredCo.dat"))
    inpfile = file.path(wd, paste0("NaplsData_", nvar, "_", bvar, "_SyChangePredCo.inp"))
    # write output file
    write.table(df,datafile,quote = F,row.names = F, col.names = F,append = F, sep = "\t")
    # define the MPLUS input file containing the model
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
  
    Analysis:   ESTIMATOR = MLR;
    
    Model:
    
    ",
    
    bvar, "_CH on ", nvar_mplus, "_CH
               age female site02 site03 site04
               site05 site06 site07 site08
               site09;
    
    ",
    
    nvar_mplus, "_CH WITH age female site02 site03 site04
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
    # specify the output file
    outfile = file.path(wd, paste0("naplsdata_", tolower(nvar), "_", tolower(bvar), "_sychangepredco.out"))
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
    write.csv(df_res, file = paste0("SingleIndicatorSyChangePredCont_PathsUnstandardized_",bvar,"_",nvar, ".csv"), row.names = FALSE)
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
    write.csv(df_res, file = paste0("SingleIndicatorSyChangePredCont_PathsStandardized_",bvar,"_",nvar, ".csv"), row.names = FALSE)
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
    write.csv(df_model, file = paste0("SingleIndicatorSyChangePredCont_Model_", bvar, "_", nvar, ".csv"), row.names = FALSE)
    
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

# LATENT CHANGE PREDICTION: SUMMARY OF RESULTS WITH LATENT CLINICAL CHANGE VARIABLES -----------------------------------------
##### Unstandardized results: Correction for multiple comparisons #####
# loop to extract the results for each MRI variable, perform a p-value correction across the six tests, and merge to one data frame 
df_si_list = list()
for (nvar in neurovars) {
  # shorten variable names of MRI variables for MPLUS and convert to upper case
  nvar_mplus = nvar %>% str_replace('fc','') %>% toupper()
  # specify paths to csv files with path statistics
  tmp_paths = list.files(path = paste0(out_dir,"MplusOutputs"),
                         pattern = sprintf("^SingleIndicatorSyChangePred_PathsUnstandardized_.*_%s\\.csv$", nvar),
                         full.names = TRUE)
  # create empty list
  df_list = list()
  # import CSV files in loop
  for (path in tmp_paths) {
    df = read.csv(path)
    df_filtered = subset(df,param==paste0(nvar_mplus,'_CH') & paramHeader != 'Variances' & paramHeader != 'Means')
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
write.csv(df_si_stats_full,paste0(out_dir,"Napls_LatentVariableSyChangePred_UnstandardizedResultsSummary_MainEffBonf_FCHip.csv"), row.names = FALSE)
write_xlsx(df_si_stats_full,path = paste0(out_dir,"Napls_LatentVariableSyChangePred_UnstandardizedResultsSummary_MainEffBonf_FCHip.xlsx"))

##### Unstandardized results: Summary of full results #####
# list all csv containing full test stats of hippocampal FC
path_files_hip = list.files(path = paste0(out_dir,"MplusOutputs"),
                            pattern = "^SingleIndicatorSyChangePred_PathsUnstandardized_.*fchi\\.csv$",
                            full.names = TRUE)
# read them all and bind rowwise
df_path_hip = path_files_hip %>% lapply(read.csv) %>% bind_rows() %>% as.data.frame()
# replace 999 with NA
df_path_hip[df_path_hip == 999] = NA
# export data
write.csv(df_path_hip,file = paste0(out_dir,"Napls_LatentVariableSyChangePred_UnstandardizedResultsSummary_FullStats_FCHip.csv"), row.names = FALSE)
write_xlsx(df_path_hip,path = paste0(out_dir,"Napls_LatentVariableSyChangePred_UnstandardizedResultsSummary_FullStats_FCHip.xlsx"))


##### Standardized results: Correction for multiple comparisons #####
# loop to extract the results for each MRI variable, perform a p-value correction across the six tests, and merge to one data frame 
df_si_list = list()
for (nvar in neurovars) {
  # shorten variable names of MRI variables for MPLUS and convert to upper case
  nvar_mplus = nvar %>% str_replace('fc','') %>% toupper()
  # specify paths to csv files with path statistics
  tmp_paths = list.files(path = paste0(out_dir,"MplusOutputs"),
                         pattern = sprintf("^SingleIndicatorSyChangePred_PathsStandardized_.*_%s\\.csv$", nvar),
                         full.names = TRUE)
  # create empty list
  df_list = list()
  # import CSV files in loop
  for (path in tmp_paths) {
    df = read.csv(path)
    df_filtered = subset(df,param==paste0(nvar_mplus,"_CH") & paramHeader != 'Variances' & paramHeader != 'Means')
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
write.csv(df_si_stats_full,paste0(out_dir,"Napls_LatentVariableSyChangePred_StandardizedResultsSummary_MainEffBonf_FCHip.csv"), row.names = FALSE)
write_xlsx(df_si_stats_full,path = paste0(out_dir,"Napls_LatentVariableSyChangePred_StandardizedResultsSummary_MainEffBonf_FCHip.xlsx"))

##### Standardized results: Summary of full results #####
# list all csv containing full test stats of hippocampal FC
path_files_hip = list.files(path = paste0(out_dir,"MplusOutputs"),
                            pattern = "^SingleIndicatorSyChangePred_PathsStandardized_.*fchi\\.csv$",
                            full.names = TRUE)
# read them all and bind rowwise
df_path_hip = path_files_hip %>% lapply(read.csv) %>% bind_rows() %>% as.data.frame()
# replace 999 with NA
df_path_hip[df_path_hip == 999] = NA
# export data
write.csv(df_path_hip,file = paste0(out_dir,"Napls_LatentVariableSyChangePred_StandardizedResultsSummary_FullStats_FCHip.csv"), row.names = FALSE)
write_xlsx(df_path_hip,path = paste0(out_dir,"Napls_LatentVariableSyChangePred_StandardizedResultsSummary_FullStats_FCHip.xlsx"))


##### Summary of model fit results #####
# run loop to extract their full test statistics on model fit results
df_list = list()
for (bvar in behavvars) {
  for (nvar in neurovars) {
    p = list.files(path = paste0(out_dir,"MplusOutputs"),
                   pattern = sprintf("^SingleIndicatorSyChangePred_Model_%s_%s\\.csv$",bvar,nvar),
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
write.csv(df_modelfit,file = paste0(out_dir,"Napls_LatentVariableSyChangePred_ResultsSummary_ModelFit_FCHip.csv"), row.names = FALSE)
write_xlsx(df_modelfit, path = paste0(out_dir,"Napls_LatentVariableSyChangePred_ResultsSummary_ModelFit_FCHip.xlsx"))


##### VISUALIZATION #####
# import data
df_neg = read.csv(paste0(out_dir,"Napls_EarlySlopeFC_LateSlopeClin_fchi_NegS.csv"))
# select columns
df_neg = df_neg[c("subj","fchi_slope_first","NegS_slope_post")]
# rename columns
colnames(df_neg) = c("subj","slope_fchi_neg","slope_neg")

# plot
p = ggplot(df_neg, aes(x = slope_fchi_neg, y = slope_neg)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(alpha = 0.2, size = 2.5) +
  geom_smooth(method = "lm",color="red4") +
  scale_color_manual(values = c("grey30", "grey65")) +
  scale_x_continuous(limits = c(-0.3, 0.3),breaks = c(-0.2,-0.1, 0, 0.1,0.2)) +
  scale_y_continuous(limits = c(-6.5, 6.5),breaks = c(-5,-2.5, 0,2.5, 5)) +
  labs(x = "Early Change in Hippocampal FC",y = "Late Change in Negative Symp") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 14, color = "black", face = "bold"),
        axis.text.x  = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        axis.text.y  = element_text(size = 12, color = "black"),
        legend.position = "none")
# export plot as jpeg
ggsave(paste0(out_dir,"Napls_LatentVariablePredict_Cases_ScatterPlot_HippoFC_NegS.jpeg"),
       plot=p,height=100,width=100,unit="mm",dpi=500)



# LATENT CHANGE PREDICTION: PREDICTING FC CHANGE BASED ON PRIOR CLINICAL CHANGE  ---------------------
# run slope estimations, imputations and latent variable regression in one loop
for (bvar in behavvars) {
  for (nvar in neurovars) {
    ##### Slope estimations based on single subject regressions #####
    # select variables of interest
    df = df_napls_mlr_imp[c("subj","CHR","age","female","site02","site03","site04","site05","site06","site07","site08","site09",
                            paste0("nt_",bvar),paste0("nt_",nvar),
                            paste0(bvar,"_BL"),paste0(bvar,"_M2"),paste0(bvar,"_M4"),paste0(bvar,"_M6"),paste0(bvar,"_M8"),
                            paste0(nvar,"_BL"),paste0(nvar,"_M2"),paste0(nvar,"_M4"),paste0(nvar,"_M6"),paste0(nvar,"_M8"))]
    # select only subjects with at least 2 measurement timepoints in FC and clinical variable and baseline value in clinical
    df = df[!is.na(df[[paste0(bvar, "_BL")]]) & df[[paste0("nt_",bvar)]] > 1 & df[[paste0("nt_",nvar)]] > 1, ]
    # reset rownames
    rownames(df) = NULL
    # define follow-up timepoints and the corresponding months
    timepoints = c("M2", "M4", "M6", "M8")
    months = c(2, 4, 6, 8)
    # compute slope per month from baseline for each follow-up timepoint
    baseline_col = paste0(bvar, "_BL") # save baseline column
    for (i in seq_along(timepoints)) {
      tp = timepoints[i] # define timepoint
      month = months[i] # define month
      followup_col = paste0(bvar, "_", tp) # define follow-up column
      slope_col = paste0(bvar, "_slope_", tp)  # new column for slope
      df[[slope_col]] = (df[[followup_col]] - df[[baseline_col]]) / month # compute slope
    }
    # initialize new columns
    df[[paste0(bvar, "_slope_first")]] = NA
    df[[paste0(bvar, "_slope_first_end")]]  = NA
    df[[paste0(nvar, "_slope_post")]] = NA
    df[[paste0(nvar, "_slope_post_se")]] = NA
    df[[paste0(nvar, "_slope_post_var_month")]] = NA
    df[[paste0(nvar, "_slope_post_var_value")]] = NA
    df[[paste0(nvar, "_slope_post_nt")]] = NA
    df[[paste0(nvar, "_slope_post_start")]] = NA
    df[[paste0(nvar, "_slope_post_end")]]   = NA
    # extract earliest slope for FC and compute slope of clinical variable for the subsequent time period
    for (i in 1:nrow(df)) {
      # extract slopes for clinical variable
      slopes = sapply(timepoints, function(tp) df[[paste0(bvar, "_slope_", tp)]][i])
      # get earliest non-missing slope for FC variable
      first_slope_idx = which(!is.na(slopes))[1]
      if (!is.na(first_slope_idx)) {
        df[[paste0(bvar, "_slope_first")]][i] = slopes[first_slope_idx]
        slope_month = months[first_slope_idx]
        df[[paste0(bvar, "_slope_first_end")]][i] = slope_month
        
        # get values of FC variable starting at or after slope_month
        nvar_values = sapply(timepoints, function(tp) df[[paste0(nvar, "_", tp)]][i])
        valid_idx = which(!is.na(nvar_values) & months >= slope_month)
        
        # run slope calculation for FC variable if there is more than two timepoints after earliest slope timepoint of FC variable
        # and if there is a valid FC value at the earliest slope timepoint of FC variable
        if (length(valid_idx) >= 2 && months[valid_idx[1]] == slope_month) {
          # extract timepoints and values of FC variable after earliest slope timepoint for clinical
          t_vals = months[valid_idx]
          v_vals = nvar_values[valid_idx]
          # compute regression slope for FC variable
          fit = lm(v_vals ~ t_vals)
          coefs = coef(summary(fit)) # get coefficients
          est_slope = coefs["t_vals", "Estimate"] # extract slope
          se_slope  = coefs["t_vals", "Std. Error"] # extract standard error of slope
          var_month = var(t_vals)  # compute variance of timepoints
          var_value = var(v_vals)  # compute variance of value
          df[[paste0(nvar, "_slope_post")]][i] = est_slope
          df[[paste0(nvar, "_slope_post_se")]][i] = se_slope
          df[[paste0(nvar, "_slope_post_var_month")]][i] = var_month
          df[[paste0(nvar, "_slope_post_var_value")]][i] = var_value
          df[[paste0(nvar, "_slope_post_nt")]][i] = length(valid_idx)
          df[[paste0(nvar, "_slope_post_start")]][i] = min(t_vals)
          df[[paste0(nvar, "_slope_post_end")]][i]   = max(t_vals)
        }
      }
    }
    # replace NaN with NA
    df[df == "NaN"] = NA
    # remove subjects without slope in clinical variable
    df = df[!is.na(df[paste0(nvar, "_slope_post")]),]
    # export data
    write.csv(df, file = paste0(out_dir,"Napls_EarlySlopeClin_LateSlopeFC_", bvar, ".csv"), row.names = FALSE)
    
    ##### Imputation of standard errors #####
    # select variables of interest
    df = df[c("subj","CHR","age","female","site02","site03","site04","site05","site06","site07","site08","site09",
              paste0(bvar,"_slope_first"),paste0(nvar,"_slope_post"),
              paste0(nvar, "_slope_post_se"),paste0(nvar, "_slope_post_var_month"),paste0(nvar, "_slope_post_var_value"),
              paste0(nvar, "_slope_post_nt"))]
    # create formula
    formula_str = paste(paste0(nvar, "_slope_post_se"),"~",
                        paste0(nvar, "_slope_post"),"+",
                        paste0(nvar, "_slope_post_nt"),"+",
                        paste0(nvar, "_slope_post_var_month"),"+",
                        paste0(nvar, "_slope_post_var_value"))
    # run regression model
    res = lm(as.formula(formula_str), data = df)
    # extract model-level statistics
    model_summary = summary(res)
    r2 = model_summary$r.squared # r2
    adj_r2 = model_summary$adj.r.squared # adjusted r2
    fstat = model_summary$fstatistic # F stats
    model_p = pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE) # model p-value
    # count missing SEs and compute imputation stats
    outcome_se_var = paste0(nvar, "_slope_post_se")
    outcome_nt_var = paste0(nvar, "_slope_post_nt")
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
    write.csv(model_stats_df, file = paste0(out_dir,"ImputationSE_model_FCChangePred_", bvar, ".csv"), row.names = FALSE)
    
    ##### Latent variable regression #####
    # re-specify working directory which must not be too long for MPLUS
    wd = "/home/lukas/mplus"
    setwd(wd)
    # print process
    print(paste0("Running Single Indicator Approach for ",bvar, " and ",nvar,"..."))
    # shorten variable names of MRI variables for MPLUS
    nvar_mplus = nvar %>% str_replace("fc","")
    # create data frame with variables of interest
    df = df[c(paste0(nvar,"_slope_post"),paste0(bvar,"_slope_first"),paste0(nvar,"_slope_post_se"),
              "age","female","site02","site03","site04","site05","site06","site07","site08","site09")]
    # get average error variances across group from standard errors of both variables
    se_nvar = mean(df[[paste0(nvar,"_slope_post_se")]],na.rm=T)^2
    # replace NA with -99
    df[is.na(df)] = -99
    # drop standard error columns again
    df = df[ , !grepl("_slope_post_se$", names(df)) ]
    # define the input and output files
    datafile = file.path(wd, paste0("NaplsData_", nvar, "_", bvar, "_FCChangePred.dat"))
    inpfile = file.path(wd, paste0("NaplsData_", nvar, "_", bvar, "_FCChangePred.inp"))
    # write output file
    write.table(df,datafile,quote = F,row.names = F, col.names = F,append = F, sep = "\t")
    # define the MPLUS input file containing the model
    inp_text = paste0(
      "Title:      LATENT VARIABLE REGRESSION

    Data:       File is ", datafile, ";
    
    Variable:
    names are
    ", nvar_mplus, "_CH
    ", bvar, "_CH
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
    ", nvar_mplus, "_CH
    ", bvar, "_CH
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
    
    LV_", nvar_mplus, " on ", bvar, "_CH 
               age female site02 site03 site04
               site05 site06 site07 site08
               site09;
    
    ",
      
      bvar, "_CH WITH age female site02 site03 site04
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
    
    SAVEDATA: ESTIMATES = power_",nvar,"_",bvar,"_fcchangepred.dat;
    Output:     tech1 tech4 stdyx sampstat;
    ")
    # write as input file for MPLUS
    writeLines(inp_text, con = inpfile)
    # run MPLUS
    system(paste("wine Mplus-8.6.exe", basename(inpfile)))
    # specify the output file
    outfile = file.path(wd, paste0("naplsdata_", tolower(nvar), "_", tolower(bvar), "_fcchangepred.out"))
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
    write.csv(df_res, file = paste0("SingleIndicatorFCChangePred_PathsUnstandardized_",bvar,"_",nvar, ".csv"), row.names = FALSE)
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
    write.csv(df_res, file = paste0("SingleIndicatorFCChangePred_PathsStandardized_",bvar,"_",nvar, ".csv"), row.names = FALSE)
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
    write.csv(df_model, file = paste0("SingleIndicatorFCChangePred_Model_", bvar, "_", nvar, ".csv"), row.names = FALSE)
    
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
for (bvar in behavvars) {
  # specify paths to csv files with path statistics
  tmp_paths = list.files(path = paste0(out_dir,"MplusOutputs"),
                         pattern = sprintf("^SingleIndicatorFCChangePred_PathsUnstandardized_.*_%s\\.csv$", nvar),
                         full.names = TRUE)
  # create empty list
  df_list = list()
  # import CSV files in loop
  for (path in tmp_paths) {
    df = read.csv(path)
    df_filtered = subset(df,param==toupper(paste0(bvar,'_CH')) & paramHeader != 'Variances' & paramHeader != 'Means')
    df_list[[length(df_list) + 1]] = df_filtered
  }
  # create one data frame
  df_si_stats = do.call(rbind, df_list)
  # add data frame to list
  df_si_list[[length(df_si_list) + 1]] = df_si_stats
}

##### Unstandardized results: Summary of main effect results #####
# create one data frame
df_si_stats_full = as.data.frame(do.call(rbind, df_si_list))
# insert adjusted p-values and mark significance using different methods
df_si_stats_full$p_fdr = p.adjust(df_si_stats_full$pval_exact,method="fdr")
df_si_stats_full$p_fdr_sig = df_si_stats_full$p_fdr < 0.05
df_si_stats_full$p_bonf = p.adjust(df_si_stats_full$pval_exact,method="bonferroni")
df_si_stats_full$p_bonf_sig = df_si_stats_full$p_bonf < 0.05
# export data
write.csv(df_si_stats_full,paste0(out_dir,"Napls_LatentVariableFCChangePred_UnstandardizedResultsSummary_MainEffBonf_FCHip.csv"), row.names = FALSE)
write_xlsx(df_si_stats_full,path = paste0(out_dir,"Napls_LatentVariableFCChangePred_UnstandardizedResultsSummary_MainEffBonf_FCHip.xlsx"))

##### Unstandardized results: Summary of full results #####
# list all csv containing full test stats of hippocampal FC
path_files_hip = list.files(path = paste0(out_dir,"MplusOutputs"),
                            pattern = "^SingleIndicatorFCChangePred_PathsUnstandardized_.*fchi\\.csv$",
                            full.names = TRUE)
# read them all and bind rowwise
df_path_hip = path_files_hip %>% lapply(read.csv) %>% bind_rows() %>% as.data.frame()
# replace 999 with NA
df_path_hip[df_path_hip == 999] = NA
# export data
write.csv(df_path_hip,file = paste0(out_dir,"Napls_LatentVariableFCChangePred_UnstandardizedResultsSummary_FullStats_FCHip.csv"), row.names = FALSE)
write_xlsx(df_path_hip,path = paste0(out_dir,"Napls_LatentVariableFCChangePred_UnstandardizedResultsSummary_FullStats_FCHip.xlsx"))


##### Standardized results: Correction for multiple comparisons #####
# loop to extract the results for each MRI variable, perform a p-value correction across the six tests, and merge to one data frame 
df_si_list = list()
for (bvar in behavvars) {
  # specify paths to csv files with path statistics
  tmp_paths = list.files(path = paste0(out_dir,"MplusOutputs"),
                         pattern = sprintf("^SingleIndicatorFCChangePred_PathsStandardized_.*_%s\\.csv$", nvar),
                         full.names = TRUE)
  # create empty list
  df_list = list()
  # import CSV files in loop
  for (path in tmp_paths) {
    df = read.csv(path)
    df_filtered = subset(df,param==toupper(paste0(bvar,'_CH')) & paramHeader != 'Variances' & paramHeader != 'Means')
    df_list[[length(df_list) + 1]] = df_filtered
  }
  # create one data frame
  df_si_stats = do.call(rbind, df_list)
  # add data frame to list
  df_si_list[[length(df_si_list) + 1]] = df_si_stats
}

##### Standardized results: Summary of main effect results #####
# create one data frame
df_si_stats_full = as.data.frame(do.call(rbind, df_si_list))
# insert adjusted p-values and mark significance using different methods
df_si_stats_full$p_fdr = p.adjust(df_si_stats_full$pval_exact,method="fdr")
df_si_stats_full$p_fdr_sig = df_si_stats_full$p_fdr < 0.05
df_si_stats_full$p_bonf = p.adjust(df_si_stats_full$pval_exact,method="bonferroni")
df_si_stats_full$p_bonf_sig = df_si_stats_full$p_bonf < 0.05
# export data
write.csv(df_si_stats_full,paste0(out_dir,"Napls_LatentVariableFCChangePred_StandardizedResultsSummary_MainEffBonf_FCHip.csv"), row.names = FALSE)
write_xlsx(df_si_stats_full,path = paste0(out_dir,"Napls_LatentVariableFCChangePred_StandardizedResultsSummary_MainEffBonf_FCHip.xlsx"))

##### Standardized results: Summary of full results #####
# list all csv containing full test stats of hippocampal FC
path_files_hip = list.files(path = paste0(out_dir,"MplusOutputs"),
                            pattern = "^SingleIndicatorFCChangePred_PathsStandardized_.*fchi\\.csv$",
                            full.names = TRUE)
# read them all and bind rowwise
df_path_hip = path_files_hip %>% lapply(read.csv) %>% bind_rows() %>% as.data.frame()
# replace 999 with NA
df_path_hip[df_path_hip == 999] = NA
# export data
write.csv(df_path_hip,file = paste0(out_dir,"Napls_LatentVariableFCChangePred_StandardizedResultsSummary_FullStats_FCHip.csv"), row.names = FALSE)
write_xlsx(df_path_hip,path = paste0(out_dir,"Napls_LatentVariableFCChangePred_StandardizedResultsSummary_FullStats_FCHip.xlsx"))

##### Summary of model fit results #####
# run loop to extract their full test statistics on model fit results
df_list = list()
for (bvar in behavvars) {
  for (nvar in neurovars) {
    p = list.files(path = paste0(out_dir,"MplusOutputs"),
                   pattern = sprintf("^SingleIndicatorFCChangePred_Model_%s_%s\\.csv$",bvar,nvar),
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
write.csv(df_modelfit,file = paste0(out_dir,"Napls_LatentVariableFCChangePred_ResultsSummary_ModelFit_FCHip.csv"), row.names = FALSE)
write_xlsx(df_modelfit, path = paste0(out_dir,"Napls_LatentVariableFCChangePred_ResultsSummary_ModelFit_FCHip.xlsx"))
