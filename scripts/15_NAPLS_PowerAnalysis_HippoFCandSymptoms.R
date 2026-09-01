# POWER ANALYSIS: CHANGE OF HIPPCAMPUS FC AND CHANGE OF CLINICAL VARIABLES ---------------------------------------------
# Description:  This script runs the power analysis analysis on the link between changes in hippocampal FC and clinical changes.
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
in_dir = "results/napls/07_NAPLS_SingleSubjectRegressions/"
# define output directory
out_dir = "results/napls/15_NAPLS_PowerAnalysis_HippoFCandSymptoms/"
# create output directory for this script
dir.create(out_dir, recursive = TRUE,showWarnings=F)
# read data including behavioral and FC data and all slopes with imputed SEs
df_napls_mlr_imp = read.csv(paste0(in_dir,"Napls_BehavFCandSlopes_impSE.csv"))
# remove subjects that have non-enhanced psychosis risk, as only 1 of these subjects has valid data for FC at two time points
df_napls_mlr_imp = subset(df_napls_mlr_imp,group !=2)
# rename group column to CHR
df_napls_mlr_imp = df_napls_mlr_imp %>% rename(CHR = group)
# create vectors with behavioral variables
behavvars = c("PosS","NegS","DepS","gaf","VeMem","SyCod")
# create vector for hippocampal FC
neurovars = c("fchi")


# POWER ANALYSIS OF LATENT VARIABLE REGRESSION: CASE-CONTROL COMPARISON ----------------------------------------------------------------------------------------------------
# Note: This code runs the power analysis for the main latent variable regression models in MPLUS comparing the link between hippocampal
# FC change and clinical severity change over time between CHR and healthy controls.

##### Power analysis loop running MPLUS #####
# Note: in this analysis cases and controls are used that have a slope in the respective behavioral AND the MRI variable.
# re-define and set working directory in which the MPLUS .exe file is stored -> necessary because the path of the working directory
# must not be too long for MPLUS
wd = "/home/lukas/mplus"
setwd(wd)
# run latent variable regression and power analysis in loop
for (bvar in behavvars) {
  for (nvar in neurovars) {
    print(paste0("Running Power Analysis for ",bvar, " and ",nvar,"..."))
    # shorten variable names of MRI variables for MPLUS
    nvar_mplus = nvar %>% str_replace("fc","")
    # create data frame with variables of interest
    df = df_napls_mlr_imp[c(paste0(bvar,"_CH"),paste0(nvar,"_CH"),"CHR",paste0(bvar,"_SE"),paste0(nvar,"_SE"),
                            "age","female","site02","site03","site04","site05","site06","site07","site08","site09")]
    # remove subjects with any missings
    df = df[complete.cases(df[c(paste0(bvar, "_CH"), paste0(nvar, "_CH"))]), ]
    # get unstandardized variance of the slopes
    var_slope_bvar = var(df[c(paste0(bvar, "_CH"))])
    var_slope_nvar = var(df[c(paste0(nvar, "_CH"))])
    # get average error variances across group from standard errors of both variables
    se_bvar = mean(df[[paste0(bvar,"_SE")]],na.rm=T)^2
    se_nvar = mean(df[[paste0(nvar,"_SE")]],na.rm=T)^2
    # compute reliability of unstandardized slopes
    reliab_bvar = (var_slope_bvar-se_bvar)/var_slope_bvar
    reliab_nvar = (var_slope_nvar-se_nvar)/var_slope_nvar
    # compute standardized error variance
    se_bvar = 1-reliab_bvar
    se_nvar = 1-reliab_nvar
    # z-standardize slope columns
    df[c(paste0(bvar, "_CH"), paste0(nvar, "_CH"))] = scale(df[c(paste0(bvar, "_CH"), paste0(nvar, "_CH"))])
    # drop standard error columns again
    df = df[ , !grepl("_SE$", names(df)) ]
    # define sample size
    n = nrow(df)
    # define the input and output files
    datafile = file.path(wd, paste0("NaplsData_", nvar, "_", bvar, ".dat"))
    inpfile = file.path(wd, paste0("NaplsData_", nvar, "_", bvar, ".inp"))
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
        CHR
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
        CHR
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
      
        ! auxiliary is (m) c1 c2 c3;
        
        Analysis:   type = random;
                    Estimator = MLR;
                    ALGORITHM=INTEGRATION;
                    STARTS=20;
                    STITERATION=50000;
        
        Model:
        
        LV_", nvar_mplus, " by ", nvar_mplus, "_CH;
        ", nvar_mplus, "_CH@", format(se_nvar, scientific = FALSE), ";
        
        LV_", bvar, " by ", bvar, "_CH;
        ", bvar, "_CH@", format(se_bvar, scientific = FALSE), ";
        
        ", nvar_mplus, "xCHR | LV_", nvar_mplus, " xwith CHR;
        
        LV_", bvar, " on LV_", nvar_mplus, " CHR ", nvar_mplus, "xCHR age
                   female site02 site03 site04
                   site05 site06 site07 site08
                   site09", ";
        
        LV_", nvar_mplus, " WITH CHR age female site02 site03 site04
                           site05 site06 site07 site08
                           site09", ";
        
        CHR WITH age female site02 site03 site04
                         site05 site06 site07 site08
                         site09", ";
        
        age WITH female site02 site03 site04 site05
                        site06 site07 site08
                        site09", ";
        
        female WITH site02 site03 site04 site05
                        site06 site07 site08
                        site09", ";
        
        site02 WITH site03 site04 site05
                        site06 site07 site08
                        site09", ";
        
        site03 WITH site04 site05 site06 site07 site08
                        site09", ";
        
        site04 WITH site05 site06 site07 site08
                        site09", ";
        
        site05 WITH site06 site07 site08
                        site09", ";
        
        site06 WITH site07 site08 site09", ";
        
        site07 WITH site08 site09", ";
        
        site08 WITH site09", ";
        
        SAVEDATA: ESTIMATES = power_",nvar,"_",bvar,".dat;
        Output:     tech1 tech4 stdyx sampstat;
        "
    )
    # write as input file for MPLUS
    writeLines(inp_text, con = inpfile)
    # run MPLUS
    system(paste("wine Mplus-8.6.exe", basename(inpfile)))
    # specify the output file
    outfile = file.path(wd, paste0("naplsdata_", tolower(nvar), "_", tolower(bvar), ".out"))
    
    # define .inp file for power analysis
    inpfile = file.path(wd, paste0("Napls_", nvar, "_", bvar, "_PA.inp"))
    # create the Mplus syntax dynamically
    inp_text =  paste0("
TITLE: Post-hoc Monte Carlo power analysis 
                  
MONTECARLO:
NAMES = ", bvar, "_CH ", nvar_mplus, "_CH CHR age female site02 site03 site04 site05 site06 
site07 site08 site09 ;
                    
NOBSERVATIONS = ", n, ";
NREPS = 10000;
SEED = 10000;
POPULATION = power_",nvar,"_", bvar, ".dat;
COVERAGE   = power_",nvar,"_", bvar, ".dat;
                
ANALYSIS:
TYPE = RANDOM;
ALGORITHM = INTEGRATION;
PROCESSORS = 4; 
                
MODEL POPULATION:
LV_", nvar_mplus, " by ", nvar_mplus, "_CH;
", nvar_mplus, "_CH@", format(se_nvar, scientific = FALSE), ";
                
LV_", bvar, " by ", bvar, "_CH;
", bvar, "_CH@", format(se_bvar, scientific = FALSE), ";
                
", nvar_mplus, "xCHR | LV_", nvar_mplus, " xwith CHR;
                    
LV_", bvar, " on LV_", nvar_mplus, " CHR ", nvar_mplus, "xCHR age female site02 site03 site04 site05 site06 
site07 site08 site09;
                    
LV_", nvar_mplus, " WITH CHR age female site02 site03 site04 site05 site06 
site07 site08 site09;
                    
CHR WITH age female site02 site03 site04 site05 site06 
site07 site08 site09;
                    
age WITH female site02 site03 site04 site05 site06 
site07 site08 site09;
                    
female WITH site02 site03 site04 site05 site06 
site07 site08 site09;
                    
site02 WITH site03 site04 site05 site06 
site07 site08 site09;
                    
site03 WITH site04 site05 site06 
site07 site08 site09;
                    
site04 WITH site05 site06 
site07 site08 site09;
                    
site05 WITH site06 
site07 site08 site09;
                    
site06 WITH site07 site08 site09;
                    
site07 WITH site08 site09;
                    
site08 WITH site09;
                    
MODEL:
LV_", nvar_mplus, " by ", nvar_mplus, "_CH;
", nvar_mplus, "_CH@", format(se_nvar, scientific = FALSE), ";
                
LV_", bvar, " by ", bvar, "_CH;
", bvar, "_CH@", format(se_bvar, scientific = FALSE), ";
                
", nvar_mplus, "xCHR | LV_", nvar_mplus, " xwith CHR;
                    
LV_", bvar, " on LV_", nvar_mplus, " CHR ", nvar_mplus, "xCHR age female site02 site03 site04 site05 site06 
site07 site08 site09;
                    
LV_", nvar_mplus, " WITH CHR age female site02 site03 site04 site05 site06 
site07 site08 site09;
                    
CHR WITH age female site02 site03 site04 site05 site06 
site07 site08 site09;
                    
age WITH female site02 site03 site04 site05 site06 
site07 site08 site09;
                    
female WITH site02 site03 site04 site05 site06 
site07 site08 site09;
                    
site02 WITH site03 site04 site05 site06 
site07 site08 site09;
                    
site03 WITH site04 site05 site06 
site07 site08 site09;
                    
site04 WITH site05 site06 
site07 site08 site09;
                    
site05 WITH site06 
site07 site08 site09;
                    
site06 WITH site07 site08 site09;
                    
site07 WITH site08 site09;
                    
site08 WITH site09;

OUTPUT: TECH9;")
    
    # export to file
    writeLines(inp_text, con = inpfile)
    # run MPLUS
    system(paste("wine Mplus-8.6.exe", basename(inpfile)))
    # specify the output file
    outfile = file.path(wd, paste0("napls_", tolower(nvar), "_", tolower(bvar), "_pa.out"))
    # read the output
    model_out = readModels(outfile)
    # extract unstandardized estimates
    df_res = as.data.frame(model_out$parameters$unstandardized)
    # export model stats
    write.csv(df_res, file = paste0("SingleIndicator_Case-Control_PathsUnstandardized_",bvar,"_",nvar, "_PowerAnalysis.csv"), row.names = FALSE)
    
  }
}
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


# SUMMARY OF RESULTS ----------------------------------------------------------------------------------------------------
##### Summary of full results #####
# list all csv containing full test stats of hippocampal FC
path_files_hip = list.files(path = paste0(out_dir,"MplusOutputs"),
                            pattern = "^SingleIndicator_Case-Control_PathsUnstandardized_.*_PowerAnalysis\\.csv$",
                            full.names = TRUE)
# read them all and bind rowwise
df_path_hip = path_files_hip %>% lapply(read.csv) %>% bind_rows() %>% as.data.frame()
# replace 999 with NA
df_path_hip[df_path_hip == 999] = NA
# export data
write.csv(df_path_hip,file = paste0(out_dir,"Napls_LatentVariableRegressPower_Case-Control_UnstandardizedResultsSummary_FullStats_FCHip.csv"), row.names = FALSE)
write_xlsx(df_path_hip,path = paste0(out_dir,"Napls_LatentVariableRegressPower_Case-Control_UnstandardizedResultsSummary_FullStats_FCHip.xlsx"))


# POWER ANALYSIS OF LATENT VARIABLE REGRESSION: PROSEPCTIVE ASSOCIATION ----------------------------------------------------------------------------------------------------
# Note: This code runs the power analysis for thelatent variable regression models in MPLUS including the prospective association.

##### Power analysis loop running MPLUS #####
# get only CHR-P
df_napls_mlr_imp = subset(df_napls_mlr_imp,CHR == 1)
# define variables
behavvars = c("NegS","DepS","gaf")
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
    # get unstandardized variance of the slopes
    var_slope_bvar = var(df[c(paste0(bvar, "_slope_post"))])
    # get average error variances across group from standard errors of both variables
    se_bvar = mean(df[[paste0(bvar,"_slope_post_se")]],na.rm=T)^2
    # compute reliability of unstandardized slopes
    reliab_bvar = (var_slope_bvar-se_bvar)/var_slope_bvar
    # compute standardized error variance
    se_bvar = 1-reliab_bvar
    # z-standardize slope columns
    df[c(paste0(bvar, "_slope_post"), paste0(nvar, "_slope_first"))] = scale(df[c(paste0(bvar, "_slope_post"), paste0(nvar, "_slope_first"))])
    # drop standard error columns again
    df = df[ , !grepl("_slope_post_se$", names(df)) ]
    # define sample size
    n = nrow(df)
    # replace NA with -99
    df[is.na(df)] = -99
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
    
    # define .inp file for power analysis
    inpfile = file.path(wd, paste0("Napls_", nvar, "_", bvar, "_SyChangePred_PA.inp"))
    inp_text = paste0("
    TITLE: Post-hoc Monte Carlo power analysis 
                  
    MONTECARLO:
    NAMES = ", bvar, "_CH ", nvar_mplus, "_CH age female site02 site03 site04 site05 site06 
    site07 site08 site09 ;
                        
    NOBSERVATIONS = ", n, ";
    NREPS = 10000;
    SEED = 10000;
    POPULATION = power_",nvar,"_", bvar, "_sychangepred.dat;
    COVERAGE   = power_",nvar,"_", bvar, "_sychangepred.dat;
                    
    ANALYSIS:
    TYPE = RANDOM;
    ALGORITHM = INTEGRATION;
    PROCESSORS = 4; 
                    
    MODEL POPULATION:
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
    
    MODEL:
    
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

    Output:     tech9;
    "
    )
    # export to file
    writeLines(inp_text, con = inpfile)
    # run MPLUS
    system(paste("wine Mplus-8.6.exe", basename(inpfile)))
    # specify the output file
    outfile = file.path(wd, paste0("napls_", tolower(nvar), "_", tolower(bvar), "_sychangepred_pa.out"))
    # read the output
    model_out = readModels(outfile)
    # extract unstandardized estimates
    df_res = as.data.frame(model_out$parameters$unstandardized)
    # export model stats
    write.csv(df_res, file = paste0("SingleIndicator_OnlyPatients_PathsUnstandardized_",bvar,"_",nvar, "_PowerAnalysis_Prediction.csv"), row.names = FALSE)
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


# SUMMARY OF RESULTS ----------------------------------------------------------------------------------------------------
##### Summary of full results #####
# list all csv containing full test stats of hippocampal FC
path_files_hip = list.files(path = paste0(out_dir,"MplusOutputs"),
                            pattern = "^SingleIndicator_OnlyPatients_PathsUnstandardized_.*_PowerAnalysis_Prediction\\.csv$",
                            full.names = TRUE)
# read them all and bind rowwise
df_path_hip = path_files_hip %>% lapply(read.csv) %>% bind_rows() %>% as.data.frame()
# replace 999 with NA
df_path_hip[df_path_hip == 999] = NA
# export data
write.csv(df_path_hip,file = paste0(out_dir,"Napls_LatentVariablePredictionPower_OnlyPatients_UnstandardizedResultsSummary_FullStats_FCHip.csv"), row.names = FALSE)
write_xlsx(df_path_hip,path = paste0(out_dir,"Napls_LatentVariablePredictionPower_OnlyPatients_UnstandardizedResultsSummary_FullStats_FCHip.xlsx"))





















