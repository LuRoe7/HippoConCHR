# ROBUSTNESS CHECK: CHANGE OF HIPPOCAMPAL FC AND CLINICAL SCORES ONLY IN PATIENTS ---------------------------------------------
# Description:  This script runs a robustness check on the link between changes in hippocampal FC and clinical changes only in patients.
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
df_napls_mlr_imp = read.csv(paste0(in_dir,"Napls_BehavFCandSlopes_impSE.csv"))
# keep only CHR individuals
df_napls_mlr_imp = subset(df_napls_mlr_imp,group ==1)
# create vectors with behavioral variables
behavvars = c("PosS","NegS","DepS","gaf","SyCod", "VeMem")
# create vector for hippocampal FC
neurovars = c("fchi")


# LATENT VARIABLE REGRESSION: ONLY IN PATIENTS  ---------------------
# Note: This code runs the latent variable regression models in Mplus studying the link between hippocampal FC change and
# symptom severity change over time only in CHR-P individuals, while controlling for medication change.

##### Latent variable regression loop running MPLUS #####
# Note: in this analysis cases are used that have a slope in the respective behavioral AND the MRI variable.
# re-define and set working directory in which the MPLUS .exe file is stored -> necessary because the path of the working directory
# must not be too long for MPLUS
wd = "PATH/TO/MPLUS/DIRECTORY"
setwd(wd)
# run latent variable regression in loop
for (bvar in behavvars) {
  for (nvar in neurovars) {
    print(paste0("Running Single Indicator Approach for ",bvar, " and ",nvar,"..."))
    # shorten variable names of MRI variables for MPLUS
    nvar_mplus = nvar %>% str_replace("fc","")
    # create data frame with variables of interest
    df = df_napls_mlr_imp[c(paste0(bvar,"_CH"),paste0(nvar,"_CH"),paste0(bvar,"_SE"),paste0(nvar,"_SE"),
                            "age","female","site02","site03","site04","site05","site06","site07","site08","site09",
                            "cpz_CH","cpz_SE")]
    # remove subjects with any missings
    df = df[complete.cases(df[c(paste0(bvar, "_CH"), paste0(nvar, "_CH"))]), ]
    # replace NA with -99
    df[is.na(df)] = -99
    # get average error variances across group from standard errors of both variables
    se_bvar = mean(df[[paste0(bvar,"_SE")]],na.rm=T)^2
    se_nvar = mean(df[[paste0(nvar,"_SE")]],na.rm=T)^2
    se_cpz = mean(df$cpz_SE, na.rm = TRUE)^2
    # drop standard error columns again
    df = df[ , !grepl("_SE$", names(df)) ]
    # define the input and output files
    datafile = file.path(wd, paste0("NaplsData_", nvar, "_", bvar, "_OnlyCHR.dat"))
    inpfile = file.path(wd, paste0("NaplsData_", nvar, "_", bvar, "_OnlyCHR.inp"))
    # write output file
    write.table(df,datafile,quote = F,row.names = F, col.names = F,append = F, sep = "\t")
    # define the MPLUS input file containing the model
    inp_text = paste0(
        "Title:      LATENT VARIABLE REGRESSION

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
        site09
        cpz_CH;
        
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
        site09
        cpz_CH;
        
        missing are all (-99);
      
        ! auxiliary is (m) c1 c2 c3;
        
        Analysis:   type = random;
                    ALGORITHM=INTEGRATION;
                    STARTS=20;
                    STITERATION=50000;
        
        Model:
        
        LV_", nvar_mplus, " by ", nvar_mplus, "_CH;
        ", nvar_mplus, "_CH@", format(se_nvar, scientific = FALSE), ";
        
        LV_", bvar, " by ", bvar, "_CH;
        ", bvar, "_CH@", format(se_bvar, scientific = FALSE), ";
        
        LV_CPZ by cpz_CH;
        cpz_CH@", format(se_cpz, scientific = FALSE), ";
        
        LV_", bvar, " on LV_", nvar_mplus, " age
                   female site02 site03 site04
                   site05 site06 site07 site08
                   site09 LV_CPZ", ";
        
        LV_", nvar_mplus, " WITH age female site02 site03 site04
                           site05 site06 site07 site08
                           site09 LV_CPZ", ";
        
        age WITH female site02 site03 site04
                         site05 site06 site07 site08
                         site09 LV_CPZ", ";
        
        female WITH site02 site03 site04 site05
                        site06 site07 site08
                        site09 LV_CPZ", ";
        
        site02 WITH site03 site04 site05
                        site06 site07 site08
                        site09 LV_CPZ", ";
        
        site03 WITH site04 site05 site06 site07 site08
                        site09 LV_CPZ", ";
        
        site04 WITH site05 site06 site07 site08
                        site09 LV_CPZ", ";
        
        site05 WITH site06 site07 site08
                        site09 LV_CPZ", ";
        
        site06 WITH site07 site08 site09 LV_CPZ", ";
        
        site07 WITH site08 site09 LV_CPZ", ";
        
        site08 WITH site09 LV_CPZ", ";
        
        site09 WITH LV_CPZ", ";
        
        SAVEDATA: ESTIMATES = power_",nvar,"_",bvar,"_OnlyCHR.dat;
        Output:     tech1 tech4 stdyx sampstat;
        "
    )
    # write as input file for MPLUS
    writeLines(inp_text, con = inpfile)
    # run MPLUS
    system(paste("wine Mplus-8.6.exe", basename(inpfile)))
    # specify the output file
    outfile = file.path(wd, paste0("naplsdata_", tolower(nvar), "_", tolower(bvar), "_onlychr.out"))
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
    write.csv(df_res, file = paste0("SingleIndicator_OnlyCHR_PathsUnstandardized_",bvar,"_",nvar, ".csv"), row.names = FALSE)
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
    write.csv(df_res, file = paste0("SingleIndicator_OnlyCHR_PathsStandardized_",bvar,"_",nvar, ".csv"), row.names = FALSE)
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
    write.csv(df_model, file = paste0("SingleIndicator_OnlyCHR_Model_", bvar, "_", nvar, ".csv"), row.names = FALSE)
  }
}
# define and set old working directory
old_wd = "PATH/TO/DIRECTORY/"
setwd(old_wd)
# create sub-directory in output folder
dir.create(paste0(out_dir,"MplusOutputs"), recursive = TRUE,showWarnings=F)
# specify paths to output files
res_files = list.files(path = wd,pattern = "\\.(csv|dat|inp|out)$",recursive = TRUE, full.names = TRUE)
# copy to output sub-directory in output directory
copied = file.copy(res_files,paste0(out_dir,"MplusOutputs"), overwrite = TRUE)
# remove copied files from intermediate working directory
file.remove(res_files[copied])


# SUMMARY OF UNSTANDARDIZED RESULTS ----------------------------------------------------------------------------------------------------

##### Correction for multiple comparisons #####
# loop to extract the results for each MRI variable, perform a p-value correction across the six tests, and merge to one data frame 
df_si_list = list()
for (nvar in neurovars) {
  # shorten variable names of MRI variables for MPLUS and convert to upper case
  nvar_mplus = nvar %>% str_replace('fc','') %>% toupper()
  # specify paths to csv files with path statistics
  tmp_paths = list.files(path = paste0(out_dir,"MplusOutputs"),
                         pattern = sprintf("^SingleIndicator_OnlyCHR_PathsUnstandardized_.*_%s\\.csv$", nvar),
                         full.names = TRUE)
  # create empty list
  df_list = list()
  # import CSV files in loop
  for (path in tmp_paths) {
    df = read.csv(path)
    df_filtered = subset(df,param==paste0('LV_',nvar_mplus) & paramHeader != 'Variances')
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

##### Summary of main effect results #####
# create one data frame
df_si_stats_full = as.data.frame(do.call(rbind, df_si_list))
# export data
write.csv(df_si_stats_full,paste0(out_dir,"Napls_LatentVariableRegress_OnlyCHR_UnstandardizedResultsSummary_MainEffBonf_FCHip.csv"), row.names = FALSE)
write_xlsx(df_si_stats_full,path = paste0(out_dir,"Napls_LatentVariableRegress_OnlyCHR_UnstandardizedResultsSummary_MainEffBonf_FCHip.xlsx"))

##### Summary of full results #####
# list all csv containing full test stats of hippocampal FC
path_files_hip = list.files(path = paste0(out_dir,"MplusOutputs"),
                            pattern = "^SingleIndicator_OnlyCHR_PathsUnstandardized_.*fchi\\.csv$",
                            full.names = TRUE)
# read them all and bind rowwise
df_path_hip = path_files_hip %>% lapply(read.csv) %>% bind_rows() %>% as.data.frame()
# replace 999 with NA
df_path_hip[df_path_hip == 999] = NA
# export data
write.csv(df_path_hip,file = paste0(out_dir,"Napls_LatentVariableRegress_OnlyCHR_UnstandardizedResultsSummary_FullStats_FCHip.csv"), row.names = FALSE)
write_xlsx(df_path_hip,path = paste0(out_dir,"Napls_LatentVariableRegress_OnlyCHR_UnstandardizedResultsSummary_FullStats_FCHip.xlsx"))


# SUMMARY OF STANDARDIZED RESULTS ----------------------------------------------------------------------------------------------------

##### Correction for multiple comparisons #####
# loop to extract the results for each MRI variable, perform a p-value correction across the six tests, and merge to one data frame 
df_si_list = list()
for (nvar in neurovars) {
  # shorten variable names of MRI variables for MPLUS and convert to upper case
  nvar_mplus = nvar %>% str_replace('fc','') %>% toupper()
  # specify paths to csv files with path statistics
  tmp_paths = list.files(path = paste0(out_dir,"MplusOutputs"),
                         pattern = sprintf("^SingleIndicator_OnlyCHR_PathsStandardized_.*_%s\\.csv$", nvar),
                         full.names = TRUE)
  # create empty list
  df_list = list()
  # import CSV files in loop
  for (path in tmp_paths) {
    df = read.csv(path)
    df_filtered = subset(df,param==paste0('LV_',nvar_mplus) & paramHeader != 'Variances')
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

##### Summary of main effect results #####
# create one data frame
df_si_stats_full = as.data.frame(do.call(rbind, df_si_list))
# export data
write.csv(df_si_stats_full,paste0(out_dir,"Napls_LatentVariableRegress_OnlyCHR_StandardizedResultsSummary_MainEffBonf_FCHip.csv"), row.names = FALSE)
write_xlsx(df_si_stats_full,path = paste0(out_dir,"Napls_LatentVariableRegress_OnlyCHR_StandardizedResultsSummary_MainEffBonf_FCHip.xlsx"))

##### Summary of full results #####
# list all csv containing full test stats of hippocampal FC
path_files_hip = list.files(path = paste0(out_dir,"MplusOutputs"),
                            pattern = "^SingleIndicator_OnlyCHR_PathsStandardized_.*fchi\\.csv$",
                            full.names = TRUE)
# read them all and bind rowwise
df_path_hip = path_files_hip %>% lapply(read.csv) %>% bind_rows() %>% as.data.frame()
# replace 999 with NA
df_path_hip[df_path_hip == 999] = NA
# export data
write.csv(df_path_hip,file = paste0(out_dir,"Napls_LatentVariableRegress_OnlyCHR_StandardizedResultsSummary_FullStats_FCHip.csv"), row.names = FALSE)
write_xlsx(df_path_hip,path = paste0(out_dir,"Napls_LatentVariableRegress_OnlyCHR_StandardizedResultsSummary_FullStats_FCHip.xlsx"))


# SUMMARY OF MODEL FITS --------------------------------------------------------------------------------------------------
# run loop to extract their full test statistics on model fit results
df_list = list()
for (bvar in behavvars) {
  for (nvar in neurovars) {
    p = list.files(path = paste0(out_dir,"MplusOutputs"),
                   pattern = sprintf("^SingleIndicator_OnlyCHR_Model_%s_%s\\.csv$",bvar,nvar),
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
write.csv(df_modelfit,file = paste0(out_dir,"Napls_LatentVariableRegress_Case-Control_ResultsSummary_ModelFit_FCHip.csv"), row.names = FALSE)
write_xlsx(df_modelfit, path = paste0(out_dir,"Napls_LatentVariableRegress_Case-Control_ResultsSummary_ModelFit_FCHip.xlsx"))

