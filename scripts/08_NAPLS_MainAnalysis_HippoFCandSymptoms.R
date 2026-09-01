# MAIN STATISTICAL ANALYSIS: CHANGE OF HIPPCAMPUS FC AND CHANGE OF CLINICAL VARIABLES ---------------------------------------------
# Description:  This script runs the main statistical analysis on the link between changes in hippocampal FC and clinical changes.
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
# remove subjects that have non-enhanced psychosis risk, as only 1 of these subjects has valid data for FC at two time points
df_napls_mlr_imp = subset(df_napls_mlr_imp,group !=2)
# rename group column to CHR
df_napls_mlr_imp = df_napls_mlr_imp %>% rename(CHR = group)
# create vectors with behavioral variables
behavvars = c("PosS","NegS","DepS","gaf","SyCod", "VeMem")
# create vector for hippocampal FC
neurovars = c("fchi")

# LATENT VARIABLE REGRESSION: CASE-CONTROL COMPARISON ----------------------------------------------------------------------------------------------------
# Note: This code runs the latent variable regression models in MPLUS comparing the link between hippocampal FC
# change and clinical severity change over time between CHR and healthy controls. The latent variable regression is used to account for
# the measurement error from the change estimates. The accuracy of the estimation of linear change over time per subject in each outcome
# differs depending on the available time points per subject in each outcome.

##### Latent variable regression loop running MPLUS #####
# Note: in this analysis cases and controls are used that have a slope in the respective behavioral AND the MRI variable.
# re-define and set working directory in which the MPLUS .exe file is stored -> necessary because the path of the working directory
# must not be too long for MPLUS
wd = "PATH/TO/MPLUS/DIRECTORY/"
setwd(wd)
# run latent variable regression in loop
for (bvar in behavvars) {
  for (nvar in neurovars) {
    print(paste0("Running Single Indicator Approach for ",bvar, " and ",nvar,"..."))
    # shorten variable names of MRI variables for MPLUS
    nvar_mplus = nvar %>% str_replace("fc","")
    # create data frame with variables of interest
    df = df_napls_mlr_imp[c(paste0(bvar,"_CH"),paste0(nvar,"_CH"),"CHR",paste0(bvar,"_SE"),paste0(nvar,"_SE"),
                            "age","female","site02","site03","site04","site05","site06","site07","site08","site09")]
    # remove subjects with any missings
    df = df[complete.cases(df[c(paste0(bvar, "_CH"), paste0(nvar, "_CH"))]), ]
    # get average error variances across group from standard errors of both variables
    se_bvar = mean(df[[paste0(bvar,"_SE")]],na.rm=T)^2
    se_nvar = mean(df[[paste0(nvar,"_SE")]],na.rm=T)^2
    # drop standard error columns again
    df = df[ , !grepl("_SE$", names(df)) ]
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
  
    Analysis:   type = random;
                ESTIMATOR = MLR;
                ALGORITHM=INTEGRATION;
                STARTS=20;
                STITERATION=50000;
    
    Model:
    
    LV_", nvar_mplus, " by ", nvar_mplus, "_CH;
    ", nvar_mplus, "_CH@", format(se_nvar, scientific = FALSE), ";
    
    LV_", bvar, " by ", bvar, "_CH;
    ", bvar, "_CH@", format(se_bvar, scientific = FALSE), ";
    
    ", nvar_mplus, "xCHR | LV_", nvar_mplus, " xwith CHR;
    
    LV_", bvar, " on LV_", nvar_mplus, " (b1)
               CHR
               ", nvar_mplus, "xCHR (b2)
               age female site02 site03 site04
               site05 site06 site07 site08
               site09;
    
    LV_", nvar_mplus, " WITH CHR age female site02 site03 site04
                       site05 site06 site07 site08
                       site09;
    
    CHR WITH age female site02 site03 site04
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
    
    Model Constraint:
        NEW(effect_ctrl effect_pat);
        effect_ctrl = b1;
        effect_pat  = b1 + b2;
    
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
    # read the output
    model_out = readModels(outfile)
    # extract unstandardized estimates
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
    write.csv(df_res, file = paste0("SingleIndicator_Case-Control_PathsUnstandardized_",bvar,"_",nvar, ".csv"), row.names = FALSE)
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
    write.csv(df_res, file = paste0("SingleIndicator_Case-Control_PathsStandardized_",bvar,"_",nvar, ".csv"), row.names = FALSE)
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
    write.csv(df_model, file = paste0("SingleIndicator_Case-Control_Model_", bvar, "_", nvar, ".csv"), row.names = FALSE)
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
df_res_list = list()
for (nvar in neurovars) {
  # specify paths to csv files with path statistics
  tmp_paths = list.files(path = paste0(out_dir,"MplusOutputs"),
                         pattern = sprintf("^SingleIndicator_Case-Control_PathsUnstandardized_.*_%s\\.csv$", nvar),
                         full.names = TRUE)
  # create empty list
  df_list = list()
  # import CSV files in loop
  for (path in tmp_paths) {
    df = read.csv(path)
    df_filtered = df[grepl("XCHR", df$param), ]
    df_list[[length(df_list) + 1]] = df_filtered
  }
  # create one data frame
  df_res_stats = do.call(rbind, df_list)
  # insert adjusted p-values and mark significance using different methods
  df_res_stats$p_fdr = p.adjust(df_res_stats$pval_exact,method="fdr")
  df_res_stats$p_fdr_sig = df_res_stats$p_fdr < 0.05
  df_res_stats$p_bonf = p.adjust(df_res_stats$pval_exact,method="bonferroni")
  df_res_stats$p_bonf_sig = df_res_stats$p_bonf < 0.05
  # add data frame to list
  df_res_list[[length(df_res_list) + 1]] = df_res_stats
}

##### Summary of interaction results #####
# create one data frame
df_res_stats_full = data.frame(do.call(rbind, df_res_list))
# export data
write.csv(df_res_stats_full,paste0(out_dir,"Napls_LatentVariableRegress_Case-Control_UnstandardizedResultsSummary_InteractionsBonf_FCHip.csv"), row.names = FALSE)
write_xlsx(df_res_stats_full,path = paste0(out_dir,"Napls_LatentVariableRegress_Case-Control_UnstandardizedResultsSummary_InteractionsBonf_FCHip.xlsx"))

##### Summary of full results #####
# list all csv containing full test stats of hippocampal FC
path_files_hip = list.files(path = paste0(out_dir,"MplusOutputs"),
                            pattern = "^SingleIndicator_Case-Control_PathsUnstandardized_.*fchi\\.csv$",
                            full.names = TRUE)
# read them all and bind rowwise
df_path_hip = path_files_hip %>% lapply(read.csv) %>% bind_rows() %>% as.data.frame()
# replace 999 with NA
df_path_hip[df_path_hip == 999] = NA
# export data
write.csv(df_path_hip,file = paste0(out_dir,"Napls_LatentVariableRegress_Case-Control_UnstandardizedResultsSummary_FullStats_FCHip.csv"), row.names = FALSE)
write_xlsx(df_path_hip,path = paste0(out_dir,"Napls_LatentVariableRegress_Case-Control_UnstandardizedResultsSummary_FullStats_FCHip.xlsx"))

##### Summary of effects within each group #####
# read full test stats
df_path_hip = read.csv(paste0(out_dir,"Napls_LatentVariableRegress_Case-Control_UnstandardizedResultsSummary_FullStats_FCHip.csv"))
# define behavioral variables and how they're named in MPLUS
behavvars_mplus = c("LV_POSS.ON","LV_NEGS.ON","LV_DEPS.ON","LV_GAF.ON","LV_SYCOD.ON","LV_VEMEM.ON")
# create empty list
df_estimates_list = list()
# run loop
for (bvar in behavvars_mplus) {
  # extract estimate of predictor LV_HI which reflects the effect of the controls
  b_fc = subset(df_path_hip,paramHeader==bvar&param=="LV_HI")$est
  # extract the interaction effects
  b_interaction = subset(df_path_hip,paramHeader==bvar&param=="HIXCHR")$est
  # extract the effect in patients which is the sum of the effect in controls and the interaction effect
  b_patients = b_fc + b_interaction
  # create data frame with these estimates
  df_est = data.frame(dv = bvar,est_hc = b_fc,est_chr = b_patients,est_interaction = b_interaction)
  # add to list
  df_estimates_list[[bvar]] = df_est
}
# merge data frames
df_estimates = as.data.frame(do.call(rbind, df_estimates_list))
# export data
write.csv(df_estimates,paste0(out_dir,"Napls_LatentVariableRegress_Case-Control_UnstandardizedResultsSummary_EstimatesPerGroup_FCHip.csv"), row.names = FALSE)
write_xlsx(df_estimates,path = paste0(out_dir,"Napls_LatentVariableRegress_Case-Control_UnstandardizedResultsSummary_EstimatesPerGroup_FCHip.xlsx"))


# SUMMARY OF STANDARDIZED RESULTS ----------------------------------------------------------------------------------------------------

##### Correction for multiple comparisons #####
# loop to extract the results for each MRI variable, perform a p-value correction across the six tests, and merge to one data frame 
df_res_list = list()
for (nvar in neurovars) {
  # specify paths to csv files with path statistics
  tmp_paths = list.files(path = paste0(out_dir,"MplusOutputs"),
                         pattern = sprintf("^SingleIndicator_Case-Control_PathsStandardized_.*_%s\\.csv$", nvar),
                         full.names = TRUE)
  # create empty list
  df_list = list()
  # import CSV files in loop
  for (path in tmp_paths) {
    df = read.csv(path)
    df_filtered = df[grepl("XCHR", df$param), ]
    df_list[[length(df_list) + 1]] = df_filtered
  }
  # create one data frame
  df_res_stats = do.call(rbind, df_list)
  # insert adjusted p-values and mark significance using different methods
  df_res_stats$p_fdr = p.adjust(df_res_stats$pval_exact,method="fdr")
  df_res_stats$p_fdr_sig = df_res_stats$p_fdr < 0.05
  df_res_stats$p_bonf = p.adjust(df_res_stats$pval_exact,method="bonferroni")
  df_res_stats$p_bonf_sig = df_res_stats$p_bonf < 0.05
  # add data frame to list
  df_res_list[[length(df_res_list) + 1]] = df_res_stats
}

##### Summary of interaction results #####
# create one data frame
df_res_stats_full = as.data.frame(do.call(rbind, df_res_list))
# export data
write.csv(df_res_stats_full,paste0(out_dir,"Napls_LatentVariableRegress_Case-Control_StandardizedResultsSummary_InteractionsBonf_FCHip.csv"), row.names = FALSE)
write_xlsx(df_res_stats_full,path = paste0(out_dir,"Napls_LatentVariableRegress_Case-Control_StandardizedResultsSummary_InteractionsBonf_FCHip.xlsx"))

##### Summary of full results #####
# list all csv containing full test stats of hippocampal FC
path_files_hip = list.files(path = paste0(out_dir,"MplusOutputs"),
                            pattern = "^SingleIndicator_Case-Control_PathsStandardized_.*fchi\\.csv$",
                            full.names = TRUE)
# read them all and bind rowwise
df_path_hip = path_files_hip %>% lapply(read.csv) %>% bind_rows() %>% as.data.frame()
# replace 999 with NA
df_path_hip[df_path_hip == 999] = NA
# export data
write.csv(df_path_hip,file = paste0(out_dir,"Napls_LatentVariableRegress_Case-Control_StandardizedResultsSummary_FullStats_FCHip.csv"), row.names = FALSE)
write_xlsx(df_path_hip,path = paste0(out_dir,"Napls_LatentVariableRegress_Case-Control_StandardizedResultsSummary_FullStats_FCHip.xlsx"))

##### Summary of effects within each group #####
# read full test stats
df_path_hip = read.csv(paste0(out_dir,"Napls_LatentVariableRegress_Case-Control_StandardizedResultsSummary_FullStats_FCHip.csv"))
# define behavioral variables and how they're named in MPLUS
behavvars_mplus = c("LV_POSS.ON","LV_NEGS.ON","LV_DEPS.ON","LV_GAF.ON","LV_SYCOD.ON","LV_VEMEM.ON")
# create empty list
df_estimates_list = list()
# run loop
for (bvar in behavvars_mplus) {
  # extract estimate of predictor LV_HI which reflects the effect of the controls
  b_fc = subset(df_path_hip,paramHeader==bvar&param=="LV_HI")$est
  # extract the interaction effects
  b_interaction = subset(df_path_hip,paramHeader==bvar&param=="HIXCHR")$est
  # extract the effect in patients which is the sum of the effect in controls and the interaction effect
  b_patients = b_fc + b_interaction
  # create data frame with these estimates
  df_est = data.frame(dv = bvar,est_hc = b_fc,est_chr = b_patients,est_interaction = b_interaction)
  # add to list
  df_estimates_list[[bvar]] = df_est
}
# merge data frames
df_estimates = as.data.frame(do.call(rbind, df_estimates_list))
# export data
write.csv(df_estimates,paste0(out_dir,"Napls_LatentVariableRegress_Case-Control_StandardizedResultsSummary_EstimatesPerGroup_FCHip.csv"), row.names = FALSE)
write_xlsx(df_estimates,path = paste0(out_dir,"Napls_LatentVariableRegress_Case-Control_StandardizedResultsSummary_EstimatesPerGroup_FCHip.xlsx"))


# SUMMARY OF MODEL FITS ----------------------------------------------------------------------------------------------------
# run loop to extract their full test statistics on model fit results
df_list = list()
for (bvar in behavvars) {
  for (nvar in neurovars) {
    p = list.files(path = paste0(out_dir,"MplusOutputs"),
                   pattern = sprintf("^SingleIndicator_Case-Control_Model_%s_%s\\.csv$",bvar,nvar),
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


# VISUALIZATION OF RESULTS ----------------------------------------------------------------------------------------

##### Descriptive scatter plots #####
# select variables of interest
df_hippo = df_napls_mlr_imp[c("NegS_CH","DepS_CH","gaf_CH","fchi_CH","CHR")]
# rename CHR column
df_hippo = df_hippo %>% rename(group = CHR)
# convert CHR column
df_hippo$group = as.character(df_hippo$group) %>% str_replace("1","CHR-P") %>% str_replace("0","HC")
# create list of columns to be plotted
vars_plot_clinical = c("NegS_CH","DepS_CH","gaf_CH")
# plot
plot_list = list()
for (i in seq_along(vars_plot_clinical)) {
  clin = vars_plot_clinical[i]
  plot_list[[clin]] = local({
    clin = clin
    clin_axis = clin %>%
      str_replace("PosS_CH", "Change in Attenuated Positive Symp") %>%
      str_replace("NegS_CH", "Change in Negative Symp") %>%
      str_replace("VeMem_CH", "Change in Verbal Memory") %>%
      str_replace("SyCod_CH", "Change in Symbol Coding") %>%
      str_replace("DepS_CH", "Change in Depressive Symp") %>%
      str_replace("gaf_CH", "Change in Psych Func")
    if (clin == "NegS_CH" | clin == "DepS_CH") {
      col_new = "red4"
    } else {
      col_new = "dodgerblue4"
    }
    p = ggplot(df_hippo, aes(x = unlist(df_hippo[,"fchi_CH"]), y = unlist(df_hippo[, clin]))) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
      geom_point(color = "grey60", alpha = 0.2, size = 2.5) +
      geom_smooth(aes(color=group),method = "lm") +
      scale_color_manual(values = c(col_new, "grey40")) +
      scale_x_continuous(limits = c(-0.2, 0.2),breaks = c(-0.1, 0, 0.1)) +
      scale_y_continuous(limits = c(-10, 10),breaks = c(-10, -5, 0, 5, 10)) +
      labs(x = "Change in Hippocampal FC",y = clin_axis) +
      theme_classic() +
      facet_wrap(~group) +
      theme(axis.title.x = element_text(size = 14, color = "black", face = "bold"),
            axis.text.x  = element_text(size = 12, color = "black"),
            axis.title.y = element_text(size = 14, color = "black", face = "bold"),
            axis.text.y  = element_text(size = 12, color = "black"),
            strip.text   = element_text(size = 14, color = "black", face = "bold"),
            legend.position = "none")
  })
}
# arrange plots
plots = grid.arrange(grobs = plot_list, nrow = 1)
# export plot as jpeg
ggsave(paste0(out_dir,"Napls_LatentVariableRegress_Case-Control_ScatterPlot_HippoFC.jpeg"),
       plot=plots,height=100,width=300,unit="mm",dpi=500)


##### Forest Plots #####
# get test statistics with interactions
df_betas = read.csv(paste0(out_dir,"Napls_LatentVariableRegress_Case-Control_StandardizedResultsSummary_InteractionsBonf_FCHip.csv"))
# rename clinical variables
df_betas$paramHeader = df_betas$paramHeader %>% str_replace("LV_","") %>% str_replace(".ON","") %>%
  str_replace("DEPS","Depressive Symp Change") %>% str_replace("NEGS","Negative Symp Change") %>%
  str_replace("POSS","Positive Symp Change") %>% str_replace("GAF","Psych Func Change") %>%
  str_replace("SYCOD","Symbol Coding Change") %>% str_replace("VEMEM","Verbal Memory Change")
# change factor levels
df_betas$paramHeader = factor(df_betas$paramHeader, levels = c("Symbol Coding Change","Verbal Memory Change",
                                                               "Positive Symp Change","Depressive Symp Change",
                                                               "Psych Func Change","Negative Symp Change"))
# add column indicating strength of significance after bonferroni correction
df_betas = df_betas %>%
  mutate(stars = case_when(
    p_bonf < 0.001 ~ "***",
    p_bonf < 0.01  ~ "**",
    p_bonf < 0.05  ~ "*",
    TRUE           ~ ""
  ),
  label = ifelse(p_bonf_sig, paste0(round(est,digits=2), stars), ""))
# forest plot
p = ggplot(df_betas, aes(x = est, y = paramHeader)) +
  geom_errorbarh(aes(xmin = est - se, xmax = est + se, color = est), height = 0.15, linewidth = 0.9) +
  geom_point(aes(color = est), size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_gradient2(low = "red4", mid = "white", high = "dodgerblue4", midpoint = 0, na.value = "grey80") +
  geom_text(aes(label = label), vjust = -1, size = 2, fontface = "bold", color = "black", na.rm = TRUE) +
  scale_x_continuous(limits = c(-0.25, 0.25), breaks = c(-0.2, -0.1, 0, 0.1, 0.2)) +
  theme_classic() +
  labs(x = "Beta", y = "") +
  theme(
    axis.text.x = element_text(size = 6, color = "black"),
    axis.text.y = element_text(size = 7, color = "black", face = "bold"),
    axis.title.x = element_text(size = 7, face = "bold"),
    legend.position = "none")
# export plot as jpeg
ggsave(paste0(out_dir,"Napls_LatentVariableRegress_Case-Control_ForestPlot_HippoFC_v6.jpeg"),
       plot=p,height=50,width=140,unit="mm",dpi=500)


# SLOPE CHARACTERISTICS PER GROUP -----------------------------------------------------------------------
library(e1071)

# function to compute summary stats including skewness and kurtosis
get_stats = function(df, var_col, bvar_name, var_label) {
  stats_df = df |>
    dplyr::group_by(CHR) |>
    dplyr::summarise(
      bvar = bvar_name,
      variable = var_label,
      n = dplyr::n(),
      mean = mean(.data[[var_col]], na.rm = TRUE),
      sd = stats::sd(.data[[var_col]], na.rm = TRUE),
      median = median(.data[[var_col]], na.rm = TRUE),
      min = min(.data[[var_col]], na.rm = TRUE),
      max = max(.data[[var_col]], na.rm = TRUE),
      skewness = e1071::skewness(.data[[var_col]], type = 2, na.rm = TRUE),
      kurtosis = e1071::kurtosis(.data[[var_col]], type = 2, na.rm = TRUE),
      .groups = "drop"
    )
  return(stats_df)
}

bvars = c("NegS","PosS","DepS","gaf","VeMem","SyCod")
nvar = "fchi"

results = data.frame()

for (bvar in bvars) {
  df = df_napls_mlr_imp[
    c(paste0(bvar,"_CH"),
      paste0(nvar,"_CH"),
      "CHR",
      "age","female",
      "site02","site03","site04","site05",
      "site06","site07","site08","site09")
  ]
  
  # keep only complete cases for the two main variables
  df = df[complete.cases(df[c(paste0(bvar, "_CH"),
                              paste0(nvar, "_CH"))]), ]
  
  # stats for behavioral variable
  bvar_stats = get_stats(df,
                         paste0(bvar, "_CH"),
                         bvar,
                         paste0(bvar, "_CH"))
  
  results = dplyr::bind_rows(results, bvar_stats)
  
  # for NegS also compute fchi stats on same sample
  if (bvar == "NegS") {
    fchi_stats = get_stats(df,
                           paste0(nvar, "_CH"),
                           bvar,
                           paste0(nvar, "_CH"))
    
    results = dplyr::bind_rows(results, fchi_stats)
  }
}

# export
write.csv(results,
          paste0(out_dir, "SlopeStatsByGroup.csv"),
          row.names = FALSE)



plot_data = data.frame()  # empty container to accumulate raw data for plotting

for (bvar in bvars) {  # loop through each bvar
  # build data frame with variables of interest for this bvar
  df = df_napls_mlr_imp[c(paste0(bvar,"_CH"),paste0(nvar,"_CH"),"CHR",
                          "age","female","site02","site03","site04","site05","site06","site07","site08","site09")]
  # keep only subjects with complete data on the bvar change score and the fchi change score
  df = df[complete.cases(df[c(paste0(bvar, "_CH"), paste0(nvar, "_CH"))]), ]
  
  # pull out the bvar change score column and tag it with variable/bvar labels for later facetting
  bvar_long = data.frame(
    CHR = df$CHR,  # group variable for x-axis
    value = df[[paste0(bvar, "_CH")]],  # the change score values for y-axis
    variable = paste0(bvar, "_CH"),  # label identifying which change score this is
    bvar = bvar  # label identifying which loop iteration this came from
  )
  plot_data = dplyr::bind_rows(plot_data, bvar_long)  # add to master plotting data frame
  
  # for NegS only, also pull out the fchi change score (same sample)
  if (bvar == "NegS") {
    fchi_long = data.frame(
      CHR = df$CHR,  # group variable for x-axis
      value = df[[paste0(nvar, "_CH")]],  # the fchi change score values for y-axis
      variable = paste0(nvar, "_CH"),  # label identifying this as the fchi change score
      bvar = bvar  # label identifying which loop iteration this came from
    )
    plot_data = dplyr::bind_rows(plot_data, fchi_long)  # add to master plotting data frame
  }
}

# make the plot: boxplot per CHR group, faceted by variable, with individual points jittered on top
p = ggplot(plot_data, ggplot2::aes(x = as.factor(CHR), y = value, color = as.factor(CHR))) +  # x = group, y = change score, color by group
  geom_boxplot(outlier.shape = NA) +  # boxplot without showing outlier points twice (jitter will show all points instead)
  geom_jitter(width = 0.2, alpha = 0.5) +  # add individual data points, spread horizontally so they don't overlap
  facet_wrap(~ variable, scales = "free_y") +  # one panel per variable, each with its own y-axis scale
  theme_classic() +  # clean plot theme with no gridlines
  labs(x = "CHR group", y = "Change score", color = "CHR group")  # axis and legend labels

ggsave(
  filename = paste0(out_dir,"SlopeStatsByGroup_Boxplot.png"),  # name of the output file
  plot = p,  # the plot object to export
  width = 10,  # width of the image in inches
  height = 7,  # height of the image in inches
  dpi = 500  # resolution of the image
)

# make the plot: histogram of change scores, rows = CHR group (HC vs patient), columns = variable
p_hist = ggplot(plot_data, ggplot2::aes(x = value, fill = as.factor(CHR))) +  # x = change score, fill color by group
  geom_histogram(position = "identity", alpha = 0.4, bins = 50) +  # overlapping bars, semi-transparent, 20 bins
  facet_grid(CHR ~ variable, scales = "free") +  # rows split by CHR group, columns split by variable, free scales on both axes
  theme_classic() +  # clean plot theme with no gridlines
  labs(x = "Change score", y = "Count", fill = "CHR group")  # axis and legend labels

print(p_hist)  # display the histogram

ggsave(
  filename = paste0(out_dir,"SlopeStatsByGroup_Histogram.png"),  # name of the output file
  plot = p_hist,  # the plot object to export
  width = 10,  # width of the image in inches
  height = 5,  # height of the image in inches
  dpi = 500  # resolution of the image
)


# PROPORTION OF SLOPE COMBINATIONS ----------------------------------------------------------------------------------------------------
# get group values
chrs <- sort(unique(df_napls_mlr_imp$CHR))
# 1) NegS_CH vs fchi_CH
df1 <- df_napls_mlr_imp %>%
  filter(!is.na(NegS_CH) & !is.na(fchi_CH)) %>%
  mutate(
    clinical = ifelse(NegS_CH > 0, "pos", "neg"),
    fc       = ifelse(fchi_CH > 0, "pos", "neg")
  ) %>%
  group_by(CHR, clinical, fc) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(CHR, clinical, fc) %>%
  mutate(pair = "NegS_CH_vs_fchi_CH")

# 2) gaf_CH vs fchi_CH
df2 <- df_napls_mlr_imp %>%
  filter(!is.na(gaf_CH) & !is.na(fchi_CH)) %>%
  mutate(
    clinical = ifelse(gaf_CH > 0, "pos", "neg"),
    fc       = ifelse(fchi_CH > 0, "pos", "neg")
  ) %>%
  group_by(CHR, clinical, fc) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(CHR, clinical, fc) %>%
  mutate(pair = "gaf_CH_vs_fchi_CH")

# 3) DepS_CH vs fchi_CH
df3 <- df_napls_mlr_imp %>%
  filter(!is.na(DepS_CH) & !is.na(fchi_CH)) %>%
  mutate(
    clinical = ifelse(DepS_CH > 0, "pos", "neg"),
    fc       = ifelse(fchi_CH > 0, "pos", "neg")
  ) %>%
  group_by(CHR, clinical, fc) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(CHR, clinical, fc) %>%
  mutate(pair = "DepS_CH_vs_fchi_CH")
# Combine all
df_out = bind_rows(df1, df2, df3)
df_out <- df_out %>% filter(!is.na(clinical) & !is.na(fc))
# Export
write.csv(df_out, file = paste0(out_dir,"SlopeProportionsPerGroup.csv"), row.names = FALSE)










