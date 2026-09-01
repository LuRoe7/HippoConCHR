# ROBUSTNESS CHECK: CHANGE OF FC ACROSS BRAIN REGIONS AND CHANGE OF CLINICAL VARIABLES ---------------------------------------------
# Description:  This script runs the robustness check on the link between changes in FC across the brain and clinical changes.
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
in_dir = "results/napls/P1_HippoFC/07_NAPLS_SingleSubjectRegressions/"
# define output directory
out_dir = "results/napls/P1_HippoFC/10_NAPLS_RobustnessCheck_OtherBrainAreas/"
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
neurovars = c("fcsal","fcdmn","fcfpn","fclim","fcvis","fcsmn","fcamy","fcbg","fctha","fcpal","fcput","fccau","fcacu")


# LATENT VARIABLE REGRESSION: CASE-CONTROL COMPARISON ----------------------------------------------------------------------------------------------------
# Note: This code runs the latent variable regression models in MPLUS comparing the link between hippocampal FC
# change and clinical severity change over time between CHR and healthy controls. The latent variable regression is used to account for
# the measurement error from the change estimates. The accuracy of the estimation of linear change over time per subject in each outcome
# differs depending on the available time points per subject in each outcome.

##### Latent variable regression loop running MPLUS #####
# Note: in this analysis cases and controls are used that have a slope in the respective behavioral AND the MRI variable.
# re-define and set working directory in which the MPLUS .exe file is stored -> necessary because the path of the working directory
# must not be too long for MPLUS
wd = "/home/lukas/mplus"
setwd(wd)
# initialize error log
error_log_file = file.path(wd, "error_log.txt")
writeLines("Error log for latent variable regression\n", con = error_log_file)
# run latent variable regression in loop
for (bvar in behavvars) {
  for (nvar in neurovars) {
    tryCatch({
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
      
      # define the MPLUS input file
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
      writeLines(inp_text, con = inpfile)
      
      # run MPLUS
      system(paste("wine Mplus-8.6.exe", basename(inpfile)))
      
      # read the output
      outfile = file.path(wd, paste0("naplsdata_", tolower(nvar), "_", tolower(bvar), ".out"))
      model_out = readModels(outfile)
      
      # extract unstandardized estimates
      df_res = as.data.frame(model_out$parameters$unstandardized)
      df_res$zval = NA
      df_res$pval_exact = NA
      valid = df_res$se > 0 & df_res$pval != 999
      df_res$zval[valid] = with(df_res[valid, ], est / se)
      df_res$pval_exact[valid] = with(df_res[valid, ], 2 * (1 - pnorm(abs(zval))))
      df_res$p_sig = df_res$pval_exact < 0.05
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
        } else { NA }
      )
      write.csv(df_model, file = paste0("SingleIndicator_Case-Control_Model_", bvar, "_", nvar, ".csv"), row.names = FALSE)
      
    }, error = function(e) {
      # log the error
      msg = paste(Sys.time(), "Error for", bvar, "and", nvar, ":", e$message, "\n")
      cat(msg, file = error_log_file, append = TRUE)
      message("Error occurred for ", bvar, " and ", nvar, ". Check error_log.txt.")
    })
  }
}
# define and set old working directory
old_wd = "/home/lukas/Desktop/LukasLinux/Projects/LongitMechanisms/"
setwd(old_wd)
# create sub-directory in output folder
dir.create(paste0(out_dir,"MplusOutputs"), recursive = TRUE,showWarnings=F)
# specify paths to output files
res_files = list.files(path = wd,pattern = "\\.(csv|dat|inp|out|txt)$",recursive = TRUE, full.names = TRUE)
# copy to output sub-directory in output directory
copied = file.copy(res_files,paste0(out_dir,"MplusOutputs"), overwrite = TRUE)
# remove copied files from intermediate working directory
file.remove(res_files[copied])


# SUMMARY OF UNSTANDARDIZED RESULTS ----------------------------------------------------------------------------------------------------

##### Correction for multiple comparisons #####
# loop to extract the results for each MRI variable, perform a p-value correction across the six tests, and merge to one data frame 
df_si_list = list()
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
  df_si_stats = do.call(rbind, df_list)
  # insert adjusted p-values and mark significance using different methods
  df_si_stats$p_fdr = p.adjust(df_si_stats$pval_exact,method="fdr")
  df_si_stats$p_fdr_sig = df_si_stats$p_fdr < 0.05
  df_si_stats$p_bonf = p.adjust(df_si_stats$pval_exact,method="bonferroni")
  df_si_stats$p_bonf_sig = df_si_stats$p_bonf < 0.05
  # add data frame to list
  df_si_list[[length(df_si_list) + 1]] = df_si_stats
}

##### Summary of interaction results #####
# create one data frame
df_si_stats_full = as.data.frame(do.call(rbind, df_si_list))
# export data
write.csv(df_si_stats_full,paste0(out_dir,"Napls_LatentVariableRegress_Case-Control_UnstandardizedResultsSummary_InteractionsBonf_AllNetw.csv"), row.names = FALSE)
write_xlsx(df_si_stats_full,path = paste0(out_dir,"Napls_LatentVariableRegress_Case-Control_UnstandardizedResultsSummary_InteractionsBonf_AllNetw.xlsx"))

##### Summary of full results #####
## FC in netwpocampus (main analysis)
# list all csv containing full test stats of netwpocampal FC
path_files_netw = list.files(path = paste0(out_dir,"MplusOutputs"),
                            pattern = "^SingleIndicator_Case-Control_PathsUnstandardized_.*\\.csv$",
                            full.names = TRUE)
# read them all and bind rowwise
df_path_netw = path_files_netw %>% lapply(read.csv) %>% bind_rows() %>% as.data.frame()
# replace 999 with NA
df_path_netw[df_path_netw == 999] = NA
# export data
write.csv(df_path_netw,file = paste0(out_dir,"Napls_LatentVariableRegress_Case-Control_UnstandardizedResultsSummary_FullStats_AllNetw.csv"), row.names = FALSE)
write_xlsx(df_path_netw,path = paste0(out_dir,"Napls_LatentVariableRegress_Case-Control_UnstandardizedResultsSummary_FullStats_AllNetw.xlsx"))


# SUMMARY OF STANDARDIZED RESULTS ----------------------------------------------------------------------------------------------------

##### Correction for multiple comparisons #####
# loop to extract the results for each MRI variable, perform a p-value correction across the six tests, and merge to one data frame 
df_si_list = list()
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
  df_si_stats = do.call(rbind, df_list)
  # insert adjusted p-values and mark significance using different methods
  df_si_stats$p_fdr = p.adjust(df_si_stats$pval_exact,method="fdr")
  df_si_stats$p_fdr_sig = df_si_stats$p_fdr < 0.05
  df_si_stats$p_bonf = p.adjust(df_si_stats$pval_exact,method="bonferroni")
  df_si_stats$p_bonf_sig = df_si_stats$p_bonf < 0.05
  # add data frame to list
  df_si_list[[length(df_si_list) + 1]] = df_si_stats
}

##### Summary of interaction results #####
# create one data frame
df_si_stats_full = as.data.frame(do.call(rbind, df_si_list))
# export data
write.csv(df_si_stats_full,paste0(out_dir,"Napls_LatentVariableRegress_Case-Control_StandardizedResultsSummary_InteractionsBonf_AllNetw.csv"), row.names = FALSE)
write_xlsx(df_si_stats_full,path = paste0(out_dir,"Napls_LatentVariableRegress_Case-Control_StandardizedResultsSummary_InteractionsBonf_AllNetw.xlsx"))

##### Summary of full results #####
## FC in netwpocampus (main analysis)
# list all csv containing full test stats of netwpocampal FC
path_files_netw = list.files(path = paste0(out_dir,"MplusOutputs"),
                            pattern = "^SingleIndicator_Case-Control_PathsStandardized_.*\\.csv$",
                            full.names = TRUE)
# read them all and bind rowwise
df_path_netw = path_files_netw %>% lapply(read.csv) %>% bind_rows() %>% as.data.frame()
# replace 999 with NA
df_path_netw[df_path_netw == 999] = NA
# export data
write.csv(df_path_netw,file = paste0(out_dir,"Napls_LatentVariableRegress_Case-Control_StandardizedResultsSummary_FullStats_AllNetw.csv"), row.names = FALSE)
write_xlsx(df_path_netw,path = paste0(out_dir,"Napls_LatentVariableRegress_Case-Control_StandardizedResultsSummary_FullStats_AllNetw.xlsx"))


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
write.csv(df_modelfit,file = paste0(out_dir,"Napls_LatentVariableRegress_Case-Control_ResultsSummary_ModelFit_AllNetw.csv"), row.names = FALSE)
write_xlsx(df_modelfit, path = paste0(out_dir,"Napls_LatentVariableRegress_Case-Control_ResultsSummary_ModelFit_AllNetw.xlsx"))


# VISUALIZATION OF RESULTS ----------------------------------------------------------------------------------------
# read test stats on group * FC interactions for each network
df_si_stats_netw = read.csv(paste0(out_dir,"Napls_LatentVariableRegress_Case-Control_StandardizedResultsSummary_InteractionsBonf_AllNetw.csv"))
# define directory and read stats of main analysis on hippocampus
out_dir_hip = "results/napls/P1_HippoFC/08_NAPLS_MainAnalysis_HippoFCandSymptoms/"
df_si_stats_hip = read.csv(paste0(out_dir_hip,"Napls_LatentVariableRegress_Case-Control_StandardizedResultsSummary_InteractionsBonf_FCHip.csv"))
# define directory and read stats of sensitivity analysis on FC patterns within hippocampus
out_dir_hipw = "results/napls/P1_HippoFC/09_NAPLS_SensitivityAnalysis_FCwithinHippo/"
df_si_stats_hipw = read.csv(paste0(out_dir_hipw,"Napls_LatentVariableRegress_Case-Control_StandardizedResultsSummary_InteractionsBonf_WithinHip.csv"))
# merge stat data frames row-wise
df_betas = as.data.frame(rbind(df_si_stats_hip,df_si_stats_hipw,df_si_stats_netw))
# adjust variable names
df_betas$param = df_betas$param %>%
  str_replace("XCHR", "") %>% str_replace("HIH", "H Hippocampal FC Change") %>% str_replace("HIL", "L Hippocampal FC Change") %>%
  str_replace("HIR", "R Hippocampal FC Change") %>% str_replace("HI", "Hippocampal FC Change") %>% 
  str_replace("SAL", "SAL FC Change") %>% str_replace("DMN", "DMN FC Change") %>% str_replace("FPN", "FPN FC Change") %>%
  str_replace("LIM", "LIM FC Change") %>% str_replace("SMN", "SMN FC Change") %>% str_replace("VIS", "VIS FC Change") %>%
  str_replace("THA", "Thalamus FC Change") %>% str_replace("AMY", "Amygdala FC Change") %>%
  str_replace("BG",  "Basal Ganglia FC Change") %>% str_replace("PUT", "Putamen FC Change") %>%
  str_replace("PAL", "Pallidum FC Change") %>% str_replace("CAU", "Caudate FC Change") %>%
  str_replace("ACU", "Accumbens FC Change")
df_betas$paramHeader = df_betas$paramHeader %>% str_replace("LV_","") %>% str_replace(".ON","") %>%
  str_replace("DEPS","Depressive Symp Change") %>% str_replace("NEGS","Negative Symp Change") %>%
  str_replace("POSS","Positive Symp Change") %>% str_replace("GAF","Psychosocial Functioning Change") %>%
  str_replace("SYCOD","Symbol Coding Change") %>% str_replace("VEMEM","Verbal Memory Change")
# change factor levels
df_betas$param = factor(df_betas$param, levels = c("Hippocampal FC Change","R Hippocampal FC Change","L Hippocampal FC Change",
                                                   "H Hippocampal FC Change","SAL FC Change","DMN FC Change","FPN FC Change",
                                                   "LIM FC Change","SMN FC Change","VIS FC Change","Thalamus FC Change",
                                                   "Amygdala FC Change", "Basal Ganglia FC Change","Putamen FC Change",
                                                   "Pallidum FC Change", "Caudate FC Change","Accumbens FC Change"))
df_betas$paramHeader = factor(df_betas$paramHeader, levels = c("Symbol Coding Change","Verbal Memory Change","Positive Symp Change",
                                                               "Depressive Symp Change","Psychosocial Functioning Change",
                                                               "Negative Symp Change"))
# add column indicating strength of significance after bonferroni correction
df_betas = df_betas %>% mutate(stars = case_when(p_bonf < 0.001 ~ "***",p_bonf < 0.01  ~ "**",p_bonf < 0.05  ~ "*",TRUE ~ ""),
                               label = ifelse(p_bonf_sig, paste0(round(est,digits = 2), stars), ""))
# plot heatmap
p = ggplot(df_betas, aes(x = param, y = paramHeader)) +
  geom_tile(aes(fill = est), color = "white") +
  geom_tile(data = df_betas %>% filter(p_bonf_sig),fill = NA, color = "black",linewidth = 1) +
  scale_fill_gradient2(low = "red4", mid = "white", high = "dodgerblue4", midpoint = 0,na.value = "grey90",name = "Std. Beta",
                       guide = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = 3, barheight = 15),
                       limits = c(-0.25,0.25),breaks = c(-0.2,-0.1, 0, 0.1,0.2)) +
  geom_text(aes(label = label), size = 5, show.legend = FALSE,color="white",fontface="bold") +
  scale_color_identity() +
  theme_classic() +
  labs(x = "", y = "") +
  theme(
    axis.text.x = element_text(size = 14, color = 'black', face = 'bold', angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, color = 'black', face = 'bold'),
    legend.title = element_text(size = 16, face = 'bold'),
    legend.text = element_text(size = 16)
  )
p
# export plot as jpeg
ggsave(paste0(out_dir,"Napls_LatentVariableRegress_Case-Control_HeatmapSpecificity.jpeg"),
       plot=p,height=200,width=450,unit="mm",dpi=500)


### some tests with distribution plots
library(ggridges)
# base density ridgeline plot
p_density = ggplot(df_betas, aes(x = est, y = paramHeader, fill = paramHeader)) +
  geom_density_ridges(alpha = 0.7, scale = 1.2, rel_min_height = 0.01,
                      color = "white") +
  geom_text(
    data = df_betas %>% filter(p_bonf_sig),
    aes(x = est, y = paramHeader, label = param),
    color = "red4", size = 3, fontface = "bold", nudge_y = 0.2,angle = 90,
  ) +
  scale_fill_cyclical(values = c("grey80", "grey60")) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red4")) +
  theme_classic() +
  labs(x = "Std. Beta", y = "Clinical Variable") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold")
  )

# save
ggsave(paste0(out_dir,"Napls_LatentVariableRegress_DensityDistributions.jpeg"),
       plot = p_density, height = 200, width = 300, unit = "mm", dpi = 500)




















