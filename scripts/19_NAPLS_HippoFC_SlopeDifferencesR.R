# ADD-ON ANALYSIS: GROUP DIFFERENCES OF CHANGE OF HIPPCAMPUS FC ---------------------------------------------
# Description:  This script compares hippocampal FC slopes between CHR and HCs
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
# create object for hippocampal FC
nvar = "fchi"


# LATENT VARIABLE REGRESSION: CASE-CONTROL COMPARISON OF HIPPOCAMPAL SLOPE ----------------------------------------------------------------------------------------------------
# Note: This code runs the latent variable regression models in MPLUS comparing the hippocampal FC
# change between CHR and healthy controls. The latent variable regression is used to account for
# the measurement error from the change estimates. The accuracy of the estimation of linear change over time per subject in each outcome
# differs depending on the available time points per subject in each outcome.
# re-define and set working directory in which the MPLUS .exe file is stored -> necessary because the path of the working directory
# must not be too long for MPLUS
wd = "PATH/TO/MPLUS/DIRECTORY"
setwd(wd)
# run latent variable regression in loop
# shorten variable names of MRI variables for MPLUS
nvar_mplus = nvar %>% str_replace("fc","")
# create data frame with variables of interest
df = df_napls_mlr_imp[c(paste0(nvar,"_CH"),"CHR",paste0(nvar,"_SE"),"NegS_CH",
                        "age","female","site02","site03","site04","site05","site06","site07","site08","site09")]
# remove subjects with any missings
df = df[complete.cases(df[c(paste0(nvar, "_CH"),"NegS_CH")]), ]
# get average error variances across group from standard errors of both variables
se_nvar = mean(df[[paste0(nvar,"_SE")]],na.rm=T)^2
# drop standard error columns again
df = df[ , !grepl("_SE$", names(df)) ]
# drop Negative Symptom change column again
df = df[ , !grepl("NegS_CH", names(df)) ]
# define the input and output files
datafile = file.path(wd, paste0("NaplsData_", nvar, ".dat"))
inpfile = file.path(wd, paste0("NaplsData_", nvar, ".inp"))
# write output file
write.table(df,datafile,quote = F,row.names = F, col.names = F,append = F, sep = "\t")
# define the MPLUS input file containing the model
inp_text = paste0("
Title:      LATENT VARIABLE REGRESSION

Data:       File is ", datafile, ";

Variable:
names are
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
            ALGORITHM=INTEGRATION;
            STARTS=20;
            STITERATION=50000;

Model:

LV_", nvar_mplus, " by ", nvar_mplus, "_CH;
", nvar_mplus, "_CH@", format(se_nvar, scientific = FALSE), ";

LV_", nvar_mplus, " on CHR age female site02 site03 site04
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

SAVEDATA: ESTIMATES = power_",nvar,".dat;
Output:     tech1 tech4 stdyx sampstat;
"
)
# write as input file for MPLUS
writeLines(inp_text, con = inpfile)
# run MPLUS
system(paste("wine Mplus-8.6.exe", basename(inpfile)))
# specify the output file
outfile = file.path(wd, paste0("naplsdata_", tolower(nvar), ".out"))
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
write.csv(df_res, file = paste0("SingleIndicator_Case-Control_PathsUnstandardized_",nvar, ".csv"), row.names = FALSE)
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
write.csv(df_res, file = paste0("SingleIndicator_Case-Control_PathsStandardized_",nvar, ".csv"), row.names = FALSE)
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
write.csv(df_model, file = paste0("SingleIndicator_Case-Control_Model_", nvar, ".csv"), row.names = FALSE)

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

