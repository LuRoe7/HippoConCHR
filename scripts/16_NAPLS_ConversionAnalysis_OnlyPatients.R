# ADD-ON ANALYSIS: EFFECT OF HIPPOCAMPAL FC CHANGE ON CONVERSION ---------------------------------------------
# Description:  This script explores the effect of intra-hippocampal FC on conversion rates
#               
# Author:       Roell, Lukas      
# Created:      2025/08/28
# License:      Creative Common: CC-BY

##### Load packages #####
# load necessary packages
library(pacman)
p_load(readr,tidyr,dplyr,stringr,ggplot2,ggrepel,gridExtra,ggExtra,readxl,writexl,lme4,broom,lavaan,viridis,glue,purrr,reshape2,
       MplusAutomation,tibble,lavaan,semTools)


# BASIC PREPARATIONS -------------------------------------------------------------------------------------------------------------
# show working directory and ensure that the folder results including the data is located there
getwd()
# define input directory for behavioral and FC data
in_dir = "results/napls/P1_HippoFC/07_NAPLS_SingleSubjectRegressions/"
# define output directory
out_dir = "results/napls/P1_HippoFC/16_NAPLS_ConversionAnalysis_OnlyPatients/"
# create output directory for this script
dir.create(out_dir, recursive = TRUE,showWarnings=F)
# read data including behavioral and FC data and all slopes with imputed SEs
df_napls_mlr_imp = read.csv(paste0(in_dir,"Napls_BehavFCandSlopes_impSE.csv"))
# keep only CHR individuals
df_napls_mlr_imp = subset(df_napls_mlr_imp,group ==1)
# define input directory for conversion data
in_dir_conv = "data/napls/behavioral/"
# read data including behavioral and FC data and all slopes with imputed SEs
df_napls_conv = read.csv(paste0(in_dir_conv,"napls_behav_conversion.csv"))
# change rowname of subject ID column
df_napls_conv = df_napls_conv %>% rename(subj = mri_id) %>% select(-diaggroup)
# remove s from subject
df_napls_conv$subj = df_napls_conv$subj %>% str_replace("s","")
# merge data frames
df_napls_mlr_imp = merge(df_napls_mlr_imp,df_napls_conv,by="subj",all.x=T)
# create data frame with variables of interest
df = df_napls_mlr_imp[c("subj","converted","fchi_CH","fchi_SE",
                        "age","female","site02","site03","site04","site05","site06","site07","site08","site09",
                        "cpz_CH","cpz_SE","NegS_CH")]
# remove subjects with any missings
df = df[complete.cases(df[c("converted","fchi_CH","NegS_CH")]), ]
# rename conversion variable
df = df %>% rename(conv = converted)
# define FC variable
nvar = "fchi"

# CONVERSION COUNT -------------------------------------------------------------------------------------------------------------
# count total subjects, number of converters, and compute percentage
conversion_summary = df %>%
  summarise(
    n_total      = n(),                                    # total number of subjects in df
    n_converted  = sum(conv == 1, na.rm = TRUE),      # number of subjects who converted
    pct_converted = round(n_converted / n_total * 100, 1)  # percentage of converters rounded to 1 decimal
  )

# print result in a readable format
cat("Total subjects:", conversion_summary$n_total, "\n",       # print total N
    "Converters:   ", conversion_summary$n_converted, "\n",    # print number of converters
    "Percentage:   ", conversion_summary$pct_converted, "%\n") # print percentage


# LATENT VARIABLE REGRESSION: CONVERSION PREDICTION ----------------------------------------------------------------------------------------------------
# Note: This code runs the latent variable regression models in MPLUS predictint conversion based on hippocampal slope
# The latent variable regression is used to account for
# the measurement error from the change estimates. The accuracy of the estimation of linear change over time per subject in each outcome
# differs depending on the available time points per subject in each outcome.
# re-define and set working directory in which the MPLUS .exe file is stored -> necessary because the path of the working directory
# must not be too long for MPLUS
wd = "/home/lukas/mplus"
setwd(wd)
# run latent variable regression in loop
# shorten variable names of MRI variables for MPLUS
nvar_mplus = nvar %>% str_replace("fc","")
# create data frame with variables of interest
df = df[,c("conv",paste0(nvar,"_CH"),paste0(nvar,"_SE"),
           "age","female","site02","site03","site04","site05","site06","site07","site08","site09")]
# get average error variances across group from standard errors of both variables
se_nvar = mean(df[[paste0(nvar,"_SE")]],na.rm=T)^2
# drop standard error columns again
df = df[ , !grepl("_SE$", names(df)) ]
# define the input and output files
datafile = file.path(wd, paste0("NaplsData_", nvar, ".dat"))
inpfile = file.path(wd, paste0("NaplsData_", nvar, ".inp"))
# write output file
write.table(df,datafile,quote = F,row.names = F, col.names = F,append = F, sep = "\t")
# define the MPLUS input file containing the model
inp_text = paste0("
Title:      LATENT VARIABLE LOGISTIC REGRESSION

Data:       File is ", datafile, ";

Variable:
names are
conv
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
conv
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

categorical = conv;

missing are all (-99);

Analysis:   type = random;
            ALGORITHM=INTEGRATION;
            STARTS=20;
            STITERATION=50000;

Model:

LV_", nvar_mplus, " by ", nvar_mplus, "_CH;
", nvar_mplus, "_CH@", format(se_nvar, scientific = FALSE), ";

conv on LV_", nvar_mplus, " age female site02 site03 site04
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
