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
in_dir = "PATH/TO/DIRECTORY/"
# define output directory
out_dir = "PATH/TO/DIRECTORY/"
# create output directory for this script
dir.create(out_dir, recursive = TRUE,showWarnings=F)
# read fMRI data (saved as .RData file due to its size)
load(paste0(in_dir,"napls_fMRI_BN_RegionToRegionFC.RData")) # this file was created previously and contains the Region-To-Region FC for each subject
# make copy of FC data
df_fc = df_raw
# rename sessions
df_fc$session = df_fc$session %>% str_replace('2m','M2') %>% str_replace('4m','M4') %>% str_replace('6m','M6') %>% 
  str_replace('8m','M8') %>% factor(levels=c('BL','M2','M4','M6','M8'))
# adjust subject ID column by adding an s at the beginning
df_fc$subject_id = paste0('s',df_fc$subject_id)
# read quality control data
df_fd = read.csv(paste0(in_dir,"napls_fmri_qc.csv"))
# merge with FC data while keeping only subjects with valid FC data
df_fc = merge(df_fc,df_fd,by=c("subject_id","session"),all.x = T)
# read text file containing region names
regions = read.csv(paste0(in_dir,"BN_Atlas_246.txt"),header = F)$V1
# remove dots from column names to make it fit to region names from text file
colnames(df_fc) = colnames(df_fc) %>% str_replace_all("/","") %>% str_replace_all("[.]","") %>% str_replace_all("[+]","")


# COMPUTE HIPPOCAMPAL FUNCTIONAL CONNECTIVITY ---------------------------------------------------------------------------------------------------
# define frontal regions in atlas using the index
hip_all = regions[c(215:218)] %>% str_replace_all("/","") %>% str_replace_all("[.]","") %>% str_replace_all("[+]","")
# create combinations
hip_combi = combn(hip_all, 2)
# collapse each combination with "_and_"
hip_conns = apply(hip_combi, 2, function(x) paste(x, collapse = "_and_"))
# compute mean connectivity
df_fc$fc_hip = rowMeans(df_fc[hip_conns],na.rm=T)


# COMPUTE RIGHT, LEFT AND HOMOTOPIC HIPPOCAMPAL FUNCTIONAL CONNECTIVITY ---------------------------------------------------------------------------------------------------
# insert FC for right hemisphere
df_fc$fc_hipr = df_fc$rHipp_R_and_cHipp_R
# insert FC for left hemisphere
df_fc$fc_hipl = df_fc$rHipp_L_and_cHipp_L
# compute mean connectivity
df_fc$fc_hiph = rowMeans(df_fc[c("rHipp_L_and_rHipp_R","rHipp_L_and_cHipp_R",
                                                 "cHipp_L_and_cHipp_R","cHipp_L_and_rHipp_R")],na.rm=T)

# COMPUTE AVERAGE NETWORK CONNECTIVITY SCORES ---------------------------------------------------------------------------------------------------
# Note: The regions for each network are defined based on visual overlap between Brainnetome atlas regions and Laird networks

##### Salience Network #####
# define insular regions in atlas using the index 
sal_ins = regions[c(165:174)] %>% str_replace_all("/","") %>% str_replace_all("[.]","") %>% str_replace_all("[+]","")
# define ACC regions in atlas using the index 
sal_acc = regions[c(177:180,183,184)] %>% str_replace_all("/","") %>% str_replace_all("[.]","") %>% str_replace_all("[+]","")
# combine in one vector
sal_all = c(sal_ins,sal_acc)
# create combinations
sal_combi = combn(sal_all, 2)
# collapse each combination with "_and_"
sal_conns = apply(sal_combi, 2, function(x) paste(x, collapse = "_and_"))
# compute mean connectivity
df_fc$fc_sal = rowMeans(df_fc[sal_conns],na.rm=T)

##### Default-mode Network #####
# define frontal regions in atlas using the index 
dmn_front = regions[c(11:14)] %>% str_replace_all("/","") %>% str_replace_all("[.]","") %>% str_replace_all("[+]","")
# define temporal regions in atlas using the index 
dmn_temp = regions[c(85,86,97,98)] %>% str_replace_all("/","") %>% str_replace_all("[.]","") %>% str_replace_all("[+]","")
# define parietal regions in atlas using the index 
dmn_pariet = regions[c(151:154)] %>% str_replace_all("/","") %>% str_replace_all("[.]","") %>% str_replace_all("[+]","")
# define posterior cingulate regions in atlas using the index 
dmn_pcc = regions[c(175,176)] %>% str_replace_all("/","") %>% str_replace_all("[.]","") %>% str_replace_all("[+]","")
# combine in one vector
dmn_all = c(dmn_front,dmn_temp,dmn_pariet,dmn_pcc)
# create combinations
dmn_combi = combn(dmn_all, 2)
# collapse each combination with "_and_"
dmn_conns = apply(dmn_combi, 2, function(x) paste(x, collapse = "_and_"))
# compute mean connectivity
df_fc$fc_dmn = rowMeans(df_fc[dmn_conns],na.rm=T)

##### Fronto-parietal Network ####
# define frontal regions in atlas using the index 
fpn_front = regions[c(15,16,21,22)] %>% str_replace_all("/","") %>% str_replace_all("[.]","") %>% str_replace_all("[+]","")
# define patietal regions in atlas using the index 
fpn_pariet = regions[c(137:144)] %>% str_replace_all("/","") %>% str_replace_all("[.]","") %>% str_replace_all("[+]","")
# combine in one vector
fpn_all = c(fpn_front,fpn_pariet)
# create combinations
fpn_combi = combn(fpn_all, 2)
# collapse each combination with "_and_"
fpn_conns = apply(fpn_combi, 2, function(x) paste(x, collapse = "_and_"))
# compute mean connectivity
df_fc$fc_fpn = rowMeans(df_fc[fpn_conns],na.rm=T)

##### Limbic Network ####
# define cingulate regions in atlas using the index 
lim_cc = regions[c(175:188)] %>% str_replace_all("/","") %>% str_replace_all("[.]","") %>% str_replace_all("[+]","")
# define parahippocampal regions in atlas using the index
lim_phip = regions[c(115,116)] %>% str_replace_all("/","") %>% str_replace_all("[.]","") %>% str_replace_all("[+]","")
# define subcortical regions in atlas using the index
lim_subcort = regions[c(211:218,223,224)] %>% str_replace_all("/","") %>% str_replace_all("[.]","") %>% str_replace_all("[+]","")
# combine in one vector
lim_all = c(lim_cc,lim_phip,lim_subcort)
# create combinations
lim_combi = combn(lim_all, 2)
# collapse each combination with "_and_"
lim_conns = apply(lim_combi, 2, function(x) paste(x, collapse = "_and_"))
# compute mean connectivity
df_fc$fc_lim = rowMeans(df_fc[lim_conns],na.rm=T)

##### Sensorimotor Network #####
# define frontal regions in atlas using the index 
smn_all = regions[c(53:68,155:162)] %>% str_replace_all("/","") %>% str_replace_all("[.]","") %>% str_replace_all("[+]","")
# create combinations
smn_combi = combn(smn_all, 2)
# collapse each combination with "_and_"
smn_conns = apply(smn_combi, 2, function(x) paste(x, collapse = "_and_"))
# compute mean connectivity
df_fc$fc_smn = rowMeans(df_fc[smn_conns],na.rm=T)

##### Visual Network #####
# define frontal regions in atlas using the index 
vis_all = regions[c(189:210)] %>% str_replace_all("/","") %>% str_replace_all("[.]","") %>% str_replace_all("[+]","")
# create combinations
vis_combi = combn(vis_all, 2)
# collapse each combination with "_and_"
vis_conns = apply(vis_combi, 2, function(x) paste(x, collapse = "_and_"))
# compute mean connectivity
df_fc$fc_vis = rowMeans(df_fc[vis_conns],na.rm=T)

##### Basal Ganglia #####
# define regions in atlas using the index
bg_all = regions[c(219:230)] %>% str_replace_all("/","") %>% str_replace_all("[.]","") %>% str_replace_all("[+]","")
# create combinations
bg_combi = combn(bg_all, 2)
# collapse each combination with "_and_"
bg_conns = apply(bg_combi, 2, function(x) paste(x, collapse = "_and_"))
# compute mean connectivity
df_fc$fc_bg = rowMeans(df_fc[bg_conns],na.rm=T)

##### Amygdala #####
# define frontal regions in atlas using the index
amyg_all = regions[c(211:214)] %>% str_replace_all("/","") %>% str_replace_all("[.]","") %>% str_replace_all("[+]","")
# create combinations
amyg_combi = combn(amyg_all, 2)
# collapse each combination with "_and_"
amyg_conns = apply(amyg_combi, 2, function(x) paste(x, collapse = "_and_"))
# compute mean connectivity
df_fc$fc_amyg = rowMeans(df_fc[amyg_conns],na.rm=T)

##### Thalamus #####
# define regions in atlas using the index
tha_all = regions[c(231:246)] %>% str_replace_all("/","") %>% str_replace_all("[.]","") %>% str_replace_all("[+]","")
# create combinations
tha_combi = combn(tha_all, 2)
# collapse each combination with "_and_"
tha_conns = apply(tha_combi, 2, function(x) paste(x, collapse = "_and_"))
# compute mean connectivity
df_fc$fc_tha = rowMeans(df_fc[tha_conns],na.rm=T)

##### Putamen #####
# define regions in atlas using the index
put_all = regions[c(225,226,229,230)] %>% str_replace_all("/","") %>% str_replace_all("[.]","") %>% str_replace_all("[+]","")
# create combinations
put_combi = combn(put_all, 2)
# collapse each combination with "_and_"
put_conns = apply(put_combi, 2, function(x) paste(x, collapse = "_and_"))
# compute mean connectivity
df_fc$fc_put = rowMeans(df_fc[put_conns],na.rm=T)

##### Pallidum #####
# define regions in atlas using the index
pal_all = regions[c(221,222)] %>% str_replace_all("/","") %>% str_replace_all("[.]","") %>% str_replace_all("[+]","")
# create combinations
pal_combi = combn(pal_all, 2)
# collapse each combination with "_and_"
pal_conns = apply(pal_combi, 2, function(x) paste(x, collapse = "_and_"))
# compute mean connectivity
df_fc$fc_pal = rowMeans(df_fc[pal_conns],na.rm=T)

##### Nucleus accumbens #####
# define regions in atlas using the index
acu_all = regions[c(223,224)] %>% str_replace_all("/","") %>% str_replace_all("[.]","") %>% str_replace_all("[+]","")
# create combinations
acu_combi = combn(acu_all, 2)
# collapse each combination with "_and_"
acu_conns = apply(acu_combi, 2, function(x) paste(x, collapse = "_and_"))
# compute mean connectivity
df_fc$fc_acu = rowMeans(df_fc[acu_conns],na.rm=T)

##### Nucleus Caudatus #####
# define regions in atlas using the index
cau_all = regions[c(219,220,227,228)] %>% str_replace_all("/","") %>% str_replace_all("[.]","") %>% str_replace_all("[+]","")
# create combinations
cau_combi = combn(cau_all, 2)
# collapse each combination with "_and_"
cau_conns = apply(cau_combi, 2, function(x) paste(x, collapse = "_and_"))
# compute mean connectivity
df_fc$fc_cau = rowMeans(df_fc[cau_conns],na.rm=T)


# DATA FRAME ADJUSTMENTS ---------------------------------------------------------------------------------------------------
# select variables of interest
vars_interest = c("fc_hip","fc_hiph","fc_hipl","fc_hipr","fc_sal","fc_dmn","fc_fpn","fc_lim","fc_smn","fc_vis","fc_bg","fc_amyg",
                  "fc_tha","fc_put","fc_pal","fc_acu","fc_cau")
# select variables of interest for the current work
df_fc = df_fc[c("subject_id","session","fc_FD",vars_interest)]
# add rows with missing values for sessions that do not exist
required_sessions = c("BL", "M2", "M4", "M6", "M8")
df_fc = df_fc %>%
  distinct(subject_id) %>%
  tidyr::crossing(session = required_sessions) %>%
  left_join(df_fc, by = c("subject_id", "session"))


# REMOVE FC VALUES WITH HIGH HEAD MOTION ---------------------------------------------------------------------------------------------------
# get scans with FD > 0.4
df_noisy = df_fc %>% filter(fc_FD > 0.4)
# export data
write_csv(df_noisy,paste0(out_dir,"napls_rsfMRI_FD_HigherThan04mm.csv"))
# remove rows with FD values greater than 0.4
df_fc = df_fc %>% filter(fc_FD <= 0.4)
# remove FD from data frame
df_fc = df_fc %>% select(-fc_FD)

# VISUALIZATION ---------------------------------------------------------------------------------------------------

##### Plot functional connectivity by session #####
# plot variables of interest
plot_list = list()
for (variab in vars_interest) {
  plot_list[[variab]] = local({
    variab = variab
    p = ggplot(df_fc,aes(x=unlist(df_fc[,"session"]),y=unlist(df_fc[,variab]))) +
      geom_violin(fill = "lightblue", alpha = 0.5, trim = FALSE) +
      geom_boxplot(fill="dimgrey",outlier.shape = NA,width=0.4) +
      geom_jitter(size=1,color="grey",width = 0.1,alpha=0.3) +
      labs(x = "", y = variab) +
      theme_classic() +
      theme(axis.title.x = element_text(size=12,color="black",face="bold"),
            axis.text.x = element_text(size=10,color="black"),
            axis.title.y = element_text(size=12,color="black",face="bold"),
            axis.text.y = element_text(size=10,color="black"))
  })
}
# arrange plots with clinical scores
plots = grid.arrange(grobs=plot_list,nrow=4)
# export plot
ggsave(paste0(out_dir,"napls_rsfMRI_FunctionalConnectivity_cleaned.jpeg"),plot=plots,height=300,width=300,unit="mm",dpi=500)


##### Plot functional connectivity by site #####
# select columns to plot, rename them and pivot to long format
df_fc_site = df_fc %>%
  mutate(site = paste0("s", sprintf("%01d", as.integer(substr(subject_id, 2,2))))) %>% 
  select(subject_id, session,site, fc_hip, fc_sal, fc_dmn, fc_fpn, fc_smn, fc_vis, fc_bg) %>%
  rename(subj = subject_id,Session = session, Site = site,`Hippocampal FC` = fc_hip,`SAL FC` = fc_sal,`DMN FC` = fc_dmn,`FPN FC` = fc_fpn,
         `SMN FC` = fc_smn,`VIS FC` = fc_vis,`Basal Ganglia FC` = fc_bg) %>%
  pivot_longer(cols = -c(subj, Session,Site),names_to = "Network",values_to = "FC")
# change factor levels
df_fc_site$Network = factor(df_fc_site$Network,levels=c("Hippocampal FC","Basal Ganglia FC","SAL FC","DMN FC","FPN FC","SMN FC","VIS FC"))
# plot site effects
p = ggplot(df_fc_site,aes(x=Site,y=FC)) +
  geom_violin(aes(fill=Site), alpha = 0.7, trim = FALSE) +
  geom_boxplot(fill="dimgrey",outlier.shape = NA,width=0.4) +
  geom_jitter(size=1,width = 0.1,color="grey",alpha=0.3) +
  scale_y_continuous(labels = label_comma()) +
  scale_fill_viridis(discrete = TRUE,option = "plasma") +
  labs(x = "", y = "Average Functional Connectivity per Region/Network") +
  facet_grid(Network~Session) +
  theme_classic() +
  theme(strip.text = element_text(size=20,color="black",face="bold"),
        axis.text.x = element_text(size=18,color="black"),
        axis.title.y = element_text(size=20,color="black",face="bold"),
        axis.text.y = element_text(size=18,color="black"),
        legend.position = "none")
# export plot
ggsave(paste0(out_dir,"napls_rsfMRI_SelectedRegions_cleaned_SessionBySite.jpeg"),
       plot=p,height=500,width=450,unit="mm",dpi=500)


# DATA EXPORT ---------------------------------------------------------------------------------------------------
# export data frame
write_csv(df_fc,paste0(out_dir,"napls_rsfMRI_FunctionalConnectivity_cleaned.csv"))















