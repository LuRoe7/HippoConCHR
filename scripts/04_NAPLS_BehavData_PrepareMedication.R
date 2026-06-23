# PREPARATION OF MEDICATION DATA ---------------------------------------------
# Description:  This script prepares and cleans the medication data from the NAPLS3 study.
#               
# Author:       Roell, Lukas      
# Created:      2025/08/28
# License:      Creative Common: CC-BY

##### Load packages #####
# load necessary packages
library(pacman)
p_load(readr,tidyr,dplyr,stringr,ggplot2,ggrepel,gridExtra,readxl,chlorpromazineR)


# BASIC PREPARATIONS ---------------------------------------------------------------------------
# set path to input folder in which medication data is stored
in_path = "PATH/TO/DIRECTORY/"
# define output folder
out_path = "PATH/TO/DIRECTORY/"
# define input directory for behavioral and functional connectivity data
in_dir = "PATH/TO/DIRECTORY/"
# define output directory
out_dir = "PATH/TO/DIRECTORY/"
# create output directory for this script
dir.create(out_dir, recursive = TRUE,showWarnings=F)
# read merged NAPLS data
df_napls_wide = read.csv(paste0(in_dir,"napls_MergedData_VarsOfInterest.csv"))
# read medication data
df_med = read_excel(paste0(in_path, "N3_MedLog2.xlsx"))
# filter by VisitLabel
df_med = subset(df_med,VisitLabel == "BL" | VisitLabel == "M2" | VisitLabel == "M4" | VisitLabel == "M6" | VisitLabel == "M8")
# adjust subject column according to naming scheme used in the other data frames
df_med$subject_id = paste0("s",as.character(df_med$SiteNumber),"0",df_med$SubjectNumber)
# order data frame by subject ID and visit
df_med = df_med[order(df_med$subject_id,df_med$VisitNumber),]
# select columns of interest
df_med = df_med[c("subject_id","SubjectType","VisitLabel","DataCollectedDate","medlist_med_code","medlist_med_name","medlist_start_date",    
                  "medlist_stop_date","medlist_dosage_mg","medlist_dosage_ua","medlist_dosage_mg_prn","medlist_dosage_ua_prn", 
                  "medlist_frequency","medlist_source_of_data","medlist_compliance")]
# convert names to start with lower case
df_med$medlist_med_name = tolower(df_med$medlist_med_name)


# PREPARE ANTIPSYCHOTICS ---------------------------------------------------------------------------------
# create a list with antipsychotics that are part of the data frame. Note that the same substance has different names sometimes
antipsychotics = c("olanzapine","zyprexa","olazepine","olonziprine",
                  "clozapine",
                  "quetiapine","seroquel","quetrapine","zoloft","seroquel xr",
                  "chlorpromazine","vraylar","antipsychotic",
                  "haloperidol","haldol",
                  "aripriprazole","abilify", "ability", "aripiprazole","arpiprazole", "abilify/aripiprazole", "aripiprazol",
                  "risperidone","risperdal","risperidol","risperdone","risperidone ic", "risperidon", "resperidone","rispendal","risperidal","risperdol",
                  "rexulti", # Brexipiprazole
                  "lurasidone","latuda",
                  "ziprasidone","geodon","geodone", 
                  "paliperidone","invega",
                  "saphris", # asenapine
                  "trilafon", # Perphenazine
                  "navane" # Thiothixene
                  )
# filter data frame selecting only rows that contain these antipsychotic
df_antipsy = df_med %>% filter(medlist_med_name %in% antipsychotics)
# rename entries of medication names to make them consistent
df_antipsy$antipsychotic = df_antipsy$medlist_med_name %>% str_replace("zyprexa","olanzapine") %>% str_replace("olazepine","olanzapine") %>%
  str_replace("olonziprine","olanzapine") %>% str_replace("seroquel xr","quetiapine") %>% str_replace("seroquel","quetiapine") %>%
  str_replace("quetrapine","quetiapine") %>% str_replace("zoloft","quetiapine") %>%
  str_replace("vraylar","chlorpromazine") %>% str_replace("antipsychotic","chlorpromazine") %>% str_replace("haloperidol","haloperidole") %>%
  str_replace("haldol","haloperidol") %>% str_replace("abilify/aripiprazole","aripiprazole") %>% str_replace("arpiprazole","aripiprazole") %>%
  str_replace("aripiprazol","aripiprazole") %>% str_replace("abilify","aripiprazole") %>% str_replace("ability","aripiprazole") %>% 
  str_replace("aripriprazole","aripiprazole") %>% str_replace("rexulti","brexipiprazole") %>% str_replace("latuda","lurasidone") %>%
  str_replace("risperdal","risperidone") %>% str_replace("risperidol","risperidone") %>% str_replace("risperdone","risperidone") %>%
  str_replace("risperidone ic","risperidone") %>% str_replace("resperidone","risperidone") %>% str_replace("resperidone","risperidone") %>%
  str_replace("rispendal","risperidone") %>% str_replace("risperdol","risperidone") %>%
  str_replace("risperida","risperidone") %>% str_replace("risperidon","risperidone") %>%
  str_replace("geodon","ziprasidone") %>% str_replace("geodone","ziprasidone") %>% str_replace("invega","paliperidone") %>%
  str_replace("saphris","asenapine") %>% str_replace("trilafon","perphenazine") %>% str_replace("navane","tiotixene") %>%
  str_replace("eel","e") %>% str_replace("ee","e")
# select columns of interest
df_antipsy = df_antipsy[c("subject_id","SubjectType","VisitLabel","antipsychotic","medlist_dosage_mg","DataCollectedDate",
                          "medlist_start_date","medlist_stop_date")]
# rename columns
colnames(df_antipsy) = c("subject_id","subjecttype","session","antipsychotic","dosage","CollectionDate","StartDate","StopDate")
# order data frame by subject ID
df_antipsy = df_antipsy[order(df_antipsy$subject_id,df_antipsy$session,df_antipsy$StartDate),]


# COMPUTE CPZ EQUIVALENTS ---------------------------------------------------------------------------------
# specify defined oral daily doses according to Leucht et al., 2016
ddd_quetiapine = 400
ddd_aripiprazole = 15
ddd_risperidone = 5
ddd_clozapine = 300
ddd_lurasidone = 60
ddd_ziprasidone = 80
ddd_chlorpromazine = 300
ddd_olanzapine = 30
ddd_paliperidone = 6
ddd_asenapine = 20
ddd_perphenazine = 30
ddd_brexipiprazole = 3
ddd_haloperidol = 8
ddd_tiotixene = 30
# add empty cpz column to data frame
df_antipsy$cpz = NA
# compute chlorpromazine equivalents for each antipsychotic
for (i in 1:nrow(df_antipsy)) {
  apsy = df_antipsy[i,"antipsychotic"]
  apsy = apsy$antipsychotic
  dose = df_antipsy[i,"dosage"]
  dose = dose$dosage
  if (apsy == "quetiapine") {
    cpz = dose*(ddd_chlorpromazine/ddd_quetiapine)
    df_antipsy[i,"cpz"] = cpz
  }
  else if (apsy == "aripiprazole") {
    cpz = dose*(ddd_chlorpromazine/ddd_aripiprazole)
    df_antipsy[i,"cpz"] = cpz
  }
  else if (apsy == "risperidone") {
    cpz = dose*(ddd_chlorpromazine/ddd_risperidone)
    df_antipsy[i,"cpz"] = cpz
  }
  else if (apsy == "clozapine") {
    cpz = dose*(ddd_chlorpromazine/ddd_clozapine)
    df_antipsy[i,"cpz"] = cpz
  }
  else if (apsy == "lurasidone") {
    cpz = dose*(ddd_chlorpromazine/ddd_lurasidone)
    df_antipsy[i,"cpz"] = cpz
  }
  else if (apsy == "ziprasidone") {
    cpz = dose*(ddd_chlorpromazine/ddd_ziprasidone)
    df_antipsy[i,"cpz"] = cpz
  }
  else if (apsy == "chlorpromazine") {
    cpz = dose*(ddd_chlorpromazine/ddd_chlorpromazine)
    df_antipsy[i,"cpz"] = cpz
  }
  else if (apsy == "olanzapine") {
    cpz = dose*(ddd_chlorpromazine/ddd_olanzapine)
    df_antipsy[i,"cpz"] = cpz
  }
  else if (apsy == "paliperidone") {
    cpz = dose*(ddd_chlorpromazine/ddd_paliperidone)
    df_antipsy[i,"cpz"] = cpz
  }
  else if (apsy == "asenapine") {
    cpz = dose*(ddd_chlorpromazine/ddd_asenapine)
    df_antipsy[i,"cpz"] = cpz
  }
  else if (apsy == "perphenazine") {
    cpz = dose*(ddd_chlorpromazine/ddd_perphenazine)
    df_antipsy[i,"cpz"] = cpz
  }
  else if (apsy == "brexipiprazole") {
    cpz = dose*(ddd_chlorpromazine/ddd_brexipiprazole)
    df_antipsy[i,"cpz"] = cpz
  }
  else if (apsy == "haloperidol") {
    cpz = dose*(ddd_chlorpromazine/ddd_haloperidol)
    df_antipsy[i,"cpz"] = cpz
  }
  else if (apsy == "tiotixene") {
    cpz = dose*(ddd_chlorpromazine/ddd_tiotixene)
    df_antipsy[i,"cpz"] = cpz
  }
}
# select only those dosages that were taken within the collection date
df_antipsy_long = df_antipsy %>%
  group_by(subject_id, session, antipsychotic) %>%
  slice_tail(n = 1) %>%
  ungroup() %>%
  filter(CollectionDate >= StartDate & CollectionDate <= StopDate)
# compute subject- and session-specific CPZ sum across antipsychotics
df_antipsy_long_cpz = df_antipsy_long %>%
  group_by(subject_id, session) %>%
  summarise(cpz = sum(cpz, na.rm = TRUE), .groups = "drop")
# convert to wide format
df_antipsy_wide_cpz = df_antipsy_long_cpz %>% pivot_wider(id_cols = c(subject_id),names_from = session,values_from = cpz)
# reorder columns
df_antipsy_wide_cpz = df_antipsy_wide_cpz[c("subject_id","BL","M2","M4","M6","M8")]
# adjust column names
colnames(df_antipsy_wide_cpz) = c("subject_id","cpz__BL","cpz__M2","cpz__M4","cpz__M6","cpz__M8")
# add antipsychotics to merged data frame
df_napls_wide = merge(df_napls_wide,df_antipsy_wide_cpz,by="subject_id",all.x = T)
# replace NAs in CPZ columns by zero if the respective subject had either a GAF or SOPS value at this session. This means that the CPZ value is not missing. 
df_napls_wide = df_napls_wide %>%
  mutate(across(starts_with("cpz__"), ~ ifelse(
    is.na(.) & (!is.na(get(sub("cpz__", "gafrat__", cur_column()))) |
                  !is.na(get(sub("cpz__", "sops_positive__", cur_column())))), 0, .)))


# VISUALIZE RAW CPZ EQUIVALENTS ---------------------------------------------------------------------------------
# select CPZ values for plotting
df_cpz_plot = df_napls_wide[c("subject_id","cpz__BL","cpz__M2","cpz__M4","cpz__M6","cpz__M8")]
# bring to long format for plotting
df_cpz_plot = df_cpz_plot %>% pivot_longer(cols = all_of(c("cpz__BL","cpz__M2","cpz__M4","cpz__M6","cpz__M8")),
                                           names_to = c("outcome","session"),
                                           names_sep = "__",
                                           values_to = "value")
# insert numeric column for session
df_cpz_plot$session_num = df_cpz_plot$session %>% str_replace("BL","1") %>% str_replace("M2","2") %>% str_replace("M4","3") %>%
  str_replace("M6","4") %>% str_replace("M8","5") %>% as.numeric()
# plot to check for implausible CPZ values and strange courses
thresh = 500
p = ggplot(df_cpz_plot, aes(x = session_num, y = value)) +
  geom_point(aes(color = value > thresh, size = value > thresh)) +
  geom_line(aes(group = subject_id)) +
  geom_text_repel(aes(label = ifelse(value > thresh, subject_id, "")), 
                  size = 5, max.overlaps = Inf, box.padding = 0.5,segment.color = NA) +
  theme_classic() +
  labs(x = "", y = "CPZ") +
  scale_x_continuous(labels=c("BL","M2","M4","M6","M8")) +
  theme(axis.title = element_text(color = "black", size = 16),
        axis.text = element_text(color = "black", size = 14),
        legend.position = "none") +
  scale_color_manual(values = c("black", "red")) +
  scale_size_manual(values = c(2, 4))
# export plot
ggsave(paste0(out_dir,"napls_cpz_course.jpeg"),plot=p,height=200,width=200,unit="mm",dpi=500)


# CLEANING OF IMPLAUSIBLE MEDICATION DATA ----------------------------------------------------------------------------
# based on visual inspection the subjects s10333 (session M2) and s60353 (session M8) require further examination
# s10333
#View(subset(df_napls_wide[c(1,2,3,4,grep("sops_p",colnames(df_napls_wide)),grep("gafra",colnames(df_napls_wide)))],
            #subject_id == "s10333")) # this subject is a 16.5 years female with high symptom ratings
#View(subset(df_med,subject_id == "s10333")) # but has 25mg risperidone which is way too much and must be an error
df_napls_wide$cpz__M2[df_napls_wide$cpz__M2 == 1500] = NA # remove this value
# s60353
#View(subset(df_napls_wide[c(1,2,3,4,grep("sops_p",colnames(df_napls_wide)),grep("gafra",colnames(df_napls_wide)))],
            #subject_id == "s60353")) # this subject is a 16.1 years female with high symptom ratings
#View(subset(df_med,subject_id == "s60353")) # takes 150 mg aripiprazole at M8 which is way too much
df_napls_wide$cpz__M8[df_napls_wide$cpz__M8 == 3000] = 300 # replace this value with 300 given that she had 15mg aripiprazole (= 300 CPZ) in the other sessions


# EXPORT DATA --------------------------------------------------------------------------------------------------------
# export data frame
write.csv(df_napls_wide,paste0(out_dir,"napls_MergedData_IncludingMed.csv"),row.names = F)

