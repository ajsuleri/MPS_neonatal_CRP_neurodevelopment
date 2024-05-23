#########################################################################
######################PROJECT: Neonatal inflammation####################
#########################################################################

#Author: Nicole Creasey
#Project: Mapping prenatal predictors and neurobehavioral outcomes of an epigenetic marker of neonatal inflammation - a longitudinal population-based study 

#SET UP WORKSPACE 
rm(list = ls()) #clears the environment

libraries <- c('foreign', 'haven', 'tidyverse', 'broom.mixed', 'xlsx', 'corrplot', 'stringi', 'miceadds', 'mitools', 'lme4', 'lmerTest', 'rlang', 'lattice', 'ggplot2', 'ggeffects', 'splines', 'reshape2', 'gridExtra', 'grid', 'ggpubr','CBPS', 'survey', 'survival', 'ggseg', 'fmsb', 'ggcorrplot', 'Hmisc', 'pROC', 'readxl', 'pROC', 'RColorBrewer', 'DescTools') 

invisible(lapply(libraries, require, character.only = T))

#set directory 
setwd("set_path_to_data")

#load raw data 
load("/home/n.creasey/episcore_CRP_project/shadow2/mergeddata/CRP_mps_merged_rawdata.Rdata")

#### CREATE THE FINAL DATASET ####

### EXCLUSIONS: handle any exclusions ###

#select only those with mps data 
df <- df[!is.na(df$array),] 

#remove ppts that withdrew data
df <- df[is.na(df$RemoveData),]  #no ppts to exclude 

#remove twin births
df <- df[df$TWIN=='No',]  # no twin births 

# get indicator for siblings -> first sort on missingness
df$missing <- rowSums(is.na(df)) ## calculate missingness
df$missing[is.na(df$Sample_ID)] <- 1000

df <- df[order(df$missing),] 

df$first_sibling_most_data <- !duplicated(df$MOTHER) # get indicator for first child (so no siblings)

df <- df[df$first_sibling_most_data == TRUE,] #   2511 -> 2394 


#mark cases that should NOT be included in the brain imaging analyses 
df$exc_smri_f09 <- as.factor(ifelse(df$mri_consent_f09 == 'yes' & df$t1_has_nii_f09 == 'yes' & df$t1_asset_has_nii_f09 != 'exclude' & df$has_braces_mri_f09 == 'no' & df$exclude_incidental_f09 == 'include' & df$freesurfer_qc_f09 == 'usable', 'include', 'exclude')) 

df$exc_smri_f13 <- as.factor(ifelse(df$mri_consent_f13 == 'yes' & df$t1_has_nii_f13 == 'yes' & df$has_braces_mri_f13 == 'no' & df$exclude_incidental_f13 == 'include' & df$freesurfer_qc_f13 == 'usable', 'include', 'exclude')) 

df$exc_dti_f09 <- as.factor(ifelse(df$mri_consent_f09 == 'yes' & df$has_braces_mri_f09 == 'no' & df$dti_has_nii_f09 == 'yes' & df$dwi_nvols_f09 == 38 & df$exclude_incidental_f09 == 'include' & df$software_ver_short_f09 == 24 & is.na(df$dti_overall_qc_f09), 'include', 'exclude'))

df$exc_dti_f13 <- as.factor(ifelse(df$mri_consent_f13 == 'yes' & df$has_braces_mri_f13 == 'no' & df$dti_has_nii_f13 == 'yes' & df$dwi_nvols_f13 == 38 & df$exclude_incidental_f13 == 'include' & is.na(df$dti_overall_qc_f13), 'include', 'exclude'))

# Set brain imaging data that should be excluded to NA

smri9_vars <- c('Left_Cerebellum_Cortex_vol_f09', 'Left_Hippocampus_vol_f09', 'Left_Amygdala_vol_f09', 'Left_Lateral_Ventricle_vol_f09', 'Brain_Stem_vol_f09','Right_Cerebellum_Cortex_vol_f09', 'Right_Amygdala_vol_f09', 'Right_Hippocampus_vol_f09', 'Right_Lateral_Ventricle_vol_f09','TotalGrayVol_f09', 'genr_tbv_f09', 'eTIV_f09', 'lhCerebralWhiteMatterVol_f09', 'rhCerebralWhiteMatterVol_f09') 


smri13_vars <-c('Left_Cerebellum_Cortex_vol_f13', 'Left_Hippocampus_vol_f13', 'Left_Amygdala_vol_f13', 'Left_Lateral_Ventricle_vol_f13', 'Brain_Stem_vol_f13','Right_Cerebellum_Cortex_vol_f13', 'Right_Amygdala_vol_f13', 'Right_Hippocampus_vol_f13', 'Right_Lateral_Ventricle_vol_f13', 'TotalGrayVol_f13', 'genr_tbv_f13', 'eTIV_f13', 'lhCerebralWhiteMatterVol_f13', 'rhCerebralWhiteMatterVol_f13')

dti9_vars <- c('mean_FA_genr_f09', 'mean_MD_genr_f09')

dti13_vars<- c('mean_FA_genr_f13', 'mean_MD_genr_f13')


df[,smri9_vars] <- if_else(df$exc_smri_f09 == "exclude", NA, df[,smri9_vars])
df[,smri13_vars] <- if_else(df$exc_smri_f13 == "exclude", NA, df[,smri13_vars])
df[,dti9_vars] <- if_else(df$exc_dti_f09 == "exclude", NA, df[,dti9_vars])
df[,dti13_vars] <- if_else(df$exc_dti_f13 == "exclude", NA, df[,dti13_vars])


### Composites for brain data ###
#(do this before recoding other variables otherwise when remove brain variables starting 'rh' also removes 'rhinitis' variables)

## Averaging across hemispheres (because no hypothesis on lateralized effects)
#Loop over all structures and average columns of the same structure in both hemispheres (at the same measurement occasion)
aseg_cols <- colnames(dplyr::select(df, starts_with("Left_"), starts_with("Right_"))) #select subcortical structures of which we have left and right hemisphere data
tbv_cols <- colnames(dplyr::select(df, starts_with("lh"), starts_with('rh')))  #select cortical structures of which we have left and right hemisphere data

aseg_structs1 <-  stri_remove_empty(unlist(strsplit(aseg_cols, "Left_"))) #unlist to make vector
aseg_structs2 <-  stri_remove_empty(unlist(strsplit(aseg_structs1, "Right_")))
aseg_structs3 <- stri_remove_empty(unlist(strsplit(aseg_structs2, "_f09")))
aseg_structs4 <- stri_remove_empty(unlist(strsplit(aseg_structs3, "_f13")))

tbv_structs1 <- stri_remove_empty(unlist(strsplit(tbv_cols, "lh")))
tbv_structs2 <- stri_remove_empty(unlist(strsplit(tbv_structs1, "rh")))
tbv_structs3 <- stri_remove_empty(unlist(strsplit(tbv_structs2, "_f09")))
tbv_structs4 <- stri_remove_empty(unlist(strsplit(tbv_structs3, "_f13")))

structs <- c(aseg_structs4, tbv_structs4) #combine cortical and subcortical structure names in vector 

#Go over all structures
 for(x in structs){ 
    #Go over all columns to identify the first column (left hemisphere)
    for(y in 1:ncol(df)){
      #Store the column name to find the counterpart in the other (right) hemisphere
      a <- colnames(df[y])
      #Go over all columns to identify the second column (right hemisphere)
      for(z in 1:ncol(df)){
        #Store the column name to identify the same structures in each hemisphere
        b <- colnames(df[z])
        #If we have the same structure for both hemispheres, calculate the mean of those
        ##MRI T2
        #Subcortical volumes MRI F9
        if(a == paste0("Left_",x,"_f09") & b == paste0("Right_",x,"_f09")){
          print(x)
          newcolname <- paste0(x,"_subcortical_f09")
          df[,newcolname] <- (df[,y] + df[,z])/2
        }
        #Cortical volumes MRI F9
        if(a == paste0("lh",x,"_f09") & b == paste0("rh",x,"_f09")){
          print(x)
          newcolname <- paste0(x,"_cortical_f09")
          df[,newcolname] <- (df[,y] + df[,z])/2
        }
        ##MRI T3
        #Subcortical volumes MRI F13
        if(a == paste0("Left_",x,"_f13") & b == paste0("Right_",x,"_f13")){
          print(x)
          newcolname <- paste0(x,"_subcortical_f13")
          df[,newcolname] <- (df[,y] + df[,z])/2
        }
        #Cortical volumes MRI F13
        if(a == paste0("lh",x,"_f13") & b == paste0("rh",x,"_f13")){
          print(x)
          newcolname <- paste0(x,"_cortical_f13")
          df[,newcolname] <- (df[,y] + df[,z])/2
        }
      }
    }
  }

df <- dplyr::select(df, -c(starts_with("Left_"), starts_with("Right_"), starts_with("lh"), starts_with("rh"))) #Now tidy up the dataset and remove all original hemisphere specific variables


### RECODE: recoding of variables ###

#recode the education variable to three levels 
#df$EDUCM_3l <- unclass(df$EDUCM) 
#df$EDUCM_3l <- factor(recode(df$EDUCM_3l, `1` = 1, `2` = 1, `3` = 2, '4'= 2, '5'= 3, '6' = 3))


## PRENATAL INFECTION SCORE (PIS) - make the preinfection variables dichotomous

df$fever_tri1 <- ifelse(df$c1100101 == 'Yes', 1, ifelse(df$c1100101 == 'No', 0, NA))
df$flu_tri1 <- ifelse(df$c1100301 == 'Yes', 1, ifelse(df$c1100301 == 'No', 0, NA))
df$pharyngitis_tri1<- ifelse(df$c1100501 == 'Yes', 1, ifelse(df$c1100501 == 'No', 0, NA))
df$rhinitis_tri1<- ifelse(df$c1100701 == 'Yes', 1, ifelse(df$c1100701 == 'No', 0, NA))
df$sinusitis_tri1<- ifelse(df$c1100901 == 'Yes', 1, ifelse(df$c1100901 == 'No', 0, NA))
df$earinf_tri1<- ifelse(df$c1101101 == 'Yes', 1, ifelse(df$c1101101 == 'No', 0, NA))
df$pneumonia_tri1<- ifelse(df$c1101301 == 'Yes', 1, ifelse(df$c1101301 == 'No', 0, NA))
df$dermatitis_tri1<- ifelse(df$c1102301 == 'Yes', 1, ifelse(df$c1102301 == 'No', 0, NA))
df$herpeszoster_tri1<- ifelse(df$c1102901 == 'Yes', 1, ifelse(df$c1102901 == 'No', 0, NA))
df$enteritis_tri1<- ifelse(df$c1103101 == 'Yes', 1, ifelse(df$c1103101 == 'No', 0, NA))
df$cystitis_tri1<- ifelse(df$c1103301 == 'Yes', 1, ifelse(df$c1103301 == 'No', 0, NA))
df$jaundice_tri1<- ifelse(df$c1103501 == 'Yes', 1, ifelse(df$c1103501 == 'No', 0, NA))
df$otherinf_tri1 <- ifelse(df$c1103701 == 'Yes', 1, ifelse(df$c1103701 == 'No', 0, NA))
df$STD_tri1 <- ifelse(df$Chlamydia_or_other_STD == 'Yes', 1, ifelse(df$Chlamydia_or_other_STD == 'No', 0, NA))
df$lower_resp_inf_tri1 <- ifelse(df$Lower_respiratory_infection == 'Yes', 1, ifelse(df$Lower_respiratory_infection == 'No', 0, NA))
df$upper_resp_inf_tri1 <- ifelse(df$Upper_respiratory_infection == 'Yes', 1, ifelse(df$Upper_respiratory_infection == 'No', 0, NA))
df$GI_inf_tri1 <- ifelse(df$Gastro_intestinal_infection == 'Yes', 1, ifelse(df$Gastro_intestinal_infection == 'No', 0, NA))
df$UWI_tri1 <- ifelse(df$UWI == 'Yes', 1, ifelse(df$UWI == 'No', 0, NA))

df$fever_tri2 <- ifelse(df$a0900103 == 'Yes', 1, ifelse(df$a0900103 == 'No', 0, NA))
df$flu_tri2 <- ifelse(df$a0900303 == 'Yes', 1, ifelse(df$a0900303 == 'No', 0, NA))
df$pharyngitis_tri2<- ifelse(df$a0900503 == 'Yes', 1, ifelse(df$a0900503 == 'No', 0, NA))
df$rhinitis_tri2 <- ifelse(df$a0900703 == 'Yes', 1, ifelse(df$a0900703 == 'No', 0, NA))
df$sinusitis_tri2<- ifelse(df$a0900903 == 'Yes', 1, ifelse(df$a0900903 == 'No', 0, NA))
df$earinf_tri2<- ifelse(df$a0901103 == 'Yes', 1, ifelse(df$a0901103 == 'No', 0, NA))
df$pneumonia_tri2<- ifelse(df$a0901303 == 'Yes', 1, ifelse(df$a0901303 == 'No', 0, NA))
df$dermatitis_tri2<- ifelse(df$a0902303 == 'Yes', 1, ifelse(df$a0902303 == 'No', 0, NA))
df$herpeszoster_tri2<- ifelse(df$a0902903 == 'Yes', 1, ifelse(df$a0902903 == 'No', 0, NA))
df$enteritis_tri2<- ifelse(df$a0903103 == 'Yes', 1, ifelse(df$a0903103 == 'No', 0, NA))
df$cystitis_tri2<- ifelse(df$a0903303 == 'Yes', 1, ifelse(df$a0903303 == 'No', 0, NA))
df$jaundice_tri2<- ifelse(df$a0903503 == 'Yes', 1, ifelse(df$a0903503 == 'No', 0, NA))
df$otherinf_tri2 <- ifelse(df$a0903703 == 'Yes', 1, ifelse(df$a0903703 == 'No', 0, NA))
df$STD_tri2 <- ifelse(df$GR1003_Chlamydia_or_other_STD == 'Yes', 1, ifelse(df$GR1003_Chlamydia_or_other_STD == 'No', 0, NA))
df$lower_resp_inf_tri2 <- ifelse(df$GR1003_Lower_respiratory_infection == 'Yes', 1, ifelse(df$GR1003_Lower_respiratory_infection == 'No', 0, NA))
df$upper_resp_inf_tri2 <- ifelse(df$GR1003_Upper_respiratory_infection == 'Yes', 1, ifelse(df$GR1003_Upper_respiratory_infection == 'No', 0, NA))
df$GI_inf_tri2 <- ifelse(df$GR1003_Gastro_intestinal_infection == 'Yes', 1, ifelse(df$GR1003_Gastro_intestinal_infection == 'No', 0, NA))
df$UWI_tri2 <- ifelse(df$GR1003_UWI == 'Yes', 1, ifelse(df$GR1003_UWI == 'No', 0, NA))

df$fever_tri3 <- ifelse(df$b0900105 == 'Yes', 1, ifelse(df$b0900105 == 'No', 0, NA))
df$flu_tri3 <- ifelse(df$b0900305 == 'Yes', 1, ifelse(df$b0900305 == 'No', 0, NA))
df$pharyngitis_tri3<- ifelse(df$b0900505 == 'Yes', 1, ifelse(df$b0900505 == 'No', 0, NA))
df$rhinitis_tri3<- ifelse(df$b0900705 == 'Yes', 1, ifelse(df$b0900705 == 'No', 0, NA))
df$sinusitis_tri3<- ifelse(df$b0900905 == 'Yes', 1, ifelse(df$b0900905 == 'No', 0, NA))
df$earinf_tri3<- ifelse(df$b0901105 == 'Yes', 1, ifelse(df$b0901105 == 'No', 0, NA))
df$pneumonia_tri3<- ifelse(df$b0901305 == 'Yes', 1, ifelse(df$b0901305 == 'No', 0, NA))
df$dermatitis_tri3<- ifelse(df$b0902305 == 'Yes', 1, ifelse(df$b0902305 == 'No', 0, NA))
df$herpeszoster_tri3<- ifelse(df$b0902905 == 'Yes', 1, ifelse(df$b0902905 == 'No', 0, NA))
df$enteritis_tri3<- ifelse(df$b0903105 == 'Yes', 1, ifelse(df$b0903105 == 'No', 0, NA))
df$cystitis_tri3<- ifelse(df$b0903305 == 'Yes', 1, ifelse(df$b0903305 == 'No', 0, NA))
df$jaundice_tri3<- ifelse(df$b0903505 == 'Yes', 1, ifelse(df$b0903505 == 'No', 0, NA))
df$otherinf_tri3 <- ifelse(df$b0903705 == 'Yes', 1, ifelse(df$b0903705 == 'No', 0, NA))
df$STD_tri3 <- ifelse(df$GR1005_Chlamydia_or_other_STD == 'Yes', 1, ifelse(df$GR1005_Chlamydia_or_other_STD == 'No', 0, NA))
df$lower_resp_inf_tri3 <- ifelse(df$GR1005_Lower_respiratory_infection == 'Yes', 1, ifelse(df$GR1005_Lower_respiratory_infection == 'No', 0, NA))
df$upper_resp_inf_tri3 <- ifelse(df$GR1005_Upper_respiratory_infection == 'Yes', 1, ifelse(df$GR1005_Upper_respiratory_infection == 'No', 0, NA))
df$GI_inf_tri3 <- ifelse(df$GR1005_Gastro_intestinal_infection == 'Yes', 1, ifelse(df$GR1005_Gastro_intestinal_infection == 'No', 0, NA))
df$UWI_tri3 <- ifelse(df$GR1005_UWI == 'Yes', 1, ifelse(df$GR1005_UWI == 'No', 0, NA))

#remove original variables 

removevars <- c("c1100101",                           
"c1100301",
"c1100501",                           
"c1100701",
"c1100901",                           
"c1101101",
"c1101301",                           
"c1102301",
"c1102901",                          
"c1103101",
"c1103301",                          
"c1103501",
"c1103701",                          
"a0900103",
"a0900303",                          
"a0900503",
"a0900703",                          
"a0900903",
"a0901103",                          
"a0901303",
"a0902303",                          
"a0902903",
"a0903103",                          
"a0903303",
"a0903503",                          
"a0903703",
"b0900105",                          
"b0900305",
"b0900505",                          
"b0900705",
"b0900905",                          
"b0901105",
"b0901305",                          
"b0902305",
"b0902905",                          
"b0903105",
"b0903305",                          
"b0903505",
"b0903705",
"Chlamydia_or_other_STD",
"Lower_respiratory_infection",
"Upper_respiratory_infection",
"Gastro_intestinal_infection",
"UWI",
"GR1001_Chlamydia_or_other_STD",
"GR1001_Lower_respiratory_infection",
"GR1001_Upper_respiratory_infection",
"GR1001_Gastro_intestinal_infection",
"GR1001_UWI",
"GR1001_Other",
"GR1003_Chlamydia_or_other_STD",
"GR1003_Lower_respiratory_infection",
"GR1003_Upper_respiratory_infection",
"GR1003_Gastro_intestinal_infection",
"GR1003_UWI",
"GR1003_Other",
"GR1003_Any_other_infection",
"GR1005_Chlamydia_or_other_STD",
"GR1005_Lower_respiratory_infection",
"GR1005_Upper_respiratory_infection",
"GR1005_Gastro_intestinal_infection",
"GR1005_UWI",
"GR1005_Other",
"GR1005_Any_other_infection")

df <- dplyr::select(df,!(all_of(removevars)))

## LIFESTYLE PROINFLAMMATORY FACTORS SCORE (LPFS) - make variables binary 

#medications
df$SSRI <- ifelse(df$SSRITOT == 'During pregnancy-first-second-third trimester' | df$SSRITOT == 'Early pregnancy only', 1, ifelse(df$SSRITOT == 'no use' | df$SSRITOT == 'only before pregnancy', 0, NA)) 
df$TRIP <- ifelse(df$TRIPTOT == 'During pregnancy-first-second-third trimester' | df$TRIPTOT == 'Early pregnancy only', 1, ifelse(df$TRIPTOT == 'no use' | df$TRIPTOT == 'only before pregnancy', 0, NA)) 
df$PSY <- ifelse(df$PSYTOT == 'During pregnancy-first-second-third trimester' | df$PSYTOT == 'Early pregnancy only', 1, ifelse(df$PSYTOT == 'no use' | df$PSYTOT == 'only before pregnancy', 0, NA)) 
df$TCA <- ifelse(df$TCATOT == 'During pregnancy-first-second-third trimester' | df$TCATOT == 'Early pregnancy only', 1, ifelse(df$TCATOT == 'no use' | df$TCATOT == 'only before pregnancy', 0, NA)) 
df$NSAID <- ifelse(df$NSAIDTOT == 'During pregnancy-first-second-third trimester' | df$NSAIDTOT == 'Early pregnancy only', 1, ifelse(df$NSAIDTOT == 'no use' | df$NSAIDTOT == 'only before pregnancy', 0, NA)) 
df$ABIO <- ifelse(df$ABIOTOT == 'During pregnancy-first-second-third trimester' | df$ABIOTOT == 'Early pregnancy only', 1, ifelse(df$ABIOTOT == 'no use' | df$ABIOTOT == 'only before pregnancy', 0, NA)) 
df$PMOL <- ifelse(df$PMOLTOT == 'During pregnancy-first-second-third trimester' | df$PMOLTOT == 'Early pregnancy only', 1, ifelse(df$PMOLTOT == 'no use' | df$PMOLTOT == 'only before pregnancy', 0, NA)) 
df$CORT <- ifelse(df$CORTTOT == 'During pregnancy-first-second-third trimester' | df$CORTTOT == 'Early pregnancy only', 1, ifelse(df$CORTTOT == 'no use' | df$CORTTOT == 'only before pregnancy', 0, NA)) 
df$MUCO <- ifelse(df$MUCOTOT == 'During pregnancy-first-second-third trimester' | df$MUCOTOT == 'Early pregnancy only', 1, ifelse(df$MUCOTOT == 'no use' | df$MUCOTOT == 'only before pregnancy', 0, NA)) 
df$COUGH <- ifelse(df$COUGHTOT == 'During pregnancy-first-second-third trimester' | df$COUGHTOT == 'Early pregnancy only', 1, ifelse(df$COUGHTOT == 'no use' | df$COUGHTOT == 'only before pregnancy', 0, NA)) 

removevars2 <- c("SSRITOT",
"TRIPTOT",
"PSYTOT",
"TCATOT",
"NSAIDTOT",
"ABIOTOT",
"PMOLTOT",
"CORTTOT",
"MUCOTOT",
"COUGHTOT")

df <- dplyr::select(df,!(all_of(removevars2)))

#other lifestyle factors

df$smoke_score <- ifelse(df$SMOKE_ALL == 'smoked until pregnancy was known' | df$SMOKE_ALL == 'continued smoking in pregnancy', 1, ifelse(df$SMOKE_ALL == 'never smoked during pregnancy', 0, NA)) #sum score smoking in which we combine pre-preg and prenatal smoking 

df$bmi_score <- ifelse(df$BMI_0 > 30, 1, ifelse(df$BMI_0 == "NA's", NA, 0)) #bmi score based on bmi=30 cut off 

df$diet_score <- ifelse(df$DietScore_pregnancy < 6.423, 1, ifelse(df$DietScore_pregnancy == "NA's", NA, 0)) #diet score, 0-15 points, higher is better, so cutoff < 1st quartile 

### PREGNANCY-RElATED INFLAMMATORY CONDITIONS SCORE (PICS) - make variables binary ###

df$HELLP <- ifelse(df$PE_total == 'HELLP' | df$PE_total == 'PE en HELLP' , 1, ifelse(df$PE_total == "NA's", NA, 0))
df$gest_diab <- ifelse(df$DIAB_GRA == 'Yes', 1, ifelse(df$DIAB_GRA == 'No', 0, NA))
df$preg_hypt <- ifelse(df$PIH_v1 == 'Yes', 1, ifelse(df$PIH_v1 == 'No', 0, NA))
df$preclampsia <- ifelse(df$PE == 'Yes', 1, ifelse(df$PE == 'No', 0, NA))
df$promm <- ifelse(df$PROMM == 'Yes', 1, ifelse(df$PROMM == 'No', 0, NA))
df$sectio <- ifelse(df$UITDRIJF == 'prim. sectio' | df$UITDRIJF == 'sec. sectio' | df$UITDRIJF == 'sectio (onbekend prim./sec.)', 1, ifelse(df$UITDRIJF == "NA's", NA, 0))

removevars4 <- c("PE_total", "DIAB_GRA", "PIH_v1", "PE", "PROMM", "UITDRIJF") 
df <- dplyr::select(df,!(all_of(removevars4))) 

### MEDICAL INFLAMMATORY CONIDTIONS SCORE (MICS) - make variables binary ###
df$intestinal <- ifelse(df$d2100101 == 'Yes', 1, ifelse(df$d2100101 == 'No',0,NA))
df$sle <- ifelse(df$d2300101 == 'Yes', 1, ifelse(df$d2300101 == 'No',0,NA))
df$arthritis <- ifelse(df$d2400101 == 'Yes', 1, ifelse(df$d2400101 == 'No',0,NA))
df$ms <- ifelse(df$d2500101 == 'Yes', 1, ifelse(df$d2500101 == 'No',0,NA))
df$thyroid <- ifelse(df$d2600101 == 'Yes', 1, ifelse(df$d2600101 == 'No',0,NA))
df$diabetes <- ifelse(df$d1100101 == 'Yes', 1, ifelse(df$d1100101 == 'No',0,NA))

removevars5 <- c("d2100101", "d2300101", "d2400101", "d2500101", "d2600101", "d1100101")
df <- dplyr::select(df,!(all_of(removevars5))) 


##make changes to data type where needed 
df$HsCRPmgL_g1 <- as.numeric(df$HsCRPmgL_g1) 
df$HsCRPmgL_g2 <- as.numeric(df$HsCRPmgL_g2) 

##remove other variables that are not needed 
names(df)
removevars6 <- c("suffix_HsCRPmgL_g1","suffix_HsCRPmgL_g2","SERUMCHILD5","MOTHER","plasmaCordblood", 
"suffix_CRP_birth", "RemoveData","TWIN",#"EDUCM", 
"visitChildF5", "first_sibling_most_data","DietScore_pregnancy",
"dwi_nvols_f09","dwi_nvols_f13", "dti_overall_qc_f09", "dti_overall_qc_f13", "t1_asset_has_nii_f09",
"missing", "exc_smri_f09", "exc_smri_f13", "exc_dti_f09", "exc_dti_f13", "t1_has_nii_f13", "dti_has_nii_f13", "mri_consent_f09")


df <- dplyr::select(df,!(all_of(removevars6))) 


#### Make the scale scores for those with full data ####

## save items needed to string 
PSS_LE <- c("family_member_died", "friend_relative_died", "family_member_ill_pregnancy", "admitted_to_hospital", "health", "unemployed", 
            "work_study_problems", "moved_house", "blood_loss", "examination", "baby_worried", "pregnancy_worried", "obstetric_care", 
            "pregnancy_planned", "victim_robbery")
PSS_CR <- c("financial_problems","trouble_pay_pregnancy","income_reduced","housing_defects","housing_adequacy","housing_basic_living","m_education_pregnancy", "p_education_pregnancy")
PSS_PR <- c("early_pregnancy","m_depression_pregnancy","m_anxiety_pregnancy","m_interp_sensitivity_pregnancy","p_depression_pregnancy","p_anxiety_pregnancy","p_interp_sensitivity_pregnancy",
            "m_violence_people","m_violence_property","m_criminal_record","p_criminal_record")
PSS_IR <- c("difficulties_contacts", "difficulties_partner", "difficulties_family_friend", "marital_status_pregnancy", "divorce_pregnancy", "family_support", 
            "family_acceptance", "family_affection", "family_acception","family_trust", "family_painful_feelings", "family_decisions", "family_conflict", "family_decisions_problems", 
            "family_plans", "family_talk_sadness", "family_talk_worries", "family_size_pregnancy")

PSS_domains <- c('PSS_LE_score', 'PSS_CR_score', 'PSS_PR_score', 'PSS_IR_score')

## for PRENATAL INFECTION SCORE (PIS)

PIS_tri1 <- c("fever_tri1","flu_tri1","pharyngitis_tri1","rhinitis_tri1","sinusitis_tri1","earinf_tri1","pneumonia_tri1","dermatitis_tri1",
#"herpeszoster_tri1",removed as all zero in dataset
"enteritis_tri1",
"cystitis_tri1",
#"jaundice_tri1", removed as all zero in dataset
"otherinf_tri1","STD_tri1","lower_resp_inf_tri1","upper_resp_inf_tri1","GI_inf_tri1","UWI_tri1")

PIS_tri2 <- c("fever_tri2",
"flu_tri2",
"pharyngitis_tri2",
"rhinitis_tri2",
"sinusitis_tri2",
"earinf_tri2",
"pneumonia_tri2",
"dermatitis_tri2",
"herpeszoster_tri2",
"enteritis_tri2",
"cystitis_tri2",
#"jaundice_tri2", removed as all zero in dataset
"otherinf_tri2",
"STD_tri2",
"lower_resp_inf_tri2",
"upper_resp_inf_tri2",
"GI_inf_tri2",
"UWI_tri2")

PIS_tri3 <- c("fever_tri3",
"flu_tri3",
"pharyngitis_tri3",
"rhinitis_tri3",
"sinusitis_tri3",
"earinf_tri3",
"pneumonia_tri3",
"dermatitis_tri3",
"herpeszoster_tri3",
"enteritis_tri3",
"cystitis_tri3",
#"jaundice_tri3", removed as all zero in dataset
"otherinf_tri3",
"STD_tri3",
"lower_resp_inf_tri3",
"upper_resp_inf_tri3",
"GI_inf_tri3",
"UWI_tri3")

PIS_T1 <- c("URI_T1", "LRI_T1", "UTI_T1", "GII_T1", "flu_tri1", "dermatitis_tri1", "STD_tri1", "fever_tri1")

PIS_T2 <- c("URI_T2", "LRI_T2", "UTI_T2", "GII_T2", "flu_tri2", "dermatitis_tri2", "herpeszoster_tri2", "STD_tri2", "fever_tri2")

PIS_T3 <- c("URI_T3", "LRI_T3", "UTI_T3", "GII_T3", "flu_tri3", "dermatitis_tri3", "herpeszoster_tri3", "STD_tri3", "fever_tri3")

PIS_domains <- c("PIS_T1_score", "PIS_T2_score", "PIS_T3_score")

## for LIFESTYLE PROINFLAMMATORY FACTORS SCORE (LPFS)

LPFS_meds <- c("SSRI",
"TRIP",
"PSY",
"TCA",
"NSAID",
"ABIO",
"PMOL",
"CORT",
"MUCO",
"COUGH")

LPFS_med_domains <- c("psymeds", "CORT", "inflams")

LPFS_lifestyle <- c("smoke_score", "bmi_score", "diet_score")

LPFS_domains <- c("smoke_score", "bmi_score", "diet_score", "psymeds", "CORT", "inflams")

## for PREGNANCY-RElATED INFLAMMATORY CONDITIONS SCORE (PICS) 

PICS_all <- c("HELLP", "gest_diab", "preg_hypt", "preclampsia", "promm", "sectio")

PICS_score <- c("PICS_score")

## for MEDICAL INFLAMMATORY CONIDTIONS SCORE (MICS) 

MICS_all <- c("intestinal", "arthritis", "ms", "thyroid", "diabetes") #"sle", removed as all zero in dataset

MICS_score <- c("MICS_score")



## Compute scores on UNIMPUTED data - does not calculate for those with NA ##

# domain score function
domainscore <- function(df, method){
  
  # a data frame consisting of all the included items is constructed
  df <- data.frame(df)
  
  # check if all variables in df are dichotomized
  for (i in 1:ncol(df)){
    if (range(df[,i], na.rm = T)[1] == 1 & range(df[,i], na.rm = T)[2] == 2){
      df[,i] <- df[,i] - 1}
    else {
      if (range(df[,i], na.rm = T)[1] != 0 | range(df[,i], na.rm = T)[2] != 1){
        stop('Items are not dichotomized')}
    }
  }
  
  # calculate the mean number of events (weighted)
  if (method == 'mean') temp <- rowMeans(df, na.rm = F)
  
  # calculate the sum of events (weighted)
  if (method == 'cumulative') temp <- rowSums(df, na.rm = F)
  
  return(temp)
  
}

## PSS

df$PSS_LE_score <- domainscore(df[,PSS_LE], method = 'mean')

df$PSS_CR_score <- domainscore(df[,PSS_CR], method = 'mean')

df$PSS_PR_score <- domainscore(df[,c(PSS_PR, "m_age")], method = 'mean') #nb. m_age is kept sep otherwise collinearity when using PSS_PR imp 

df$PSS_IR_score <- domainscore(df[,PSS_IR], method = 'mean')

df$PSS_TOT_score <- rowSums(df[,PSS_domains], na.rm = F) #cumulative sumscore 

summary(df[,PSS_domains])
summary(df$PSS_TOT_score)
hist(df$PSS_TOT_score)

## PIS

##first make clusters
#Upper respiratory tract infection
df$URI_T1 <- ifelse(df$upper_resp_inf_tri1 == 1 | df$rhinitis_tri1 == 1 | df$pharyngitis_tri1 == 1 | df$earinf_tri1 == 1 | df$sinusitis_tri1 == 1, 1, (ifelse(df$upper_resp_inf_tri1 == 0 & df$rhinitis_tri1 == 0 & df$pharyngitis_tri1 == 0 & df$earinf_tri1 == 0 & df$sinusitis_tri1 == 0, 0, NA)))

df$URI_T2 <- ifelse(df$upper_resp_inf_tri2 == 1 | df$rhinitis_tri2 == 1 | df$pharyngitis_tri2 == 1 | df$earinf_tri2 == 1 | df$sinusitis_tri2 == 1, 1, (ifelse(df$upper_resp_inf_tri2 == 0 & df$rhinitis_tri2 == 0 & df$pharyngitis_tri2 == 0 & df$earinf_tri2 == 0 & df$sinusitis_tri2 == 0, 0, NA)))

df$URI_T3 <- ifelse(df$upper_resp_inf_tri3 == 1 | df$rhinitis_tri3 == 1 | df$pharyngitis_tri3 == 1 | df$earinf_tri3 == 1 | df$sinusitis_tri3 == 1, 1, (ifelse(df$upper_resp_inf_tri3 == 0 & df$rhinitis_tri3 == 0 & df$pharyngitis_tri3 == 0 & df$earinf_tri3 == 0 & df$sinusitis_tri3 == 0, 0, NA)))


#Lower respiratory tract infection
df$LRI_T1 <- ifelse(df$lower_resp_inf_tri1 == 1 | df$pneumonia_tri1 == 1, 1, (ifelse(df$lower_resp_inf_tri1 == 0 & df$pneumonia_tri1 == 0, 0, NA)))
df$LRI_T2 <- ifelse(df$lower_resp_inf_tri2 == 1 | df$pneumonia_tri2 == 1, 1, (ifelse(df$lower_resp_inf_tri2 == 0 & df$pneumonia_tri2 == 0, 0, NA)))
df$LRI_T3 <- ifelse(df$lower_resp_inf_tri3 == 1 | df$pneumonia_tri3 == 1, 1, (ifelse(df$lower_resp_inf_tri3 == 0 & df$pneumonia_tri3 == 0, 0, NA)))

#Cleaning cystitis
df$UTI_T1 <- ifelse(df$UWI_tri1 == 1 | df$cystitis_tri1 == 1, 1, (ifelse(df$UWI_tri1 == 0 & df$cystitis_tri1 == 0, 0, NA)))
df$UTI_T2 <- ifelse(df$UWI_tri2 == 1 | df$cystitis_tri2 == 1, 1, (ifelse(df$UWI_tri2 == 0 & df$cystitis_tri2 == 0, 0, NA)))
df$UTI_T3 <- ifelse(df$UWI_tri3 == 1 | df$cystitis_tri3 == 1, 1, (ifelse(df$UWI_tri3 == 0 & df$cystitis_tri3 == 0, 0, NA)))

#Gastrointestinal tract infection 
df$GII_T1 <- ifelse(df$GI_inf_tri1 ==1 | df$enteritis_tri1 ==1, 1, (ifelse(df$GI_inf_tri1 ==0 & df$enteritis_tri1 ==0,0, NA)))
df$GII_T2 <- ifelse(df$GI_inf_tri2 ==1 | df$enteritis_tri2 ==1, 1, (ifelse(df$GI_inf_tri2 ==0 & df$enteritis_tri2 ==0,0, NA)))
df$GII_T3 <- ifelse(df$GI_inf_tri3 ==1 | df$enteritis_tri3 ==1, 1, (ifelse(df$GI_inf_tri3 ==0 & df$enteritis_tri3 ==0,0, NA)))


#then the trimester infection scores
df$PIS_T1_score <- domainscore(df[,PIS_T1], method = 'mean')
df$PIS_T2_score <- domainscore(df[,PIS_T2], method = 'mean')
df$PIS_T3_score <- domainscore(df[,PIS_T3], method = 'mean')


#then the cumulative infection score 
df$PIS_TOT_score <- rowSums(df[,PIS_domains], na.rm = F) #cumulative sumscore 

## LPFS
#first make the medicine domains 
df$psymeds <- ifelse(df$SSRI == 1 | df$TRIP == 1 | df$PSY == 1 | df$TCA == 1, 1, (ifelse(df$SSRI == 0 & df$TRIP == 0 & df$PSY == 0 & df$TCA == 0, 0, NA)))
 
df$inflams <- ifelse(df$NSAID == 1 | df$ABIO == 1 | df$PMOL == 1 | df$MUCO == 1 | df$COUGH == 1, 1,
(ifelse(df$NSAID == 0 & df$ABIO == 0 & df$PMOL == 0 & df$MUCO == 0 & df$COUGH == 0, 0, NA))) #CORT removed so that it is not included in more than one subscale 

#make the LPFS score 

df$LPFS_score <- domainscore(df[,c(LPFS_lifestyle,LPFS_med_domains)], method = 'cumulative')


## PICS
df$PICS_score <- domainscore(df[,PICS_all], method = 'cumulative')

## MICS
df$MICS_score <- domainscore(df[,MICS_all], method = 'cumulative')

## check the scores before imputation and save summary ##
myscores <- c("PSS_TOT_score","PIS_TOT_score","LPFS_score","PICS_score","MICS_score")
summary(df[,myscores])
summary_stats <- file("/home/n.creasey/episcore_CRP_project/shadow2/output/preimp_scores_summary.txt")
sink(summary_stats, append=TRUE, type = "output")
summary(df[,myscores])
unlink(summary_stats)
sink()

#### save the final data 
save(df, file="path_to_store_results")

# END OF THIS SCRIPT, GO TO PART C
