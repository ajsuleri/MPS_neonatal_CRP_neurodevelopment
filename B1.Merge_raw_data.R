#########################################################################
######################PROJECT: Neonatal inflammation####################
#########################################################################

#Project: Mapping prenatal predictors and neurobehavioral outcomes of an epigenetic marker of neonatal inflammation - a longitudinal population-based study 

#SET UP WORKSPACE 
rm(list = ls()) #clears the environment

libraries <- c('foreign', 'haven', 'tidyverse', 'broom.mixed', 'xlsx', 'corrplot', 'stringi', 'miceadds', 'mitools', 'lme4', 'lmerTest', 'rlang', 'lattice', 'ggplot2', 'ggeffects', 'splines', 'reshape2', 'gridExtra', 'grid', 'ggpubr','CBPS', 'survey', 'survival', 'ggseg', 'fmsb', 'ggcorrplot', 'Hmisc', 'pROC', 'readxl', 'pROC', 'RColorBrewer', 'DescTools') 

invisible(lapply(libraries, require, character.only = T))

#set directory 
setwd("path_to_data")

#### EXTRACT AND MERGE THE RAW DATA ####

#load the datasets that will be needed 
allmps <- read.csv("/home/n.creasey/episcore_CRP_project/shadow2/computed_scores/allMPS_CRP.csv", header=TRUE, sep = ",")

general <- read.spss("/home/n.creasey/episcore_CRP_project/shadow2/rawdata/CHILD-ALLGENERALDATA_24102022.sav", use.value.labels = T, to.data.frame = T)

smoking <- read.spss("/home/n.creasey/episcore_CRP_project/shadow2/rawdata/MATERNALSMOKING_22112016.sav", use.value.labels = T, to.data.frame = T)

epicselec <- read.spss("/home/n.creasey/episcore_CRP_project/shadow2/rawdata/Selection_GENR_MethylEPIC_release1_birth_20230717.sav", use.value.labels = T, to.data.frame = T)

selec450 <- read.spss("/home/n.creasey/episcore_CRP_project/shadow2/rawdata/Selection_GENR_450kmeth_release3_birth_20180125.sav", use.value.labels = T, to.data.frame = T)

load.Rdata("/home/n.creasey/gc_epgs_EPIC/GENR_EPICv1METH_birth_CellTypes_combined.RData", "cellcounts_EPIC")
cellcounts_450 <- read.csv("/home/n.creasey/episcore_CRP_project/shadow2/rawdata/GENR_450kmeth_release3_birth_Salas.csv")

gpcv3 <- read.spss("/home/n.creasey/episcore_CRP_project/shadow2/rawdata/PCA_Selection GWAv3_revised def_European-October2022.sav", use.value.labels = T, to.data.frame = T)

gpcv4 <- read.spss("/home/n.creasey/episcore_CRP_project/shadow2/rawdata/PCA_Selection GWAv4_revised def_European-October2022.sav", use.value.labels = T, to.data.frame = T)

idc_idm <- read.spss("/home/n.creasey/episcore_CRP_project/shadow2/rawdata/IDC-IDM-MOTHER.sav", use.value.labels = T, to.data.frame = T)

crp_motherpreg <- read.spss("/home/n.creasey/episcore_CRP_project/shadow2/rawdata/MOTHERPREGNANCY-CRP_06062023.sav", use.value.labels = T, to.data.frame = T)
crp_cordblood <- read.spss("/home/n.creasey/episcore_CRP_project/Score_val/CHILDCORDBLOOD-CRP_23092013.sav", use.value.labels = T, to.data.frame = T)
crp_child5 <- read.spss("/home/n.creasey/episcore_CRP_project/Score_val/CHILDCRP5_24092013.sav", use.value.labels = T, to.data.frame = T)

inflam_somatic_df <- read.spss("/home/n.creasey/episcore_CRP_project/shadow2/rawdata/GR1001-D1-37_16072021.sav", use.value.labels = T, to.data.frame = T)
preg_inflam_df1 <- read.spss("/home/n.creasey/episcore_CRP_project/shadow2/rawdata/MATERNALCOMPLICATIONS_22112016.sav", use.value.labels = T, to.data.frame = T)
preg_inflam_df2 <- read.spss("/home/n.creasey/episcore_CRP_project/shadow2/rawdata/PARTUS_17072015.sav", use.value.labels = T, to.data.frame = T)

prenatal_medication <- read.spss("/home/n.creasey/episcore_CRP_project/shadow2/rawdata/MEDICATIONSELFREPORTPREGNANCY_30112017.sav", use.value.labels = T, to.data.frame = T)
motherdiet <- read.spss("/home/n.creasey/episcore_CRP_project/shadow2/rawdata/MOTHERDietScore_pregnancy_18052018.sav", use.value.labels = T, to.data.frame = T)

preinfect1 <- read.spss("/home/n.creasey/episcore_CRP_project/shadow2/rawdata/GR1001-C2-13_22112016.sav", use.value.labels = T, to.data.frame = T)

preinfect2 <- read.spss("/home/n.creasey/episcore_CRP_project/shadow2/rawdata/GR1003-A8-11_02092013.sav", use.value.labels = T, to.data.frame = T)

preinfect3 <- read.spss("/home/n.creasey/episcore_CRP_project/shadow2/rawdata/GR1005-B8-11_22112016.sav", use.value.labels = T, to.data.frame = T)

preinfect4 <- read.spss("/home/n.creasey/episcore_CRP_project/shadow2/rawdata/20220412_Anna-GR1001_C11_Otherinfections.sav", use.value.labels = T, to.data.frame = T)

preinfect5<- read.spss("/home/n.creasey/episcore_CRP_project/shadow2/rawdata/PrenatalQuest_Otherinfections_25042022.sav", use.value.labels = T, to.data.frame = T)

prenatal_stress <- read.spss("/home/n.creasey/gc_epgs/dataset_NicoleCreasey_updatedGWA_20221201.sav", use.value.labels = T, to.data.frame = T)

cbcl1 <- read.spss("/home/n.creasey/episcore_CRP_project/shadow2/rawdata/20090326_GR1029-CBCL A1 (18m).sav", use.value.labels = T, to.data.frame = T)
cbcl3 <- read.spss("/home/n.creasey/episcore_CRP_project/shadow2/rawdata/CBCL_3_incl_Tscores__GR1065E2_GR1066A1_20201111.sav", use.value.labels = T, to.data.frame = T)
cbcl6 <- read.spss("/home/n.creasey/episcore_CRP_project/shadow2/rawdata/CHILDCBCL_6_incl_Tscores_20201111.sav", use.value.labels = T, to.data.frame = T)
cbcl9 <- read.spss("/home/n.creasey/episcore_CRP_project/shadow2/rawdata/CHILDCBCL9_incl_Tscores_20201111.sav", use.value.labels = T, to.data.frame = T)
cbcl14 <- read.spss("/home/n.creasey/episcore_CRP_project/shadow2/rawdata/GR1093-E1_CBCL_18062020.sav", use.value.labels = T, to.data.frame = T)


#prepare the datasets 

#arrange the epic technical covs and IDs
epicselec <- merge(epicselec, cellcounts_EPIC, by= 'SampleID', all=TRUE)
epicselec$EwasChildMethylEpic <- NULL
colnames(epicselec)[1] ="Sample_ID"
epicselec$Sample_Plate <- as.character(as.numeric(epicselec$Sample_Plate)+28)


#arrange the 450k technical covs and IDs

selec450$Sample_ID <- gsub(" ","", selec450$Sample_ID, fixed = T) #fix the issue with 450k sample names
colnames(selec450)[1] ="Sample_ID"
colnames(cellcounts_450)[1] ="Sample_ID"
selec450 <- merge(selec450, cellcounts_450, by='Sample_ID', all=TRUE)
keepvars<- colnames(epicselec)
selec450 <- selec450[,keepvars]
selec450 <- selec450[rowSums(is.na(selec450)) == 0,]

#merge the technical cov and ID info (different ppts so rbind)
identical(colnames(epicselec),colnames(selec450)) #check column names match
techcovs <- rbind(epicselec, selec450)

# add the genetic pcs 
keepvars <- c("IDC","C1", "C2", "C3","C4")
gpcv3 <- gpcv3[,keepvars]
gpcv4 <- gpcv4[,keepvars]
identical(colnames(gpcv3),colnames(gpcv4)) #check column names match
gpcs <- rbind(gpcv3, gpcv4)
techcovs <- merge(techcovs, gpcs, by='IDC', all=TRUE)

#now merge the mps and techcov data
df <- merge(allmps,techcovs,by="Sample_ID", all=TRUE)

#get the CRP data ready for merging and merge 
crp_motherpreg <- merge(crp_motherpreg,idc_idm,by="IDM", all=TRUE)
allcrp <- merge(crp_motherpreg,crp_cordblood,by="IDC", all=TRUE)
allcrp <- merge(allcrp, crp_child5,by="IDC", all=TRUE)
allcrp$GESTBIR <- NULL
df <- merge(df, allcrp,by="IDC", all=TRUE)

#get the other covariates (child sex, maternal age at birth, maternaL EDUCATION, HOUSEHOLD INCOME, MATERNAL SMOKING, PARITY, CHILD AGE, plus BMI0 for inflam score, plus twin and removedata for exclusion)
covs <- dplyr::select(general, 'IDC', 'RemoveData', 'TWIN', 'EDUCM', 'AGE_M_v2', 'GESTBIR', "WEIGHT", 'PARITY', 'GENDER','INCOME','BMI_0')

df <- merge(df, covs,by="IDC", all=TRUE)

smoking <- dplyr::select(smoking, 'idm','SMOKE_ALL')
colnames(smoking)[1] <- "IDM"
df <- merge(df, smoking,by="IDM", all=TRUE)

#get the child behavior scores and merge 
cbcl1 <- dplyr::select(cbcl1, 'IDC', 'AGE18M', 'sum_int', 'sum_ext', 'cbcl_sum')
colnames(cbcl1)[3] <- "sum_int_18m"
colnames(cbcl1)[4] <- "sum_ext_18m"
colnames(cbcl1)[5] <- "cbcl_sum_18m"
cbcl3 <- dplyr::select(cbcl3, 'IDC', 'sum_int_36m', 'sum_ext_36m', 'cbcl_sum_36m', 'age_GR1065', 'age_GR1066')
cbcl6 <- dplyr::select(cbcl6, 'IDC', 'agechild_GR1075', 'sum_int_5', 'sum_ext_5', 'cbcl_sum_5')
cbcl9 <- dplyr::select(cbcl9, 'IDC', 'AgeChild_CBCL9m', 'sum_int_9m', 'sum_ext_9m', 'cbcl_sum_9m')
cbcl14 <- dplyr::select(cbcl14, 'IDC', 'AGECHILD_GR1093','sum_int_14', 'sum_ext_14', 'cbcl_sum_14')

df <- merge(df, cbcl1,by="IDC", all=TRUE)
df <- merge(df, cbcl3,by="IDC", all=TRUE)
df <- merge(df, cbcl6,by="IDC", all=TRUE)
df <- merge(df, cbcl9,by="IDC", all=TRUE)
df <- merge(df, cbcl14,by="IDC", all=TRUE)


#prenatal infection variables 
preinfect1 <- dplyr::select(preinfect1, IDM, c1100101,c1100301,c1100501,c1100701,c1100901,c1101101,c1101301,c1102301,c1102901, c1103101, c1103301, c1103501, c1103701) #trimester 1
preinfect2 <- dplyr::select(preinfect2, IDM, a0900103, a0900303, a0900503, a0900703, a0900903, a0901103, a0901303, a0902303, a0902903, a0903103, a0903303, a0903503, a0903703) #trimester 2
preinfect3 <- dplyr::select(preinfect3, IDM, b0900105, b0900305, b0900505, b0900705, b0900905, b0901105, b0901305, b0902305, b0902905, b0903105, b0903305, b0903505, b0903705) #trimester 3

df <- merge(df, preinfect1,by="IDM", all=TRUE)
df <- merge(df, preinfect2,by="IDM", all=TRUE)
df <- merge(df, preinfect3,by="IDM", all=TRUE)
df <- merge(df, preinfect4,by="IDM", all=TRUE)
df <- merge(df, preinfect5,by="IDM", all=TRUE)


#prenatal stress variables
prestress_variables <- c("IDC", "m_age", "family_member_died", "friend_relative_died", "family_member_ill_pregnancy", "admitted_to_hospital", "health", "unemployed", 
            "work_study_problems", "moved_house", "blood_loss", "examination", "baby_worried", "pregnancy_worried", "obstetric_care", 
            "pregnancy_planned", "victim_robbery","financial_problems","trouble_pay_pregnancy","income_reduced","housing_defects","housing_adequacy","housing_basic_living","m_education_pregnancy", "p_education_pregnancy", "early_pregnancy","m_depression_pregnancy","m_anxiety_pregnancy","m_interp_sensitivity_pregnancy","p_depression_pregnancy","p_anxiety_pregnancy","p_interp_sensitivity_pregnancy",
            "m_violence_people","m_violence_property","m_criminal_record","p_criminal_record","difficulties_contacts", "difficulties_partner", "difficulties_family_friend", "marital_status_pregnancy", "divorce_pregnancy", "family_support", 
            "family_acceptance", "family_affection", "family_acception","family_trust", "family_painful_feelings", "family_decisions", "family_conflict", "family_decisions_problems","family_plans", "family_talk_sadness", "family_talk_worries", "family_size_pregnancy", "m_depression","p_depression","marital_status", "m_education") #also includes 4 postnatal variables for use in imputation as aux variables 

prenatal_stress <- prenatal_stress[,prestress_variables]
df <- merge(df, prenatal_stress,by="IDC", all=TRUE)

#pro-inflamatory lifestyle score variables (BMI_0 and SMOKE_ALL already included above)
prenatal_medication <- dplyr::select(prenatal_medication, c('IDM', 'SSRITOT', 'TRIPTOT', 'PSYTOT', 'TCATOT', 'NSAIDTOT', 'ABIOTOT', 'PMOLTOT', 'CORTTOT', 'MUCOTOT' ,'COUGHTOT')) 
df <- merge(df, prenatal_medication, by="IDM", all=TRUE)

motherdiet <- dplyr::select(motherdiet, c('IDM', 'DietScore_pregnancy'))
df <- merge(df, motherdiet, by="IDM", all=TRUE)

#pregnancy-related inflammatory conditions score variables
preg_inflam_df1 <- dplyr::select(preg_inflam_df1, c('IDM', 'PE_total','DIAB_GRA','PIH_v1','PE'))
preg_inflam_df2 <- dplyr::select(preg_inflam_df2, c('IDM', 'PROMM','UITDRIJF'))

df <- merge(df, preg_inflam_df1 , by="IDM", all=TRUE)
df <- merge(df, preg_inflam_df2, by="IDM", all=TRUE)

#medical inflammatory conditions score
inflam_somatic_df <- dplyr::select(inflam_somatic_df, c('IDM', 'd2100101', 'd2300101',
'd2400101','d2500101','d2600101','d1100101'))
df <- merge(df, inflam_somatic_df, by="IDM", all=TRUE)

#get and arrange the brain image data for merging 
df_brain1 <- readRDS("/home/n.creasey/episcore_CRP_project/shadow2/rawdata/f09_freesurfer_v6_09dec2016_aseg_stats_pull06june2017_v1.rds")
df_brain1 <- dplyr::select(df_brain1, c('idc', 'Left_Cerebellum_Cortex_vol_f09', 'Left_Hippocampus_vol_f09', 'Left_Amygdala_vol_f09', 'Left_Lateral_Ventricle_vol_f09', 'Brain_Stem_vol_f09','Right_Cerebellum_Cortex_vol_f09', 'Right_Amygdala_vol_f09', 'Right_Hippocampus_vol_f09', 'Right_Lateral_Ventricle_vol_f09')) #select columns we need 

df_brain2 <- readRDS("/home/n.creasey/episcore_CRP_project/shadow2/rawdata/f09_freesurfer_v6_09dec2016_tbv_stats_pull20june2017_v2.rds")
df_brain2 <- dplyr::select(df_brain2, c('idc', 'TotalGrayVol_f09', 'genr_tbv_f09', 'eTIV_f09', 'lhCerebralWhiteMatterVol_f09', 'rhCerebralWhiteMatterVol_f09')) #select columns we need 

df_brain3 <- readRDS("/home/n.creasey/episcore_CRP_project/shadow2/rawdata/f13_freesurfer_v6_14oct2020_aseg_stats_pull23Nov2020_v1.rds")
df_brain3 <- dplyr::select(df_brain3, c('idc', 'Left_Cerebellum_Cortex_vol_f13', 'Left_Hippocampus_vol_f13', 'Left_Amygdala_vol_f13', 'Left_Lateral_Ventricle_vol_f13', 'Brain_Stem_vol_f13','Right_Cerebellum_Cortex_vol_f13', 'Right_Amygdala_vol_f13', 'Right_Hippocampus_vol_f13', 'Right_Lateral_Ventricle_vol_f13')) #select columns we need 

df_brain4 <- readRDS("/home/n.creasey/episcore_CRP_project/shadow2/rawdata/f13_freesurfer_v6_14oct2020_tbv_stats_pull23Nov2020_v2.rds")
df_brain4 <- dplyr::select(df_brain4, c('idc', 'TotalGrayVol_f13', 'genr_tbv_f13', 'eTIV_f13', 'lhCerebralWhiteMatterVol_f13', 'rhCerebralWhiteMatterVol_f13')) #select columns we need 

df_brain5 <- readRDS("/home/n.creasey/episcore_CRP_project/shadow2/rawdata/f09_GenR_MRI_eddy_dipy_wls_14Feb2022_autoPtx_dti_stats_inc_glob_measV1.rds")
df_brain5 <- dplyr::select(df_brain5, c('idc', 'mean_FA_genr_f09', 'mean_MD_genr_f09'))

df_brain6 <- readRDS("/home/n.creasey/episcore_CRP_project/shadow2/rawdata/f13_GenR_MRI_eddy_dipy_wls_14Feb2022_autoPtx_dti_stats_inc_glob_measV1.rds")
df_brain6 <- dplyr::select(df_brain6, c('idc', 'mean_FA_genr_f13', 'mean_MD_genr_f13'))

df_brain7 <- readRDS("/home/n.creasey/episcore_CRP_project/shadow2/rawdata/genr_mri_core_data_20230505.rds")
df_brain7 <- dplyr::select(df_brain7, c('idc', 'mri_consent_f09', 'mri_consent_f13', 'age_child_mri_f09', 'age_child_mri_f13', 't1_has_nii_f09', 't1_has_nii_f13', 'dti_has_nii_f09', 'dti_has_nii_f13', 'has_braces_mri_f09', 'has_braces_mri_f13', 'exclude_incidental_f09', 'exclude_incidental_f13', 'freesurfer_qc_f09', 'freesurfer_qc_f13', 'software_ver_short_f09', 'dwi_nvols_f09', 'dwi_nvols_f13', 'dti_overall_qc_f09', 'dti_overall_qc_f13', 't1_asset_has_nii_f09'))

allbrain <- merge(df_brain1, df_brain2,by="idc", all=TRUE)
allbrain <- merge(allbrain, df_brain3,by="idc", all=TRUE)
allbrain <- merge(allbrain, df_brain4,by="idc", all=TRUE)
allbrain <- merge(allbrain, df_brain5,by="idc", all=TRUE)
allbrain <- merge(allbrain, df_brain6,by="idc", all=TRUE)
allbrain <- merge(allbrain, df_brain7,by="idc", all=TRUE)

colnames(allbrain)[1] <- "IDC"
df <- merge(df, allbrain,by="IDC", all=TRUE)

#save the merged raw dataset 
save(df, file="path_to_store_results")

# END OF THIS SCRIPT, GO TO PART B2.
