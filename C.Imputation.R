#########################################################################
######################PROJECT: Neonatal inflammation####################
#########################################################################

#Project: Mapping prenatal predictors and neurobehavioral outcomes of an epigenetic marker of neonatal inflammation - a longitudinal population-based study 

#SET UP WORKSPACE 
rm(list = ls()) #clears the environment

libraries <- c('mice', 'miceadds') 

invisible(lapply(libraries, require, character.only = T))

#set directory 
setwd("path_to_data")

#load final data 
load("/home/n.creasey/episcore_CRP_project/shadow2/mergeddata/CRP_mps_finaldata_23022024.Rdata")

##### Organize variable names for into variables to make it easier later #####
## for PRENATAL STRESS SCORE (PSS)
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

## OUTCOMES 

outcomes18m <- c("sum_int_18m", "sum_ext_18m","cbcl_sum_18m")

outcomes5y <- c("sum_int_5", "sum_ext_5", "cbcl_sum_5")

outcomes9y <- c("sum_int_9m", "sum_ext_9m", "cbcl_sum_9m")

outcomes14y <- c("sum_int_14", "sum_ext_14", "cbcl_sum_14")

all_cbcl <- c("sum_int_18m", "sum_ext_18m","cbcl_sum_18m", "sum_int_5", "sum_ext_5", "cbcl_sum_5", "sum_int_9m", "sum_ext_9m", "cbcl_sum_9m", "sum_int_14", "sum_ext_14", "cbcl_sum_14")

## CRP  variables 
CRP_vars <- c("HsCRPmgL_g1", "HsCRPmgL_g2", "CRP_birth", "CRPCHILD5")

## AUXILLARY VARIABLE LIST 
aux <- c("AGE_M_v2","SMOKE_ALL","GESTBIR","WEIGHT","GENDER", "m_depression","p_depression","marital_status", "PARITY", "INCOME", "BMI_0","EDUCM", "m_education") #nb. m_education is postnatal, EDUCM is pregnancy (latter is also included in PSS so use m_education in imputation for PSS 

aux_withmissings <- c("SMOKE_ALL","WEIGHT", "m_depression","p_depression","marital_status", "PARITY", "INCOME", "BMI_0","EDUCM","m_education")

### Explore missingness
sapply(df,function(y)sum(is.na(y))) 

#missingness for each of the prenatal scores (and save IDs as variabels ready for removing later)

df$PSS_missing <- rowMeans(is.na(df[,c(PSS_CR, PSS_LE, PSS_IR,PSS_PR)]))
hist(df$PSS_missing)
length(which(df$PSS_missing > 0.5)) #234
PSS_miss_exc <- which(df$PSS_missing > 0.5) 

df$PIS_missing <- rowMeans(is.na(df[,c(PIS_tri1, PIS_tri2, PIS_tri3)]))
hist(df$PIS_missing)
length(which(df$PIS_missing > 0.5)) #18
PIS_miss_exc <- which(df$PIS_missing > 0.5)


df$LPFS_missing <- rowMeans(is.na(df[,c(LPFS_meds, LPFS_lifestyle)]))
hist(df$LPFS_missing)
length(which(df$LPFS_missing > 0.5)) #279
LPFS_miss_exc <- which(df$LPFS_missing > 0.5)

df$PICS_missing <- rowMeans(is.na(df[,c(PICS_all)]))
hist(df$PICS_missing)
length(which(df$PICS_missing > 0.5)) #4
PICS_miss_exc <- which(df$PICS_missing > 0.5)

df$MICS_missing <- rowMeans(is.na(df[,c(MICS_all)]))
hist(df$MICS_missing)
length(which(df$MICS_missing > 0.5)) #238 
MICS_miss_exc <- which(df$MICS_missing > 0.5)

#check if anyone has above 50% on all variables 

df$high_miss <- as.factor(ifelse(df$PSS_missing > 0.5 & df$PIS_missing > 0.5 & df$LPFS_missing > 0.5 & df$PICS_missing > 0.5 & df$MICS_missing > 0.5, 1, 0)) #0


#### SET UP THE IMPUTATION ####
library(mice) 

## set up run 
imp0 <- mice(df, maxit = 0, defaultMethod = c('pmm', 'pmm', 'pmm', 'pmm'))

#after the test run some changes were made 
#sle, herpeszoster_tri1, and jaundice tri1,2,3 were removed from the scoring because they are all zeros


#function to derive mean domain scores on imputed items 
passive_imp_formula <- function(domain, add = "") {
  conc <- paste(domain, collapse = " + ")
  str <- paste0("~I( (", conc, ") / ", length(domain), ")")
  if (add != "") {
    str <- paste0("~I( (", conc, " + ", add, ") / ", length(domain)+1, ")")
  }
  return(str)
}

##adjust method
meth <-imp0$method


##passive imputation of domain scores
#PSS score
meth["PSS_LE_score"] <- passive_imp_formula(PSS_LE)

meth["PSS_CR_score"] <- passive_imp_formula(PSS_CR)

meth["PSS_PR_score"] <- passive_imp_formula(PSS_PR, add="m_age")

meth["PSS_IR_score"] <- passive_imp_formula(PSS_IR)

meth['PSS_TOT_score']  <- "~I(PSS_LE_score + PSS_CR_score + PSS_PR_score + PSS_IR_score)"

#PIS score
meth["URI_T1"] <- "~I(ifelse(df$upper_resp_inf_tri1 == 1 | df$rhinitis_tri1 == 1 | df$pharyngitis_tri1 == 1 | df$earinf_tri1 == 1 | df$sinusitis_tri1 == 1, 1, (ifelse(df$upper_resp_inf_tri1 == 0 & df$rhinitis_tri1 == 0 & df$pharyngitis_tri1 == 0 & df$earinf_tri1 == 0 & df$sinusitis_tri1 == 0, 0, NA))))"

meth["URI_T2"] <- "~I(ifelse(df$upper_resp_inf_tri2 == 1 | df$rhinitis_tri2 == 1 | df$pharyngitis_tri2 == 1 | df$earinf_tri2 == 1 | df$sinusitis_tri2 == 1, 1, (ifelse(df$upper_resp_inf_tri2 == 0 & df$rhinitis_tri2 == 0 & df$pharyngitis_tri2 == 0 & df$earinf_tri2 == 0 & df$sinusitis_tri2 == 0, 0, NA))))"


meth["URI_T3"] <- "~I(df$URI_T3 <- ifelse(df$upper_resp_inf_tri3 == 1 | df$rhinitis_tri3 == 1 | df$pharyngitis_tri3 == 1 | df$earinf_tri3 == 1 | df$sinusitis_tri3 == 1, 1, (ifelse(df$upper_resp_inf_tri3 == 0 & df$rhinitis_tri3 == 0 & df$pharyngitis_tri3 == 0 & df$earinf_tri3 == 0 & df$sinusitis_tri3 == 0, 0, NA))))"

meth["UTI_T1"] <- "~I(ifelse(df$UWI_tri1 == 1 | df$cystitis_tri1 == 1, 1, (ifelse(df$UWI_tri1 == 0 & df$cystitis_tri1 == 0, 0, NA))))"
meth["UTI_T2"] <- "~I(ifelse(df$UWI_tri2 == 1 | df$cystitis_tri2 == 1, 1, (ifelse(df$UWI_tri2 == 0 & df$cystitis_tri2 == 0, 0, NA))))"
meth["UTI_T3"] <- "~I(ifelse(df$UWI_tri3 == 1 | df$cystitis_tri3 == 1, 1, (ifelse(df$UWI_tri3 == 0 & df$cystitis_tri3 == 0, 0, NA))))"
meth["UTI_T1"] <- "~I(ifelse(df$UWI_tri1 == 1 | df$cystitis_tri1 == 1, 1, (ifelse(df$UWI_tri1 == 0 & df$cystitis_tri1 == 0, 0, NA))))"
meth["UTI_T2"] <- "~I(ifelse(df$UWI_tri2 == 1 | df$cystitis_tri2 == 1, 1, (ifelse(df$UWI_tri2 == 0 & df$cystitis_tri2 == 0, 0, NA))))"
meth["UTI_T3"] <- "~I(ifelse(df$UWI_tri3 == 1 | df$cystitis_tri3 == 1, 1, (ifelse(df$UWI_tri3 == 0 & df$cystitis_tri3 == 0, 0, NA))))"

meth["GII_T1"] <- "~I(ifelse(df$GI_inf_tri1 ==1 | df$enteritis_tri1 ==1, 1, (ifelse(df$GI_inf_tri1 ==0 & df$enteritis_tri1 ==0,0, NA))))"
meth["GII_T2"] <- "~I(ifelse(df$GI_inf_tri2 ==1 | df$enteritis_tri2 ==1, 1, (ifelse(df$GI_inf_tri2 ==0 & df$enteritis_tri2 ==0,0, NA))))"
meth["GII_T3"] <- "~I(ifelse(df$GI_inf_tri3 ==1 | df$enteritis_tri3 ==1, 1, (ifelse(df$GI_inf_tri3 ==0 & df$enteritis_tri3 ==0,0, NA))))"

meth["PIS_T1_score"] <- passive_imp_formula(PIS_T1)
meth["PIS_T2_score"] <- passive_imp_formula(PIS_T2)
meth["PIS_T3_score"] <- passive_imp_formula(PIS_T2)

meth["PIS_TOT_score"] <- "~I(PIS_T1_score + PIS_T2_score + PIS_T3_score)"

#LPFS score
#first make the medicine domains 
meth["psymeds"] <- "~I(ifelse(df$SSRI == 1 | df$TRIP == 1 | df$PSY == 1 | df$TCA == 1, 1, (ifelse(df$SSRI == 0 & df$TRIP == 0 & df$PSY == 0 & df$TCA == 0, 0, NA))))"
 
meth["inflams"] <- "~I(ifelse(df$NSAID == 1 | df$ABIO == 1 | df$PMOL == 1 | df$MUCO == 1 | df$COUGH == 1, 1, (ifelse(df$NSAID == 0 & df$ABIO == 0 & df$PMOL == 0 & df$MUCO == 0 & df$COUGH == 0, 0, NA))))"

meth["LPFS_score"] <- "~I(smoke_score + bmi_score + diet_score + psymeds + inflams + CORT)"

#PICS score
meth["PICS_score"] <- "~I(HELLP + gest_diab + preg_hypt + preclampsia + promm + sectio)"

#MICS score
meth["MICS_score"] <- "~I(intestinal + arthritis + ms + thyroid + diabetes)"

#education - try different method as not imputing well
meth["EDUCM" ] <- "polr"

### Adjust imputation model ###
# change predictor matrix
predictormatrix <- imp0$predictorMatrix

# set predictor matrix to 0
predictormatrix[,] <- 0

# set matrix to impute auxiliary and outcome variables (i.e., those that don't need passive imputation at the item level; exc. brain outcomes)

predictormatrix[c(aux),
              # impute using:
              c(# auxiliary variables
                aux[!aux == 'EDUCM'], PSS_domains, PIS_domains, LPFS_domains, MICS_score, PICS_score,
                # outcome variables
                outcomes5y[!outcomes5y == 'cbcl_sum_5m']
              )] <- 1
              
              
predictormatrix[c(outcomes18m),
              # impute using:
              c(# auxiliary variables
                aux[!aux == 'EDUCM'], PSS_domains, PIS_domains, LPFS_domains, MICS_score, PICS_score,
                # outcome variables
                outcomes5y[!outcomes5y == 'cbcl_sum_5m']
              )] <- 1
              
predictormatrix[c(outcomes5y),
              # impute using:
              c(# auxiliary variables
                aux[!aux == 'EDUCM'], PSS_domains, PIS_domains, LPFS_domains, MICS_score, PICS_score,
                # outcome variables
                outcomes9y[!outcomes9y == 'cbcl_sum_9m']
              )] <- 1
              
predictormatrix[c(outcomes9y),
              # impute using:
              c(# auxiliary variables
                aux[!aux == 'EDUCM'], PSS_domains, PIS_domains, LPFS_domains, MICS_score, PICS_score,
                # outcome variables
                outcomes5y[!outcomes5y == 'cbcl_sum_5m']
              )] <- 1
              
predictormatrix[c(outcomes14y),
              # impute using:
              c(# auxiliary variables
                aux[!aux == 'EDUCM'], PSS_domains, PIS_domains, LPFS_domains, MICS_score, PICS_score,
                # outcome variables
                outcomes5y[!outcomes5y == 'cbcl_sum_5m']
              )] <- 1
              
                 
## PSS
#set matrix to impute prenatal items per domain (based on other items in domain, other domain scores, auxiliary items, and outcomes)

PSS_aux <- c("AGE_M_v2","SMOKE_ALL","GESTBIR","WEIGHT","GENDER", "m_depression","p_depression","marital_status","m_education") # we include some factors usually used for postnatal score to act as auxilary variables for imputation
             
#pre_life events 
predictormatrix[c(PSS_LE),
                
                ## impute using:
                c(#other items on domain
                  PSS_LE,
                  
                  # auxiliary variables 
                  PSS_aux,
                  
                  # domain scores (but not for current scale computing)
                  PSS_domains[!PSS_domains == 'PSS_LE_score'],
                  # outcome variables
                  outcomes5y[!outcomes5y == 'cbcl_sum_5m'] 
                )] <- 1

#pre_contextual_risk
predictormatrix[c(PSS_CR),
                
                ## impute using:
                c(#other items on domain
                  PSS_CR,
                  
                  # auxiliary variables 
                  PSS_aux,
                  
                  # domain scores (but not for current scale computing)
                  PSS_domains[!PSS_domains == 'PSS_CR_score'],
                  # outcome variables
                  outcomes5y[!outcomes5y == 'cbcl_sum_5m']
                )] <- 1

#pre_parental risk
predictormatrix[c(PSS_PR),
                
                ## impute using:
                c(#other items on domain
                  PSS_PR,
                  
                  # auxiliary variables 
                  PSS_aux,
                  
                  # domain scores (but not for current scale computing)
                  PSS_domains[!PSS_domains == 'PSS_PR_score'],
                  # outcome variables
                  outcomes5y[!outcomes5y == 'cbcl_sum_5m']
                )] <- 1

#pre_interpersonal risk
predictormatrix[c(PSS_IR),
                
                ## impute using:
                c(#other items on domain
                  PSS_IR,
                  
                  # auxiliary variables 
                  PSS_aux,
                  
                  # domain scores (but not for current scale computing)
                  PSS_domains[!PSS_domains == 'PSS_IR_score'], 
                  # outcome variables
                  outcomes5y[!outcomes5y == 'cbcl_sum_5m']
                )] <- 1
 
## PIS
## set matrix to predict PIS items(PER TRIMESTER DOMAIN) using other items in trimetser domain and the other trmester domains scores, plus aux variables and outcomes   

PIS_aux <- c("AGE_M_v2","GESTBIR","GENDER","marital_status", "PARITY") #smoke and BMI removed as included in lifestyle scores 

#PIS tri 1             
predictormatrix[c(PIS_tri1),
                
                ## impute using:
                c(#other items on domain
                  PIS_tri1,
                  
                  # auxiliary variables 
                  PIS_aux,
                
                  # domain scores (but not for current scale computing)
                  PIS_domains[!PIS_domains == 'PIS_T1_score'], 
                  
                  # outcome variables
                  outcomes5y[!outcomes5y == 'cbcl_sum_5m']
                )] <- 1

#PIS tri 2
predictormatrix[c(PIS_tri2),
                
                ## impute using:
                c(#other items on domain
                  PIS_tri2,
                  
                  # auxiliary variables 
                  PIS_aux,  
                  
                  # domain scores (but not for current scale computing)
                  PIS_domains[!PIS_domains == 'PIS_T2_score'], 
                 
                  # outcome variables
                  outcomes5y[!outcomes5y == 'cbcl_sum_5m']
                )] <- 1

#PIS tri 3
 predictormatrix[c(PIS_tri3),
                
                ## impute using:
                c(#other items on domain
                  PIS_tri3,
                  
                  # auxiliary variables 
                  PIS_aux,
                  
                  # domain scores (but not for current scale computing)
                  PIS_domains[!PIS_domains == 'PIS_T3_score'], 
                  
                  # outcome variables
                  outcomes5y[!outcomes5y == 'cbcl_sum_5m']
                )] <- 1    

## LPFS

LPFS_aux <- c("AGE_M_v2","GESTBIR","GENDER","marital_status", "PARITY", "m_depression","p_depression")#smoke and BMI removed as included in lifestyle scores 

#medicines
predictormatrix[c(LPFS_meds),
                
                ## impute using:
                c(#other items on domain
                  LPFS_meds,
                  
                  #other LPFS domain
                  LPFS_lifestyle,

                  # auxiliary variables 
                  LPFS_aux,
                  
                  # outcome variables
                  outcomes5y[!outcomes5y == 'cbcl_sum_5m']
                )] <- 1  


#lifstyle scores 
predictormatrix[c(LPFS_lifestyle),
                
                ## impute using:
                c(#other items on domain
                  LPFS_lifestyle,
                  
                  #other LPFS domain
                  LPFS_meds,
                  
                  # auxiliary variables 
                  LPFS_aux,

                  # outcome variables
                  outcomes5y[!outcomes5y == 'cbcl_sum_5m']
                )] <- 1  


## PICS

PICS_aux <- c("AGE_M_v2","GESTBIR","GENDER","marital_status", "PARITY")

predictormatrix[c(PICS_all),
                
                ## impute using:
                c(#other items on domain
                  PICS_all,
                  
                  # auxiliary variables 
                  PICS_aux,
                  
                  # outcome variables
                  outcomes5y[!outcomes5y == 'cbcl_sum_5m']
                )] <- 1  

## MICS

MICS_aux <- c("AGE_M_v2","GESTBIR","GENDER","marital_status", "PARITY")

predictormatrix[c(MICS_all),
                
                ## impute using:
                c(#other items on domain
                  MICS_all,
                  
                  # auxiliary variables 
                  MICS_aux,
                  
                  # outcome variables
                  outcomes5y[!outcomes5y == 'cbcl_sum_5m']
                )] <- 1  



##### Final predictor matrix #####
#set diagonal to zero
diag(predictormatrix) <-0

## Adapting the predictor matrix based on the logged events of test run

#predictormatrix[c('p_criminal_record'),c("m_violence_property")] <- 0

predictormatrix[c(all_cbcl, "cbcl_sum_18m", "cbcl_sum_5", "cbcl_sum_9m", "m_depression", "m_education", "marital_status"),
                
                ## do not impute using:
                c("smoke_score","INCOME")] <- 0  

predictormatrix[c("EDUCM", "INCOME", "PARITY", "WEIGHT"),
                
                ## do not impute using:
                c("smoke_score")] <- 0  

predictormatrix[c("p_depression"),
                
                ## do not impute using:
                c("smoke_score", "INCOME", "bmi_score")] <- 0  


predictormatrix[c("preclampsia"),
                
                ## do not impute using:
                c("preg_hypt")] <- 0

predictormatrix[c("preclampsia"),
                
                ## do not impute using:
                c("preg_hypt")] <- 0

## adapting the predictor matrix based on the logged events of second test run 

predictormatrix[c("BMI_0"),
                
                ## do not impute using:
                c("smoke_score")] <- 0

predictormatrix[c("diet_score"),
                
                ## do not impute using:
                c("PSY")] <- 0


#check the matrix looks legit
pheatmap::pheatmap(predictormatrix, cluster_rows = F, cluster_cols = F)

##### Set which items to impute with a WHERE matrix ####

where <- make.where(df, keyword = c("missing")) #first we set it to impute any missing items 

#create a list of items that should not be imputed 
noimp_list <- c("C1","C2","C3","C4",
"gestage_plasma_g1",
"HsCRPmgL_g1",
"gestage_plasma_g2",
"HsCRPmgL_g2",
"CRP_birth",
"ageChildMF5",
"ageChildYF5",
"CRPCHILD5",
"age_GR1065",
"age_GR1066",
"agechild_GR1075",
"AgeChild_CBCL9m",
"AGECHILD_GR1093",
"Brain_Stem_vol_f09",
"TotalGrayVol_f09",
"genr_tbv_f09",
"eTIV_f09",
"Brain_Stem_vol_f13",
"TotalGrayVol_f13",
"genr_tbv_f13",
"eTIV_f13",
"mean_FA_genr_f09",
"mean_MD_genr_f09",
"mean_FA_genr_f13",
"mean_MD_genr_f13",
"mri_consent_f13",
"age_child_mri_f09",
"age_child_mri_f13",
"t1_has_nii_f09",
"dti_has_nii_f09",
"has_braces_mri_f09",
"has_braces_mri_f13",
"exclude_incidental_f09",
"exclude_incidental_f13",
"freesurfer_qc_f09",
"freesurfer_qc_f13",
"software_ver_short_f09",
"Cerebellum_Cortex_vol_subcortical_f09",
"Cerebellum_Cortex_vol_subcortical_f13",
"Hippocampus_vol_subcortical_f09",
"Hippocampus_vol_subcortical_f13",
"Amygdala_vol_subcortical_f09",
"Amygdala_vol_subcortical_f13",
"Lateral_Ventricle_vol_subcortical_f09",
"Lateral_Ventricle_vol_subcortical_f13",
"CerebralWhiteMatterVol_cortical_f09",
"CerebralWhiteMatterVol_cortical_f13",
"jaundice_tri1",
"jaundice_tri2",
"jaundice_tri3",
"herpeszoster_tri1",
"sle")


#set the where matrix to false for these variables
where[,noimp_list] <- as.logical("FALSE")

#duplicate the where matrix so we can edit to create 'where2' and test methods to exclude ppt with too little data
where2 <- where

#set the where matrix to FALSE for ppts with too little data (<0.5 items scores for each predictor)

where2[PSS_miss_exc,c(PSS_LE, PSS_CR, PSS_PR, PSS_IR, PSS_domains, "PSS_TOT_score")] <- as.logical("FALSE")
where2[PIS_miss_exc,c(PIS_tri1, PIS_tri2, PIS_tri3, PIS_domains, "PIS_TOT_score")] <- as.logical("FALSE")
where2[PICS_miss_exc,c(PICS_all,"PICS_score")] <- as.logical("FALSE")
where2[MICS_miss_exc,c(MICS_all,"MICS_score")] <- as.logical("FALSE")
#nb. LPFS scale has no one with >0.5 missing 

#visit the sequence 
VisSeq <- imp0$visitSequence


#### Testruns ####

#with where matrix that removes variables do not want imputed only
testrun_implist <- mice(df, m = 5, maxit = 10, where = where, seed = 08121996, method = meth, visitSequence = VisSeq, predictorMatrix = predictormatrix)
summary(complete(testrun_implist, action=1))

#without where matrix at all (everything and everyone imputed)
testrun_implist2 <- mice(df, m = 5, maxit = 10, seed = 08121996, method = meth, visitSequence = VisSeq, predictorMatrix = predictormatrix)
summary(complete(testrun_implist2, action=1))

#with a where matrix that excludes both variables and ppts that should not be imputed for specific variables)
testrun_implist3 <- mice(df, m = 5, maxit = 10, where = where2, seed = 08121996, method = meth, visitSequence = VisSeq, predictorMatrix = predictormatrix)
summary(complete(testrun_implist3, action=1))

#the imputation with the first where matrix (variables excluded but NOT ppts seems to work as well as with no where matrix, but where2 does not impute as many ppts 
#we go forward with where = where version


##### Run full imputation #####
implist <- mice(df, m = 30, maxit = 60, where = where, seed = 08121996, method = meth, visitSequence = VisSeq, predictorMatrix = predictormatrix)
#save the imputed dataset
saveRDS(implist, "/home/n.creasey/episcore_CRP_project/shadow2/output/implist_23022024.rds") 

#load the final imputation
implist <- readr::read_rds("/home/n.creasey/episcore_CRP_project/shadow2/output/implist_23022024.rds")
 
head(implist$loggedEvents) #no logged events 

#summary stats first imputed set of test run (change value to see stats for different iterations)
summary(complete(implist, action=1))

#some diagnostics on test run to check running correctly 
pdf("/home/n.creasey/episcore_CRP_project/shadow2/output/finalrun_densityplots_23022024.pdf")

densityplot(implist, ~PSS_TOT_score)
densityplot(implist, ~PIS_TOT_score)
densityplot(implist, ~LPFS_score)
densityplot(implist, ~PICS_score)
densityplot(implist, ~MICS_score)
densityplot(implist, ~sum_int_18m)
densityplot(implist, ~sum_ext_18m)
densityplot(implist, ~cbcl_sum_18m)
densityplot(implist, ~sum_int_5)
densityplot(implist, ~sum_ext_5)
densityplot(implist, ~cbcl_sum_5)
densityplot(implist, ~sum_ext_9m)
densityplot(implist, ~sum_int_9m)
densityplot(implist, ~cbcl_sum_9m)
densityplot(implist, ~sum_ext_14)
densityplot(implist, ~sum_int_14)
densityplot(implist, ~cbcl_sum_14)
densityplot(implist, ~INCOME)
densityplot(implist, ~EDUCM)

dev.off()


##### Set which items should be set to NA because too little data or should not be imputed ####

library(dplyr)

#make long 
long.impdata <- complete(implist, 'long', include = TRUE)

#transform the education variable into 3 levels

long.impdata$EDUCM_3l <- unclass(long.impdata$EDUCM) 
long.impdata$EDUCM_3l <- factor(recode(long.impdata$EDUCM_3l, `1` = 1, `2` = 1, `3` = 2, '4'= 2, '5'= 3, '6' = 3))

#change those with too much missing data to NA

long.impdata$PSS_TOT_score <- if_else(long.impdata$PSS_missing > 0.5, NA, long.impdata$PSS_TOT_score)
long.impdata$PIS_TOT_score <- if_else(long.impdata$PIS_missing > 0.5, NA, long.impdata$PIS_TOT_score)
long.impdata$PICS_score <- if_else(long.impdata$PICS_missing > 0.5, NA, long.impdata$PICS_score)
long.impdata$MICS_score <- if_else(long.impdata$MICS_missing > 0.5, NA, long.impdata$PICS_score)
long.impdata$LPFS_score <- if_else(long.impdata$LPFS_missing > 0.5, NA, long.impdata$LPFS_score)

# convert to mids
implist_clean <- as.mids(long.impdata)

#compare the original and clean implist to check correct
summary(complete(implist_clean, action=1))
summary(complete(implist, action=1))

#save the clean implist
saveRDS(implist_clean, "path_to_store_results") 

# END OF THIS SCRIPT, GO TO PART D
