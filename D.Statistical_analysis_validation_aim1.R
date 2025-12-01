#########################################################################
######################PROJECT: Neonatal inflammation####################
#########################################################################

#Project: Mapping prenatal predictors and neurobehavioral outcomes of an epigenetic marker of neonatal inflammation - a longitudinal population-based study 

###Parts in this script###
#Of note, in the previous parts of the script we constructed the methylation score, created clinical inflammation scores, created a final dataframe with all relevant data, applied inclusion and exclusion criteria, imputed missingness. 

#Now in this script we have the following parts: 
#' Part 0: Loading relevant libraries and data 
#' Part 1: Flowchart, baseline table
#' Part 2: Statistical analysis aim 1 [descriptive figures, validation MPS-CRP, prenatal clinical scores]

#\

###-----------------------------------------------###
#-----------------------PART 0----------------------#
###----------------------------------------------###

# Clear environment
rm(list = ls())

# Set seed 
set.seed(2024)

# Load relevant libraries 
libraries <- c('foreign', 'haven', 'tidyverse', 'broom.mixed', 'writexl', 'readxl', 'corrplot', 'stringi', 'miceadds', 'mitools', 'lme4', 'lmerTest', 'rlang', 'lattice', 'ggplot2', 'ggeffects', 'splines', 'reshape2', 'gridExtra', 'grid', 'ggpubr', 'survey', 'survival', 'ggseg', 'fmsb', 'ggcorrplot', 'Hmisc', 'RColorBrewer', 'car', 'DescTools', 'lavaan', 'semTools', 'sjPlot', 'pROC', 'metafor', 'pbmcapply', 'variancePartition', 'viridis')

invisible(lapply(libraries, require, character.only = T)) 

# Set working directory
setwd('set_path_to_data')

# Load final imputed data
final_imputed_data_mids <- readRDS('FINAL_implist_23022024.rds')

#\

###-----------------------------------------------###
#-----------------------PART 1----------------------#
###----------------------------------------------###

## Create flowchart (we start with 2511 with available epigenetic good QC data from the original 9901 children)
flowchart_df <- complete(final_imputed_data_mids, 30)
mother_ID <- read.spss('IDC-IDM-MOTHER.sav', to.data.frame = T)

flowchart_df_combined <- merge(flowchart_df, mother_ID, by = c('IDC', 'IDM'), all.x = T)

# only select singletons (n = 2511 > 2394)
flowchart_df_combined$na_count <- apply(flowchart_df_combined, 1, function(x) sum(is.na(x)))
flowchart_df_combined2 <- flowchart_df_combined[order(flowchart_df_combined$na_count),]
flowchart_df_combined3 <- flowchart_df_combined2[!duplicated(flowchart_df_combined2$MOTHER, fromLast = T),] 

# select sample size aim 2 (2394 > smri = 1439 / dti = 1479)
sMRI_df <- subset(flowchart_df_combined3, complete.cases(genr_tbv_f09) | genr_tbv_f13)
dti_df <- subset(flowchart_df_combined3, complete.cases(mean_MD_genr_f09) | mean_MD_genr_f13)

# select sample size aim 3 (2394 > 2394)
cbcl_df <- subset(flowchart_df_combined3, complete.cases(cbcl_sum_18m) | complete.cases(cbcl_sum_36m) | complete.cases(cbcl_sum_5) | complete.cases(cbcl_sum_9m) | complete.cases(cbcl_sum_14))

## Create baseline table for aim 1 sample (on last imputed dataset)
# Recode vars
flowchart_df_combined3$INCOME <- as.factor(ifelse(flowchart_df_combined3$INCOME == 'less than 450' | flowchart_df_combined3$INCOME == '450-600 euro' | flowchart_df_combined3$INCOME == '600-700 euro' | flowchart_df_combined3$INCOME == '700-800 euro' | flowchart_df_combined3$INCOME == '800-900 euro' | flowchart_df_combined3$INCOME == '900-1200 euro' | flowchart_df_combined3$INCOME == '1200-1400 euro' | flowchart_df_combined3$INCOME == '1400-1600 euro' | flowchart_df_combined3$INCOME == '1600-1800 euro' | flowchart_df_combined3$INCOME == '1800-2000 euro' | flowchart_df_combined3$INCOME == '2000-2200 euro', '<2000', ifelse(is.na(flowchart_df_combined3$INCOME), NA, '>2000')))

# Loop over baseline vars for sum stats
baselinevars <- c('AGE_M_v2', 'BMI_0', 'EDUCM_3l', 'INCOME', 'SMOKE_ALL', 'GENDER', 'GESTBIR')

for(i in baselinevars){
  #x = i vars that are columns in the dataframe df
  x <- flowchart_df_combined3[, i] 
  #show column name as heading per output 
  message(i) 
  #function for continuous variables
  summary_continuous <- function(x){
    standev <- sd(x, na.rm = T)
    meanvar <- mean(x, na.rm = T)
    print(paste0(round(meanvar, 1), '(', round(standev, 1), ')'))
  }
  #function for categorical variables 
  summary_categorical <- function(x){
    tab1 <- prop.table(table(x, useNA = 'always'))
    tab2 <- table(x, useNA = "always")
    print(paste(round(tab1 * 100, 1), '%', names(tab1), collapse = ','))
    print(paste(tab2, names(tab2)))
  }
  #if else to apply correct function for vars type 
  if (class(x) == 'numeric') {
    output <- summary_continuous(x)
  } 
  else 
  {
    output <- summary_categorical(x) 
  }
}

## Manipulations to imputed dataset before continuing to analyses
imp.test_long <- complete(final_imputed_data_mids, include = T, action = "long")

#Exclude siblings from imputed dataset before continuing to analyses
imp.test_long2 <- imp.test_long[(imp.test_long$IDC %in% flowchart_df_combined3$IDC),]

# Recode income
imp.test_long2$INCOME <- as.factor(ifelse(imp.test_long2$INCOME == 'less than 450' | imp.test_long2$INCOME == '450-600 euro' | imp.test_long2$INCOME == '600-700 euro' | imp.test_long2$INCOME == '700-800 euro' | imp.test_long2$INCOME == '800-900 euro' | imp.test_long2$INCOME == '900-1200 euro' | imp.test_long2$INCOME == '1200-1400 euro' | imp.test_long2$INCOME == '1400-1600 euro' | imp.test_long2$INCOME == '1600-1800 euro' | imp.test_long2$INCOME == '1800-2000 euro' | imp.test_long2$INCOME == '2000-2200 euro', '<2000', ifelse(is.na(imp.test_long2$INCOME), NA, '>2000')))

# Winsorize CRP
serum_crp_vars <- c("HsCRPmgL_g1", "HsCRPmgL_g2", "CRP_birth", "CRPCHILD5")

for (i in serum_crp_vars) {
  imp.test_long2[[paste0(i, "_win")]] <- Winsorize(imp.test_long2[[i]], probs = c(0.05, 0.95), na.rm = TRUE)
}

# Convert back to mids object 
imp.data.mids <- as.mids(imp.test_long2)
saveRDS(imp.data.mids,'imp.data.mids.rds')

#\

###-----------------------------------------------###
#-----------------------PART 2----------------------#
###----------------------------------------------###

# Select dataframe for figures
figure_df <- complete(imp.data.mids, 30)

### Step 1: make descriptive figures 

## Data collection time point
# Prenatal 
prenatal_vars <- c('Maternal age', 'Gestational age at birth', 'Maternal tobacco use', 'Maternal education', 'Household income', 'Child sex', 'Parity', 'Prenatal stress score', 'Lifestyle pro-inflammatory factors', 'Pregnancy related inflammatory conditions', 'Inflammatory medical conditions', 'C-reactive protein serum levels', 'Methylation profile score of CRP (Ligthart/Wielscher)') 

# Create dummy variables of mean timepoint per var
timepoint1 <- c(rep(40, 11), 13, 40) 
timepoint2 <- c(rep(NA, 11), 20, NA)
timepoint3 <- c(rep(NA, 11), 40, NA) 

# Create id vector of length of prenatal vars
id <- seq(1:length(prenatal_vars)) 

# Combine in one df 
prenatal_df <- data.frame(id, prenatal_vars, timepoint1, timepoint2, timepoint3) 

# Wide to long format 
prenatal_df_long <- reshape(prenatal_df, idvar = 'id', varying = c("timepoint1", "timepoint2", "timepoint3"), v.names = c('time'), direction ='long') 

# Order by ID 
prenatal_df_long <- prenatal_df_long[order(prenatal_df_long$id),] 

# Create factor of measurement number 
prenatal_df_long[, "Measurement"] <- rep(1:3, length.out = nrow(prenatal_df_long)) %>%
  factor() 

# Create prenatal plot 
prenatal_timepoint_plot <- ggplot(prenatal_df_long, aes(prenatal_vars, time)) + 
  geom_point(aes(group = id,
                 color = Measurement,
                 shape = Measurement,
                 size = 2),
             na.rm = TRUE) + 
  coord_flip() + 
  labs(y= 'Time (weeks gestation)', 
       x = "", 
       title = 'Prenatal period') + 
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        plot.title = element_text(size=13, face = 'bold'),
        legend.title = element_text(size = 11),
        panel.background = element_rect(fill = 'white')) +
  scale_size(guide = 'none')

# Postnatal
postnatal_vars <- c('CBCL total behavioral symptoms score', 'CBCL internalizing symptoms score', 'CBCL externalizing symptoms score', 'Total brain volume', 'Total grey matter volume', 'Total white matter volume', 'Cerebellum volume', 'Hippocampus volume', 'Amygdala volume', 'Brain stem volume', 'Lateral ventricles volume', 'Global fractional anisotropy', 'Global mean diffusivity', 'C-reactive protein serum levels')

# Create dummy variables of mean timepoint per var
timepoint1 <- c(rep(1.5, 3), rep(10, 10), 5) 
timepoint2 <- c(rep(3, 3), rep(14, 10), NA)
timepoint3 <- c(rep(6, 3), rep(NA, 11))
timepoint4 <- c(rep(10, 3), rep(NA, 11))
timepoint5 <- c(rep(14, 3), rep(NA, 11))

# Create id vector of length of prenatal vars
id2 <- seq(1:length(postnatal_vars)) 

# Combine in one df 
postnatal_df <- data.frame(id2, postnatal_vars, timepoint1, timepoint2, timepoint3, timepoint4, timepoint5) 

# Wide to long format 
postnatal_df_long <- reshape(postnatal_df, idvar = 'id2', varying = c("timepoint1", "timepoint2", "timepoint3", "timepoint4", "timepoint5"), v.names = c('time'), direction ='long') 

# Order by ID 
postnatal_df_long <- postnatal_df_long[order(postnatal_df_long$id2),] 

# Create factor of measurement number 
postnatal_df_long[, "Measurement"] <- rep(1:5, length.out = nrow(postnatal_df_long)) %>%
  factor() 

# Create postnatal plot 
postnatal_timepoint_plot <- ggplot(postnatal_df_long, aes(postnatal_vars, time)) + 
  geom_point(aes(group = id2,
                 color = Measurement,
                 shape = Measurement,
                 size = 2),
             na.rm = TRUE) + 
  coord_flip() + 
  labs(y= 'Time (child age in years)', 
       x = "", 
       title = 'Postnatal period') + 
  theme_bw() +
  theme(axis.title = element_text(size = 12,),
        plot.title = element_text(size=13, face = 'bold'),
        legend.title = element_text(size = 11),
        panel.background = element_rect(fill = 'white')) +
  scale_size(guide = 'none')

# Combine prenatal and postnatal plots
ggarrange(prenatal_timepoint_plot, postnatal_timepoint_plot, labels = 'AUTO', nrow = 2, align = 'v')

#\

## Correlation plot
# First residualize batch effects out of exposure variable for correlation plot, because it adds to so much noise 
resid_lig <- lm(mpsz_lig ~ Sample_Plate, data = figure_df)
resid_wiel <- lm(mpsz_wiel ~ Sample_Plate, data = figure_df)
figure_df$resid_mpsz_lig <- residuals(resid_lig)
figure_df$resid_mpsz_wiel <- residuals(resid_wiel)

# Subset cols for the correlation plot
correlation_df <- dplyr::select(figure_df, c('resid_mpsz_lig', 'resid_mpsz_wiel', 'Bcell', 'CD4T', 'CD8T', 'Gran', 'Mono', 'NK', 'nRBC', 'HsCRPmgL_g1', 'HsCRPmgL_g1', 'HsCRPmgL_g2' ,'CRP_birth', 'CRPCHILD5', 'AGE_M_v2', 'GESTBIR', 'PARITY', 'GENDER', 'INCOME', 'BMI_0', 'SMOKE_ALL', 'EDUCM_3l', 'cbcl_sum_14', 'genr_tbv_f13', 'mean_MD_genr_f13', 'mean_FA_genr_f13', 'PSS_TOT_score', 'PIS_TOT_score', 'MICS_score', 'PICS_score', 'LPFS_score'))

# Converting factor to numeric only for correlation plot purpose
factors <- c("GENDER", "INCOME", "SMOKE_ALL", "EDUCM_3l")
correlation_df[factors] <- lapply(correlation_df[factors], as.numeric)

# Rename variables
correlation_df2 <- rename(correlation_df, '1. MPS by Ligthart' = resid_mpsz_lig, '2. MPS by Wielscher' = resid_mpsz_wiel, '15. B-cells' = Bcell, '16. CD4 T-cells' = CD4T, '17. CD8 T-cells' = CD8T, '18. Granulocytes' = Gran, '19. Monocytes' = Mono, '20. NK-cells' = NK, '21. nucleated red blood cells' = nRBC, '3. serum CRP (trimester 1)' = HsCRPmgL_g1, '4. serum CRP (trimester 2)' = HsCRPmgL_g2, '5. serum CRP (cord blood)' = CRP_birth, '6. serum CRP (child age 5)' = CRPCHILD5, '7. Maternal age' = AGE_M_v2, '8. Pre-pregnancy BMI' = BMI_0, '9. Maternal tobacco use' = SMOKE_ALL, '10. Maternal education' = EDUCM_3l, '22. CBCL total score (child age 14)' = cbcl_sum_14, '23. Total brain volume (child age 14)' = genr_tbv_f13, '24. Global mean diffusivity (child age 14)' = mean_MD_genr_f13, '25. Global fractional anisotropy (child age 14)' = mean_FA_genr_f13, '26. Prenatal stress score' = PSS_TOT_score, '27. Prenatal infection sum score' = PIS_TOT_score, '28. Inflammatory medical conditions score' = MICS_score, '29. Pregnancy related inflammatory factors score' = PICS_score, '30. Lifestyle pro-inflammatory factors score' = LPFS_score, '11. Gestational age at birth' = GESTBIR, '12. Parity' = PARITY, '13. Child sex' = GENDER, '14. Household income' = INCOME)

# Specify the desired order of variables
corr_order <- c("1. MPS by Ligthart", "2. MPS by Wielscher", "3. serum CRP (trimester 1)", "4. serum CRP (trimester 2)", "5. serum CRP (cord blood)", "6. serum CRP (child age 5)", "7. Maternal age", "8. Pre-pregnancy BMI", "9. Maternal tobacco use", "10. Maternal education", "11. Gestational age at birth", "12. Parity", "13. Child sex", "14. Household income", "15. B-cells", "16. CD4 T-cells", "17. CD8 T-cells", "18. Granulocytes", "19. Monocytes", "20. NK-cells", "21. nucleated red blood cells", "22. CBCL total score (child age 14)", "23. Total brain volume (child age 14)", "24. Global mean diffusivity (child age 14)", "25. Global fractional anisotropy (child age 14)", "26. Prenatal stress score", "27. Prenatal infection sum score", "28. Inflammatory medical conditions score", "29. Pregnancy related inflammatory factors score", "30. Lifestyle pro-inflammatory factors score")

# Reorder the correlation matrix based on the desired order
correlation <- cor(correlation_df2, use="pairwise.complete.obs")
correlation_ordered <- correlation[corr_order, corr_order]

# Plot the correlation matrix with the specified order
corrplot(correlation_ordered, method = 'color', addCoef.col = "black", number.cex = 0.4, type = 'lower', diag = FALSE, tl.col = 'black', tl.cex = 0.5, col = colorRampPalette(c("midnightblue", "white", "darkred"))(100), colnames = colnames(correlation_ordered))

#\

## Histogram and box plot MPS ligthart/wielscher
figure_450k_df <- subset(figure_df, array == '450k')
figure_epic_df <- subset(figure_df, array == 'epic')

hist_ligthart_450k <- ggplot(figure_450k_df, aes(x = mpsz_lig)) + geom_histogram(bins = 50, color = 'white', fill = '#58508d') + theme_bw() + xlab('Ligthart (450k assay)') + ylab('Count')
boxplot_ligthart_450k <- ggplot(figure_450k_df, aes(x = "", y = mpsz_lig)) + geom_boxplot(color = 'black', fill = '#58508d') + theme_bw() + ylab('Ligthart (450k assay)') + xlab("") 
hist_ligthart_epic <- ggplot(figure_epic_df, aes(x = mpsz_lig)) + geom_histogram(bins = 50, color = 'white', fill = '#58508d') + theme_bw() + xlab('Ligthart (epic assay)') + ylab('Count')
boxplot_ligthart_epic <- ggplot(figure_epic_df, aes(x = "", y = mpsz_lig)) + geom_boxplot(color = 'black', fill = '#58508d') + theme_bw() + ylab('Ligthart (epic assay)') + xlab("") 
hist_wielscher_450k <- ggplot(figure_450k_df, aes(x = mpsz_wiel)) + geom_histogram(bins = 50, color = 'white', fill = '#ffa600') + theme_bw() + xlab('Wielscher (450k assay)') + ylab('Count')
boxplot_wielscher_450k <- ggplot(figure_450k_df, aes(x = "", y = mpsz_wiel)) + geom_boxplot(color = 'black', fill = '#ffa600') + theme_bw() + ylab('Wielscher (450k assay)') + xlab("") 
hist_wielscher_epic <- ggplot(figure_epic_df, aes(x = mpsz_wiel)) + geom_histogram(bins = 50, color = 'white', fill = '#ffa600') + theme_bw() + xlab('Wielscher (epic assay)') + ylab('Count')
boxplot_wielscher_epic <- ggplot(figure_epic_df, aes(x = "", y = mpsz_wiel)) + geom_boxplot(color = 'black', fill = '#ffa600') + theme_bw() + ylab('Wielscher (epic assay)') + xlab("") 

ggarrange(hist_ligthart_450k, boxplot_ligthart_450k,hist_ligthart_epic,boxplot_ligthart_epic,hist_wielscher_450k,boxplot_wielscher_450k,hist_wielscher_epic,boxplot_wielscher_epic, labels = 'AUTO', ncol = 2, nrow = 4) 

#\

## ROC/AUC MPSes
# Residualize scores for batch effects 
resid_ligthart_450k <- lm(mpsz_lig ~ Sample_Plate, data = figure_450k_df) 
resid_wielscher_450k <- lm(mpsz_wiel ~ Sample_Plate, data = figure_450k_df) 

figure_450k_df$resid_mpsz_lig <- residuals(resid_ligthart_450k)
figure_450k_df$resid_mpsz_wiel <- residuals(resid_wielscher_450k)

resid_ligthart_epic <- lm(mpsz_lig ~ Sample_Plate, data = figure_epic_df) 
resid_wielscher_epic <- lm(mpsz_wiel ~ Sample_Plate, data = figure_epic_df) 

figure_epic_df$resid_mpsz_lig <- residuals(resid_ligthart_epic)
figure_epic_df$resid_mpsz_wiel <- residuals(resid_wielscher_epic)

# Create cutoffs for crp based on clinical recommendations 
figure_450k_df$CRP_birth_categorical <- as.factor(ifelse(figure_450k_df$CRP_birth > 1, 1, 0)) 
figure_epic_df$CRP_birth_categorical <- as.factor(ifelse(figure_epic_df$CRP_birth > 1, 1, 0))

# Calculate AUC
roc1 <- roc(CRP_birth_categorical ~ resid_mpsz_lig, data = figure_450k_df) #auc = 0.76
roc2 <- roc(CRP_birth_categorical ~ resid_mpsz_lig, data = figure_epic_df) #auc = 0.55
roc3 <- roc(CRP_birth_categorical ~ resid_mpsz_wiel, data = figure_450k_df) #auc = 0.68
roc4 <- roc(CRP_birth_categorical ~ resid_mpsz_wiel, data = figure_epic_df) #auc = 0.53

# Create plots 
par <- par(mfrow=c(2,2))

par <- plot(roc1)
text(x = 0.3, y = 0.1, 'ROC curve for Ligthart 450k MPS score (AUC = 0.76)', font = 0.8)

par <- plot(roc2)
text(x = 0.25, y = 0.1, 'ROC curve for Wielscher 450k MPS score (AUC = 0.55)', font = 0.8)

par <- plot(roc3)
text(x = 0.25, y = 0.1, 'ROC curve for Ligthart epic MPS score (AUC = 0.68)', font = 0.8)

par <- plot(roc4)
text(x = 0.25, y = 0.1, 'ROC curve for Wielscher epic MPS score (AUC = 0.53)', font = 0.8)

#\

## Distribution plots prenatal scores
hist_pis <- ggplot(figure_df, aes(x = PIS_TOT_score)) + geom_histogram(fill = 'rosybrown1') + theme_bw() + labs(x = 'Prenatal infection sum score', y = 'Frequency')
hist_pss <- ggplot(figure_df, aes(x = PSS_TOT_score)) + geom_histogram(fill = 'palegreen2') + theme_bw() + labs(x = 'Prenatal stress score', y = 'Frequency')
hist_lpfs <- ggplot(figure_df, aes(x = LPFS_score)) + geom_histogram(fill = 'paleturquoise3') + theme_bw() + labs(x = 'Lifestyle pro-inflammatory factors score', y = 'Frequency')
hist_pics <- ggplot(figure_df, aes(x = PICS_score)) + geom_histogram(fill = 'palegoldenrod', stat = 'count') + theme_bw() + labs(x = 'Pregnancy related inflammatory conditions score', y = 'Frequency')
hist_mics <- ggplot(figure_df, aes(x = MICS_score)) + geom_histogram(fill = 'darksalmon', stat = 'count') + theme_bw() + labs(x = 'Inflammatory medical conditions score', y = 'Frequency')

ggarrange(hist_pis, hist_pss, hist_lpfs, hist_pics, hist_mics, labels = 'AUTO')

#\

## Radar chart prenatal scores
# Prenatal stress score 
pss_df <- dplyr::select(figure_df, c('PSS_LE_score', 'PSS_CR_score', 'PSS_PR_score', 'PSS_IR_score'))

col_max <- apply(pss_df, 2, max, na.rm = T)
col_min <- apply(pss_df, 2, min, na.rm = T)
col_mean <- apply(pss_df, 2, mean, na.rm = T)
col_summary <- t(data.frame(Max = col_max, Min = col_min, Average = col_mean))

df_scaled2 <- as.data.frame(rbind(col_summary, pss_df))
df_average <- df_scaled2[c("Max", "Min", "Average"),]
df_average2 <- rename(df_average, "Life events" = PSS_LE_score, "Contextual risk" = PSS_CR_score, "Parental risk" = PSS_PR_score,  "Interpersonal stress" = PSS_IR_score)

radarchart(df_average2, axistype = 1, pcol= rgb(0.9,0.3,0.4,1), pfcol= rgb(0.8,0.3,0.4,0.3),  plwd=3, cglcol="grey", cglty=1,  axislabcol="grey",  caxislabels=seq(0,20,5), cglwd=0.8, vlcex=1.0, title = 'Prenatal stress score')

# Prenatal infection score (redo!)
figure_df$Fever <- figure_df$fever_tri1 + figure_df$fever_tri2 + figure_df$fever_tri3
figure_df$'Upper respiratory infections' <- figure_df$URI_T1 + figure_df$URI_T2 + figure_df$URI_T3
figure_df$'Lower respiratory infections' <- figure_df$LRI_T1 + figure_df$LRI_T2 + figure_df$LRI_T3  
figure_df$'Urinary tract infections' <- figure_df$UTI_T1 + figure_df$UTI_T2 + figure_df$UTI_T3
figure_df$'Gastrointestinal infections' <- figure_df$GII_T1 + figure_df$GII_T2 + figure_df$GII_T3
figure_df$Flu <- figure_df$flu_tri1 + figure_df$flu_tri2 + figure_df$flu_tri3
figure_df$Dermatitis <- figure_df$dermatitis_tri1 + figure_df$dermatitis_tri2 + figure_df$dermatitis_tri3
figure_df$'Sexually transmitted diseases' <- figure_df$STD_tri1 + figure_df$STD_tri2 + figure_df$STD_tri3
figure_df$'Herpes zoster' <- figure_df$herpeszoster_tri1 + figure_df$herpeszoster_tri2 + figure_df$herpeszoster_tri3
figure_df$'Eye infections' <- figure_df$jaundice_tri1 + figure_df$jaundice_tri2 + figure_df$jaundice_tri3

piss_df <- dplyr::select(figure_df, c('Fever', 'Upper respiratory infections', 'Lower respiratory infections', 'Urinary tract infections', 'Gastrointestinal infections', 'Flu', 'Dermatitis', 'Sexually transmitted diseases', 'Herpes zoster', 'Eye infections'))

col_max <- apply(piss_df, 2, max, na.rm = T)
col_min <- apply(piss_df, 2, min, na.rm = T)
col_mean <- apply(piss_df, 2, mean, na.rm = T)
col_summary <- t(data.frame(Max = col_max, Min = col_min, Average = col_mean))

df_scaled2 <- as.data.frame(rbind(col_summary, piss_df))
df_average <- df_scaled2[c("Max", "Min", "Average"),]

radarchart(df_average, axistype = 1, pcol= rgb(0.3,0.3,0.4,1), pfcol= rgb(0.3,0.3,0.4,0.3),  plwd=3, cglcol="grey", cglty=1,  axislabcol="grey",  caxislabels=seq(0,20,5), cglwd=0.8, vlcex=1.0, title = 'Prenatal infection score')

# Lifestyle pro-inflammatory factors score
lpfs_df <- dplyr::select(figure_df, c("smoke_score", "bmi_score", "diet_score", "psymeds", "CORT", "inflams"))

col_max <- apply(lpfs_df, 2, max, na.rm = T)
col_min <- apply(lpfs_df, 2, min, na.rm = T)
col_mean <- apply(lpfs_df, 2, mean, na.rm = T)
col_summary <- t(data.frame(Max = col_max, Min = col_min, Average = col_mean))

df_scaled2 <- as.data.frame(rbind(col_summary, lpfs_df))
df_average <- df_scaled2[c("Max", "Min", "Average"),]
df_average2 <- rename(df_average, 'Tobacco use' = smoke_score, 'Pre-pregnancy obesity' = bmi_score, 'Corticosteroid use' = CORT, 'Poor maternal diet' = diet_score, 'Maternal psychotropic medication use' = psymeds, 'Maternal anti-inflammatory or infection medication use' = inflams)

radarchart(df_average2, axistype = 1, pcol= rgb(0.1,0.6,0.4,1), pfcol= rgb(0.1,0.6,0.4,0.3),  plwd=3, cglcol="grey", cglty=1,  axislabcol="grey",  caxislabels=seq(0,20,5), cglwd=0.8, vlcex=1.0, title = 'Lifestyle pro-inflammatory factors score')

# Pregnancy related inflammatory conditions score
pics_df <- dplyr::select(figure_df, c("HELLP", "gest_diab", "preg_hypt", "preclampsia", "promm", "sectio"))

col_max <- apply(pics_df, 2, max, na.rm = T)
col_min <- apply(pics_df, 2, min, na.rm = T)
col_mean <- apply(pics_df, 2, mean, na.rm = T)
col_summary <- t(data.frame(Max = col_max, Min = col_min, Average = col_mean))

df_scaled2 <- as.data.frame(rbind(col_summary, pics_df))
df_average <- df_scaled2[c("Max", "Min", "Average"),]
df_average2 <- rename(df_average, 'HELLP' = HELLP, 'Gestational diabetes' = gest_diab, 'Pregnancy induced hypertension' = preg_hypt, 'Preeclampsia' = preclampsia, 'Caesarian delivery' = sectio, 'Premature ruptured membranes' = promm)

radarchart(df_average2, axistype = 1, pcol= rgb(0.7,0.7,0.4,1), pfcol= rgb(0.7,0.7,0.4,0.3),  plwd=3, cglcol="grey", cglty=1,  axislabcol="grey",  caxislabels=seq(0,20,5), cglwd=0.8, vlcex=1.0, title = 'Pregnancy related inflammatory conditions score')

# Inflammatory medical conditions score
mics_df <- dplyr::select(figure_df, c("intestinal", "arthritis", "ms", "thyroid", "diabetes"))

col_max <- apply(mics_df, 2, max, na.rm = T)
col_min <- apply(mics_df, 2, min, na.rm = T)
col_mean <- apply(mics_df, 2, mean, na.rm = T)
col_summary <- t(data.frame(Max = col_max, Min = col_min, Average = col_mean))

df_scaled2 <- as.data.frame(rbind(col_summary, mics_df))
df_average <- df_scaled2[c("Max", "Min", "Average"),]
df_average2 <- rename(df_average, 'Intestinal disorders' = intestinal, 'Arthritis' = arthritis, 'Multiple sclerosis' = ms, 'Thyroid disorder' = thyroid, 'Diabetes' = diabetes)

radarchart(df_average2, axistype = 1, pcol= rgb(0.8,0.5,0.1,1), pfcol= rgb(0.8,0.5,0.1,0.3),  plwd=3, cglcol="grey", cglty=1,  axislabcol="grey",  caxislabels=seq(0,20,5), cglwd=0.8, vlcex=1.0, title = 'Inflammatory medical conditions score')

#\

## Histogram CRP trimester 1/2, cord blood and child age 5 (median/IQR)
# Untransformed 
crp1 <- ggplot(figure_df, aes(x = HsCRPmgL_g1)) + geom_histogram(fill = 'sienna2', bins = 50) + theme_bw() + labs(x = 'C-reactive protein (trimester 1)', y = 'Frequency') #+ geom_label(x = 200, y = 2050, label = '4.4 (2.3 - 7.8) [mg/L]', color = 'black', size = 4) 
crp2 <- ggplot(figure_df, aes(x = HsCRPmgL_g2)) + geom_histogram(fill = 'violetred3', bins = 50) + theme_bw() + labs(x = 'C-reactive protein (trimester 2)', y = 'Frequency') #+ geom_label(x = 150, y = 3050, label = '4.2 (2.4 - 7.0) [mg/L]', color = 'black', size = 4)
crp3 <- ggplot(figure_df, aes(x = CRP_birth)) + geom_histogram(fill = 'turquoise4', bins = 50) + theme_bw() + labs(x = 'C-reactive protein (cord blood)', y = 'Frequency') #+ geom_label(x = 37, y = 3700, label = '0.2 (0.2 - 0.2) [mg/L]', color = 'black', size = 4)
crp4 <- ggplot(figure_df, aes(x = CRPCHILD5)) + geom_histogram(fill = 'slateblue3', bins = 50) + theme_bw() + labs(x = 'C-reactive protein (child age 5 years)', y = 'Frequency') #+ geom_label(x = 73, y = 2700, label = '0.3 (0.1 - 0.9) [mg/L]', color = 'black', size = 4)

# Log transformed
crp5 <- ggplot(figure_df, aes(x = log(HsCRPmgL_g1))) + geom_histogram(fill = 'sienna', bins = 50) + theme_bw() + labs(x = 'Log transformed C-reactive protein (trimester 1)', y = 'Frequency')
crp6 <- ggplot(figure_df, aes(x = log(HsCRPmgL_g2))) + geom_histogram(fill = 'violetred4', bins = 50) + theme_bw() + labs(x = 'Log transformed C-reactive protein (trimester 2)', y = 'Frequency')
crp7 <- ggplot(figure_df, aes(x = log(CRP_birth))) + geom_histogram(fill = 'turquoise', bins = 50) + theme_bw() + labs(x = 'Log transformed C-reactive protein (cord blood)', y = 'Frequency')
crp8 <- ggplot(figure_df, aes(x = log(CRPCHILD5))) + geom_histogram(fill = 'slateblue4', bins = 50) + theme_bw() + labs(x = 'Log transformed C-reactive protein (child age 5 years)', y = 'Frequency')

ggarrange(crp1, crp5, crp2, crp6, crp3, crp7, crp4, crp8, labels = 'AUTO', nrow = 4, ncol = 2) 

#\

## Scatter plot MPSes vs observed serum CRP in cord blood
figure_450k_df <- subset(figure_df, array = '450k')
figure_epic_df <- subset(figure_df, array = 'epic')

# 450k
scatter_450k <- ggplot(figure_450k_df) +
  # Group 1
  geom_point(aes(x = resid_mpsz_lig, y = CRP_birth, color = 'Group 1'), 
             size = 1, alpha = 0.7, position = position_jitter(width = 0.2, height = 0.2)) +
  stat_smooth(aes(x = resid_mpsz_lig, y = CRP_birth), 
              method = "lm", formula = y ~ x, 
              geom = "smooth", color = "#58508d", size = 1, se = FALSE) +
  # Group 2
  geom_point(aes(x = resid_mpsz_wiel, y = CRP_birth, color = 'Group 2'), 
             size = 1, alpha = 0.7, position = position_jitter(width = 0.2, height = 0.2)) +
  stat_smooth(aes(x = resid_mpsz_wiel, y = CRP_birth), 
              method = "lm", formula = y ~ x, 
              geom = "smooth", color = "#ffa600", size = 1, se = FALSE) +
  # Overall styling
  theme_bw() +
  scale_color_manual(values = c("Group 1" = "#58508d", "Group 2" = "#ffa600"), 
                     labels = c("MPS by Ligthart - 450k", "MPS by Wielscher - 450k")) +
  labs(x = 'MPS score - 450k (scaled and residualized for batch effects)', 
       y = 'CRP serum levels at birth', 
       color = 'MPS score') +
  theme(legend.position = 'top', legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4))) 

print(scatter_450k)

# Epic
scatter_epic <- ggplot(figure_epic_df) +
  # Group 1
  geom_point(aes(x = resid_mpsz_lig, y = CRP_birth, color = 'Group 1'), 
             size = 1, alpha = 0.7, position = position_jitter(width = 0.2, height = 0.2)) +
  stat_smooth(aes(x = resid_mpsz_lig, y = CRP_birth), 
              method = "lm", formula = y ~ x, 
              geom = "smooth", color = "#58508d", size = 1, se = FALSE) +
  # Group 2
  geom_point(aes(x = resid_mpsz_wiel, y = CRP_birth, color = 'Group 2'), 
             size = 1, alpha = 0.7, position = position_jitter(width = 0.2, height = 0.2)) +
  stat_smooth(aes(x = resid_mpsz_wiel, y = CRP_birth), 
              method = "lm", formula = y ~ x, 
              geom = "smooth", color = "#ffa600", size = 1, se = FALSE) +
  # Overall styling
  theme_bw() +
  scale_color_manual(values = c("Group 1" = "#58508d", "Group 2" = "#ffa600"), 
                     labels = c("MPS by Ligthart - epic", "MPS by Wielscher - epic")) +
  labs(x = 'MPS score - epic (scaled and residualized for batch effects)', 
       y = 'CRP serum levels at birth', 
       color = 'MPS score') +
  theme(legend.position = 'top', legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4))) 

print(scatter_epic)

# Combine plots
ggarrange(scatter_450k, scatter_epic, labels = 'AUTO', nrow = 2, common.legend = T)

#\

## Plot prenatal scores vs MPS ligthart/wielscher 
# Ligthart 
ligthart <- ggplot(figure_df) +
  stat_smooth(aes(x = scale(PSS_TOT_score), 
                  y = resid_mpsz_lig, color = "Group 1"), 
               method = "gam", 
              formula = y ~ x, 
              geom = "smooth", 
              size = 1, 
              alpha = 0.1) +
  stat_smooth(aes(x = scale(PIS_TOT_score), 
                  y = resid_mpsz_lig, color = "Group 2"), 
               method = "gam", 
              formula = y ~ x, 
              geom = "smooth", 
              size = 1, 
              alpha = 0.1)+
  stat_smooth(aes(x = scale(PICS_score),
                  y = resid_mpsz_lig, 
                  color = "Group 3"), 
              method = "gam", 
              formula = y ~ x, 
              geom = "smooth", 
              size = 1, 
              alpha = 0.1) +
  stat_smooth(aes(x = scale(MICS_score), 
                  y = resid_mpsz_lig, 
                  color = "Group 4"), 
              method = "gam", 
              formula = y ~ x, 
              geom = "smooth", 
              size = 1, 
              alpha = 0.1)+
  stat_smooth(aes(x = scale(LPFS_score), 
                  y = resid_mpsz_lig, 
                  color = "Group 5"), 
              method = "gam", 
              formula = y ~ x, 
              geom = "smooth", 
              size = 1, 
              alpha = 0.1) +
  scale_color_manual(values = c("Group 1" = "coral3", 
                                "Group 2" = "chartreuse3", 
                                "Group 3" = "dodgerblue4", 
                                "Group 4" = "darkorchid3", 
                                "Group 5" = "goldenrod3"),
                     labels = c("Prenatal stress", "Prenatal infections", 
                                "Pregnancy related inflammatory
                                conditions", 
                                "Inflammatory medical conditions", 
                                "Lifestyle pro-inflammatory factors")) +
  facet_wrap(~ array, scales = "free_y", ncol = 2) +
  theme_bw() +
  labs(x = "Prenatal clinical inflammatory scores", 
       y = "Methylation profile score - Ligthart",  
       color = 'Prenatal predictor') + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        strip.background = element_rect(fill = 'moccasin'),
        strip.text = element_text(colour = 'black', 
                                  face = 'bold')) +
  geom_point(aes(x = scale(PSS_TOT_score), 
                 y = resid_mpsz_lig, 
                 color = "Group 1"), 
             size = 0.01, 
             alpha = 0.35,
             position = position_jitter(width = 0.2, 
                                        height = 0.2)) +
  geom_point(aes(x = scale(PIS_TOT_score),
                 y = resid_mpsz_lig,
                 color = "Group 2"),
             size = 0.01, 
             alpha = 0.35,
             position = position_jitter(width = 0.2, 
                                        height = 0.2)) +
  geom_point(aes(x = scale(PICS_score), 
                 y = resid_mpsz_lig, 
                 color = "Group 3"), 
             size = 0.01, 
             alpha = 0.35,
             position = position_jitter(width = 0.2, 
                                        height = 0.2)) +
  geom_point(aes(x = scale(MICS_score), 
                 y = resid_mpsz_lig, 
                 color = "Group 4"), 
             size = 0.01, 
             alpha = 0.35,
             position = position_jitter(width = 0.2, 
                                        height = 0.2)) +
  geom_point(aes(x = scale(LPFS_score), 
                 y = resid_mpsz_lig, color = "Group 5"), 
             size = 0.01, 
             alpha = 0.35,
             position = position_jitter(width = 0.2, 
                                        height = 0.2))

print(ligthart)

# Wielscher
wielscher <- ggplot(figure_df) +
  stat_smooth(aes(x = scale(PSS_TOT_score), 
                  y = resid_mpsz_wiel, color = "Group 1"), 
              method = "gam", 
              formula = y ~ x, 
              geom = "smooth", 
              size = 1, 
              alpha = 0.1) +
  stat_smooth(aes(x = scale(PIS_TOT_score), 
                  y = resid_mpsz_wiel, color = "Group 2"), 
              method = "gam", 
              formula = y ~ x, 
              geom = "smooth", 
              size = 1, 
              alpha = 0.1)+
  stat_smooth(aes(x = scale(PICS_score),
                  y = resid_mpsz_wiel, 
                  color = "Group 3"), 
              method = "gam", 
              formula = y ~ x, 
              geom = "smooth", 
              size = 1, 
              alpha = 0.1) +
  stat_smooth(aes(x = scale(MICS_score), 
                  y = resid_mpsz_wiel, 
                  color = "Group 4"), 
              method = "gam", 
              formula = y ~ x, 
              geom = "smooth", 
              size = 1, 
              alpha = 0.1)+
  stat_smooth(aes(x = scale(LPFS_score), 
                  y = resid_mpsz_wiel, 
                  color = "Group 5"), 
              method = "gam", 
              formula = y ~ x, 
              geom = "smooth", 
              size = 1, 
              alpha = 0.1) +
  scale_color_manual(values = c("Group 1" = "coral3", 
                                "Group 2" = "chartreuse3", 
                                "Group 3" = "dodgerblue4", 
                                "Group 4" = "darkorchid3", 
                                "Group 5" = "goldenrod3"),
                     labels = c("Prenatal stress", "Prenatal infections", 
                                "Pregnancy related inflammatory 
                                conditions", 
                                "Inflammatory medical conditions", 
                                "Lifestyle pro-inflammatory factors")) +
  facet_wrap(~ array, scales = "free_y", ncol = 2) +
  theme_bw() +
  labs(x = "Prenatal clinical inflammatory scores", 
       y = "Methylation profile score - Ligthart",  
       color = 'Prenatal predictor') + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        strip.background = element_rect(fill = 'moccasin'),
        strip.text = element_text(colour = 'black', 
                                  face = 'bold')) +
  geom_point(aes(x = scale(PSS_TOT_score), 
                 y = resid_mpsz_wiel, 
                 color = "Group 1"), 
             size = 0.01, 
             alpha = 0.35,
             position = position_jitter(width = 0.2, 
                                        height = 0.2)) +
  geom_point(aes(x = scale(PIS_TOT_score),
                 y = resid_mpsz_wiel,
                 color = "Group 2"),
             size = 0.01, 
             alpha = 0.35,
             position = position_jitter(width = 0.2, 
                                        height = 0.2)) +
  geom_point(aes(x = scale(PICS_score), 
                 y = resid_mpsz_wiel, 
                 color = "Group 3"), 
             size = 0.01, 
             alpha = 0.35,
             position = position_jitter(width = 0.2, 
                                        height = 0.2)) +
  geom_point(aes(x = scale(MICS_score), 
                 y = resid_mpsz_wiel, 
                 color = "Group 4"), 
             size = 0.01, 
             alpha = 0.35,
             position = position_jitter(width = 0.2, 
                                        height = 0.2)) +
  geom_point(aes(x = scale(LPFS_score), 
                 y = resid_mpsz_wiel, color = "Group 5"), 
             size = 0.01, 
             alpha = 0.35,
             position = position_jitter(width = 0.2, 
                                        height = 0.2))

print(wielscher)

# Combine plots 
ggarrange(ligthart, wielscher, labels = 'AUTO', common.legend = T, legend = 'bottom', ncol = 1, nrow = 2) 

#\

## Plot variance partitioning for MPS-CRP (with covariates and vars of aim 1) (MPSes combined)
# 450k 
df <- dplyr::select(figure_450k_df, c('IDC', 'Sample_Plate', 'LPFS_score', 'PICS_score', 'MICS_score', 'HsCRPmgL_g1', 'HsCRPmgL_g2', 'CRP_birth','PSS_TOT_score', 'PIS_TOT_score', 'GESTBIR', 'Bcell', 'CD4T', 'CD8T', 'Gran', 'Mono', 'NK', 'nRBC', 'EDUCM', 'AGE_M_v2', 'PARITY', 'SMOKE_ALL', 'INCOME'))

df2 <- dplyr::rename(df, 
             'Sample plate' = Sample_Plate, 
             'Lifestyle pro-inflammatory factors score' = LPFS_score,
             'Prenatal related inflammatory conditions score' = PICS_score,
             'Inflammatory medical conditions score ' = MICS_score,
             'CRP trimester 1' = HsCRPmgL_g1,
             'CRP trimester 2' = HsCRPmgL_g2,
             'CRP in cord blood' = CRP_birth,
             'Prenatal stress score' = PSS_TOT_score,
             'Prenatal infection sum score' = PIS_TOT_score,
             'Gestational age at birth' = GESTBIR,
             'B-cells' = Bcell,
             'CD4 T-cells' = CD4T,
             'CD8 T-cells' = CD8T,
             'Granulocytes' = Gran,
             'Monocytes' = Mono,
             'Natural killer cells' = NK,
             'Nucleated red blood cells' = nRBC,
             'Maternal education' = EDUCM,
             'Maternal age' = AGE_M_v2,
             'Parity' = PARITY,
             'Maternal tobacco use' = SMOKE_ALL,
             'Household income' = INCOME) 

formula <- ~ `Lifestyle pro-inflammatory factors score` + `Prenatal related inflammatory conditions score` + `Prenatal stress score` + `Prenatal infection sum score` +  `Maternal education` + `Household income` + `B-cells` + `CD4 T-cells` + `CD8 T-cells` + `Granulocytes` + `Monocytes` + `Natural killer cells` + `Nucleated red blood cells` + `Sample plate`

mpsz <- dplyr::select(figure_450k_df, 'mpsz_lig', 'mpsz_wiel')
mpsz <- dplyr::rename(mpsz, 'MPS - Ligthart' = mpsz_lig, 'MPS - Wielscher' = mpsz_wiel)
mpsz_t <- t(mpsz)

array_450k <- fitExtractVarPartModel(mpsz_t, formula, df2)
vp1 <- plotPercentBars(array_450k,col = magma(15))

# epic
df3 <- dplyr::select(figure_epic_df, c('IDC', 'Sample_Plate', 'LPFS_score', 'PICS_score', 'MICS_score', 'HsCRPmgL_g1', 'HsCRPmgL_g2', 'CRP_birth','PSS_TOT_score', 'PIS_TOT_score', 'GESTBIR', 'Bcell', 'CD4T', 'CD8T', 'Gran', 'Mono', 'NK', 'nRBC', 'EDUCM', 'AGE_M_v2', 'PARITY', 'SMOKE_ALL', 'INCOME'))

df4 <- dplyr::rename(df3, 
                     'Sample plate' = Sample_Plate, 
                     'Lifestyle pro-inflammatory factors score' = LPFS_score,
                     'Prenatal related inflammatory conditions score' = 
                       PICS_score,
                     'Inflammatory medical conditions score ' = MICS_score,
                     'CRP trimester 1' = HsCRPmgL_g1,
                     'CRP trimester 2' = HsCRPmgL_g2,
                     'CRP in cord blood' = CRP_birth,
                     'Prenatal stress score' = PSS_TOT_score,
                     'Prenatal infection sum score' = PIS_TOT_score,
                     'Gestational age at birth' = GESTBIR,
                     'B-cells' = Bcell,
                     'CD4 T-cells' = CD4T,
                     'CD8 T-cells' = CD8T,
                     'Granulocytes' = Gran,
                     'Monocytes' = Mono,
                     'Natural killer cells' = NK,
                     'Nucleated red blood cells' = nRBC,
                     'Maternal education' = EDUCM,
                     'Maternal age' = AGE_M_v2,
                     'Parity' = PARITY,
                     'Maternal tobacco use' = SMOKE_ALL,
                     'Household income' = INCOME) 

mpsz_epic <- dplyr::select(figure_epic_df, 'mpsz_lig', 'mpsz_wiel')
mpsz_epic <- dplyr::rename(mpsz_epic, 'MPS - Ligthart' = mpsz_lig, 'MPS - Wielscher' = mpsz_wiel)
mpsz_epic_t <- t(mpsz_epic)

array_epic <- fitExtractVarPartModel(mpsz_epic_t, formula, df4)

vp2 <- plotPercentBars(array_epic,col = magma(15))

# combine 450k and epic plots
ggarrange(vp1, vp2, ncol = 1, nrow = 2, labels = 'AUTO', common.legend = T, legend = 'bottom')

## Plot variance partitioning for MPS-CRP (with covariates and vars of aim 1) (MPS ligthart only)
# 450k 
library(viridis)
df <- dplyr::select(figure_450k_df, c('IDC', 'Sample_Plate', 'LPFS_score', 'PICS_score', 'MICS_score', 'HsCRPmgL_g1', 'HsCRPmgL_g2', 'CRP_birth','PSS_TOT_score', 'PIS_TOT_score', 'GESTBIR', 'Bcell', 'CD4T', 'CD8T', 'Gran', 'Mono', 'NK', 'nRBC', 'EDUCM', 'AGE_M_v2', 'PARITY', 'SMOKE_ALL', 'INCOME'))

df2 <- dplyr::rename(df, 
                     'Sample plate' = Sample_Plate, 
                     'Lifestyle pro-inflammatory factors score' = LPFS_score,
                     'Prenatal related inflammatory conditions score' = PICS_score,
                     'Inflammatory medical conditions score ' = MICS_score,
                     'CRP trimester 1' = HsCRPmgL_g1,
                     'CRP trimester 2' = HsCRPmgL_g2,
                     'CRP in cord blood' = CRP_birth,
                     'Prenatal stress score' = PSS_TOT_score,
                     'Prenatal infection sum score' = PIS_TOT_score,
                     'Gestational age at birth' = GESTBIR,
                     'B-cells' = Bcell,
                     'CD4 T-cells' = CD4T,
                     'CD8 T-cells' = CD8T,
                     'Granulocytes' = Gran,
                     'Monocytes' = Mono,
                     'Natural killer cells' = NK,
                     'Nucleated red blood cells' = nRBC,
                     'Maternal education' = EDUCM,
                     'Maternal age' = AGE_M_v2,
                     'Parity' = PARITY,
                     'Maternal tobacco use' = SMOKE_ALL,
                     'Household income' = INCOME) 

formula <- ~ `Lifestyle pro-inflammatory factors score` + `Prenatal related inflammatory conditions score` + `Prenatal stress score` + `Prenatal infection sum score` +  `Maternal education` + `Household income` + `B-cells` + `CD4 T-cells` + `CD8 T-cells` + `Granulocytes` + `Monocytes` + `Natural killer cells` + `Nucleated red blood cells` + `Sample plate`

mpsz <- dplyr::select(figure_450k_df, 'mpsz_lig')
mpsz <- dplyr::rename(mpsz, 'MPS - Ligthart' = mpsz_lig)
mpsz_t <- t(mpsz)

array_450k <- fitExtractVarPartModel(mpsz_t, formula, df2)
vp1 <- plotPercentBars(array_450k, col = magma(15))

# epic
df3 <- dplyr::select(figure_epic_df, c('IDC', 'Sample_Plate', 'LPFS_score', 'PICS_score', 'MICS_score', 'HsCRPmgL_g1', 'HsCRPmgL_g2', 'CRP_birth','PSS_TOT_score', 'PIS_TOT_score', 'GESTBIR', 'Bcell', 'CD4T', 'CD8T', 'Gran', 'Mono', 'NK', 'nRBC', 'EDUCM', 'AGE_M_v2', 'PARITY', 'SMOKE_ALL', 'INCOME'))

df4 <- dplyr::rename(df3, 
                     'Sample plate' = Sample_Plate, 
                     'Lifestyle pro-inflammatory factors score' = LPFS_score,
                     'Prenatal related inflammatory conditions score' = 
                       PICS_score,
                     'Inflammatory medical conditions score ' = MICS_score,
                     'CRP trimester 1' = HsCRPmgL_g1,
                     'CRP trimester 2' = HsCRPmgL_g2,
                     'CRP in cord blood' = CRP_birth,
                     'Prenatal stress score' = PSS_TOT_score,
                     'Prenatal infection sum score' = PIS_TOT_score,
                     'Gestational age at birth' = GESTBIR,
                     'B-cells' = Bcell,
                     'CD4 T-cells' = CD4T,
                     'CD8 T-cells' = CD8T,
                     'Granulocytes' = Gran,
                     'Monocytes' = Mono,
                     'Natural killer cells' = NK,
                     'Nucleated red blood cells' = nRBC,
                     'Maternal education' = EDUCM,
                     'Maternal age' = AGE_M_v2,
                     'Parity' = PARITY,
                     'Maternal tobacco use' = SMOKE_ALL,
                     'Household income' = INCOME) 

mpsz_epic <- dplyr::select(figure_epic_df, 'mpsz_lig')
mpsz_epic <- dplyr::rename(mpsz_epic, 'MPS - Ligthart' = mpsz_lig)
mpsz_epic_t <- t(mpsz_epic)

array_epic <- fitExtractVarPartModel(mpsz_epic_t, formula, df4)

vp2 <- plotPercentBars(array_epic,col = magma(15))

# combine 450k and epic plots
ggarrange(vp1, vp2, ncol = 1, nrow = 2, labels = 'AUTO', common.legend = T, legend = 'bottom')

#\

#====================================================

### Step 2: associate MPS ligthart/wielscher to serum CRP variables [MPS validation]

## Linear regression 
# Create 450k / epic specific dataframes
imp.data_long <- complete(imp.data.mids, include = T, action = "long")
imp.data.long_450k <- subset(imp.data_long, array == '450k')
imp.data.long_epic <- subset(imp.data_long, array == 'epic')

imp.data.mids_450k <- as.mids(imp.data.long_450k)
imp.data.mids_epic <- as.mids(imp.data.long_epic)

saveRDS(imp.data.mids_450k, 'imp.data.mids_450k.rds')
saveRDS(imp.data.mids_epic, 'imp.data.mids_epic.rds')

# Specify exposure and outcome vars 
crp_vars <- c("HsCRPmgL_g1", "HsCRPmgL_g2", "CRP_birth", "CRPCHILD5")
mps_vars <- c('mpsz_lig', 'mpsz_wiel')

# Create empty dataframes 
results_450k_univar_model <- data.frame()
results_450k_batch_model <- data.frame()
results_450k_batch_cell_model <- data.frame()

results_epic_univar_model <- data.frame()
results_epic_batch_model <- data.frame()
results_epic_batch_cell_model <- data.frame()

# Specify loop 
for (outcome in crp_vars) {
  
  for (exposure in mps_vars) {
    
    # Specify model formulas 
    univar_model <- paste0("scale(",outcome,") ~",exposure) 
    batch_model <- paste0("scale(",outcome,") ~ ", exposure, "+ Sample_Plate")
    batch_cell_model <- paste0("scale(",outcome,") ~", exposure, "+ Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC")
    
    # Extract coefficients of interest for exposure
    univar_model_output_450k <- summary(pool(with(imp.data.mids_450k, lm(as.formula(univar_model)))), conf.int = TRUE)[2, c(2, 3, 6, 7, 8)]
    univar_model_output_epic <- summary(pool(with(imp.data.mids_epic, lm(as.formula(univar_model)))), conf.int = TRUE)[2, c(2, 3, 6, 7, 8)]
    
    batch_model_output_450k <- summary(pool(with(imp.data.mids_450k, lm(as.formula(batch_model)))), conf.int = TRUE)[2, c(2, 3, 6, 7, 8)]
    batch_model_output_epic <- summary(pool(with(imp.data.mids_epic, lm(as.formula(batch_model)))), conf.int = TRUE)[2, c(2, 3, 6, 7, 8)]
    
    batch_cell_model_output_450k <- summary(pool(with(imp.data.mids_450k, lm(as.formula(batch_cell_model)))), conf.int = TRUE)[2, c(2, 3, 6, 7, 8)]
    batch_cell_model_output_epic <- summary(pool(with(imp.data.mids_epic, lm(as.formula(batch_cell_model)))), conf.int = TRUE)[2, c(2, 3, 6, 7, 8)]
    
    # Adding outcome and exposure information
    univar_model_output_450k$outcome <- outcome
    univar_model_output_450k$exposure <- exposure
    univar_model_output_epic$outcome <- outcome
    univar_model_output_epic$exposure <- exposure
    
    batch_model_output_450k$outcome <- outcome
    batch_model_output_450k$exposure <- exposure
    batch_model_output_epic$exposure <- exposure
    batch_model_output_epic$outcome <- outcome
    
    batch_cell_model_output_450k$outcome <- outcome
    batch_cell_model_output_450k$exposure <- exposure
    batch_cell_model_output_epic$exposure <- exposure
    batch_cell_model_output_epic$outcome <- outcome
  
    # Adding results to data frames
    results_450k_univar_model <- rbind(results_450k_univar_model,
                                       univar_model_output_450k)
    results_epic_univar_model <- rbind(results_epic_univar_model,
                                       univar_model_output_epic)
    
    results_450k_batch_model <- rbind(results_450k_batch_model,
                                      batch_model_output_450k)
    results_epic_batch_model <- rbind(results_epic_batch_model,
                                      batch_model_output_epic)
    
    results_450k_batch_cell_model <- rbind(results_450k_batch_cell_model,
                                           batch_cell_model_output_450k)
    results_epic_batch_cell_model <- rbind(results_epic_batch_cell_model,
                                           batch_cell_model_output_epic)
  }
 }

## Meta-analysis
# Specify function
meta_analyze <- function(data) {
  rma(yi = data$estimate, sei = data$std.error, method = "FE")
}

# Create input files 
combined_results_univar_model_ligthart <- bind_rows(
  subset(results_450k_univar_model, exposure == 'mpsz_lig'),
  subset(results_epic_univar_model, exposure == 'mpsz_lig')
)

combined_results_univar_model_wielscher <- bind_rows(
  subset(results_450k_univar_model, exposure == 'mpsz_wiel'),
  subset(results_epic_univar_model, exposure == 'mpsz_wiel')
)

combined_results_batch_model_ligthart <- bind_rows(
  subset(results_450k_batch_model, exposure == 'mpsz_lig'),
  subset(results_epic_batch_model, exposure == 'mpsz_lig')
)

combined_results_batch_model_wielscher <- bind_rows(
  subset(results_450k_batch_model, exposure == 'mpsz_wiel'),
  subset(results_epic_batch_model, exposure == 'mpsz_wiel')
)

combined_results_batch_cell_model_ligthart <- bind_rows(
  subset(results_450k_batch_cell_model, exposure == 'mpsz_lig'),
  subset(results_epic_batch_cell_model, exposure == 'mpsz_lig')
)

combined_results_batch_cell_model_wielscher <- bind_rows(
  subset(results_450k_batch_cell_model, exposure == 'mpsz_wiel'),
  subset(results_epic_batch_cell_model, exposure == 'mpsz_wiel')
)

# Specify  unique outcome
unique_outcomes_univar_ligthart <- unique(combined_results_univar_model_ligthart$outcome)
unique_outcomes_univar_wielscher <- unique(combined_results_univar_model_wielscher$outcome)

unique_outcomes_batch_ligthart <- unique(combined_results_batch_model_ligthart$outcome)
unique_outcomes_batch_wielscher <- unique(combined_results_batch_model_wielscher$outcome)

unique_outcomes_batch_cell_ligthart <- unique(combined_results_batch_cell_model_ligthart$outcome)
unique_outcomes_batch_cell_wielscher <- unique(combined_results_batch_cell_model_wielscher$outcome)

# Apply meta-analysis function to each unique outcome
meta_analyze_outcomes <- function(unique_outcomes, combined_results_model) {
  lapply(unique_outcomes, function(outcome) {
    subset_data <- combined_results_model[combined_results_model$outcome == outcome, ]
    meta_analyze(subset_data)
  })
}

meta_results_univar_ligthart <- meta_analyze_outcomes(unique_outcomes_univar_ligthart, combined_results_univar_model_ligthart)
meta_results_univar_wielscher <- meta_analyze_outcomes(unique_outcomes_univar_wielscher, combined_results_univar_model_wielscher)

meta_results_batch_ligthart <- meta_analyze_outcomes(unique_outcomes_batch_ligthart, combined_results_batch_model_ligthart)
meta_results_batch_wielscher <- meta_analyze_outcomes(unique_outcomes_batch_wielscher, combined_results_batch_model_wielscher)

meta_results_batch_cell_ligthart <- meta_analyze_outcomes(unique_outcomes_batch_cell_ligthart, combined_results_batch_cell_model_ligthart)
meta_results_batch_cell_wielscher <- meta_analyze_outcomes(unique_outcomes_batch_cell_wielscher, combined_results_batch_cell_model_wielscher)

#\

#====================================================

### Step 3: associate prenatal predictors to MPS ligthart/wielscher [aim 1]

## Linear regression

# First check assumptions
check_model_assumptions <- function(model) {
  # Distribution of residuals
  res <- resid(model)
  par(mfrow = c(2, 2))
  qqnorm(res)
  qqline(res)
  plot(density(res), main = "Density Plot of Residuals")
  # Heteroscedasticity assumption
  plot(fitted(model), res, main = "Residuals vs. Fitted Values")
  abline(h = 0, lty = 2, col = "red")
}

data_450k <- complete(imp.data.mids_450k, 30)
data_epic <- complete(imp.data.mids_epic, 30)

check_model_assumptions(lm(mpsz_wiel ~ EDUCM_3l + INCOME + PIS_TOT_score + PSS_TOT_score + LPFS_score + PICS_score + MICS_score + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data = data_450k))
check_model_assumptions(lm(mpsz_wiel ~ EDUCM_3l + INCOME + PIS_TOT_score + PSS_TOT_score + LPFS_score + PICS_score + MICS_score + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC, data = data_epic))

# Create empty dataframes 
results_full_model_aim1_450k <- data.frame()
results_nocell_model_aim1_450k <- data.frame()

results_full_model_aim1_epic <- data.frame()
results_nocell_model_aim1_epic <- data.frame()

# Specify loop for multiple regressions (for all predictors individually)
prenatal_scores_vars <- c("PSS_TOT_score", "PIS_TOT_score", "LPFS_score", "PICS_score", "MICS_score")

for (outcome in mps_vars) {
  for (predictor in prenatal_scores_vars) {
    
    # Specify model formulas
    full_model <- paste0(outcome, "~ scale(",predictor,") + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + GENDER +  EDUCM_3l + INCOME")

    nocell_model <- paste0(outcome, "~ scale(",predictor,") + Sample_Plate + GENDER +  EDUCM_3l + INCOME")

    # Extract coefficients of interest for exposure
    full_model_450k <- summary(pool(with(imp.data.mids_450k, lm(as.formula(full_model)))), conf.int = TRUE)[2, c(2, 3, 6, 7, 8)]
    nocell_model_450k <- summary(pool(with(imp.data.mids_450k, lm(as.formula(nocell_model)))), conf.int = TRUE)[2, c(2, 3, 6, 7, 8)]
  
    full_model_epic <- summary(pool(with(imp.data.mids_epic, lm(as.formula(full_model)))), conf.int = TRUE)[2, c(2, 3, 6, 7, 8)]
    nocell_model_epic <- summary(pool(with(imp.data.mids_epic, lm(as.formula(nocell_model)))), conf.int = TRUE)[2, c(2, 3, 6, 7, 8)]
    
    # Adding outcome and exposure information
    full_model_450k$outcome <- outcome
    full_model_450k$predictor <- predictor
    nocell_model_450k$outcome <- outcome
    nocell_model_450k$predictor <- predictor
    
    full_model_epic$outcome <- outcome
    full_model_epic$predictor <- predictor
    nocell_model_epic$outcome <- outcome
    nocell_model_epic$predictor <- predictor
    
    # Adding results to data frames
    results_full_model_aim1_450k <- rbind(results_full_model_aim1_450k,
                                          full_model_450k)
    results_nocell_model_aim1_450k <- rbind(results_nocell_model_aim1_450k, 
                                            nocell_model_450k)
    
    results_full_model_aim1_epic <- rbind(results_full_model_aim1_epic,
                                          full_model_epic)
    results_nocell_model_aim1_epic <- rbind(results_nocell_model_aim1_epic, 
                                            nocell_model_epic)
  }
}

## Meta-analysis
# Create input files 
combined_results_full_model_ligthart <- bind_rows(
  subset(results_full_model_aim1_450k, outcome == 'mpsz_lig'),
  subset(results_full_model_aim1_epic, outcome == 'mpsz_lig'),
)
combined_results_full_model_wielscher <- bind_rows(
  subset(results_full_model_aim1_450k, outcome == 'mpsz_wiel'),
  subset(results_full_model_aim1_epic, outcome == 'mpsz_wiel'),
)

combined_results_nocell_model_ligthart <- bind_rows(
  subset(results_nocell_model_aim1_450k, outcome == 'mpsz_lig'),
  subset(results_nocell_model_aim1_epic, outcome == 'mpsz_lig'),
)
combined_results_nocell_model_wielscher <- bind_rows(
  subset(results_nocell_model_aim1_450k, outcome == 'mpsz_wiel'),
  subset(results_nocell_model_aim1_epic, outcome == 'mpsz_wiel'),
)

# Specify unique predictor
unique_predictors <- unique(combined_results_full_model_ligthart$predictor)

# Apply meta-analysis function to each unique predictor
meta_analyze_predictors <- function(unique_predictors, combined_results_model) {
  lapply(unique_predictors, function(predictor) {
    subset_data <- combined_results_model[combined_results_model$predictor == predictor, ]
    meta_analyze(subset_data)
  })
}

meta_results_full_model_ligthart <- meta_analyze_predictors(unique_predictors, combined_results_full_model_ligthart)
meta_results_full_model_wielscher <- meta_analyze_predictors(unique_predictors, combined_results_full_model_wielscher)

meta_results_nocell_model_ligthart <- meta_analyze_predictors(unique_predictors, combined_results_nocell_model_ligthart)
meta_results_nocell_model_wielscher <- meta_analyze_predictors(unique_predictors, combined_results_nocell_model_wielscher)

#\

## Specify loop for linear multivariate regression witih all predictors in one model
# Create empty dataframes 
results_full_model_aim1_450k_multivariate <- data.frame()
results_nocell_model_aim1_450k_multivariate <- data.frame()

results_full_model_aim1_epic_multivariate <- data.frame()
results_nocell_model_aim1_epic_multivariate <- data.frame()

# Specify loop for multiple regressions (for all predictors individually)
#' given the high colinearity between PICS and MICS, we will exclude MICS in this model, otherwise the model will not run 

for (outcome in mps_vars) {
    # Specify model formulas
    full_model <- paste0(outcome, "~ scale(PSS_TOT_score) + scale(PIS_TOT_score) + scale(LPFS_score) + scale(PICS_score) + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + GENDER +  EDUCM_3l + INCOME")
    
    nocell_model <- paste0(outcome, "~ scale(PSS_TOT_score) + scale(PIS_TOT_score) + scale(LPFS_score) + scale(PICS_score) + Sample_Plate + GENDER +  EDUCM_3l + INCOME")
    
    # Extract coefficients of interest for exposure
    full_model_450k <- summary(pool(with(imp.data.mids_450k, lm(as.formula(full_model)))), conf.int = TRUE)[2:5, c(2, 3, 6, 7, 8)]
    nocell_model_450k <- summary(pool(with(imp.data.mids_450k, lm(as.formula(nocell_model)))), conf.int = TRUE)[2:5, c(2, 3, 6, 7, 8)]
    
    full_model_epic <- summary(pool(with(imp.data.mids_epic, lm(as.formula(full_model)))), conf.int = TRUE)[2:5, c(2, 3, 6, 7, 8)]
    nocell_model_epic <- summary(pool(with(imp.data.mids_epic, lm(as.formula(nocell_model)))), conf.int = TRUE)[2:5, c(2, 3, 6, 7, 8)]
    
    # Adding outcome and exposure information
    full_model_450k$outcome <- outcome
    nocell_model_450k$outcome <- outcome
    
    full_model_epic$outcome <- outcome
    nocell_model_epic$outcome <- outcome
    
    # Adding results to data frames
    results_full_model_aim1_450k_multivariate <- rbind(results_full_model_aim1_450k_multivariate, full_model_450k)
    results_nocell_model_aim1_450k_multivariate <- rbind(results_nocell_model_aim1_450k_multivariate, nocell_model_450k)
    
    results_full_model_aim1_epic_multivariate <- rbind (results_full_model_aim1_epic_multivariate, full_model_epic)
    results_nocell_model_aim1_epic_multivariate <- rbind(results_nocell_model_aim1_epic_multivariate, nocell_model_epic)
}

# Add predictor names to dataframe
results_full_model_aim1_450k_multivariate$predictor <- c('scale(PSS_TOT_score)', 'scale(PIS_TOT_score)', 'scale(LPFS_score)', 'scale(PICS_score)','scale(PSS_TOT_score)', 'scale(PIS_TOT_score)', 'scale(LPFS_score)', 'scale(PICS_score)')
results_nocell_model_aim1_450k_multivariate$predictor <- c('scale(PSS_TOT_score)', 'scale(PIS_TOT_score)', 'scale(LPFS_score)', 'scale(PICS_score)','scale(PSS_TOT_score)', 'scale(PIS_TOT_score)', 'scale(LPFS_score)', 'scale(PICS_score)')
results_full_model_aim1_epic_multivariate$predictor <- c('scale(PSS_TOT_score)', 'scale(PIS_TOT_score)', 'scale(LPFS_score)', 'scale(PICS_score)','scale(PSS_TOT_score)', 'scale(PIS_TOT_score)', 'scale(LPFS_score)', 'scale(PICS_score)')
results_nocell_model_aim1_epic_multivariate$predictor <- c('scale(PSS_TOT_score)', 'scale(PIS_TOT_score)', 'scale(LPFS_score)', 'scale(PICS_score)','scale(PSS_TOT_score)', 'scale(PIS_TOT_score)', 'scale(LPFS_score)', 'scale(PICS_score)')

## Meta-analysis
# Create input files 
combined_results_full_model_ligthart_multivariate <- bind_rows(
  subset(results_full_model_aim1_450k_multivariate, outcome == 'mpsz_lig'),
  subset(results_full_model_aim1_epic_multivariate, outcome == 'mpsz_lig'),
)
combined_results_full_model_wielscher_multivariate <- bind_rows(
  subset(results_full_model_aim1_450k_multivariate, outcome == 'mpsz_wiel'),
  subset(results_full_model_aim1_epic_multivariate, outcome == 'mpsz_wiel'),
)

combined_results_nocell_model_ligthart_multivariate <- bind_rows(
  subset(results_nocell_model_aim1_450k_multivariate, outcome == 'mpsz_lig'),
  subset(results_nocell_model_aim1_epic_multivariate, outcome == 'mpsz_lig'),
)
combined_results_nocell_model_wielscher_multivariate <- bind_rows(
  subset(results_nocell_model_aim1_450k_multivariate, outcome == 'mpsz_wiel'),
  subset(results_nocell_model_aim1_epic_multivariate, outcome == 'mpsz_wiel'),
)

# Specify unique predictor
unique_predictors <- unique(combined_results_full_model_ligthart_multivariate$predictor)

# Apply meta-analysis function to each unique predictor
meta_analyze_predictors <- function(unique_predictors, combined_results_model) {
  lapply(unique_predictors, function(predictor) {
    subset_data <- combined_results_model[combined_results_model$predictor == predictor, ]
    meta_analyze(subset_data)
  })
}

meta_results_full_model_ligthart_multivariate <- meta_analyze_predictors(unique_predictors, combined_results_full_model_ligthart_multivariate)
meta_results_full_model_wielscher_multivariate <- meta_analyze_predictors(unique_predictors, combined_results_full_model_wielscher_multivariate)

meta_results_nocell_model_ligthart_multivariate <- meta_analyze_predictors(unique_predictors, combined_results_nocell_model_ligthart_multivariate)
meta_results_nocell_model_wielscher_multivariate <- meta_analyze_predictors(unique_predictors, combined_results_nocell_model_wielscher_multivariate)

#\

#====================================================

### Step 4: associate child prs of CRP to MPS ligthart/wielscher [aim 1 extra]

# Load PRS-CRP data 
load('final_crp_pgs_15032024.rdata')
prs_crp_df <- crp_pgs

### First perform validation checks with serum crp 
# Check distribution - nicely normally distributed 
pgs_crp_hist <- ggplot(prs_crp_df, aes(x = ZR_crp_pgs)) +
  geom_histogram(binwidth = 1, fill = "purple", color = "white", alpha = 0.4) +
  labs(x = "PGS of CRP", y = "Frequency", title = "Histogram of PGS of CRP") +
  theme_minimal()

# Extract one imputed dataframe and merge dataframes
imp.data_long <- complete(imp.data.mids, include = T, action = "long")
single.imp.data <- subset(imp.data_long, .imp == 30)
single.imp.data_prs <- merge(single.imp.data, prs_crp_df, all.x = T, by = 'IDC')

# Check distribution with untransformed serum CRP
scatter1 <- ggplot(single.imp.data_prs, aes(x = CRP_birth, y = ZR_crp_pgs)) +
  geom_point(color = 'red', alpha = 0.2, size = 2) + labs(x = "CRP in cord blood", y = "PGS of CRP", title = "Serum crp in cord blood vs pgs of CRP") + theme_minimal() +
  geom_smooth(method = "gam", se = T, color = "black", linetype = 'dotted')

scatter2 <-ggplot(single.imp.data_prs, aes(x = CRPCHILD5, y = ZR_crp_pgs)) +
  geom_point(color = 'orange', alpha = 0.2, size = 2) + labs(x = "CRP in child at age 5", y = "PGS of CRP", title = "Serum crp in whole blood (child age 5) vs pgs of CRP") + theme_minimal()+
  geom_smooth(method = "gam", se = T, color = "black", linetype = 'dotted') 

scatter3 <- ggplot(single.imp.data_prs, aes(x = HsCRPmgL_g1, y = ZR_crp_pgs)) +
  geom_point(color = 'blue', alpha = 0.2, size = 2) + labs(x = "CRP in trimester 1", y = "PGS of CRP", title = "Serum crp in whole blood (trimester 1) vs pgs of CRP") + theme_minimal() +
  geom_smooth(method = "gam", se = T, color = "black", linetype = 'dotted') 

scatter4 <-ggplot(single.imp.data_prs, aes(x = HsCRPmgL_g2, y = ZR_crp_pgs)) +
  geom_point(color = 'darkgreen', alpha = 0.2, size = 2) + labs(x = "CRP in trimester 2", y = "PGS of CRP", title = "Serum crp in whole blood (trimester 2) vs pgs of CRP") + theme_minimal() + geom_smooth(method = "gam", se = T, color = "black", linetype = 'dotted') 

ggarrange(pgs_crp_hist, scatter1, scatter2, scatter3, scatter4, labels = 'auto')

# Check distribution with log2 transformed serum CRP
scatter1log <- ggplot(single.imp.data_prs, aes(x = log2(CRP_birth), y = ZR_crp_pgs)) +
  geom_point(color = 'red', alpha = 0.2, size = 2) + labs(x = "CRP in cord blood", y = "PGS of CRP", title = "Log2 serum crp in cord blood vs pgs of CRP") + theme_minimal() +
  geom_smooth(method = "gam", se = T, color = "black", linetype = 'dotted')

scatter2log <-ggplot(single.imp.data_prs, aes(x = log2(CRPCHILD5), y = ZR_crp_pgs)) +
  geom_point(color = 'orange', alpha = 0.2, size = 2) + labs(x = "CRP in child at age 5", y = "PGS of CRP", title = "Log2 serum crp in whole blood (child age 5) vs pgs of CRP") + theme_minimal()+
  geom_smooth(method = "gam", se = T, color = "black", linetype = 'dotted') 

scatter3log <- ggplot(single.imp.data_prs, aes(x = log2(HsCRPmgL_g1), y = ZR_crp_pgs)) +
  geom_point(color = 'blue', alpha = 0.2, size = 2) + labs(x = "CRP in trimester 1", y = "PGS of CRP", title = "Log2 serum crp in whole blood (trimester 1) vs pgs of CRP") + theme_minimal() +
  geom_smooth(method = "gam", se = T, color = "black", linetype = 'dotted') 

scatter4log <-ggplot(single.imp.data_prs, aes(x = log2(HsCRPmgL_g2), y = ZR_crp_pgs)) +
  geom_point(color = 'darkgreen', alpha = 0.2, size = 2) + labs(x = "CRP in trimester 2", y = "PGS of CRP", title = "Log2 serum crp in whole blood (trimester 2) vs pgs of CRP") + theme_minimal() + geom_smooth(method = "gam", se = T, color = "black", linetype = 'dotted') 

ggarrange(pgs_crp_hist, scatter1log, scatter2log, scatter3log, scatter4log, labels = 'auto')

# Check correlation with serum CRP
cor(single.imp.data_prs$HsCRPmgL_g2, single.imp.data_prs$ZR_crp_pgs, use = 'complete.obs', method = 'spearman') #0.21
cor(single.imp.data_prs$HsCRPmgL_g1, single.imp.data_prs$ZR_crp_pgs, use = 'complete.obs',method = 'spearman') #0.22
cor(single.imp.data_prs$CRPCHILD5, single.imp.data_prs$ZR_crp_pgs, use = 'complete.obs',method = 'spearman') #0.23
cor(single.imp.data_prs$CRP_birth, single.imp.data_prs$ZR_crp_pgs, use = 'complete.obs', method = 'spearman') #0.04

pgs_crp_validation <- data.frame(
  Outcome = c("CRP_tri1", "CRP_tri2", "CRP_cord_blood", "CRP_child_age5"),
  Spearman_correlation = c(0.2208082, 0.2085726, 0.0435209, 0.2329319)
)

# Check association with serum CRP
summary(lm(scale(HsCRPmgL_g1) ~ ZR_crp_pgs, data = single.imp.data_prs), conf.int = T)
summary(lm(scale(HsCRPmgL_g2) ~ ZR_crp_pgs, data = single.imp.data_prs), conf.int = T)
summary(lm(scale(CRP_birth) ~ ZR_crp_pgs, data = single.imp.data_prs), conf.int = T)
summary(lm(scale(CRPCHILD5) ~ ZR_crp_pgs, data = single.imp.data_prs), conf.int = T)

pgs_crp_validation$univar_lm_beta <- c(0.16853, 0.17344, 0.016136, 0.0808956)
pgs_crp_validation$univar_lm_pval <- c(2.95e-13, 4.69e-16, 0.478, 0.00244)

summary(lm(scale(HsCRPmgL_g1) ~ ZR_crp_pgs + GENDER +  EDUCM_3l + INCOME, data = single.imp.data_prs), conf.int = T)
summary(lm(scale(HsCRPmgL_g2) ~ ZR_crp_pgs + GENDER +  EDUCM_3l + INCOME, data = single.imp.data_prs), conf.int = T)
summary(lm(scale(CRP_birth) ~ ZR_crp_pgs + GENDER +  EDUCM_3l + INCOME, data = single.imp.data_prs), conf.int = T)
summary(lm(scale(CRPCHILD5) ~ ZR_crp_pgs + GENDER +  EDUCM_3l + INCOME, data = single.imp.data_prs), conf.int = T)

pgs_crp_validation$covar_adjusted_lm_beta <- c(0.154630, 0.16247, 0.01332, 0.07708)
pgs_crp_validation$covar_adjusted_lm_pval <- c(1.43e-11, 1.48e-14, 0.5586, 0.00387)

write_xlsx(pgs_crp_validation, 'pgs-crp_serum-crp_validation.xlsx')

# Check incremental R2 with serum child CRP
summary(lm(scale(CRPCHILD5) ~ GENDER +  EDUCM_3l + INCOME, data = single.imp.data_prs), conf.int = T) #r2 = 0.01073

summary(lm(scale(CRPCHILD5) ~ ZR_crp_pgs + GENDER +  EDUCM_3l + INCOME, data = single.imp.data_prs), conf.int = T)  #r2 = 0.01677

# Incremental r2
0.01677 - 0.01073 # 0.00604

#\

### Check association with MPS of CRP
## Merge dataframes 
imp.data_long <- complete(imp.data.mids, include = T, action = "long")

# Initialize a list to store the results
imp_data_prs_list <- list()

# Loop through the imputed values
for (i in 0:30) {
  # Subset the data for each imputation
  single_imp_data <- subset(imp.data_long, .imp == i)
  # Merge with prs_crp_df
  single_imp_data_prs <- merge(single_imp_data, prs_crp_df, all.x = TRUE, by = 'IDC')
  # Store the result in the list
  imp_data_prs_list[[i+1]] <- single_imp_data_prs
}

imp.data_mids_prs <- datalist2mids(imp_data_prs_list)

# Split into 450k and epic and save as mids 
imp.data_mids_prs_long <- complete(imp.data_mids_prs, include = T, action = "long")

imp.data_mids_prs_450k_long <- subset(imp.data_mids_prs_long, array == '450k')
imp.data_mids_prs_epic_long <- subset(imp.data_mids_prs_long, array == 'epic')

imp.data_mids_prs_450k <- as.mids(imp.data_mids_prs_450k_long)
imp.data_mids_prs_epic <- as.mids(imp.data_mids_prs_epic_long)

## Run regression for epic and 450k 
# Create empty dataframes 
results_full_model_prs_450k_multivariate <- data.frame()
results_nocell_model_prs_450k_multivariate <- data.frame()

results_full_model_prs_epic_multivariate <- data.frame()
results_nocell_model_prs_epic_multivariate <- data.frame()

# Specify loop for multiple regressions (for all predictors individually)
for (outcome in mps_vars) {
  # Specify model formulas
  full_model <- paste0(outcome, "~  ZR_crp_pgs + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + GENDER +  EDUCM_3l + INCOME")
  
  nocell_model <- paste0(outcome, "~ ZR_crp_pgs + Sample_Plate + GENDER +  EDUCM_3l + INCOME")
  
  # Extract coefficients of interest for exposure
  full_model_450k <- summary(pool(with(imp.data_mids_prs_450k, lm(as.formula(full_model)))), conf.int = TRUE)[2, c(2, 3, 6, 7, 8)]
  nocell_model_450k <- summary(pool(with(imp.data_mids_prs_450k, lm(as.formula(nocell_model)))), conf.int = TRUE)[2, c(2, 3, 6, 7, 8)]
  
  full_model_epic <- summary(pool(with(imp.data_mids_prs_epic, lm(as.formula(full_model)))), conf.int = TRUE)[2, c(2, 3, 6, 7, 8)]
  nocell_model_epic <- summary(pool(with(imp.data_mids_prs_epic, lm(as.formula(nocell_model)))), conf.int = TRUE)[2, c(2, 3, 6, 7, 8)]
  
  # Adding outcome and exposure information
  full_model_450k$outcome <- outcome
  nocell_model_450k$outcome <- outcome
  
  full_model_epic$outcome <- outcome
  nocell_model_epic$outcome <- outcome
  
  # Adding results to data frames
  results_full_model_prs_450k_multivariate <- rbind(results_full_model_prs_450k_multivariate, full_model_450k)
  results_nocell_model_prs_450k_multivariate <- rbind(results_nocell_model_prs_450k_multivariate, nocell_model_450k)
  
  results_full_model_prs_epic_multivariate <- rbind (results_full_model_prs_epic_multivariate, full_model_epic)
  results_nocell_model_prs_epic_multivariate <- rbind(results_nocell_model_prs_epic_multivariate, nocell_model_epic)
}

## Meta analyze results
# Create input files 
combined_results_full_model_ligthart_multivariate <- bind_rows(
  subset(results_full_model_prs_450k_multivariate, outcome == 'mpsz_lig'),
  subset(results_full_model_prs_epic_multivariate, outcome == 'mpsz_lig'),
)
combined_results_full_model_ligthart_multivariate$exposure <- c('prs_crp', 'prs_crp')

combined_results_full_model_wielscher_multivariate <- bind_rows(
  subset(results_full_model_prs_450k_multivariate, outcome == 'mpsz_wiel'),
  subset(results_full_model_prs_epic_multivariate, outcome == 'mpsz_wiel'),
)
combined_results_full_model_wielscher_multivariate$exposure <- c('prs_crp', 'prs_crp')

combined_results_nocell_model_ligthart_multivariate <- bind_rows(
  subset(results_nocell_model_prs_450k_multivariate, outcome == 'mpsz_lig'),
  subset(results_nocell_model_prs_epic_multivariate, outcome == 'mpsz_lig'),
)
combined_results_nocell_model_ligthart_multivariate$exposure <- c('prs_crp', 'prs_crp')

combined_results_nocell_model_wielscher_multivariate <- bind_rows(
  subset(results_nocell_model_prs_450k_multivariate, outcome == 'mpsz_wiel'),
  subset(results_nocell_model_prs_epic_multivariate, outcome == 'mpsz_wiel'),
)
combined_results_nocell_model_wielscher_multivariate$exposure <- c('prs_crp', 'prs_crp')

# Apply meta-analysis
rma(yi = estimate, sei = std.error, data = combined_results_full_model_ligthart_multivariate)
rma(yi = estimate, sei = std.error, data = combined_results_full_model_wielscher_multivariate)
rma(yi = estimate, sei = std.error, data = combined_results_nocell_model_ligthart_multivariate)
rma(yi = estimate, sei = std.error, data = combined_results_nocell_model_wielscher_multivariate)

## Repeat variance partitioning plot (both MPSes)
figure_450k_df2 <- complete(imp.data_mids_prs_450k, 30)
figure_epic_df2 <- complete(imp.data_mids_prs_epic, 30)

# 450k 
df3 <- dplyr::select(figure_450k_df2, c('IDC', 'Sample_Plate', 'LPFS_score', 'PICS_score', 'MICS_score', 'HsCRPmgL_g1', 'HsCRPmgL_g2', 'CRP_birth','PSS_TOT_score', 'PIS_TOT_score', 'GESTBIR', 'Bcell', 'CD4T', 'CD8T', 'Gran', 'Mono', 'NK', 'nRBC', 'EDUCM', 'AGE_M_v2', 'PARITY', 'SMOKE_ALL', 'INCOME', 'ZR_crp_pgs'))

df4 <- dplyr::rename(df3, 
                     'Sample plate' = Sample_Plate, 
                     'Lifestyle pro-inflammatory factors score' = LPFS_score,
                     'Prenatal related inflammatory conditions score' = PICS_score,
                     'Inflammatory medical conditions score ' = MICS_score,
                     'CRP trimester 1' = HsCRPmgL_g1,
                     'CRP trimester 2' = HsCRPmgL_g2,
                     'CRP in cord blood' = CRP_birth,
                     'Prenatal stress score' = PSS_TOT_score,
                     'Prenatal infection sum score' = PIS_TOT_score,
                     'Gestational age at birth' = GESTBIR,
                     'B-cells' = Bcell,
                     'CD4 T-cells' = CD4T,
                     'CD8 T-cells' = CD8T,
                     'Granulocytes' = Gran,
                     'Monocytes' = Mono,
                     'Natural killer cells' = NK,
                     'Nucleated red blood cells' = nRBC,
                     'Maternal education' = EDUCM,
                     'Maternal age' = AGE_M_v2,
                     'Parity' = PARITY,
                     'Maternal tobacco use' = SMOKE_ALL,
                     'Household income' = INCOME,
                      'PRS of CRP (child)' = ZR_crp_pgs) 

formula <- ~ `PRS of CRP (child)` + `Lifestyle pro-inflammatory factors score` + `Prenatal related inflammatory conditions score` + `Prenatal stress score` + `Prenatal infection sum score` +  `Maternal education` + `Household income` + `B-cells` + `CD4 T-cells` + `CD8 T-cells` + `Granulocytes` + `Monocytes` + `Natural killer cells` + `Nucleated red blood cells` + `Sample plate`

num_colors <- 16
palette <- viridis_pal(option = "A")(num_colors)

mpsz <- dplyr::select(figure_450k_df2, 'mpsz_lig', 'mpsz_wiel')
mpsz <- dplyr::rename(mpsz, 'MPS - Ligthart' = mpsz_lig, 'MPS - Wielscher' = mpsz_wiel)
mpsz_t <- t(mpsz)

array_450k <- fitExtractVarPartModel(mpsz_t, formula, df4)

a <- plotPercentBars(array_450k, col = palette)
b <- plotVarPart(array_450k, col = palette)

# epic
df5 <- dplyr::select(figure_epic_df2, c('IDC', 'Sample_Plate', 'LPFS_score', 'PICS_score', 'MICS_score', 'HsCRPmgL_g1', 'HsCRPmgL_g2', 'CRP_birth','PSS_TOT_score', 'PIS_TOT_score', 'GESTBIR', 'Bcell', 'CD4T', 'CD8T', 'Gran', 'Mono', 'NK', 'nRBC', 'EDUCM', 'AGE_M_v2', 'PARITY', 'SMOKE_ALL', 'INCOME', 'ZR_crp_pgs'))

df6 <- dplyr::rename(df5, 
                     'Sample plate' = Sample_Plate, 
                     'Lifestyle pro-inflammatory factors score' = LPFS_score,
                     'Prenatal related inflammatory conditions score' = 
                       PICS_score,
                     'Inflammatory medical conditions score ' = MICS_score,
                     'CRP trimester 1' = HsCRPmgL_g1,
                     'CRP trimester 2' = HsCRPmgL_g2,
                     'CRP in cord blood' = CRP_birth,
                     'Prenatal stress score' = PSS_TOT_score,
                     'Prenatal infection sum score' = PIS_TOT_score,
                     'Gestational age at birth' = GESTBIR,
                     'B-cells' = Bcell,
                     'CD4 T-cells' = CD4T,
                     'CD8 T-cells' = CD8T,
                     'Granulocytes' = Gran,
                     'Monocytes' = Mono,
                     'Natural killer cells' = NK,
                     'Nucleated red blood cells' = nRBC,
                     'Maternal education' = EDUCM,
                     'Maternal age' = AGE_M_v2,
                     'Parity' = PARITY,
                     'Maternal tobacco use' = SMOKE_ALL,
                     'Household income' = INCOME,
                     'PRS of CRP (child)' = ZR_crp_pgs) 

mpsz_epic <- dplyr::select(figure_epic_df2, 'mpsz_lig', 'mpsz_wiel')
mpsz_epic <- dplyr::rename(mpsz_epic, 'MPS - Ligthart' = mpsz_lig, 'MPS - Wielscher' = mpsz_wiel)
mpsz_epic_t <- t(mpsz_epic)

array_epic <- fitExtractVarPartModel(mpsz_epic_t, formula, df6)

c <- plotPercentBars(array_epic, col = palette)
d <- plotVarPart(array_epic, col = palette)

# combine 450k and epic plots
ggarrange(a,c, ncol = 1, nrow = 2, labels = 'AUTO', common.legend = T, legend = 'bottom')

## Repeat variance partitioning plot (Ligthart only)

# 450k 
df3 <- dplyr::select(figure_450k_df2, c('IDC', 'Sample_Plate', 'LPFS_score', 'PICS_score', 'MICS_score', 'HsCRPmgL_g1', 'HsCRPmgL_g2', 'CRP_birth','PSS_TOT_score', 'PIS_TOT_score', 'GESTBIR', 'Bcell', 'CD4T', 'CD8T', 'Gran', 'Mono', 'NK', 'nRBC', 'EDUCM', 'AGE_M_v2', 'PARITY', 'SMOKE_ALL', 'INCOME', 'ZR_crp_pgs'))

df4 <- dplyr::rename(df3, 
                     'Sample plate' = Sample_Plate, 
                     'Lifestyle pro-inflammatory factors score' = LPFS_score,
                     'Prenatal related inflammatory conditions score' = PICS_score,
                     'Inflammatory medical conditions score ' = MICS_score,
                     'CRP trimester 1' = HsCRPmgL_g1,
                     'CRP trimester 2' = HsCRPmgL_g2,
                     'CRP in cord blood' = CRP_birth,
                     'Prenatal stress score' = PSS_TOT_score,
                     'Prenatal infection sum score' = PIS_TOT_score,
                     'Gestational age at birth' = GESTBIR,
                     'B-cells' = Bcell,
                     'CD4 T-cells' = CD4T,
                     'CD8 T-cells' = CD8T,
                     'Granulocytes' = Gran,
                     'Monocytes' = Mono,
                     'Natural killer cells' = NK,
                     'Nucleated red blood cells' = nRBC,
                     'Maternal education' = EDUCM,
                     'Maternal age' = AGE_M_v2,
                     'Parity' = PARITY,
                     'Maternal tobacco use' = SMOKE_ALL,
                     'Household income' = INCOME,
                     'PRS of CRP (child)' = ZR_crp_pgs) 

formula <- ~ `PRS of CRP (child)` + `Lifestyle pro-inflammatory factors score` + `Prenatal related inflammatory conditions score` + `Prenatal stress score` + `Prenatal infection sum score` +  `Maternal education` + `Household income` + `B-cells` + `CD4 T-cells` + `CD8 T-cells` + `Granulocytes` + `Monocytes` + `Natural killer cells` + `Nucleated red blood cells` + `Sample plate`

num_colors <- 16
palette <- viridis_pal(option = "A")(num_colors)

mpsz <- dplyr::select(figure_450k_df2, 'mpsz_lig')
mpsz <- dplyr::rename(mpsz, 'MPS - Ligthart' = mpsz_lig)
mpsz_t <- t(mpsz)

array_450k <- fitExtractVarPartModel(mpsz_t, formula, df4)

a <- plotPercentBars(array_450k, col = palette)

# epic
df5 <- dplyr::select(figure_epic_df2, c('IDC', 'Sample_Plate', 'LPFS_score', 'PICS_score', 'MICS_score', 'HsCRPmgL_g1', 'HsCRPmgL_g2', 'CRP_birth','PSS_TOT_score', 'PIS_TOT_score', 'GESTBIR', 'Bcell', 'CD4T', 'CD8T', 'Gran', 'Mono', 'NK', 'nRBC', 'EDUCM', 'AGE_M_v2', 'PARITY', 'SMOKE_ALL', 'INCOME', 'ZR_crp_pgs'))

df6 <- dplyr::rename(df5, 
                     'Sample plate' = Sample_Plate, 
                     'Lifestyle pro-inflammatory factors score' = LPFS_score,
                     'Prenatal related inflammatory conditions score' = 
                       PICS_score,
                     'Inflammatory medical conditions score ' = MICS_score,
                     'CRP trimester 1' = HsCRPmgL_g1,
                     'CRP trimester 2' = HsCRPmgL_g2,
                     'CRP in cord blood' = CRP_birth,
                     'Prenatal stress score' = PSS_TOT_score,
                     'Prenatal infection sum score' = PIS_TOT_score,
                     'Gestational age at birth' = GESTBIR,
                     'B-cells' = Bcell,
                     'CD4 T-cells' = CD4T,
                     'CD8 T-cells' = CD8T,
                     'Granulocytes' = Gran,
                     'Monocytes' = Mono,
                     'Natural killer cells' = NK,
                     'Nucleated red blood cells' = nRBC,
                     'Maternal education' = EDUCM,
                     'Maternal age' = AGE_M_v2,
                     'Parity' = PARITY,
                     'Maternal tobacco use' = SMOKE_ALL,
                     'Household income' = INCOME,
                     'PRS of CRP (child)' = ZR_crp_pgs) 

mpsz_epic <- dplyr::select(figure_epic_df2, 'mpsz_lig')
mpsz_epic <- dplyr::rename(mpsz_epic, 'MPS - Ligthart' = mpsz_lig)
mpsz_epic_t <- t(mpsz_epic)

array_epic <- fitExtractVarPartModel(mpsz_epic_t, formula, df6)

c <- plotPercentBars(array_epic, col = palette)

# combine 450k and epic plots
ggarrange(a,c, ncol = 1, nrow = 2, labels = 'AUTO', common.legend = T, legend = 'bottom')

# END OF THIS SCRIPT, GO TO FINAL PART E
