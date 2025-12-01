#########################################################################
######################PROJECT: Neonatal inflammation####################
#########################################################################

#Project: Mapping prenatal predictors and neurobehavioral outcomes of an epigenetic marker of neonatal inflammation - a longitudinal population-based study 

###Parts in this script###
#Of note, in the previous parts of the script we constructed the methylation score, created clinical inflammation scores, created a final dataframe with all relevant data, applied inclusion and exclusion criteria, imputed missingness, created the flowchart and baseline table, created descriptive figures, run the validation of MPS analyses as well as the analyses of aim 1. 

#Now in this script we have the following parts: 
#' Part 0: Loading relevant libraries and data & pivotting data from wide to long
#' Part 1: Assumption checks for models in aim 2 and 3
#' Part 2: Statistical analysis aim 2 [brain development]
#' Part 3: Statistical analysis aim 3 [behavior development]
#' Part 4: Apply FDR correction for all results for aims 1-3
#' Part 5: Create a forest plot for main model results for aims 1-3

#\

###-----------------------------------------------###
#-----------------------PART 0----------------------#
###----------------------------------------------###

#====================================================
## Preparation script 
# Clear environment
rm(list = ls())

# Set seed 
set.seed(2024)

# Load relevant libraries 
libraries <- c('foreign', 'haven', 'tidyverse', 'broom.mixed', 'writexl', 'readxl', 'corrplot', 'stringi', 'miceadds', 'mitools', 'lme4', 'lmerTest', 'rlang', 'lattice', 'ggplot2', 'ggeffects', 'splines', 'reshape2', 'gridExtra', 'grid', 'ggpubr', 'survey', 'survival', 'ggseg', 'fmsb', 'ggcorrplot', 'Hmisc', 'RColorBrewer', 'car', 'DescTools', 'lavaan', 'semTools', 'sjPlot', 'pROC', 'metafor', 'pbmcapply', 'ggseg', 'RColorBrewer', 'Matrix', 'readxl')

invisible(lapply(libraries, require, character.only = T)) 

# Set working directory
setwd('path_to_data')

#====================================================
## Transformations dataset 
# Load final imputed data
imp.data.mids <- readRDS('imp.data.mids.rds')

# General manipulations (convert CBCL age to years if it is in months)
imp.data_long <- complete(imp.data.mids, include = T, action = "long")
imp.data_long$AGE18M_y <- imp.data_long$AGE18M / 12
imp.data_long$age_GR1065_y <- imp.data_long$age_GR1065 / 12
imp.data_long$agechild_GR1075_y <- imp.data_long$agechild_GR1075 / 12
imp.data.mids2 <- as.mids(imp.data_long)

# Transform to long datasets for MRI (smri + dti)
imp.data_mri_combined <- list()

for (i in 0:30) {
  imp.data <- complete(imp.data.mids2, i)
  reshape_imp.data <- reshape(imp.data, idvar = "IDC", varying = list(
    c("age_child_mri_f09", "age_child_mri_f13"), 
    c("Brain_Stem_vol_f09", "Brain_Stem_vol_f13"),
    c("TotalGrayVol_f09", "TotalGrayVol_f13"),
    c("genr_tbv_f09", "genr_tbv_f13"), 
    c("eTIV_f09", "eTIV_f13"),
    c("mean_FA_genr_f09", "mean_FA_genr_f13"),
    c("mean_MD_genr_f09", "mean_MD_genr_f13"),
    c("Cerebellum_Cortex_vol_subcortical_f09", "Cerebellum_Cortex_vol_subcortical_f13"
      ),
    c("Hippocampus_vol_subcortical_f09", "Hippocampus_vol_subcortical_f13"),
    c("Amygdala_vol_subcortical_f09", "Amygdala_vol_subcortical_f13"),
    c("Lateral_Ventricle_vol_subcortical_f09", "Lateral_Ventricle_vol_subcortical_f13"
      ),
    c("CerebralWhiteMatterVol_cortical_f09", "CerebralWhiteMatterVol_cortical_f13")),
    v.names = c("time", 'brain_stem', 'gray_matter', 'total_brain_volume', 
                'intracranial_volume','global_mean_FA','global_mean_MD',
                'cerebellum_cortex', 'hippocampus', 'amygdala', 'lateral_ventricles'
                ,'white_matter'),
    direction = 'long')
  imp.data_mri_combined[[i+1]] <- reshape_imp.data
}

reshape_imp.data_mri <- datalist2mids(imp.data_mri_combined)

# Transform to long datasets for CBCL
imp.data_cbcl_combined <- list()

for (i in 0:30) {
  imp.data <- complete(imp.data.mids2, i)
  reshape_imp.data <- reshape(imp.data, idvar = "IDC", varying = list(
    c("AGE18M_y", "age_GR1065_y", "agechild_GR1075_y", "AgeChild_CBCL9m",
      "AGECHILD_GR1093"), 
    c("cbcl_sum_18m", "cbcl_sum_36m", "cbcl_sum_5", "cbcl_sum_9m", "cbcl_sum_14"),
    c("sum_int_18m", "sum_int_36m", "sum_int_5", "sum_int_9m", "sum_int_14"),
    c("sum_ext_18m", "sum_ext_36m", "sum_ext_5", "sum_ext_9m", "sum_ext_14")),
    v.names = c("time", "CBCL_total", "CBCL_internalizing", "CBCL_externalizing"),
    direction = 'long')
  imp.data_cbcl_combined[[i+1]] <- reshape_imp.data
}

reshape_imp.data_cbcl <- datalist2mids(imp.data_cbcl_combined)

# Check if pivotting went correctly 
reshape_imp.data_cbcl_long <- complete(reshape_imp.data_cbcl, include = T, action = "long")

reshape_imp.data_mri_long <- complete(reshape_imp.data_mri, include = T, action = "long")

# Split up on array: 450k and epic separately 
reshape_imp.data_mri_long_450k <- subset(reshape_imp.data_mri_long, array == '450k')
reshape_imp.data_mri_long_epic <- subset(reshape_imp.data_mri_long, array == 'epic')

reshape_imp.data_cbcl_long_450k <- subset(reshape_imp.data_cbcl_long, array == '450k')
reshape_imp.data_cbcl_long_epic <- subset(reshape_imp.data_cbcl_long, array == 'epic')

# Convert back to mids
long_imp.data_mri_450k <- as.mids(reshape_imp.data_mri_long_450k)
long_imp.data_mri_epic <- as.mids(reshape_imp.data_mri_long_epic)

long_imp.data_cbcl_450k <- as.mids(reshape_imp.data_cbcl_long_450k)
long_imp.data_cbcl_epic <- as.mids(reshape_imp.data_cbcl_long_epic)

saveRDS(long_imp.data_mri_450k, "long_imp.data_mri_450k.rds")
saveRDS(long_imp.data_mri_epic, "long_imp.data_mri_epic.rds")

saveRDS(long_imp.data_cbcl_450k, "long_imp.data_cbcl_450k.rds")
saveRDS(long_imp.data_cbcl_epic, "long_imp.data_cbcl_450k.rds")

#\ 

###-----------------------------------------------###
#-----------------------PART 1----------------------#
###----------------------------------------------###
setwd('path_to_results')

# Specify function for assumption checks 
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

# Specify models (global outcomes for wielscher MPS as exposure for 450k/epic separately)
models <- list(
  lm1 = lmer(scale(total_brain_volume) ~ mpsz_wiel + time + AGE_M_v2 + SMOKE_ALL + EDUCM_3l + INCOME + GENDER + PARITY + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + (1| IDC), data = single_df_mri_450k),
  lm2 = lmer(scale(global_mean_MD) ~ mpsz_wiel + time + AGE_M_v2 + SMOKE_ALL + EDUCM_3l + INCOME + GENDER + PARITY + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + (1| IDC), data = single_df_mri_450k),
  lm3 = lmer(scale(total_brain_volume) ~ mpsz_wiel + time + AGE_M_v2 + SMOKE_ALL + EDUCM_3l + INCOME + GENDER + PARITY + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + (1| IDC), data = single_df_mri_epic),
  lm4 = lmer(scale(global_mean_MD) ~ mpsz_wiel + time + AGE_M_v2 + SMOKE_ALL + EDUCM_3l + INCOME + GENDER + PARITY + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + (1| IDC), data = single_df_mri_epic),
  lm5 = lmer(scale(CBCL_total) ~ mpsz_wiel + time + AGE_M_v2 + SMOKE_ALL + EDUCM_3l + INCOME + GENDER + PARITY + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + (1| IDC), data = single_df_cbcl_450k),
  lm6 = lmer(scale(CBCL_total) ~ mpsz_wiel + time + AGE_M_v2 + SMOKE_ALL + EDUCM_3l + INCOME + GENDER + PARITY + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + (1| IDC), data = single_df_cbcl_epic),
  lm7 = lmer(scale(sqrt(CBCL_total)) ~ mpsz_wiel + time + AGE_M_v2 + SMOKE_ALL + EDUCM_3l + INCOME + GENDER + PARITY + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + (1| IDC), data = single_df_cbcl_450k),
  lm8 = lmer(scale(sqrt(CBCL_total)) ~ mpsz_wiel + time + AGE_M_v2 + SMOKE_ALL + EDUCM_3l + INCOME + GENDER + PARITY + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + (1| IDC), data = single_df_cbcl_epic)
)

# Check assumptions (adding sqrt transformation for cbcl outcomes improves assumptions)
for (lm in names(models)) {
  cat("Checking assumptions for", lm, "\n")
  check_model_assumptions(models[[lm]])
}

#\

###-----------------------------------------------###
#-----------------------PART 2----------------------#
###----------------------------------------------###

#====================================================
### Step 1: associate MPS ligthart/wielscher to sMRI & DTI variables
# Specify exposure eand outcome vars 
mps_vars <- c('mpsz_lig', 'mpsz_wiel')
brain_outcome_vars <- c('brain_stem', 'gray_matter', 'total_brain_volume', 'cerebellum_cortex', 'hippocampus', 'amygdala', 'white_matter', 'lateral_ventricles', 'global_mean_MD', 'global_mean_FA')

# Specify empty dataframes 
results_aim2_mri_main_450k <- data.frame()
results_aim2_mri_interaction_450k <- data.frame()
results_aim2_mri_gestbir_450k <- data.frame()
results_aim2_mri_icv_450k <- data.frame()
results_aim2_mri_nocell_450k <- data.frame()

results_aim2_mri_main_epic <- data.frame()
results_aim2_mri_interaction_epic <- data.frame()
results_aim2_mri_gestbir_epic <- data.frame()
results_aim2_mri_icv_epic <- data.frame()
results_aim2_mri_nocell_epic <- data.frame()

# Apply all models and loop over exposures and outcomes 
for (outcome in brain_outcome_vars) {
  for (exposure in mps_vars) {
    
    # Model formulas
    main <- paste0('scale(',outcome, ') ~', exposure, '+ time + AGE_M_v2 + SMOKE_ALL + EDUCM_3l + INCOME + GENDER + PARITY + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + (1| IDC)')  
    interaction <- paste0('scale(',outcome, ') ~', exposure, '*time + AGE_M_v2 + SMOKE_ALL + EDUCM_3l + INCOME + GENDER + PARITY + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + (1| IDC)')  
    gestbir_mod <- paste0('scale(',outcome, ') ~', exposure, '*scale(GESTBIR)*time + AGE_M_v2 + SMOKE_ALL + EDUCM_3l + INCOME + GENDER + PARITY + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + (1| IDC)')  
    main_icv <- paste0('scale(',outcome, ') ~', exposure, '+ time + AGE_M_v2 + SMOKE_ALL + EDUCM_3l + INCOME + GENDER + PARITY + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + scale(intracranial_volume) + (1| IDC)') 
    main_nocell <- paste0('scale(',outcome, ')~', exposure, '+ time + AGE_M_v2 + SMOKE_ALL + EDUCM_3l + INCOME + GENDER + PARITY + Sample_Plate + (1| IDC)')  
    
    # Calculating beta, confidence interval, and p-value for all models
    main_model_output_450k <- summary(pool(with(long_imp.data_mri_450k, lmer(as.formula(main)))), conf.int = TRUE)[2, c(2, 3, 6, 7, 8)]
    interaction_model_output_450k <- summary(pool(with(long_imp.data_mri_450k, lmer(as.formula(interaction)))), conf.int = TRUE)[46, c(2,3, 6, 7, 8)]
    gestbir_mod_output_450k <- summary(pool(with(long_imp.data_mri_450k, lmer(as.formula(gestbir_mod)))), conf.int = TRUE)[50, c(2,3, 6, 7, 8)]
    icv_adjust_output_450k <- summary(pool(with(long_imp.data_mri_450k, lmer(as.formula(main_icv)))), conf.int = TRUE)[2, c(2,3, 6, 7, 8)]
    no_cell_output_450k <- summary(pool(with(long_imp.data_mri_450k, lmer(as.formula(main_nocell)))), conf.int = TRUE)[2, c(2, 3,6, 7, 8)]
    
    main_model_output_epic <- summary(pool(with(long_imp.data_mri_epic, lmer(as.formula(main)))), conf.int = TRUE)[2, c(2, 3, 6, 7, 8)]
    interaction_model_output_epic <- summary(pool(with(long_imp.data_mri_epic, lmer(as.formula(interaction)))), conf.int = TRUE)[30, c(2,3, 6, 7, 8)]
    gestbir_mod_output_epic <- summary(pool(with(long_imp.data_mri_epic, lmer(as.formula(gestbir_mod)))), conf.int = TRUE)[34, c(2,3, 6, 7, 8)]
    icv_adjust_output_epic <- summary(pool(with(long_imp.data_mri_epic, lmer(as.formula(main_icv)))), conf.int = TRUE)[2, c(2,3, 6, 7, 8)]
    no_cell_output_epic <- summary(pool(with(long_imp.data_mri_epic, lmer(as.formula(main_nocell)))), conf.int = TRUE)[2, c(2, 3,6, 7, 8)]
    
    # Adding outcome and exposure information
    main_model_output_450k$outcome <- outcome
    main_model_output_450k$exposure <- exposure
    interaction_model_output_450k$outcome <- outcome
    interaction_model_output_450k$exposure <- exposure
    gestbir_mod_output_450k$outcome <- outcome
    gestbir_mod_output_450k$exposure <- exposure
    icv_adjust_output_450k$outcome <- outcome
    icv_adjust_output_450k$exposure <- exposure
    no_cell_output_450k$outcome <- outcome
    no_cell_output_450k$exposure <- exposure
    
    main_model_output_epic$outcome <- outcome
    main_model_output_epic$exposure <- exposure
    interaction_model_output_epic$outcome <- outcome
    interaction_model_output_epic$exposure <- exposure
    gestbir_mod_output_epic$outcome <- outcome
    gestbir_mod_output_epic$exposure <- exposure
    icv_adjust_output_epic$outcome <- outcome
    icv_adjust_output_epic$exposure <- exposure
    no_cell_output_epic$outcome <- outcome
    no_cell_output_epic$exposure <- exposure
    
    # Adding results to data frames
    results_aim2_mri_main_450k <- rbind(results_aim2_mri_main_450k, 
                                         main_model_output_450k)
    results_aim2_mri_interaction_450k <- rbind(results_aim2_mri_interaction_450k, 
                                                interaction_model_output_450k)
    results_aim2_mri_gestbir_450k <- rbind(results_aim2_mri_gestbir_450k, 
                                            gestbir_mod_output_450k)
    results_aim2_mri_icv_450k <- rbind(results_aim2_mri_icv_450k, 
                                        icv_adjust_output_450k)
    results_aim2_mri_nocell_450k <- rbind(results_aim2_mri_nocell_450k, 
                                           no_cell_output_450k)
    
    results_aim2_mri_main_epic <- rbind(results_aim2_mri_main_epic, 
                                         main_model_output_epic)
    results_aim2_mri_interaction_epic <- rbind(results_aim2_mri_interaction_epic, 
                                                interaction_model_output_epic)
    results_aim2_mri_gestbir_epic <- rbind(results_aim2_mri_gestbir_epic, 
                                            gestbir_mod_output_epic)
    results_aim2_mri_icv_epic <- rbind(results_aim2_mri_icv_epic, 
                                        icv_adjust_output_epic)
    results_aim2_mri_nocell_epic <- rbind(results_aim2_mri_nocell_epic, 
                                           no_cell_output_epic)
  }
}

#\

#====================================================
### Step 2: Meta-analyses 
# Specify function
meta_analyze <- function(data) {
  rma(yi = data$estimate, sei = data$std.error, method = "FE")
}

# Create input files 
combined_results_main_model_aim2_ligthart <- bind_rows(
  subset(results_aim2_mri_main_450k, exposure == 'mpsz_lig'),
  subset(results_aim2_mri_main_epic, exposure == 'mpsz_lig'),
)
combined_results_main_model_aim2_wielscher <- bind_rows(
  subset(results_aim2_mri_main_450k, exposure == 'mpsz_wiel'),
  subset(results_aim2_mri_main_epic, exposure == 'mpsz_wiel'),
)

combined_results_interaction_model_aim2_ligthart <- bind_rows(
  subset(results_aim2_mri_interaction_450k, exposure == 'mpsz_lig'),
  subset(results_aim2_mri_interaction_epic, exposure == 'mpsz_lig'),
)
combined_results_interaction_model_aim2_wielscher <- bind_rows(
  subset(results_aim2_mri_interaction_450k, exposure == 'mpsz_wiel'),
  subset(results_aim2_mri_interaction_epic, exposure == 'mpsz_wiel'),
)

combined_results_gestbir_model_aim2_ligthart <- bind_rows(
  subset(results_aim2_mri_gestbir_450k, exposure == 'mpsz_lig'),
  subset(results_aim2_mri_gestbir_epic, exposure == 'mpsz_lig'),
)
combined_results_gestbir_model_aim2_wielscher <- bind_rows(
  subset(results_aim2_mri_gestbir_450k, exposure == 'mpsz_wiel'),
  subset(results_aim2_mri_gestbir_epic, exposure == 'mpsz_wiel'),
)

combined_results_icv_model_aim2_ligthart <- bind_rows(
  subset(results_aim2_mri_icv_450k, exposure == 'mpsz_lig'),
  subset(results_aim2_mri_icv_epic, exposure == 'mpsz_lig'),
)
combined_results_icv_model_aim2_wielscher <- bind_rows(
  subset(results_aim2_mri_icv_450k, exposure == 'mpsz_wiel'),
  subset(results_aim2_mri_icv_epic, exposure == 'mpsz_wiel'),
)

combined_results_nocell_model_aim2_ligthart <- bind_rows(
  subset(results_aim2_mri_nocell_450k, exposure == 'mpsz_lig'),
  subset(results_aim2_mri_nocell_epic, exposure == 'mpsz_lig'),
)
combined_results_nocell_model_aim2_wielscher <- bind_rows(
  subset(results_aim2_mri_nocell_450k, exposure == 'mpsz_wiel'),
  subset(results_aim2_mri_nocell_epic, exposure == 'mpsz_wiel'),
)

# Specify unique outcome
unique_outcomes <- unique(combined_results_main_model_aim2_ligthart$outcome)

# Apply meta-analysis function to each unique outcome
meta_analyze_outcomes <- function(unique_outcomes, combined_results_model) {
  lapply(unique_outcomes, function(outcome) {
    subset_data <- combined_results_model[combined_results_model$outcome == outcome, ]
    meta_analyze(subset_data)
  })
}

meta_results_main_model_ligthart <- meta_analyze_outcomes(unique_outcomes, combined_results_main_model_aim2_ligthart)
meta_results_interaction_model_ligthart <- meta_analyze_outcomes(unique_outcomes, combined_results_interaction_model_aim2_ligthart)
meta_results_gestbir_model_ligthart <- meta_analyze_outcomes(unique_outcomes, combined_results_gestbir_model_aim2_ligthart)
meta_results_icv_model_ligthart <- meta_analyze_outcomes(unique_outcomes, combined_results_icv_model_aim2_ligthart)
meta_results_nocell_model_ligthart <- meta_analyze_outcomes(unique_outcomes, combined_results_nocell_model_aim2_ligthart)

meta_results_main_model_wielscher <- meta_analyze_outcomes(unique_outcomes, combined_results_main_model_aim2_wielscher)
meta_results_interaction_model_wielscher <- meta_analyze_outcomes(unique_outcomes, combined_results_interaction_model_aim2_wielscher)
meta_results_gestbir_model_wielscher <- meta_analyze_outcomes(unique_outcomes, combined_results_gestbir_model_aim2_wielscher)
meta_results_icv_model_wielscher <- meta_analyze_outcomes(unique_outcomes, combined_results_icv_model_aim2_wielscher)
meta_results_nocell_model_wielscher <- meta_analyze_outcomes(unique_outcomes, combined_results_nocell_model_aim2_wielscher)

#====================================================
### Step 3: Plots for suggestive findings
## Anatomy plot (interaction models)
suggestive_sign_brain_findings <- data.frame(
  region = c('amygdala', 'cerebellum cortex'),
  interaction_beta = c(-0.007, 0.016)
)

anat_plot <- suggestive_sign_brain_findings %>%
  ggplot() +
  geom_brain(atlas = aseg,
             aes(fill = interaction_beta)) +
  theme_minimal() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank() ) + 
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 8, angle =25)) +
  scale_fill_gradient(low = 'burlywood1', high = 'coral3') + 
  guides(fill = guide_colorbar(barheight = 1)) + 
  labs(fill = 'Standardized interaction beta coefficient')

## Longitudinal plots
# main - ligthart - cerebellum
# interaction - ligthart - cerebellum/amygdala/global mean FA; wielscher - amygdala/global mean FA*

# Select one of the imputed dataset only for visualization purpose 
figure_df_450k <- complete(long_imp.data_mri_450k, 30)
figure_df_epic <- complete(long_imp.data_mri_epic, 30)
figure_df_combined <- rbind(figure_df_450k, figure_df_epic)

# Specify models 450k
lm_interaction_ligthart1_450k <- lmer(scale(cerebellum_cortex) ~ mpsz_lig*time + AGE_M_v2 + SMOKE_ALL + EDUCM_3l + INCOME + GENDER + PARITY + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + (1| IDC), data = figure_df_450k)
lm_interaction_ligthart2_450k <- lmer(scale(amygdala) ~ mpsz_lig*time + AGE_M_v2 + SMOKE_ALL + EDUCM_3l + INCOME + GENDER + PARITY + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + (1| IDC), data = figure_df_450k)
lm_interaction_ligthart3_450k  <- lmer(scale(global_mean_FA) ~ mpsz_lig*time + AGE_M_v2 + SMOKE_ALL + EDUCM_3l + INCOME + GENDER + PARITY + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + (1| IDC), data = figure_df_450k) 

lm_interaction_wielscher0_450k<- lmer(scale(cerebellum_cortex) ~ mpsz_wiel*time + AGE_M_v2 + SMOKE_ALL + EDUCM_3l + INCOME + GENDER + PARITY + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + (1| IDC), data = figure_df_450k)
lm_interaction_wielscher1_450k<- lmer(scale(amygdala) ~ mpsz_wiel*time + AGE_M_v2 + SMOKE_ALL + EDUCM_3l + INCOME + GENDER + PARITY + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + (1| IDC), data = figure_df_450k)
lm_interaction_wielscher2_450k <- lmer(scale(global_mean_FA) ~ mpsz_wiel*time + AGE_M_v2 + SMOKE_ALL + EDUCM_3l + INCOME + GENDER + PARITY + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + (1| IDC), data = figure_df_450k)

# Specify models epic
lm_interaction_ligthart1_epic <- lmer(scale(cerebellum_cortex) ~ mpsz_lig*time + AGE_M_v2 + SMOKE_ALL + EDUCM_3l + INCOME + GENDER + PARITY + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + (1| IDC), data = figure_df_epic)
lm_interaction_ligthart2_epic <- lmer(scale(amygdala) ~ mpsz_lig*time + AGE_M_v2 + SMOKE_ALL + EDUCM_3l + INCOME + GENDER + PARITY + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + (1| IDC), data = figure_df_epic)
lm_interaction_ligthart3_epic  <- lmer(scale(global_mean_FA) ~ mpsz_lig*time + AGE_M_v2 + SMOKE_ALL + EDUCM_3l + INCOME + GENDER + PARITY + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + (1| IDC), data = figure_df_epic) 

lm_interaction_wielscher0_epic<- lmer(scale(cerebellum_cortex) ~ mpsz_wiel*time + AGE_M_v2 + SMOKE_ALL + EDUCM_3l + INCOME + GENDER + PARITY + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + (1| IDC), data = figure_df_epic)
lm_interaction_wielscher1_epic<- lmer(scale(amygdala) ~ mpsz_wiel*time + AGE_M_v2 + SMOKE_ALL + EDUCM_3l + INCOME + GENDER + PARITY + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + (1| IDC), data = figure_df_epic)
lm_interaction_wielscher2_epic <- lmer(scale(global_mean_FA) ~ mpsz_wiel*time + AGE_M_v2 + SMOKE_ALL + EDUCM_3l + INCOME + GENDER + PARITY + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + (1| IDC), data = figure_df_epic)

# Specify models combined
lm_interaction_ligthart1_combined <- lmer(scale(cerebellum_cortex) ~ mpsz_lig*time + AGE_M_v2 + SMOKE_ALL + EDUCM_3l + INCOME + GENDER + PARITY + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + (1| IDC), data = figure_df_combined)
lm_interaction_ligthart2_combined <- lmer(scale(amygdala) ~ mpsz_lig*time + AGE_M_v2 + SMOKE_ALL + EDUCM_3l + INCOME + GENDER + PARITY + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + (1| IDC), data = figure_df_combined)
lm_interaction_ligthart3_combined  <- lmer(scale(global_mean_FA) ~ mpsz_lig*time + AGE_M_v2 + SMOKE_ALL + EDUCM_3l + INCOME + GENDER + PARITY + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + (1| IDC), data = figure_df_combined) 

lm_interaction_wielscher0_combined<- lmer(scale(cerebellum_cortex) ~ mpsz_wiel*time + AGE_M_v2 + SMOKE_ALL + EDUCM_3l + INCOME + GENDER + PARITY + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + (1| IDC), data = figure_df_combined)
lm_interaction_wielscher1_combined<- lmer(scale(amygdala) ~ mpsz_wiel*time + AGE_M_v2 + SMOKE_ALL + EDUCM_3l + INCOME + GENDER + PARITY + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + (1| IDC), data = figure_df_combined)
lm_interaction_wielscher2_combined <- lmer(scale(global_mean_FA) ~ mpsz_wiel*time + AGE_M_v2 + SMOKE_ALL + EDUCM_3l + INCOME + GENDER + PARITY + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + (1| IDC), data = figure_df_combined)

# Create plots 450k
# Plotting 1st quartile and 3rd quartile for interactions
plot_interaction_ligthart1_450k <- plot_model(lm_interaction_ligthart1_450k, type = "pred", terms = c("time", "mpsz_lig[-0.66579, 0.66164]")) + ggplot2::xlab("Child Age (years)") +
  ggplot2::ylab("Cerebellum cortex") +
  ggplot2::ggtitle("450k array") + 
  ggplot2::labs(color = 'Standardized MPS - Ligthart') +
  scale_color_manual(values = c("-0.66579" = "deeppink3", "0.66164" = "deepskyblue4"))+ theme_bw() 

plot_interaction_ligthart2_450k <- plot_model(lm_interaction_ligthart2_450k, type = "pred", terms = c("time", "mpsz_lig[-0.66579, 0.66164]")) + ggplot2::xlab("Child Age (years)") +
  ggplot2::ylab("Amygdala") +
  ggplot2::ggtitle("450k array") + 
  ggplot2::labs(color = 'Standardized MPS - Ligthart') +
  scale_color_manual(values = c("-0.66579" = "deeppink3", "0.66164" = "deepskyblue4"))+ theme_bw() 

plot_interaction_ligthart3_450k <- plot_model(lm_interaction_ligthart3_450k, type = "pred", terms = c("time", "mpsz_lig[-0.66579, 0.66164]")) + ggplot2::xlab("Child Age (years)") +
  ggplot2::ylab("Global mean FA") +
  ggplot2::ggtitle("450k array") + 
  ggplot2::labs(color = 'Standardized MPS - Ligthart') +
  scale_color_manual(values = c("-0.66579" = "deeppink3", "0.66164" = "deepskyblue4"))+ theme_bw() 

plot_interaction_wielscher0_450k <- plot_model(lm_interaction_wielscher0_450k, type = "pred", terms = c("time", "mpsz_wiel[-0.71993, 0.66103]")) + ggplot2::xlab("Child Age (years)") +
  ggplot2::ylab("Cerebellum cortex") +
  ggplot2::ggtitle("450k array") + 
  ggplot2::labs(color = 'Standardized MPS - Wielscher') +
  scale_color_manual(values = c("-0.71993" = "deeppink3", "0.66103" = "deepskyblue4"))+ theme_bw() 

plot_interaction_wielscher1_450k <- plot_model(lm_interaction_wielscher1_450k, type = "pred", terms = c("time", "mpsz_wiel[-0.71993, 0.66103]")) + ggplot2::xlab("Child Age (years)") +
  ggplot2::ylab("Amygdala") +
  ggplot2::ggtitle("450k array") + 
  ggplot2::labs(color = 'Standardized MPS - Wielscher') +
  scale_color_manual(values = c("-0.71993" = "deeppink3", "0.66103" = "deepskyblue4"))+ theme_bw() 

plot_interaction_wielscher2_450k <- plot_model(lm_interaction_wielscher2_450k, type = "pred", terms = c("time", "mpsz_wiel[-0.71993, 0.66103]")) + ggplot2::xlab("Child Age (years)") +
  ggplot2::ylab("Global mean FA") +
  ggplot2::ggtitle("450k array") + 
  ggplot2::labs(color = 'Standardized MPS - Wielscher') +
  scale_color_manual(values = c("-0.71993" = "deeppink3", "0.66103" = "deepskyblue4"))+ theme_bw() 

# Create plots epic
plot_interaction_ligthart1_epic <- plot_model(lm_interaction_ligthart1_epic, type = "pred", terms = c("time", "mpsz_lig[-0.608189, 0.662727]")) + ggplot2::xlab("Child Age (years)") +
  ggplot2::ylab("Cerebellum cortex") +
  ggplot2::ggtitle("epic array") + 
  ggplot2::labs(color = 'Standardized MPS - Ligthart') +
  scale_color_manual(values = c("-0.608189" = "deeppink3", "0.662727" = "deepskyblue4"))+ theme_bw() 

plot_interaction_ligthart2_epic <- plot_model(lm_interaction_ligthart2_epic, type = "pred", terms = c("time", "mpsz_lig[-0.608189, 0.662727]")) + ggplot2::xlab("Child Age (years)") +
  ggplot2::ylab("Amygdala") +
  ggplot2::ggtitle("epic array") + 
  ggplot2::labs(color = 'Standardized MPS - Ligthart') +
  scale_color_manual(values = c("-0.608189" = "deeppink3", "0.662727" = "deepskyblue4"))+ theme_bw() 

plot_interaction_ligthart3_epic <- plot_model(lm_interaction_ligthart3_epic, type = "pred", terms = c("time", "mpsz_lig[-0.608189, 0.662727]")) + ggplot2::xlab("Child Age (years)") +
  ggplot2::ylab("Global mean FA") +
  ggplot2::ggtitle("epic array") + 
  ggplot2::labs(color = 'Standardized MPS - Ligthart') +
  scale_color_manual(values = c("-0.608189" = "deeppink3", "0.66164" = "deepskyblue4"))+ theme_bw() 

plot_interaction_wielscher0_epic <- plot_model(lm_interaction_wielscher0_epic, type = "pred", terms = c("time", "mpsz_wiel[-0.59280, 0.69538]")) + ggplot2::xlab("Child Age (years)") +
  ggplot2::ylab("Cerebellum cortex") +
  ggplot2::ggtitle("450k array") + 
  ggplot2::labs(color = 'Standardized MPS - Wielscher') +
  scale_color_manual(values = c("-0.5928" = "deeppink3", "0.69538" = "deepskyblue4"))+ theme_bw() 

plot_interaction_wielscher1_epic <- plot_model(lm_interaction_wielscher1_epic, type = "pred", terms = c("time", "mpsz_wiel[-0.59280, 0.69538]")) + ggplot2::xlab("Child Age (years)") +
  ggplot2::ylab("Amygdala") +
  ggplot2::ggtitle("epic array") + 
  ggplot2::labs(color = 'Standardized MPS - Wielscher') +
  scale_color_manual(values = c("-0.5928" = "deeppink3", "0.69538" = "deepskyblue4"))+ theme_bw() 

plot_interaction_wielscher2_epic <- plot_model(lm_interaction_wielscher2_epic, type = "pred", terms = c("time", "mpsz_wiel[-0.59280, 0.69538]")) + ggplot2::xlab("Child Age (years)") +
  ggplot2::ylab("Global mean FA") +
  ggplot2::ggtitle("epic array") + 
  ggplot2::labs(color = 'Standardized MPS - Wielscher') +
  scale_color_manual(values = c("-0.5928" = "deeppink3", "0.69538" = "deepskyblue4"))+ theme_bw() 

# Combine plots for assays specific
ggarrange(plot_interaction_ligthart1_450k, plot_interaction_ligthart2_450k, 
          plot_interaction_ligthart3_450k, plot_interaction_wielscher0_450k, 
          plot_interaction_wielscher1_450k, plot_interaction_wielscher2_450k,
          plot_interaction_ligthart1_epic, plot_interaction_ligthart2_epic, 
          plot_interaction_ligthart3_epic, plot_interaction_wielscher0_450k
          ,plot_interaction_wielscher1_epic, plot_interaction_wielscher2_epic
          , labels = 'AUTO', ncol = 3, nrow =5)

# Create plots combined assay
plot_interaction_ligthart1_combined <- plot_model(lm_interaction_ligthart1_combined, type = "pred", terms = c("time", "mpsz_lig[-0.650262, 0.661943]")) + ggplot2::xlab("Child Age (years)") +
  ggplot2::ylab("Cerebellum cortex volume*") +
  ggplot2::labs(color = 'Standardized MPS - Ligthart') +
  scale_color_manual(values = c("-0.650262" = "burlywood1", "0.661943" = "coral3"))+ theme_bw() + ggplot2::ggtitle("")

plot_interaction_ligthart2_combined <- plot_model(lm_interaction_ligthart2_combined, type = "pred", terms = c("time", "mpsz_lig[-0.650262, 0.661943]")) + ggplot2::xlab("Child Age (years)") +
  ggplot2::ylab("Amygdala volume*") +
  ggplot2::labs(color = 'Standardized MPS - Ligthart') +
  scale_color_manual(values = c("-0.650262" = "burlywood1", "0.661943" = "coral3"))+ theme_bw() + ggplot2::ggtitle("")

plot_interaction_ligthart3_combined <- plot_model(lm_interaction_ligthart3_combined, type = "pred", terms = c("time", "mpsz_lig[-0.650262, 0.661943]")) + ggplot2::xlab("Child Age (years)") +
  ggplot2::ylab("Global mean FA*") +
  ggplot2::labs(color = 'Standardized MPS - Ligthart') +
  scale_color_manual(values = c("-0.650262" = "burlywood1", "0.661943" = "coral3"))+ theme_bw() + ggplot2::ggtitle("")

plot_interaction_wielscher0_combined <- plot_model(lm_interaction_wielscher0_combined, type = "pred", terms = c("time", "mpsz_wiel[-0.671707, 0.682004]")) + ggplot2::xlab("Child Age (years)") +
  ggplot2::ylab("Cerebellum cortex volume") +
  ggplot2::labs(color = 'Standardized MPS - Wielscher') +
  scale_color_manual(values = c("-0.671707" = "burlywood1", "0.682004" = "coral3"))+ theme_bw() + ggplot2::ggtitle("")

plot_interaction_wielscher1_combined <- plot_model(lm_interaction_wielscher1_combined, type = "pred", terms = c("time", "mpsz_wiel[-0.671707, 0.682004]")) + ggplot2::xlab("Child Age (years)") +
  ggplot2::ylab("Amygdala volume*") +
  ggplot2::labs(color = 'Standardized MPS - Wielscher') +
  scale_color_manual(values = c("-0.671707" = "burlywood1", "0.682004" = "coral3"))+ theme_bw() + ggplot2::ggtitle("")

plot_interaction_wielscher2_combined <- plot_model(lm_interaction_wielscher2_combined, type = "pred", terms = c("time", "mpsz_wiel[-0.671707, 0.682004]")) + ggplot2::xlab("Child Age (years)") +
  ggplot2::ylab("Global mean FA**") +
  ggplot2::labs(color = 'Standardized MPS - Wielscher') +
  scale_color_manual(values = c("-0.671707" = "burlywood1", "0.682004" = "coral3"))+ theme_bw() + ggplot2::ggtitle("")

# Combine plots for combined df
ligthart <- ggarrange(plot_interaction_ligthart1_combined, plot_interaction_ligthart2_combined, plot_interaction_ligthart3_combined, common.legend = T, legend = 'bottom', ncol = 3)

wielscher <- ggarrange(plot_interaction_wielscher0_combined, plot_interaction_wielscher1_combined, plot_interaction_wielscher2_combined, common.legend = T, legend = 'bottom', ncol = 3)

ggarrange(anat_plot, ligthart,wielscher, nrow =3, labels = 'AUTO', align = 'v')
ggarrange(anat_plot, ligthart, nrow = 2, labels = 'AUTO', align = 'v')

#\

###-----------------------------------------------###
#-----------------------PART 3----------------------#
###----------------------------------------------###

#====================================================
### Step 1: associate MPS ligthart/wielscher to CBCL variables
cbcl_outcome_vars <- c('CBCL_total', 'CBCL_internalizing', 'CBCL_externalizing')

# Specify empty dataframes 
results_aim3_cbcl_main_450k <- data.frame()
results_aim3_cbcl_interaction_450k <- data.frame()
results_aim3_cbcl_gestbir_450k <- data.frame()
results_aim3_cbcl_icv_450k <- data.frame()
results_aim3_cbcl_nocell_450k <- data.frame()

results_aim3_cbcl_main_epic <- data.frame()
results_aim3_cbcl_interaction_epic <- data.frame()
results_aim3_cbcl_gestbir_epic <- data.frame()
results_aim3_cbcl_icv_epic <- data.frame()
results_aim3_cbcl_nocell_epic <- data.frame()

# Apply all models and loop over exposures and outcomes 
for (outcome in cbcl_outcome_vars) {
  for (exposure in mps_vars) {
    
    # Model formulas
    main <- paste0('scale(sqrt(',outcome, ')) ~', exposure, '+ time + AGE_M_v2 + SMOKE_ALL + EDUCM_3l + INCOME + GENDER + PARITY + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + (1| IDC)')  
    interaction <- paste0('scale(sqrt(',outcome, ')) ~', exposure, '*time + AGE_M_v2 + SMOKE_ALL + EDUCM_3l + INCOME + GENDER + PARITY + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + (1| IDC)')  
    gestbir_mod <- paste0('scale(sqrt(',outcome, ')) ~', exposure, '*scale(GESTBIR)*time + AGE_M_v2 + SMOKE_ALL + EDUCM_3l + INCOME + GENDER + PARITY + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK + nRBC + (1| IDC)')  
    main_nocell <- paste0('scale(sqrt(',outcome, '))~', exposure, '+ time + AGE_M_v2 + SMOKE_ALL + EDUCM_3l + INCOME + GENDER + PARITY + Sample_Plate + (1| IDC)')  
    
    # Calculating beta, confidence interval, and p-value for all models
    main_model_output_450k <- summary(pool(with(long_imp.data_cbcl_450k, lmer(as.formula(main)))), conf.int = TRUE)[2, c(2, 3, 6, 7, 8)]
    interaction_model_output_450k <- summary(pool(with(long_imp.data_cbcl_450k, lmer(as.formula(interaction)))), conf.int = TRUE)[46, c(2,3, 6, 7, 8)]
    gestbir_mod_output_450k <- summary(pool(with(long_imp.data_cbcl_450k, lmer(as.formula(gestbir_mod)))), conf.int = TRUE)[50, c(2,3, 6, 7, 8)]
    no_cell_output_450k <- summary(pool(with(long_imp.data_cbcl_450k, lmer(as.formula(main_nocell)))), conf.int = TRUE)[2, c(2, 3,6, 7, 8)]
    
    main_model_output_epic <- summary(pool(with(long_imp.data_cbcl_epic, lmer(as.formula(main)))), conf.int = TRUE)[2, c(2, 3, 6, 7, 8)]
    interaction_model_output_epic <- summary(pool(with(long_imp.data_cbcl_epic, lmer(as.formula(interaction)))), conf.int = TRUE)[30, c(2,3, 6, 7, 8)]
    gestbir_mod_output_epic <- summary(pool(with(long_imp.data_cbcl_epic, lmer(as.formula(gestbir_mod)))), conf.int = TRUE)[34, c(2,3, 6, 7, 8)]
    no_cell_output_epic <- summary(pool(with(long_imp.data_cbcl_epic, lmer(as.formula(main_nocell)))), conf.int = TRUE)[2, c(2, 3,6, 7, 8)]
    
    # Adding outcome and exposure information
    main_model_output_450k$outcome <- outcome
    main_model_output_450k$exposure <- exposure
    interaction_model_output_450k$outcome <- outcome
    interaction_model_output_450k$exposure <- exposure
    gestbir_mod_output_450k$outcome <- outcome
    gestbir_mod_output_450k$exposure <- exposure
    no_cell_output_450k$outcome <- outcome
    no_cell_output_450k$exposure <- exposure
    
    main_model_output_epic$outcome <- outcome
    main_model_output_epic$exposure <- exposure
    interaction_model_output_epic$outcome <- outcome
    interaction_model_output_epic$exposure <- exposure
    gestbir_mod_output_epic$outcome <- outcome
    gestbir_mod_output_epic$exposure <- exposure
    no_cell_output_epic$outcome <- outcome
    no_cell_output_epic$exposure <- exposure
    
    # Adding results to data frames
    results_aim3_cbcl_main_450k <- rbind(results_aim3_cbcl_main_450k, 
                                        main_model_output_450k)
    results_aim3_cbcl_interaction_450k <- rbind(results_aim3_cbcl_interaction_450k, 
                                               interaction_model_output_450k)
    results_aim3_cbcl_gestbir_450k <- rbind(results_aim3_cbcl_gestbir_450k, 
                                           gestbir_mod_output_450k)
    results_aim3_cbcl_nocell_450k <- rbind(results_aim3_cbcl_nocell_450k, 
                                          no_cell_output_450k)
    
    results_aim3_cbcl_main_epic <- rbind(results_aim3_cbcl_main_epic, 
                                        main_model_output_epic)
    results_aim3_cbcl_interaction_epic <- rbind(results_aim3_cbcl_interaction_epic, 
                                               interaction_model_output_epic)
    results_aim3_cbcl_gestbir_epic <- rbind(results_aim3_cbcl_gestbir_epic, 
                                           gestbir_mod_output_epic)
    results_aim3_cbcl_nocell_epic <- rbind(results_aim3_cbcl_nocell_epic, 
                                          no_cell_output_epic)
  }
}

#\

#====================================================
### Step 2: Meta-analyses 

# Create input files 
combined_results_main_model_aim3_ligthart <- bind_rows(
  subset(results_aim3_cbcl_main_450k, exposure == 'mpsz_lig'),
  subset(results_aim3_cbcl_main_epic, exposure == 'mpsz_lig'),
)
combined_results_main_model_aim3_wielscher <- bind_rows(
  subset(results_aim3_cbcl_main_450k, exposure == 'mpsz_wiel'),
  subset(results_aim3_cbcl_main_epic, exposure == 'mpsz_wiel'),
)

combined_results_interaction_model_aim3_ligthart <- bind_rows(
  subset(results_aim3_cbcl_interaction_450k, exposure == 'mpsz_lig'),
  subset(results_aim3_cbcl_interaction_epic, exposure == 'mpsz_lig'),
)
combined_results_interaction_model_aim3_wielscher <- bind_rows(
  subset(results_aim3_cbcl_interaction_450k, exposure == 'mpsz_wiel'),
  subset(results_aim3_cbcl_interaction_epic, exposure == 'mpsz_wiel'),
)

combined_results_gestbir_model_aim3_ligthart <- bind_rows(
  subset(results_aim3_cbcl_gestbir_450k, exposure == 'mpsz_lig'),
  subset(results_aim3_cbcl_gestbir_epic, exposure == 'mpsz_lig'),
)
combined_results_gestbir_model_aim3_wielscher <- bind_rows(
  subset(results_aim3_cbcl_gestbir_450k, exposure == 'mpsz_wiel'),
  subset(results_aim3_cbcl_gestbir_epic, exposure == 'mpsz_wiel'),
)

combined_results_nocell_model_aim3_ligthart <- bind_rows(
  subset(results_aim3_cbcl_nocell_450k, exposure == 'mpsz_lig'),
  subset(results_aim3_cbcl_nocell_epic, exposure == 'mpsz_lig'),
)
combined_results_nocell_model_aim3_wielscher <- bind_rows(
  subset(results_aim3_cbcl_nocell_450k, exposure == 'mpsz_wiel'),
  subset(results_aim3_cbcl_nocell_epic, exposure == 'mpsz_wiel'),
)

# Specify unique outcome
unique_outcomes_aim3 <- unique(combined_results_main_model_aim3_ligthart$outcome)

# Apply meta-analysis function to each unique outcome
meta_results_main_model_ligthart <- meta_analyze_outcomes(unique_outcomes_aim3, combined_results_main_model_aim3_ligthart)
meta_results_interaction_model_ligthart <- meta_analyze_outcomes(unique_outcomes_aim3, combined_results_interaction_model_aim3_ligthart)
meta_results_gestbir_model_ligthart <- meta_analyze_outcomes(unique_outcomes_aim3, combined_results_gestbir_model_aim3_ligthart)
meta_results_nocell_model_ligthart <- meta_analyze_outcomes(unique_outcomes_aim3, combined_results_nocell_model_aim3_ligthart)

meta_results_main_model_wielscher <- meta_analyze_outcomes(unique_outcomes_aim3, combined_results_main_model_aim3_wielscher)
meta_results_interaction_model_wielscher <- meta_analyze_outcomes(unique_outcomes_aim3, combined_results_interaction_model_aim3_wielscher)
meta_results_gestbir_model_wielscher <- meta_analyze_outcomes(unique_outcomes_aim3, combined_results_gestbir_model_aim3_wielscher)
meta_results_nocell_model_wielscher <- meta_analyze_outcomes(unique_outcomes_aim3, combined_results_nocell_model_aim3_wielscher)

#\

###-----------------------------------------------###
#-----------------------PART 4----------------------#
###----------------------------------------------###
# function for fdr correction 
add_fdr_pvalue <- function(data_frame, pvalue_column_index) {
  pval <- unlist(data_frame[, pvalue_column_index])
  data_frame$fdr_pvalue <- p.adjust(pval, method = 'fdr')
  return(data_frame)
}

# load excel files for aims 1-3
MA_results_univariate_model <- as.data.frame(read_excel("MA_results_univariate_model.xlsx"))
MA_results_batch_model <- as.data.frame(read_excel("MA_results_batch_model.xlsx"))
MA_results_batch_cell_model <- as.data.frame(read_excel("MA_results_batch_cell_model.xlsx"))

MA_results_full_model_aim1_multivariate_model <- as.data.frame(read_excel("MA_results_full_model_aim1_multivariate_model.xlsx"))
MA_results_nocell_model_aim1_multivariate_model <- as.data.frame(read_excel("MA_results_nocell_model_aim1_multivariate_model.xlsx"))

meta_results_gestbir_model_aim2 <- as.data.frame(read_excel("meta_results_gestbir_model_aim2.xlsx"))
meta_results_icv_model_aim2 <- as.data.frame(read_excel("meta_results_icv_model_aim2.xlsx"))
meta_results_interaction_model_aim2 <- as.data.frame(read_excel("meta_results_interaction_model_aim2.xlsx"))
meta_results_main_model_aim2 <- as.data.frame(read_excel("meta_results_main_model_aim2.xlsx"))
meta_results_nocell_model_aim2 <- as.data.frame(read_excel("meta_results_nocell_model_aim2.xlsx"))

meta_results_gestbir_model_aim3 <- as.data.frame(read_excel("meta_results_gestbir_model_aim3.xlsx"))
meta_results_interaction_model_aim3 <- as.data.frame(read_excel("meta_results_interaction_model_aim3.xlsx"))
meta_results_main_model_aim3 <- as.data.frame(read_excel("meta_results_main_model_aim3.xlsx"))
meta_results_nocell_model_aim3 <- as.data.frame(read_excel("meta_results_nocell_model_aim3.xlsx"))

# apply fdr
MA_results_univariate_model_fdr <- add_fdr_pvalue(MA_results_univariate_model, 4)
MA_results_batch_model_fdr <- add_fdr_pvalue(MA_results_batch_model, 4)
MA_results_batch_cell_model_fdr <- add_fdr_pvalue(MA_results_batch_cell_model, 4)

MA_results_full_model_aim1_multivariate_model_fdr <- add_fdr_pvalue(MA_results_full_model_aim1_multivariate_model, 4)
MA_results_nocell_model_aim1_multivariate_model_fdr <- add_fdr_pvalue(MA_results_nocell_model_aim1_multivariate_model, 4)

meta_results_gestbir_model_aim2_fdr <- add_fdr_pvalue(meta_results_gestbir_model_aim2, 4)
meta_results_icv_model_aim2_fdr <- add_fdr_pvalue(meta_results_icv_model_aim2, 4)
meta_results_interaction_model_aim2_fdr <- add_fdr_pvalue(meta_results_interaction_model_aim2, 4)
meta_results_main_model_aim2_fdr <- add_fdr_pvalue(meta_results_main_model_aim2, 4)
meta_results_nocell_model_aim2_fdr <- add_fdr_pvalue(meta_results_nocell_model_aim2, 4)

meta_results_gestbir_model_aim3_fdr <- add_fdr_pvalue(meta_results_gestbir_model_aim3, 4)
meta_results_interaction_model_aim3_fdr <- add_fdr_pvalue(meta_results_interaction_model_aim3, 4)
meta_results_main_model_aim3_fdr <- add_fdr_pvalue(meta_results_main_model_aim3, 4)
meta_results_nocell_model_aim3_fdr <- add_fdr_pvalue(meta_results_nocell_model_aim3, 4)

# save to excel
write_xlsx(MA_results_univariate_model_fdr, 'MA_results_univariate_model_fdr.xlsx')
write_xlsx(MA_results_batch_model_fdr, 'MA_results_batch_model_fdr.xlsx')
write_xlsx(MA_results_batch_cell_model_fdr, 'MA_results_batch_cell_model_fdr.xlsx')

write_xlsx(MA_results_full_model_aim1_multivariate_model_fdr, 'MA_results_full_model_aim1_multivariate_model_fdr.xlsx')
write_xlsx(MA_results_nocell_model_aim1_multivariate_model_fdr, 'MA_results_nocell_model_aim1_multivariate_model_fdr.xlsx')

write_xlsx(meta_results_gestbir_model_aim2_fdr, 'meta_results_gestbir_model_aim2_fdr.xlsx')
write_xlsx(meta_results_icv_model_aim2_fdr, 'meta_results_icv_model_aim2_fdr.xlsx')
write_xlsx(meta_results_interaction_model_aim2_fdr, 'meta_results_interaction_model_aim2_fdr.xlsx')
write_xlsx(meta_results_main_model_aim2_fdr, 'meta_results_main_model_aim2_fdr.xlsx')
write_xlsx(meta_results_nocell_model_aim2_fdr, 'meta_results_nocell_model_aim2_fdr.xlsx')

write_xlsx(meta_results_gestbir_model_aim3_fdr, 'meta_results_gestbir_model_aim3_fdr.xlsx')
write_xlsx(meta_results_interaction_model_aim3_fdr, 'meta_results_interaction_model_aim3_fdr.xlsx')
write_xlsx(meta_results_main_model_aim3_fdr, 'meta_results_main_model_aim3_fdr.xlsx')
write_xlsx(meta_results_nocell_model_aim3_fdr, 'meta_results_nocell_model_aim3_fdr.xlsx')

#\

###-----------------------------------------------###
#-----------------------PART 6----------------------#
###----------------------------------------------###

#====================================================
### Step 1: Create a forest plot of aim 1 results
# create dataframes with results aim 1 per mps 
df_forest_plot_ligthart_aim1 <- data.frame(
  Predictor = c('Prenatal stress score', 'Prenatal infection score', 'Lifestyle pro-inflammatory factors score', 'Prenatal inflammatory conditions score'),
  beta_estimates = c(-0.005, 0.006, 0.033, -0.009),
  se_estimates = c(0.018, 0.014, 0.015, 0.014)
)
df_forest_plot_ligthart_aim1$Outcome <- "MPS of Ligthart" 

df_forest_plot_wielscher_aim1 <- data.frame(
  Predictor = c('Prenatal stress score', 'Prenatal infection score', 'Lifestyle pro-inflammatory factors score', 'Prenatal inflammatory conditions score'),
  beta_estimates = c(-0.001, 0.006, -0.003, -0.017),
  se_estimates = c(0.014, 0.012, 0.012, 0.012)
)
df_forest_plot_wielscher_aim1$Outcome <- "MPS of Wielscher" 

# create combined data frame 
df_forest_plot_aim1 <- rbind(df_forest_plot_ligthart_aim1, df_forest_plot_wielscher_aim1)

# create plot for both MPSES
aim1 <- ggplot(df_forest_plot_aim1, aes(x = Predictor, y = beta_estimates, fill = Outcome)) +
  geom_point(position = position_dodge(width = 0.5), aes(color = Outcome), size = 4) +
  geom_errorbar(aes(ymin = beta_estimates - 1.96 * se_estimates, 
                    ymax = beta_estimates + 1.96 * se_estimates), 
                position = position_dodge(width = 0.5), width = 0.25) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "red") +
  labs(title = "Aim 1: Mapping prenatal predictors of MPS-CRP", 
       y = "Standardized beta coefficients", 
       x = "Prenatal predictor") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        axis.title.x = element_text(face = 'bold'),
        axis.title.y = element_text(face = 'bold'),
        title = element_text(face = 'bold')) +
  scale_fill_manual(values = c("MPS of Ligthart" = "#ffa600", "MPS of Wielscher" = "#58508d")) +
  scale_color_manual(values = c("#ffa600", "#58508d")) 

aim1

# create plot for only ligthart 
aim1_ligthart <- ggplot(df_forest_plot_ligthart_aim1, aes(x = Predictor, y = beta_estimates, fill = Outcome)) +
  geom_point(position = position_dodge(width = 0.5), aes(color = Outcome), size = 4) +
  geom_errorbar(aes(ymin = beta_estimates - 1.96 * se_estimates, 
                    ymax = beta_estimates + 1.96 * se_estimates), 
                position = position_dodge(width = 0.5), width = 0.25) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "red") +
  labs(title = "Aim 1: Mapping prenatal predictors of MPS-CRP", 
       y = "Standardized beta coefficients", 
       x = "Prenatal predictor") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        axis.title.x = element_text(face = 'bold'),
        axis.title.y = element_text(face = 'bold'),
        title = element_text(face = 'bold')) +
  scale_fill_manual(values = c("MPS of Ligthart" = "#ffa600")) +
  scale_color_manual(values = c("#ffa600")) 

aim1_ligthart

#====================================================
### Step 2: Create a forest plot of aim 2 results (interaction model)
# create dataframes with results aim 1 per mps 
df_forest_plot_ligthart_aim2 <- data.frame(
  Outcome = c('Brain stem volume', 'Gray matter volume', 'Total brain volume', 'Cerebellum volume', 'Hippocampus volume', 'Amygdala volume', 'White matter volume', 'Lateral ventricles volume', 'Global mean MD', 'Global mean FA'),
  beta_estimates = c(-0.001, -0.007, -0.005, -0.007, 0.001, 0.016, -0.002, 0.001, 0.004, -0.014),
  se_estimates = c(0.004, 0.004, 0.003, 0.003, 0.005, 0.006, 0.002, 0.003, 0.006, 0.007)
)
df_forest_plot_ligthart_aim2$Predictor <- "MPS of Ligthart" 

df_forest_plot_wielscher_aim2 <- data.frame(
  Outcome = c('Brain stem volume', 'Gray matter volume', 'Total brain volume', 'Cerebellum volume', 'Hippocampus volume', 'Amygdala volume', 'White matter volume', 'Lateral ventricles volume', 'Global mean MD', 'Global mean FA'),
  beta_estimates = c(-0.001, -0.005, -0.004, -0.005, 0.002, 0.012, -0.004, 0, 0.001, -0.018),
  se_estimates = c(0.004, 0.004, 0.003, 0.003, 0.004, 0.006, 0.002, 0.003, 0.006, 0.006)
)
df_forest_plot_wielscher_aim2$Predictor <- "MPS of Wielscher" 

# create combined data frame 
df_forest_plot_aim2 <- rbind(df_forest_plot_ligthart_aim2, df_forest_plot_wielscher_aim2)

# create plot
aim2 <- ggplot(df_forest_plot_aim2, aes(x = Outcome, y = beta_estimates, fill = Predictor)) +
  geom_point(position = position_dodge(width = 0.5), aes(color = Predictor), size = 4) +
  geom_errorbar(aes(ymin = beta_estimates - 1.96 * se_estimates, 
                    ymax = beta_estimates + 1.96 * se_estimates), 
                position = position_dodge(width = 0.5), width = 0.25)  +
  geom_hline(yintercept = 0, linetype = "dotted", color = "red") +
  labs(title = "Aim 2: Mapping brain outcomes of MPS-CRP (interaction model)", 
       y = "Standardized beta coefficients", 
       x = "sMRI and DTI outcome") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        axis.title.x = element_text(face = 'bold'),
        axis.title.y = element_text(face = 'bold'),
        title = element_text(face = 'bold')) +
  scale_fill_manual(values = c("MPS of Ligthart" = "#ffa600", "MPS of Wielscher" = "#58508d")) +
  scale_color_manual(values = c("#ffa600", "#58508d")) 

aim2

# create plot for only ligthart 
aim2_ligthart <- ggplot(df_forest_plot_ligthart_aim2, aes(x = Outcome, y = beta_estimates, fill = Predictor)) +
  geom_point(position = position_dodge(width = 0.5), aes(color = Predictor), size = 4) +
  geom_errorbar(aes(ymin = beta_estimates - 1.96 * se_estimates, 
                    ymax = beta_estimates + 1.96 * se_estimates), 
                position = position_dodge(width = 0.5), width = 0.25)  +
  geom_hline(yintercept = 0, linetype = "dotted", color = "red") +
  labs(title = "Aim 2: Mapping brain outcomes of MPS-CRP (interaction model)", 
       y = "Standardized beta coefficients", 
       x = "sMRI and DTI outcome") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        axis.title.x = element_text(face = 'bold'),
        axis.title.y = element_text(face = 'bold'),
        title = element_text(face = 'bold')) +
  scale_fill_manual(values = c("MPS of Ligthart" = "#ffa600")) +
  scale_color_manual(values = c("#ffa600")) 

aim2_ligthart

#\

#====================================================
### Step 3: Create a forest plot of aim 3 results (main effects model)
# create dataframes with results aim 1 per mps 
df_forest_plot_ligthart_aim3 <- data.frame(
  Outcome = c('CBCL - total behavioral problems score', 'CBCL - internalizing problems score', 'CBCL - externalizing problems score'),
  beta_estimates = c(0.023, 0.016, 0.012),
  se_estimates = c(0.03, 0.028, 0.027)
)
df_forest_plot_ligthart_aim3$Predictor <- "MPS of Ligthart" 

df_forest_plot_wielscher_aim3 <- data.frame(
  Outcome = c('CBCL - total behavioral problems score', 'CBCL - internalizing problems score', 'CBCL - externalizing problems score'),
  beta_estimates = c(-0.011, 0.013, -0.025),
  se_estimates = c(0.037, 0.035, 0.033)
)
df_forest_plot_wielscher_aim3$Predictor <- "MPS of Wielscher" 

# create combined data frame 
df_forest_plot_aim3 <- rbind(df_forest_plot_ligthart_aim3, df_forest_plot_wielscher_aim3)

# creating plot
aim3 <- ggplot(df_forest_plot_aim3, aes(x = Outcome, y = beta_estimates, fill = Predictor)) +
  geom_point(position = position_dodge(width = 0.5), aes(color = Predictor), size = 4) +
  geom_errorbar(aes(ymin = beta_estimates - 1.96 * se_estimates, 
                    ymax = beta_estimates + 1.96 * se_estimates), 
                position = position_dodge(width = 0.5), width = 0.25)  +
  geom_hline(yintercept = 0, linetype = "dotted", color = "red") +
  labs(title = "Aim 3: Mapping behavior outcomes of MPS-CRP (main model)", 
       y = "Standardized beta coefficients", 
       x = "CBCL outcome") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        axis.title.x = element_text(face = 'bold'),
        axis.title.y = element_text(face = 'bold'),
        title = element_text(face = 'bold')) +
  scale_fill_manual(values = c("MPS of Ligthart" = "#ffa600", "MPS of Wielscher" = "#58508d")) +
  scale_color_manual(values = c("#ffa600", "#58508d")) 

aim3

# create plot for only ligthart 
aim3_ligthart <- ggplot(df_forest_plot_ligthart_aim3, aes(x = Outcome, y = beta_estimates, fill = Predictor)) +
  geom_point(position = position_dodge(width = 0.5), aes(color = Predictor), size = 4) +
  geom_errorbar(aes(ymin = beta_estimates - 1.96 * se_estimates, 
                    ymax = beta_estimates + 1.96 * se_estimates), 
                position = position_dodge(width = 0.5), width = 0.25)  +
  geom_hline(yintercept = 0, linetype = "dotted", color = "red") +
  labs(title = "Aim 3: Mapping behavior outcomes of MPS-CRP (main model)", 
       y = "Standardized beta coefficients", 
       x = "CBCL outcome") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        axis.title.x = element_text(face = 'bold'),
        axis.title.y = element_text(face = 'bold'),
        title = element_text(face = 'bold')) +
  scale_fill_manual(values = c("MPS of Ligthart" = "#ffa600")) +
  scale_color_manual(values = c("#ffa600")) 

aim3_ligthart

#\ 

#====================================================
### Step 4: Combine all forest plots 
ggarrange(aim1, aim2, aim3, ncol = 1, nrow = 3, align = 'v')

# Ligthart only plot 
ggarrange(aim1_ligthart, 
          aim2_ligthart, 
          aim3_ligthart, 
          ncol = 1, 
          nrow = 3, 
          align = 'v')

# FINAL PART OF SCRIPTS. 
