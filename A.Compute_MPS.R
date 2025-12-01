#########################################################################
######################PROJECT: Neonatal inflammation####################
#########################################################################

#Project: Mapping prenatal predictors and neurobehavioral outcomes of an epigenetic marker of neonatal inflammation - a longitudinal population-based study 

#### STEP 1 - create scores ####
#set directory 
setwd("set_your_path")

#loading imputed matrix 
load("/home/n.creasey/gc_epgs/450kbetas_imputed.Rdata")

#load sumstats 
sumstats_ligthart <- read.csv("/home/n.creasey/episcore_CRP_project/EpiscoreCRP_CPGsites_ligthart.csv")
sumstats_wielscher <- read.csv("/home/n.creasey/episcore_CRP_project/EpiscoreCRP_CPGsites_wielscher.csv")

methscore <- function(betas,sumstats,filename){

  library(dplyr)

  #extract the betas
  epgs_betas <- as.data.frame(betas[rownames(betas) %in% sumstats$CpG,])
  
  #identify missings CpGs and remove from sumstats
  missing <- setdiff(sumstats$CpG, row.names(epgs_betas))
  sumstats <- sumstats %>% filter(!CpG %in% (missing))
  myfilename <- filename
  write.table(missing,file=myfilename, row.names = FALSE, col.names = FALSE) 
  
  #transpose the betas
  t_epgs_betas <- data.frame(t(epgs_betas))
  
  #prepare coefs
  coefs <- sumstats$Weights
  names(coefs) <- sumstats$CpG
  coefs <- sumstats$Weights
  names(coefs) <- sumstats$CpG
  
  #compute
  score <- (t(epgs_betas)) %*% coefs
  
  #create dataframe
  score <- cbind(rownames(score), data.frame(score, row.names=NULL))
  colnames(score)[1] <- "Sample_ID" #change this to your matching variable name
  colnames(score)[2] <- "mps"
  
  #z score
  score$zmps <- scale(score$mps,center = TRUE,scale = TRUE)
  colnames(score)[3] <- "zmps"
  
  #return 
  score
  
}

#compute 450k scores
mps_Lighart450k <- methscore(betas=x_imputedmatrix,sumstats=sumstats_ligthart,filename="lighart450k_missingcgs")
write.csv(mps_Lighart450k, "/home/n.creasey/episcore_CRP_project/shadow2/computed_scores/mps_Lighart450k.csv",row.names = FALSE)

mps_wielscher450k <- methscore(betas=x_imputedmatrix,sumstats=sumstats_wielscher,filename="wielscher450k_missingcgs")
write.csv(mps_wielscher450k, "/home/n.creasey/episcore_CRP_project/shadow2/computed_scores/mps_wielscher450k_withimputation.csv",row.names = FALSE)

#remove 450k data 
rm(x_imputedmatrix)

#load EPIC data 
library(miceadds)
load.Rdata("/home/n.creasey/GENR3/Methylation/GENR_EPICMETH_Norm/GENR_EPICv1METH_Norm_Betas_birth_ALL.RData", "y")

#compute EPIC scores
mps_LighartEPIC <- methscore(betas=y,sumstats=sumstats_ligthart,filename="lighartEPIC_missingcgs")
write.csv(mps_LighartEPIC, "/home/n.creasey/episcore_CRP_project/shadow2/computed_scores/mps_LighartEPIC.csv",row.names = FALSE)

mps_wielscherEPIC <- methscore(betas=y,sumstats=sumstats_wielscher,filename="wielscherEPIC_missingcgs")
write.csv(mps_wielscherEPIC, "/home/n.creasey/episcore_CRP_project/shadow2/computed_scores/mps_wielscherEPIC.csv",row.names = FALSE)

#merge the data to one set 
colnames(mps_Lighart450k)[2] <- "mps_lig" 
colnames(mps_Lighart450k)[3] <- "mpsz_lig"

colnames(mps_wielscher450k)[2] <- "mps_wiel" 
colnames(mps_wielscher450k)[3] <- "mpsz_wiel"

colnames(mps_LighartEPIC)[2] <- "mps_lig" 
colnames(mps_LighartEPIC)[3] <- "mpsz_lig"

colnames(mps_wielscherEPIC)[2] <- "mps_wiel" 
colnames(mps_wielscherEPIC)[3] <- "mpsz_wiel"

## merge (and add labels for array type)
allmps_450k <- merge(mps_Lighart450k,mps_wielscher450k, by="Sample_ID")
allmps_450k$array <- as.factor("450k")

allmps_EPIC <- merge(mps_LighartEPIC,mps_wielscherEPIC, by="Sample_ID")
allmps_EPIC$array <- as.factor("epic")

identical(colnames(allmps_450k),colnames(allmps_EPIC)) #check column names match

allmps <- rbind(allmps_450k, allmps_EPIC)
write.csv(allmps, "path_to_storing_results", row.names = FALSE)

# END OF THIS SCRIPT, GO TO PART B1.
