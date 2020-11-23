####### Script Information ########################
# Brandon P.M. Edwards
# NA-POPS: analysis
# 4-QPAD-offsets.R
# Created August 2020
# Last Updated August 2020

####### Import Libraries and External Files #######

library(mefa4)
library(data.table)

####### Read Data #################################

load("data/counts.rda")
load("data/samples.rda")
load("data/dist_design.rda")
load("data/time_design.rda")
load("data/covariates.rda")
rem_coef <- read.csv("../results/detection-coefs/removal_coef.csv")
dist_coef <- read.csv("../results/detection-coefs/distance_coef.csv")

####### Data Wrangling ############################
covariates$ForestOnly_5x5 <- covariates$ForestOnly_5x5 / 25

# Get max distance and durations for each sample
project_samples <- merge(project_samples,
                         dist_design[, c("Method", "Max_Distance")],
                         by.x = "Distance_Method",
                         by.y = "Method")
project_samples <- merge(project_samples,
                         time_design[, c("Method", "Max_Duration")],
                         by.x = "Time_Method",
                         by.y = "Method")
project_samples <- merge(project_samples,
                         covariates[,c("Sample_ID", "roaddist", "ForestOnly_5x5")],
                         by = "Sample_ID")

# Remove any NAs from dist_coef and rem_coef
dist_coef <- dist_coef[!is.na(dist_coef$d_int), ]
rem_coef <- rem_coef[!is.na(rem_coef$r_int), ]

####### Calculate Offsets #########################

#find species that had sufficient sample size to estimate both dist and rem coefs
species <- intersect(unique(dist_coef$Species), unique(rem_coef$Species))

#loop through each species and calculate offsets
for (i in species) 
{
  spp = project_counts[project_counts$Species==i,]
  spp = plyr::join(project_samples[,c("Sample_ID",
                                      "roaddist",
                                      "ForestOnly_5x5",
                                      "Max_Distance",
                                      "Max_Duration")], 
                   spp, 
                   by="Sample_ID",
                   type="left")
  
  #change NA values in Y column to 0
  spp$Species[is.na(spp$Species)] <- i
  spp$Y[is.na(spp$Abundance)] <- 0
  
  #add all species-specific coefficient estimates from dist_coef
  spp = plyr::join(spp, dist_coef, by="Species", type="left")
  
  #add species-specific logtau and tau from dist_coef (output from script 3)
  if (unique(spp$d_model)==1) {
    spp$tau = exp(spp$d_int)
  }
  if (unique(spp$d_model)==2) {
    spp$tau = exp(spp$d_int + spp$d_roaddist*spp$roaddist) 
  }
  if (unique(spp$d_model)==3) {
    spp$tau = exp(spp$d_int + spp$d_forest*spp$ForestOnly_5x5)
  }
  if (unique(spp$d_model)==4) {
    spp$tau = exp(spp$d_int + spp$d_roaddist*spp$roaddist + spp$d_forest*spp$ForestOnly_5x5)
  }
  if (unique(spp$d_model)==5) {
    spp$tau = exp(spp$d_int + spp$d_roaddist*spp$roaddist + spp$d_forest*spp$ForestOnly_5x5 + spp$d_roaddist*spp$roaddist*spp$d_forest*spp$ForestOnly_5x5)
  }
  spp$logtau <- log(spp$tau)
  
  #calculate survey area (A)
  spp$A = ifelse(spp$Max_Distance == "Inf", pi*spp$tau^2, pi*spp$Max_Distance^2)
  
  #add all species-specific coefficient estimates from rem_coef (output from script 2)
  spp = plyr::join(spp, rem_coef, by="Species", type="left")

  #Calculate probability of singing during survey (p) 
  #Phi calculation is based on which removal model was best for each species (model column)

  if (unique(spp$r_model)==1) {
    phi = exp(spp$r_int)
  }
  if (unique(spp$r_model)==2) {
    phi = exp(spp$r_int + spp$r_roaddist*spp$roaddist) 
  }
  if (unique(spp$r_model)==3) {
    phi = exp(spp$r_int + spp$r_forest*spp$ForestOnly_5x5)
  }
  if (unique(spp$r_model)==4) {
    phi = exp(spp$r_int + spp$r_roaddist*spp$roaddist + spp$r_forest*spp$ForestOnly_5x5)
  }
  if (unique(spp$r_model)==5) {
    phi = exp(spp$r_int + spp$r_roaddist*spp$roaddist + spp$r_forest*spp$ForestOnly_5x5 + spp$r_roaddist*spp$roaddist*spp$r_forest*spp$ForestOnly_5x5)
  }
  spp$p = 1-exp(-(spp$Max_Duration*phi))
  
  #Calculate perceptibility (q)
  spp$q = ifelse(spp$Max_Distance == "Inf", 1, (spp$tau^2/spp$Max_Distance^2)*(1-exp(-(spp$Max_Distance^2)/spp$tau^2)))
  
  #Calculate the QPAD offsets
  spp$offset = log(spp$A) + log(spp$p) + log(spp$q)
  
  write.table(x = spp,
              file = paste0("../results/offsets/",
                            i,
                            ".csv"),
              row.names = FALSE,
              sep = ",")
  
  #create species-specific object containing count data and offset calculation for each survey (PKEY)
  # assign(i,spp)
  # rm(spp,phi,i)
}

