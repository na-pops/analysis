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
load("data/dist_design.rda")
load("data/time_design.rda")
rem_coef <- read.csv("../results/detection-coefs/removal_coef.csv")
dist_coef <- read.csv("../results/detection-coefs/distance_coef.csv")

####### Calculate Offsets #########################
#find species that had sufficient sample size to estimate both dist and rem coefs
species <- intersect(unique(dist_coef$SPPCODE), unique(rem_coef$SPPCODE))

#loop through each species and calculate offsets
for (i in species) 
{
  spp = project_counts[project_counts$Species==i,]
  spp = plyr::join(PKEY, spp, by="PKEY", type="left")
  
  #change NA values in Y column to 0
  spp$Species[is.na(spp$Species)] <- i
  spp$Y[is.na(spp$Y)] <- 0
  
  #add species-specific logtau and tau from dist_coef (output from script 3)
  spp$logtau = dist_coef[which(dist_coef$SPPCODE==i),"logtau"]
  spp$tau = exp(spp$logtau)
  
  #add all species-specific coefficient estimates from rem_coef (output from script 2)
  spp = plyr::join(spp, rem_coef, by="Species", type="left")
  
  #calculate survey area (A)
  spp$A = ifelse(spp$MaxDist=="Inf", pi*spp$tau^2, pi*spp$MaxDist^2)
  
  #Calculate probability of singing during survey (p) 
  #Phi calculation is based on which removal model was best for each species (model column)
  #These are the phi models to choose from:
  #model 1 phi <- exp(sraint)
  #model 2 phi <- exp(sraint + (sratssr*TSSR))
  #model 4 phi <- exp(sraint + (srajday*JDAY) + (sratssr*TSSR))
  #model 5 phi <- exp(sraint + (srajday*JDAY) + (sratssr*TSSR) + (sratssrjday*TSSR*JDAY))
  #model 6 phi <- exp(sraint + (sratssr2*TSSR2))
  if (unique(spp$model)==1) {
    phi = exp(spp$sraint)  
  }
  if (unique(spp$model)==2) {
    phi = exp(spp$sraint + spp$sratssr*spp$TSSR) 
  }
  if (unique(spp$model)==3) {
    phi = exp(spp$sraint + spp$srajday*spp$JDAY)
  }
  if (unique(spp$model)==4) {
    phi = exp(spp$sraint + spp$sratssr*spp$TSSR + spp$srajday*spp$JDAY)
  }
  if (unique(spp$model)==5) {
    phi = exp(spp$sraint + spp$sratssr*spp$TSSR + spp$srajday*spp$JDAY + spp$sratssrjday*spp$TSSR*spp$JDAY)
  }
  if (unique(spp$model)==6) {
    phi =  exp(spp$sraint + spp$sratssr2*spp$TSSR2)
  }
  spp$p = 1-exp(-(spp$MaxDur*phi))
  
  #Calculate perceptibility (q)
  spp$q = ifelse(spp$MaxDist=="Inf", 1, (spp$tau^2/spp$MaxDist^2)*(1-exp(-(spp$MaxDist^2)/spp$tau^2)))
  
  #Calculate the QPAD offsets
  spp$offset = log(spp$A) + log(spp$p) + log(spp$q)
  
  #create species-specific object containing count data and offset calculation for each survey (PKEY)
  assign(i,spp)
  rm(spp,phi,i)
}