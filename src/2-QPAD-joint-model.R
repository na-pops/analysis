####### Script Information ########################
# Brandon P.M. Edwards
# NA-POPS: analysis
# 2-QPAD-joint-model.R
# Created February 2023
# Last Updated February 2023

####### Import Libraries and External Files #######

library(detect)
library(plyr)
library(doParallel)
library(foreach)
#library(Rcpp)

#Rcpp::sourceCpp("src/functions/nll_fun.cpp")
#source("cmulti_fit_joint_cpp.R")
source("src/functions/joint_fns.R")

####### Read Data #################################

load(file = "data/combined/joint_count_array.rda")
load(file = "data/combined/joint_count_metadata.rda")
load(file = "data/combined/temporal_covariates.rda")
load(file = "data/combined/time_count_design.rda")
load(file = "data/combined/landcover_covariates.rda")
load(file = "data/combined/dist_count_design.rda")

####### Set Constants #############################

n_cores <- 15
# For testing purposes
species <- c("SAVS", "AMRO", "VESP", "CCSP", "BLPW", "OVEN", "WTSP")

####### Data Wrangling ############################

# Drop method I
dist_count_design <- dist_count_design[-which(dist_count_design$Method == "I"), ]
joint_count_metadata <- joint_count_metadata[-which(joint_count_metadata$Distance_Method == "I"), ]

# Filter to only species that will be modelled
# This will eventually just be species with counts >75, but for now
# we will just set a small selection of species for testing
joint_count_metadata <- joint_count_metadata[which(joint_count_metadata$Species %in% species), ]

for (s in species)
{
  print(s)
  indices <- joint_count_metadata[which(joint_count_metadata$Species == s),
                                  "Index"]
  
  counts <- joint_count_array[indices,,]
  dist_matrix <- as.matrix(dist_count_design[indices, c(4:ncol(dist_count_design))])/100
  time_matrix <- as.matrix(time_count_design[indices, c(4:ncol(time_count_design))])
  
  fit <- cmulti.fit.joint(counts,
                          dist_matrix,
                          time_matrix)
  
  save(fit, file = paste0("data/joint/", s, ".rda"))  
}

