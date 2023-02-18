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

####### Read Data #################################

load(file = "data/combined/joint_count_array.rda")
load(file = "data/combined/joint_count_metadata.rda")
load(file = "data/combined/temporal_covariates.rda")
load(file = "data/combined/time_design.rda")
load(file = "data/combined/landcover_covariates.rda")
load(file = "data/combined/dist_design.rda")

####### Set Constants #############################

n_cores <- 15
# For testing purposes
species <- c("AMRO", "WTSP", "OVEN", "BLPW", "REVI")

####### Data Wrangling ############################

# Add indexing variable
joint_count_metadata$Index <- seq(1,nrow(joint_count_metadata))

# Drop method I
dist_design <- dist_design[-which(dist_design$Method == "I"), ]
joint_count_metadata <- joint_count_metadata[-which(joint_count_metadata$Distance_Method == "I"), ]

# Filter to only species that will be modelled
# This will eventually just be species with counts >75, but for now
# we will just set a small selection of species for testing
joint_count_metadata <- joint_count_metadata[which(joint_count_metadata$Species %in% species), ]

# Create time design matrix
names(time_design)[1] <- "Time_Method"
count_time_design <- plyr::join(joint_count_metadata[,c("Index", "Time_Method")], time_design,
                                by = "Time_Method", type = "left")

names(dist_design)[1] <- "Distance_Method"
count_dist_design <- plyr::join(joint_count_metadata[,c("Index", "Distance_Method")], dist_design,
                                by = "Distance_Method", type = "left")

