####### Script Information ########################
# Brandon P.M. Edwards
# NA-POPS: analysis
# 1-combine-data.R
# Created August 2020
# Last Updated March 2021

####### Import Libraries and External Files #######

library(reshape2)

source("../utilities/get-data.R")

####### Read Data #################################

time <- read.csv("../covariates/survey/time_lookup.csv")
dist <- read.csv("../covariates/survey/distance_lookup.csv")

# Get project names
project_list <- read.table("../utilities/proj-list")
n_proj <- nrow(project_list)

# Create combined count data frame
project_counts <- vector('list', n_proj)
project_samples <- vector('list', n_proj)
message("1/7 Now reading and combining project counts and samples.\n")
for (i in 1:n_proj)
{
  p <- project_list$V1[i]
  # Get counts
  data_dir <- paste0("../project-",
                     p,
                     "/output/",
                     p,
                     "_counts.rda")
  project_counts[[i]] <- data.frame(get_data(data_dir))
  
  # Get samples
  data_dir <- paste0("../project-",
                     p,
                     "/output/",
                     p,
                     "_samples.rda")
  project_samples[[i]] <- data.frame(get_data(data_dir))
}
project_counts <- do.call(rbind, project_counts)
project_samples <- do.call(rbind, project_samples)

# Create combined landcover covariate df and temporal covariate df
landcover_covariates <- vector('list', n_proj)
temporal_covariates <- vector('list', n_proj)
message("2/7 Now reading and combining covariates.\n")
for (i in 1:n_proj)
{
  p <- project_list$V1[i]
  
  # Landcover covariates
  data_dir <- paste0("../covariates/landcover/project-",
                     p,
                     ".csv")
  landcover_covariates[[i]] <- read.csv(data_dir)
  
  # Temporal covariates
  data_dir <- paste0("../covariates/temporal/project-",
                     p,
                     ".csv")
  temporal_covariates[[i]] <- read.csv(data_dir)
}
landcover_covariates <- do.call(rbind, landcover_covariates)
temporal_covariates <- do.call(rbind, temporal_covariates)

####### Wrangle Data ##############################

# Create time count matrix
message("3/7 Now creating time_count_matrix with dcast().\n")
time_only <- project_counts[, c(1:3, 7,8 )]
time_only <- time_only[-which(time_only$Time_Method == "ZZ"), ]
time_only <- time_only[-which(is.na(time_only$Time_Level)), ]
time_count_matrix <- dcast(time_only,
                           Sample_ID + Species + Time_Method ~ as.numeric(Time_Level),
                           value.var = "Abundance",
                           fun.aggregate = sum)

# Create distance count matrix
message("4/7 Now creating dist_count_matrix with dcast().\n")
dist_only <- project_counts[, c(1:5)]
dist_only <- dist_only[-which(dist_only$Distance_Method == "ZZ"), ]
dist_only <- dist_only[-which(is.na(dist_only$Distance_Level)), ]
dist_count_matrix <- dcast(dist_only,
                           Sample_ID + Species + Distance_Method ~ as.numeric(Distance_Level),
                           value.var = "Abundance",
                           fun.aggregate = sum)

# Create time removal design matrix
message("5/7 Now creating time_design with dcast().\n")
time <- time[-c(1), ]
time_design <- dcast(time, Method + Max_Duration ~ Level, value.var = "End_Duration")

# Create distance design matrix
message("6/7 Now creating dist_design with dcast().\n")
dist <- dist[-c(1), ]
dist_design <- dcast(dist, Method + Max_Distance ~ Level, value.var = "End_Distance")

####### Output Data ###############################
message("7/7 Now saving all data.\n")
save(project_counts, file = "data/counts.rda")
save(project_samples, file = "data/samples.rda")
save(time_count_matrix, file = "data/time_count_matrix.rda")
save(dist_count_matrix, file = "data/dist_count_matrix.rda")
save(landcover_covariates, file = "data/landcover_covariates.rda")
save(temporal_covariates, file = "data/temporal_covariates.rda")
save(time_design, file = "data/time_design.rda")
save(dist_design, file = "data/dist_design.rda")