####### Script Information ########################
# Brandon P.M. Edwards
# NA-POPS: analysis
# 1-combine-data.R
# Created August 2020
# Last Updated October 2020

####### Import Libraries and External Files #######

library(here)
library(reshape2)

source("../utilities/get-data.R")

####### Read Data #################################

# Get project names
project_list <- read.table(here::here("../utilities/proj-list"))
n_proj <- nrow(project_list)

# Create combined count data frame
project_counts <- vector('list', n_proj)
project_samples <- vector('list', n_proj)
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

####### Wrangle Data ##############################

# Create time count matrix
time_only <- project_counts[, c(1:3, 7,8 )]
time_only <- time_only[-which(time_only$Time_Method == "ZZ"), ]
time_count_matrix <- dcast(time_only,
                           Sample_ID + Species + Time_Method ~ Time_Level,
                           value.var = "Abundance",
                           fun.aggregate = sum)

# Create distance count matrix
dist_only <- project_counts[, c(1:5)]
dist_only <- dist_only[-which(dist_only$Distance_Method == "ZZ"), ]
dist_count_matrix <- dcast(dist_only,
                           Sample_ID + Species + Distance_Method ~ Distance_Level,
                           value.var = "Abundance",
                           fun.aggregate = sum)

# Create combined landcover covariate df and temporal covariate df
landcover_covariates <- vector('list', n_proj)
temporal_covariates <- vector('list', n_proj)

for (i in 1:n_proj)
{
  p <- project_list$V1[i]
  
  # Landcover covariates
  data_dir <- paste0("../covariates/landcover/project-",
                     p,
                     ".csv")
  landcover_covariates[[i]] <- read.csv(here::here(data_dir))
  
  # Temporal covariates
  data_dir <- paste0("../covariates/temporal/project-",
                     p,
                     ".csv")
  temporal_covariates[[i]] <- read.csv(here::here(data_dir))
}
landcover_covariates <- do.call(rbind, landcover_covariates)
temporal_covariates <- do.call(rbind, temporal_covariates)

# Create time removal design matrix
time <- read.csv(here::here("../utilities/time_lookup.csv"))
time <- time[-c(1), ]
time_design <- dcast(time, Method + Max_Duration ~ Level, value.var = "End_Duration")

# Create distance design matrix
dist <- read.csv(here::here("../utilities/distance_lookup.csv"))
dist <- dist[-c(1), ]
dist_design <- dcast(dist, Method + Max_Distance ~ Level, value.var = "End_Distance")

####### Output Data ###############################

save(project_counts, file = here::here("data/counts.rda"))
save(project_samples, file = here::here("data/samples.rda"))
save(time_count_matrix, file = here::here("data/time_count_matrix.rda"))
save(dist_count_matrix, file = here::here("data/dist_count_matrix.rda"))
save(landcover_covariates, file = here::here("data/landcover_covariates.rda"))
save(temporal_covariates, file = here::here("data/temporal_covariates.rda"))
save(time_design, file = here::here("data/time_design.rda"))
save(dist_design, file = here::here("data/dist_design.rda"))