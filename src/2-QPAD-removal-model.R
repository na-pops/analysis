####### Script Information ########################
# Brandon P.M. Edwards
# NA-POPS: analysis
# 2-QPAD-removal-model.R
# Created August 2020
# Last Updated February 2021

####### Import Libraries and External Files #######

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

library(plyr)

####### Read Data #################################

load("data/time_count_matrix.rda")
load("data/temporal_covariates.rda")
load("data/time_design.rda")
na_sp_list <- read.csv("../utilities/IBP-Alpha-Codes20.csv")

####### Data Wrangling ############################

# Drop rows that have NA in JD or TSSR
covariates_reduced <- temporal_covariates[which(!is.na(temporal_covariates$JD)), ]
covariates_reduced <- covariates_reduced[which(!is.na(covariates_reduced$TSSR)), ]

count_names <- c("Sample_ID", "Species", "Time_Method",
                 paste0(rep("Int", times = 10), 1:10))
names(time_count_matrix) <- count_names

design_names <- c("Time_Method", "Max_Duration",
                 paste0(rep("Interval", times = 10), 1:10))
names(time_design) <- design_names

# Join data
count_design <- plyr::join(time_count_matrix, time_design,
                           by = "Time_Method", type = "left")

# Filter out methods with only 1 interval
to_remove <- which(is.na(count_design$Interval2))
count_design <- count_design[-c(to_remove), ]

# Create separate data frames
counts <- count_design[,c("Sample_ID", "Species", "Time_Method", count_names[4:length(count_names)])]
design <- count_design[,c("Sample_ID", "Species", "Time_Method", design_names[3:length(design_names)])]
names(design) <- names(counts)
col_names <- names(counts)[4:length(names(counts))]

# Change 0s to NA in counts table where appropriate based on design table
for (i in col_names)
{
  indices <- which(counts[, i] == 0)
  
  counts[indices, i] <- ifelse(is.na(design[indices, i]), NA, 0)
}

# Add in covariates to counts matrix
covars <- plyr::join(counts,
                     covariates_reduced, 
                     by = "Sample_ID", 
                     type = "left", 
                     match = "all")

# Remove rows with NA covariates
to_remove <- which(is.na(covars$JD))
if (length(to_remove > 0))
{
  counts <- counts[-c(to_remove), ]
  design <- design[-c(to_remove), ]
  covars <- covars[-c(to_remove), ]  
}

to_remove <- which(is.na(covars$TSSR))
if (length(to_remove > 0))
{
  counts <- counts[-c(to_remove), ]
  design <- design[-c(to_remove), ]
  covars <- covars[-c(to_remove), ]  
}

# Build matrices by species
species_all <- sort(as.character(unique(counts$Species)))
landbirds <- na_sp_list[which(na_sp_list$LANDBIRD == "TRUE"), "SPEC"]
species <- species_all[which(species_all %in% landbirds)]

# List of input counts
Y <- vector(mode = "list", length = length(species))
names(Y) <- species

# List of input design
D <- vector(mode = "list", length = length(species))
names(D) <- species

# List of input covariates
C <- vector(mode = "list", length = length(species))
names(C) <- species

# Species indicator list
sp_list <- vector(mode = "list", length = length(species))
names(sp_list) <- species

for (s in species)
{
  n_counts <- nrow(counts[counts$Species == s, ])
  if (n_counts >= 75)
  {
    Y[[s]] <- as.matrix(counts[counts$Species==s, col_names])
    D[[s]] <- as.matrix(design[design$Species==s, col_names])
    C[[s]] <- subset(covars,
                     Species==s,
                     select = c(14:17))
    sp_list[[s]] <- data.frame(Species = rep(s, n_counts))
  }
}

Y <- Y[lengths(Y) != 0]; Y <- do.call(rbind, Y)
D <- D[lengths(D) != 0]; D <- do.call(rbind, D)
C <- C[lengths(C) != 0]; C <- do.call(rbind, C)
sp_list <- sp_list[lengths(sp_list) != 0]; sp_list <- do.call(rbind, sp_list)

#' Corresponds with "bands_per_sample" in removal.stan
time_bands_per_sample <- unname(apply(Y, 1, function(x) sum(!is.na(x))))

#' Total species abundance per sampling event.
#' I.e., this is the sum of Y_sij over j
#' Corresponds with "abund_per_sample" in removal.stan
total_abund_per_sample <- unname(apply(Y, 1, function(x) sum(x, na.rm = TRUE)))

#' Total samples per species
#' I.e., total i BY SPECIES
#' Corresponds with "samples_per_species" in removal.stan
samples_per_species <- as.vector(table(sp_list[,1]))

#' Correspinds with "X" in removal.stan
X_names <- names(C)
X <- cbind(rep(1, nrow(C)), C)
names(X) <- c("Intercept", X_names)

#' Corresponds with "abund_per_band" in removal.stan
abundance_per_band <- Y
abundance_per_band[is.na(abundance_per_band)] <- 0

#' Corresponds with "max_time" in removal.stan
max_time <- D
max_time[is.na(max_time)] <- 0

n_samples <- nrow(Y)
n_cov <- ncol(X)
n_species <- length(unique(sp_list[,1]))
max_intervals <- ncol(Y)

stan_data <- list(n_samples = n_samples,
                  n_covariates = n_cov,
                  n_species = n_species,
                  max_intervals = max_intervals,
                  samples_per_species = samples_per_species,
                  abund_per_band = abundance_per_band,
                  abund_per_sample = total_abund_per_sample,
                  bands_per_sample = time_bands_per_sample,
                  max_time = max_time,
                  X = X)

########### Modelling #############################

model <- stan_model(file = "models/removal.stan")
stime = system.time(stan_fit <- 
                      sampling(model,
                               data = stan_data,
                               verbose = TRUE,
                               chains = 3,
                               iter = 2000,
                               warmup = 1000,
                               cores = 3,
                               pars = c("gamma"),
                               control = list(adapt_delta = 0.8,
                                              max_treedepth = 15)))
