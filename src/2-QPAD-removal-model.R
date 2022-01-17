####### Script Information ########################
# Brandon P.M. Edwards
# NA-POPS: analysis
# 2-QPAD-removal-model.R
# Created August 2020
# Last Updated January 2022

####### Import Libraries and External Files #######

library(detect)
library(plyr)
library(doParallel)
library(foreach)

####### Read Data #################################

load(file = "data/combined/time_count_matrix.rda")
load(file = "data/combined/temporal_covariates.rda")
load(file = "data/combined/time_design.rda")
na_sp_list <- read.csv("../utilities/IBP-Alpha-Codes20.csv")

n_cores <- 4

####### Data Wrangling ############################

# Drop rows that have NA in OD or TSSR
covariates_reduced <- temporal_covariates[which(!is.na(temporal_covariates$OD)), ]
covariates_reduced <- covariates_reduced[which(!is.na(covariates_reduced$TSSR)), ]
covariates_reduced <- covariates_reduced[which(!is.na(covariates_reduced$BCR)), ]

# Drop rows that are before March 1 or later than August 31
covariates_reduced <- covariates_reduced[which(covariates_reduced$OD >= (61)), ]
covariates_reduced <- covariates_reduced[which(covariates_reduced$OD <= (244)), ]

# Drop rows that are earlier than 3 hours before sunrise
covariates_reduced <- covariates_reduced[which(covariates_reduced$TSSR >= (-3)), ]

max_bands <- ncol(time_design) - 2
count_names <- c("Sample_ID", "Species", "Time_Method",
                 paste0(rep("Int", times = max_bands), 1:max_bands))
names(time_count_matrix) <- count_names

design_names <- c("Time_Method", "Max_Duration",
                 paste0(rep("Interval", times = max_bands), 1:max_bands))
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
to_remove <- which(is.na(covars$OD))
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

# Save which sample IDs (and therefore which covariates) were used during removal modelling
rem_covars_used <- covars[, c("Sample_ID", "OD", "TSSR", "BCR")]
rem_covars_used <- cbind(rem_covars_used, design[, c("Time_Method")])
rem_covars_used <- rem_covars_used[!duplicated(rem_covars_used$Sample_ID), ]
names(rem_covars_used)[ncol(rem_covars_used)] <- "Time_Method"
save(rem_covars_used, file = "data/combined/rem_covars_used.rda")

# Build matrices by species
species_all <- sort(as.character(unique(counts$Species)))
landbirds <- na_sp_list[-which(na_sp_list$LANDBIRD == "FALSE" |
                                 na_sp_list$SP == "+"), "SPEC"]
species <- species_all[which(species_all %in% landbirds)]

for (s in species)
{
  if (nrow(counts[counts$Species == s, ]) >= 75)
  {
    counts_sp <- counts[which(counts$Species == s), col_names]
    design_sp <- design[which(design$Species == s), col_names]
    covars_sp <- covars[which(covars$Species == s), 14:15]
    
    covars_sp$OD <- (covars_sp$OD - median(covars_sp$OD)) / 365
    covars_sp$OD2 <- covars_sp$OD ^ 2
    
    covars_sp$TSSR <- (covars_sp$TSSR - median(covars_sp$TSSR)) / 24
    covars_sp$TSSR2 <- covars_sp$TSSR ^ 2
    
    assign(paste0("Y_",s), as.matrix(counts_sp))
    assign(paste0("D_",s), as.matrix(design_sp))
    assign(paste0("C_",s), covars_sp)
  }else
  {
    species <- species[!(species %in% s)]
  }
}

input_list <- vector(mode="list", length=length(species))
names(input_list) <- species
for (s in species)
{
  input_list[[s]] <- list(Y=get(paste0("Y_",s)), D=get(paste0("D_",s)), C=get(paste0("C_",s)))
}

########### Modelling #############################

message(paste0("Beginning parallel modelling for ", length(species), " species.\n"))

cluster <- makeCluster(n_cores, type = "PSOCK")
registerDoParallel(cluster)

foreach(sp = names(input_list), .packages = 'detect') %dopar%
  {
    x <- input_list[[sp]]

    m1 = cmulti(x$Y | x$D ~ 1, type="rem")
    m2 = cmulti(x$Y | x$D ~ 1 + x$C$OD, type="rem")
    m3 = cmulti(x$Y | x$D ~ 1 + x$C$OD + x$C$OD2, type = "rem")
    m4 = cmulti(x$Y | x$D ~ 1 + x$C$TSSR, type="rem")
    m5 = cmulti(x$Y | x$D ~ 1 + x$C$TSSR + x$C$TSSR2, type = "rem")
    m6 = cmulti(x$Y | x$D ~ 1 + x$C$TSSR + x$C$OD, type="rem")
    m7 = cmulti(x$Y | x$D ~ 1 + x$C$TSSR + x$C$OD + x$C$OD2, type="rem")
    m8 = cmulti(x$Y | x$D ~ 1 + x$C$TSSR + x$C$TSSR2 + x$C$OD, type="rem")
    m9 = cmulti(x$Y | x$D ~ 1 + x$C$TSSR + x$C$TSSR2 + x$C$OD + x$C$OD2, type="rem")

    removal_list <- list(m1, m2, m3, m4, m5, m6, m7, m8, m9)
    save(removal_list, file = paste0("data/removal/", sp, ".rda"))
  }

stopCluster(cluster)
