####### Script Information ########################
# Brandon P.M. Edwards
# NA-POPS: analysis
# 3-QPAD-distance-model.R
# Created August 2020
# Last Updated January 2022

####### Import Libraries and External Files #######

library(detect)
library(plyr)
library(doParallel)
library(foreach)

####### Read Data #################################

load(file = "data/combined/dist_count_matrix.rda")
load(file = "data/combined/landcover_covariates.rda")
load(file = "data/combined/dist_design.rda")

n_cores <- 15

####### Data Wrangling ############################

# Drop most of the landcover covariates for now
covariates_reduced <- landcover_covariates[, c("Sample_ID", "roaddist", "ForestOnly_5x5")]
covariates_reduced$ForestOnly_5x5 <- covariates_reduced$ForestOnly_5x5 / 25

# Drop invalid road distances, and create roadside variable
covariates_reduced <- covariates_reduced[which(covariates_reduced$roaddist >= 0), ]
covariates_reduced$roadside <- ifelse(covariates_reduced$roaddist < 30, 1, 0)

# Drop method I
dist_count_matrix <- dist_count_matrix[-which(dist_count_matrix$Distance_Method == "I"), ]
dist_design <- dist_design[-which(dist_design$Method == "I"), ]

max_bands <- ncol(dist_design) - 2
count_names <- c("Sample_ID", "Species", "Distance_Method",
                 paste0(rep("Int", times = max_bands), 1:max_bands))
names(dist_count_matrix) <- count_names

design_names <- c("Distance_Method", "Max_Distance",
                  paste0(rep("Interval", times = max_bands), 1:max_bands))
names(dist_design) <- design_names

# Join data
count_design <- plyr::join(dist_count_matrix, dist_design,
                           by = "Distance_Method", type = "left")

# Filter out methods with only 1 interval
to_remove <- which(is.na(count_design$Interval2))
count_design <- count_design[-c(to_remove), ]

# Create separate data frames
counts <- count_design[,c("Sample_ID", "Species", "Distance_Method", count_names[4:length(count_names)])]
design <- count_design[,c("Sample_ID", "Species", "Distance_Method", design_names[3:length(design_names)])]
names(design) <- names(counts)
col_names <- names(counts)[4:length(names(counts))]

# Change 0s to NA in counts table where appropriate based on design table
for (i in col_names)
{
  indices <- which(counts[, i] == 0)
  
  counts[indices, i] <- ifelse(is.na(design[indices, i]), NA, 0)
}

# Add in landcover covariates to counts matrix
covars <- plyr::join(counts, covariates_reduced, by = "Sample_ID", type = "left", match = "all")

# Remove rows with NA covariates
to_remove <- which(is.na(covars$roadside))
counts <- counts[-c(to_remove), ]
design <- design[-c(to_remove), ]
covars <- covars[-c(to_remove), ]

to_remove <- which(is.na(covars$ForestOnly_5x5))
counts <- counts[-c(to_remove), ]
design <- design[-c(to_remove), ]
covars <- covars[-c(to_remove), ]

# Save which sample IDs (and therefore which covariates) were used during distance modelling
dis_covars_used <- covars[, c("Sample_ID", "ForestOnly_5x5", "roadside")]
dis_covars_used <- cbind(dis_covars_used, design$Distance_Method)
names(dis_covars_used)[ncol(dis_covars_used)] <- "Distance_Method"
dis_covars_used <- dis_covars_used[!duplicated(dis_covars_used$Sample_ID), ]
save(dis_covars_used, file = "data/combined/dis_covars_used.rda")

# Build matrices by species
species <- sort(as.character(unique(counts$Species)))

dis_species_summary <- vector(mode = "list", length = length(species))
names(dis_species_summary) <- species

for (s in species)
{
  if (nrow(counts[counts$Species == s, ]) >= 75)
  {
    # Save covariates being used
    covars_sp <- covars[which(covars$Species == s), ]
    dis_species_summary[[s]] <- covars_sp[, c("ForestOnly_5x5", "roadside", "Distance_Method")]
    
    # Create matrices for modelling
    counts_sp <- counts[which(counts$Species == s), col_names]
    design_sp <- design[which(design$Species == s), col_names]
    covars_sp <- covars_sp[, 18:19]
    
    assign(paste0("Y_",s), as.matrix(counts_sp))
    assign(paste0("D_",s), as.matrix(design_sp))
    assign(paste0("C_",s), covars_sp)
  }else
  {
    species <- species[!(species %in% s)]
  }
}

# Remove all empty species 
dis_species_summary[sapply(dis_species_summary, is.null)] <- NULL
save(dis_species_summary, file = "../results/quant-summary/dis_species_summary.rda")

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
    m1 = cmulti(x$Y | x$D ~ 1, type="dis")
    m2 = cmulti(x$Y | x$D ~ 1 + x$C$roadside, type="dis")
    m3 = cmulti(x$Y | x$D ~ 1 + x$C$ForestOnly_5x5, type="dis")
    m4 = cmulti(x$Y | x$D ~ 1 + x$C$roadside + x$C$ForestOnly_5x5, type="dis")
    m5 = cmulti(x$Y | x$D ~ 1 + x$C$roadside + x$C$ForestOnly_5x5 + x$C$roadside:x$C$ForestOnly_5x5, type="dis")
    distance_list <- list(m1, m2, m3, m4, m5)
    save(distance_list, file = paste0("data/distance/", sp, ".rda"))    
  }
stopCluster(cluster)
