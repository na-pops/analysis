####### Script Information ########################
# Brandon P.M. Edwards
# NA-POPS: analysis
# 3-QPAD-distance-model.R
# Created August 2020
# Last Updated August 2020

####### Import Libraries and External Files #######

library(detect)
library(plyr)
library(parallel)

####### Read Data #################################

load(file = here::here("data/dist_count_matrix.rda"))
load(file = here::here("data/covariates.rda"))
load(file = here::here("data/dist_design.rda"))

####### Data Wrangling ############################

# Drop most of the landcover covariates for now
covariates_reduced <- covariates[, c("Sample_ID", "roaddist", "ForestOnly_5x5")]
covariates_reduced$ForestOnly_5x5 <- covariates_reduced$ForestOnly_5x5 / 25

# Drop method I
dist_count_matrix <- dist_count_matrix[-which(dist_count_matrix$Distance_Method == "I"), ]
dist_design <- dist_design[-which(dist_design$Method == "I"), ]

# WARNING, UGLY HARD CODED MESS! FIX PRIOR TO FULL PUBLICATION :)
dist_count_matrix$"7" <- rep(0, nrow(dist_count_matrix))

count_names <- c("Sample_ID", "Species", "Distance_Method",
                 paste0(rep("Int", times = 7), 1:7))
names(dist_count_matrix) <- count_names

design_names <- c("Distance_Method", "Max_Distance",
                  paste0(rep("Interval", times = 7), 1:7))
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
to_remove <- which(is.na(covars$roaddist))
counts <- counts[-c(to_remove), ]
design <- design[-c(to_remove), ]
covars <- covars[-c(to_remove), ]

to_remove <- which(is.na(covars$ForestOnly_5x5))
counts <- counts[-c(to_remove), ]
design <- design[-c(to_remove), ]
covars <- covars[-c(to_remove), ]

# Build matrices by species :D
species <- sort(as.character(unique(counts$Species)))

for (s in species)
{
  if (nrow(counts[counts$Species == s, ]) >= 50)
  {
    assign(paste0("Y_",s), as.matrix(counts[counts$Species==s, col_names]))
    assign(paste0("D_",s), as.matrix(design[design$Species==s, col_names]))
    assign(paste0("C_",s), subset(covars,
                                  Species==s,
                                  select = c(11:12)))
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

# Function to do multiple cmulti bois
multi_multi <- function(x) {
  
  require(detect)
  m1 = cmulti(x$Y | x$D ~ 1, type="dis")
  m2 = cmulti(x$Y | x$D ~ x$C$roaddist, type="dis")
  m3 = cmulti(x$Y | x$D ~ x$C$ForestOnly_5x5, type="dis")
  m4 = cmulti(x$Y | x$D ~ x$C$roaddist + x$C$ForestOnly_5x5, type="dis")
  m5 = cmulti(x$Y | x$D ~ x$C$roaddist * x$C$ForestOnly_5x5, type="dis")
  return(list(m1,m2,m3,m4,m5))
}

graceful_mm <- function(x)
{
  to_return <- NA
  tryCatch({to_return <- multi_multi(x);}, 
           error = function(e) {to_return <- NA});
  
  return(to_return)
}

cluster <- makeCluster(detectCores() - 1)
clusterEvalQ(cluster, library(detect))
clusterExport(cluster, "input_list"); clusterExport(cluster, "multi_multi"); clusterExport(cluster, "graceful_mm")
distance_output_list <- parLapply(cl = cluster,
                                  X = input_list,
                                  fun = graceful_mm)
stopCluster(cluster)

save(distance_output_list, file = "data/distance_output_list.rda")

