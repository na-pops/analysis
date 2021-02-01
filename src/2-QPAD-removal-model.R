####### Script Information ########################
# Brandon P.M. Edwards
# NA-POPS: analysis
# 2-QPAD-removal-model.R
# Created August 2020
# Last Updated February 2021

####### Import Libraries and External Files #######

library(detect)
library(plyr)
library(doParallel)
library(foreach)

####### Read Data #################################

load(file = here::here("data/time_count_matrix.rda"))
load(file = here::here("data/temporal_covariates.rda"))
load(file = here::here("data/time_design.rda"))
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

#' Abundance of species x during time band j in sampling event i
#' This is effectively our original Y matrix, flattened into one dimension,
#' with all NAs removed.
abund_per_time_band <- as.vector(t(Y))
abund_per_time_band <- abund_per_time_band[!is.na(abund_per_time_band)]

#' Total species abundance per sampling event.
#' I.e., this is the sum of Y_sij over j
total_abund_per_sample <- unname(apply(Y, 1, function(x) sum(x, na.rm = TRUE)))

#' For each sampling event, how many time bands were in that sampling event.
#' In index terms, this is J for each i
#' Length of this = TOTAL SAMPLES BEING CONSIDERED
time_bands_per_sample <- unname(apply(Y, 1, function(x) sum(!is.na(x))))

#' Total samples per species
#' I.e., total i BY SPECIES
samples_per_species <- as.vector(table(sp_list[,1]))

#' Total counts per species.
#' I.e., total i x j BY SPECIES
sp_list_count <- rep(sp_list[,1], times = time_bands_per_sample)
count_per_species <- as.vector(table(sp_list_count))

max_time <- as.vector(t(D))
max_time <- max_time[!is.na(max_time)]

X_names <- names(C)
X <- cbind(rep(1, nrow(C)), C)
names(X) <- c("Intercept", X_names)








input_list <- vector(mode="list", length=length(species))
names(input_list) <- species
for (s in species)
{
  input_list[[s]] <- list(Y=get(paste0("Y_",s)), D=get(paste0("D_",s)), C=get(paste0("C_",s)))
}

########### Modelling #############################

cluster <- makeCluster(3, type = "PSOCK")
registerDoParallel(cluster)

foreach(sp = names(input_list), .packages = 'detect') %dopar%
  {
    x <- input_list[[sp]]
    m1 = cmulti(x$Y | x$D ~ 1, type="rem")
    m2 = cmulti(x$Y | x$D ~ x$C$TSSR, type="rem")
    m3 = cmulti(x$Y | x$D ~ x$C$JD, type="rem")
    m4 = cmulti(x$Y | x$D ~ x$C$TSSR + x$C$TSSR2, type = "rem")
    m5 = cmulti(x$Y | x$D ~ x$C$JD + x$C$JD2, type = "rem")
    m6 = cmulti(x$Y | x$D ~ x$C$TSSR + x$C$JD, type="rem")
    m7 = cmulti(x$Y | x$D ~ x$C$TSSR + x$C$TSSR2 + x$C$JD, type="rem")
    m8 = cmulti(x$Y | x$D ~ x$C$TSSR + x$C$JD + x$C$JD2, type="rem")
    m9 = cmulti(x$Y | x$D ~ x$C$TSSR + x$C$TSSR2 + x$C$JD + x$C$JD2, type="rem")
    removal_list <- list(m1, m2, m3, m4, m5, m6, m7, m8, m9)
    save(removal_list, file = paste0("data/removal/", sp, ".rda"))    
  }

stopCluster(cluster)
