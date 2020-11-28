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
load(file = here::here("data/landcover_covariates.rda"))
load(file = here::here("data/dist_design.rda"))

####### Data Wrangling ############################

# Drop most of the landcover covariates for now
covariates_reduced <- landcover_covariates[, c("Sample_ID", "roaddist", "ForestOnly_5x5")]
covariates_reduced$ForestOnly_5x5 <- covariates_reduced$ForestOnly_5x5 / 25

# # Drop method I
# dist_count_matrix <- dist_count_matrix[-which(dist_count_matrix$Distance_Method == "I"), ]
# dist_design <- dist_design[-which(dist_design$Method == "I"), ]

# WARNING, UGLY HARD CODED MESS! FIX PRIOR TO FULL PUBLICATION :)
dist_count_matrix$"7" <- rep(0, nrow(dist_count_matrix))

count_names <- c("Sample_ID", "Species", "Distance_Method",
                 paste0(rep("Int", times = 12), 1:12))
names(dist_count_matrix) <- count_names

design_names <- c("Distance_Method", "Max_Distance",
                  paste0(rep("Interval", times = 12), 1:12))
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
                                  select = c(16:17)))
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

start_time <- Sys.time()
distance_output_list <- parLapply(cl = cluster,
                                  X = input_list,
                                  fun = graceful_mm)
end_time <- Sys.time()
elapsed_time <- end_time - start_time

stopCluster(cluster)

save(distance_output_list, file = "data/distance_output_list.rda")
save(elapsed_time, file = "data/elapsed_time_distance.rda")

########### Model Selection #######################

# create function to sapply BIC to a list
bic_multi <- function(x)
{
  tryCatch(sapply(x, FUN = BIC),
           error = function(e) NA)
}

#estimate BIC values for all models across all species
bic_list <- lapply(distance_output_list, FUN = bic_multi)
names(bic_list) <- species

#Create dataframe to contain coefficients
dist_coef <- data.frame(Species=species, d_model=NA, d_n=NA, d_int=NA,d_roaddist=NA,d_forest=NA,d_roaddistforest=NA, logtau = NA)

#Extracts model coefficients based on top model in BIC tables and assigns them to the correct column
for (i in species) {
  if (is.na(bic_list[[i]]))
  {
    next
  }
  model = which(bic_list[[i]]==min(bic_list[[i]])) #determine which model has lowest BIC
  coef = coef(distance_output_list[[i]][[model]]) #extract coefficients from best model
  if (is.null(coef))
  {
    next
  }
  dist_coef[which(dist_coef$Species==i),"d_model"]=model
  dist_coef[which(dist_coef$Species==i),"d_n"]=nrow(get(paste0("Y_",i)))
  if (model==1) {
    dist_coef[which(dist_coef$Species==i),"d_int"] = coef[1]
  }
  if (model==2) {
    dist_coef[which(dist_coef$Species==i),"d_int"] = coef[1]
    dist_coef[which(dist_coef$Species==i),"d_roaddist"] = coef[2]
  }
  if (model==3) {
    dist_coef[which(dist_coef$Species==i),"d_int"] = coef[1]
    dist_coef[which(dist_coef$Species==i),"d_forest"] = coef[2]
  }
  if (model==4) {
    dist_coef[which(dist_coef$Species==i),"d_int"] = coef[1]
    dist_coef[which(dist_coef$Species==i),"d_roaddist"] = coef[2]
    dist_coef[which(dist_coef$Species==i),"d_forest"] = coef[3]
  }
  if (model==5) {
    dist_coef[which(dist_coef$Species==i),"d_int"] = coef[1]
    dist_coef[which(dist_coef$Species==i),"d_roaddist"] = coef[2]
    dist_coef[which(dist_coef$Species==i),"d_forest"] = coef[3]
    dist_coef[which(dist_coef$Species==i),"d_roaddistforest"] = coef[4]
  }
  rm(coef,model)
}

write.table(dist_coef, 
            file = "../results/detection-coefs/distance_coef.csv",
            sep = ",",
            row.names = FALSE)
