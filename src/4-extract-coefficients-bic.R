####### Script Information ########################
# Brandon P.M. Edwards
# NA-POPS: analysis
# 4-extract-coefficients-bic.R
# Created November 2020
# Last Updated April 2021

####### Import Libraries and External Files #######

library(detect)

####### Set Constants #############################

# create function to sapply BIC to a list
bic_multi <- function(x)
{
  sapply(x, FUN = BIC)
}

n_rem_models <- 15
n_rem_bcr_models <- 6
n_rem_bcr_jd2_models <- 3
n_dis_models <- 5

####### Removal Model Coefficients ################

species <- substr(list.files(path = "data/removal"), 
                  start = 1, 
                  stop = 4)

rem_coef <- data.frame(Species = rep(species, each = n_rem_models),
                       n = NA, 
                       model = rep(seq(1,n_rem_models), times = length(species)), 
                       bic = NA,
                       intercept = NA,
                       tssr = NA,
                       tssr2 = NA,
                       jd = NA,
                       jd2 = NA)

# Using 40 here as it will cover all of the BCRs so far
bcr_df_names <- c(paste0(rep("BCR", times = 40), 1:40),
                  "Species", "model")
bcr_coef <- as.data.frame(matrix(nrow = length(species) * n_rem_bcr_models,
                                 ncol = length(bcr_df_names)))
bcr_jd_coef <- as.data.frame(matrix(nrow = length(species) * n_rem_bcr_models,
                                    ncol = length(bcr_df_names)))
bcr_jd2_coef <- as.data.frame(matrix(nrow = length(species) * n_rem_bcr_jd2_models,
                                     ncol = length(bcr_df_names)))

names(bcr_coef) <- names(bcr_jd_coef) <- names(bcr_jd2_coef) <- bcr_df_names
bcr_coef$Species <- bcr_jd_coef$Species <- rep(species, each = n_rem_bcr_models)
bcr_coef$model <- bcr_jd_coef$model <- rep(seq(10, 15), times = length(species))

bcr_jd2_coef$Species <- rep(species, each = n_rem_bcr_jd2_models)
bcr_jd2_coef$model <- rep(c(11, 14, 15), times = length(species))

rem_vcv_list <- vector(mode = "list", length = n_rem_models)
sp_list <- vector(mode = "list", length = length(species))
names(sp_list) <- species
for (m in 1:n_rem_models)
{
  rem_vcv_list[[m]] <- sp_list
}

rem_bic <- vector(mode = "list", length = length(species))
names(rem_bic) <- species

for (s in species)
{
  load(file = paste0("data/removal/", s, ".rda"))
  rem_coef[which(rem_coef$Species == s), "n"] <- nrow(removal_list[[1]]$Y)
  
  bic <- bic_multi(x = removal_list)
  rem_coef[which(rem_coef$Species == s), "bic"] <- bic
  
  bic_df <- data.frame(Model = seq(1, n_rem_models),
                       BIC = bic)
  bic_df <- bic_df[order(bic_df$BIC), ]
  bic_df$Delta_BIC <- bic_df$BIC - bic_df$BIC[1]
  rem_bic[[s]] <- bic_df
  
  for (m in 1:n_rem_models)
  {
    coef <- coef(removal_list[[m]])
    if (m == 1) {
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"intercept"] = coef[1]
    }
    if (m == 2) {
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"intercept"] = coef[1]
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"tssr"] = coef[2]
    }
    if (m == 3) {
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"intercept"] = coef[1]
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"jd"] = coef[2]
    }
    if (m == 4) {
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"intercept"] = coef[1]
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"tssr"] = coef[2]
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"tssr2"] = coef[3]
    }
    if (m == 5) {
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"intercept"] = coef[1]
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"jd"] = coef[2]
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"jd2"] = coef[3]
    }
    if (m == 6) {
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"intercept"] = coef[1]
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"tssr"] = coef[2]
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"jd"] = coef[3]
    }
    if (m == 7) {
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"intercept"] = coef[1]
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"tssr"] = coef[2]
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"tssr2"] = coef[3]
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"jd"] = coef[4]
    }
    if (m == 8) {
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"intercept"] = coef[1]
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"tssr"] = coef[2]
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"jd"] = coef[3]
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"jd2"] = coef[4]
    }
    if (m == 9) {
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"intercept"] = coef[1]
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"tssr"] = coef[2]
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"tssr2"] = coef[3]
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"jd"] = coef[4]
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"jd2"] = coef[5]
    }
    if (m == 10) {
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"intercept"] = coef[1]
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"jd"] = coef[2]
      
      #' Because we are modelling the BCRs as factors, the lowest BCR will get input
      #' into the model as the base factor, and we won't actually get a coefficient for
      #' it. Thus, this line will extract all BCRs input into the model, then later on
      #' we can add whichever one does not have a coefficient
      all_bcr_used <- as.numeric(as.character(levels(removal_list[[m]]$model$`x$C$BCR`)))
      
      #' First, get all the coefficients for the JD BCR interaction. Do this by
      #' finding which indices of the coefficient list have the string "JD:" in it,
      #' which only occur for the JD:BCR interactions. Then, get just those coefficients
      #' and save it to variable called bcr_jd_coef_list. To get the specific BCR
      #' which they belong to (which we will carry through the rest of the BCR parameter
      #' extractions), we use the sub() to get all that comes after "BCR", which is just
      #' the BCR number. Then, we can subset the over bcr_jd_coef dataframe by Species
      #' and model number in the row, and by BCR in the column, and slap the coefficients
      #' into the proper spot!
      bcr_jd_indices <- which(grepl("JD:", names(coef)))
      bcr_jd_coef_list <- coef[bcr_jd_indices]
      bcr_used <- as.integer(sub(".*BCR", "", names(bcr_jd_coef_list)))
      bcr_jd_coef[which(bcr_jd_coef$Species == s & bcr_jd_coef$model == m),
                  bcr_used] = coef[bcr_jd_indices]
      # Add a 0 for the base BCR, just so we can still track which BCRs were used
      leftover <- setdiff(all_bcr_used, bcr_used)
      bcr_jd_coef[which(bcr_jd_coef$Species == s & bcr_jd_coef$model == m),
                  leftover] = 0
      
      #' Now we want all parameters that are JUST the BCR, but NOT the interaction terms.
      #' So, get the indices of all those which contain BCR in the name (which will initially
      #' include the interaction terms), then get the differnce between that set and the
      #' set that only contain interactions, which will leave us with just the BCR coefficients.
      #' The rest of the operations follow similar to above.
      bcr_indices <- setdiff(which(grepl("BCR", names(coef))), bcr_jd_indices)
      bcr_coef_list <- coef[bcr_indices]
      bcr_used <- as.integer(sub(".*BCR", "", names(bcr_coef_list)))
      bcr_coef[which(bcr_coef$Species == s & bcr_coef$model == m),
               bcr_used] = coef[bcr_indices]
      leftover <- setdiff(all_bcr_used, bcr_used)
      bcr_coef[which(bcr_coef$Species == s & bcr_coef$model == m),
               leftover] = 0
    }
    if (m == 11) {
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"intercept"] = coef[1]
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"jd"] = coef[2]
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"jd2"] = coef[3]
      
      #' Similar to above, but we first will get the JD2 interactions
      all_bcr_used <- as.numeric(as.character(levels(removal_list[[m]]$model$`x$C$BCR`)))
      
      bcr_jd2_indices <- which(grepl("JD2:", names(coef)))
      bcr_jd2_coef_list <- coef[bcr_jd2_indices]
      bcr_used <- as.integer(sub(".*BCR", "", names(bcr_jd2_coef_list)))
      bcr_jd2_coef[which(bcr_jd2_coef$Species == s & bcr_jd2_coef$model == m),
                   bcr_used] = coef[bcr_jd2_indices]
      leftover <- setdiff(all_bcr_used, bcr_used)
      bcr_jd2_coef[which(bcr_jd2_coef$Species == s & bcr_jd2_coef$model == m),
                   leftover] = 0
      
      #' Now just JD interactions
      bcr_jd_indices <- which(grepl("JD:", names(coef)))
      bcr_jd_coef_list <- coef[bcr_jd_indices]
      bcr_used <- as.integer(sub(".*BCR", "", names(bcr_jd_coef_list)))
      bcr_jd_coef[which(bcr_jd_coef$Species == s & bcr_jd_coef$model == m),
                  bcr_used] = coef[bcr_jd_indices]
      leftover <- setdiff(all_bcr_used, bcr_used)
      bcr_jd_coef[which(bcr_jd_coef$Species == s & bcr_jd_coef$model == m),
                  leftover] = 0
      
      #' Now, we will do setdiff again, only instead we need to take the setdiff
      #' of the BCR set, and the union of the set of JD and JD2 indices
      bcr_indices <- setdiff(which(grepl("BCR", names(coef))),
                             union(bcr_jd_indices, bcr_jd2_indices))
      bcr_coef_list <- coef[bcr_indices]
      bcr_used <- as.integer(sub(".*BCR", "", names(bcr_coef_list)))
      bcr_coef[which(bcr_coef$Species == s & bcr_coef$model == m),
               bcr_used] = coef[bcr_indices]
      leftover <- setdiff(all_bcr_used, bcr_used)
      bcr_coef[which(bcr_coef$Species == s & bcr_coef$model == m),
               leftover] = 0
    }
    if (m == 12) {
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"intercept"] = coef[1]
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"tssr"] = coef[2]
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"jd"] = coef[3]
      
      # See explanation above of the following code 
      all_bcr_used <- as.numeric(as.character(levels(removal_list[[m]]$model$`x$C$BCR`)))
      
      bcr_jd_indices <- which(grepl("JD:", names(coef)))
      bcr_jd_coef_list <- coef[bcr_jd_indices]
      bcr_used <- as.integer(sub(".*BCR", "", names(bcr_jd_coef_list)))
      bcr_jd_coef[which(bcr_jd_coef$Species == s & bcr_jd_coef$model == m),
                  bcr_used] = coef[bcr_jd_indices]
      leftover <- setdiff(all_bcr_used, bcr_used)
      bcr_jd_coef[which(bcr_jd_coef$Species == s & bcr_jd_coef$model == m),
                  leftover] = 0
      
      bcr_indices <- setdiff(which(grepl("BCR", names(coef))), bcr_jd_indices)
      bcr_coef_list <- coef[bcr_indices]
      bcr_used <- as.integer(sub(".*BCR", "", names(bcr_coef_list)))
      bcr_coef[which(bcr_coef$Species == s & bcr_coef$model == m),
               bcr_used] = coef[bcr_indices]
      leftover <- setdiff(all_bcr_used, bcr_used)
      bcr_coef[which(bcr_coef$Species == s & bcr_coef$model == m),
               leftover] = 0
    }
    if (m == 13) {
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"intercept"] = coef[1]
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"tssr"] = coef[2]
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"tssr2"] = coef[3]
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"jd"] = coef[4]
      
      # See explanation above of the following code 
      all_bcr_used <- as.numeric(as.character(levels(removal_list[[m]]$model$`x$C$BCR`)))
      
      bcr_jd_indices <- which(grepl("JD:", names(coef)))
      bcr_jd_coef_list <- coef[bcr_jd_indices]
      bcr_used <- as.integer(sub(".*BCR", "", names(bcr_jd_coef_list)))
      bcr_jd_coef[which(bcr_jd_coef$Species == s & bcr_jd_coef$model == m),
                  bcr_used] = coef[bcr_jd_indices]
      leftover <- setdiff(all_bcr_used, bcr_used)
      bcr_jd_coef[which(bcr_jd_coef$Species == s & bcr_jd_coef$model == m),
                  leftover] = 0
      
      bcr_indices <- setdiff(which(grepl("BCR", names(coef))), bcr_jd_indices)
      bcr_coef_list <- coef[bcr_indices]
      bcr_used <- as.integer(sub(".*BCR", "", names(bcr_coef_list)))
      bcr_coef[which(bcr_coef$Species == s & bcr_coef$model == m),
               bcr_used] = coef[bcr_indices]
      leftover <- setdiff(all_bcr_used, bcr_used)
      bcr_coef[which(bcr_coef$Species == s & bcr_coef$model == m),
               leftover] = 0
    }
    if (m == 14) {
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"intercept"] = coef[1]
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"tssr"] = coef[2]
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"jd"] = coef[3]
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"jd2"] = coef[4]
      
      #' See above for explanation
      all_bcr_used <- as.numeric(as.character(levels(removal_list[[m]]$model$`x$C$BCR`)))
      
      bcr_jd2_indices <- which(grepl("JD2:", names(coef)))
      bcr_jd2_coef_list <- coef[bcr_jd2_indices]
      bcr_used <- as.integer(sub(".*BCR", "", names(bcr_jd2_coef_list)))
      bcr_jd2_coef[which(bcr_jd2_coef$Species == s & bcr_jd2_coef$model == m),
                   bcr_used] = coef[bcr_jd2_indices]
      leftover <- setdiff(all_bcr_used, bcr_used)
      bcr_jd2_coef[which(bcr_jd2_coef$Species == s & bcr_jd2_coef$model == m),
                   leftover] = 0
      
      bcr_jd_indices <- which(grepl("JD:", names(coef)))
      bcr_jd_coef_list <- coef[bcr_jd_indices]
      bcr_used <- as.integer(sub(".*BCR", "", names(bcr_jd_coef_list)))
      bcr_jd_coef[which(bcr_jd_coef$Species == s & bcr_jd_coef$model == m),
                  bcr_used] = coef[bcr_jd_indices]
      leftover <- setdiff(all_bcr_used, bcr_used)
      bcr_jd_coef[which(bcr_jd_coef$Species == s & bcr_jd_coef$model == m),
                  leftover] = 0
      
      bcr_indices <- setdiff(which(grepl("BCR", names(coef))),
                             union(bcr_jd_indices, bcr_jd2_indices))
      bcr_coef_list <- coef[bcr_indices]
      bcr_used <- as.integer(sub(".*BCR", "", names(bcr_coef_list)))
      bcr_coef[which(bcr_coef$Species == s & bcr_coef$model == m),
               bcr_used] = coef[bcr_indices]
      leftover <- setdiff(all_bcr_used, bcr_used)
      bcr_coef[which(bcr_coef$Species == s & bcr_coef$model == m),
               leftover] = 0
    }
    if (m == 15) {
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"intercept"] = coef[1]
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"tssr"] = coef[2]
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"tssr2"] = coef[3]
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"jd"] = coef[4]
      rem_coef[which(rem_coef$Species == s & rem_coef$model == m),"jd2"] = coef[5]
      
      #' See above for explanation
      all_bcr_used <- as.numeric(as.character(levels(removal_list[[m]]$model$`x$C$BCR`)))
      
      bcr_jd2_indices <- which(grepl("JD2:", names(coef)))
      bcr_jd2_coef_list <- coef[bcr_jd2_indices]
      bcr_used <- as.integer(sub(".*BCR", "", names(bcr_jd2_coef_list)))
      bcr_jd2_coef[which(bcr_jd2_coef$Species == s & bcr_jd2_coef$model == m),
                   bcr_used] = coef[bcr_jd2_indices]
      leftover <- setdiff(all_bcr_used, bcr_used)
      bcr_jd2_coef[which(bcr_jd2_coef$Species == s & bcr_jd2_coef$model == m),
                   leftover] = 0
      
      bcr_jd_indices <- which(grepl("JD:", names(coef)))
      bcr_jd_coef_list <- coef[bcr_jd_indices]
      bcr_used <- as.integer(sub(".*BCR", "", names(bcr_jd_coef_list)))
      bcr_jd_coef[which(bcr_jd_coef$Species == s & bcr_jd_coef$model == m),
                  bcr_used] = coef[bcr_jd_indices]
      leftover <- setdiff(all_bcr_used, bcr_used)
      bcr_jd_coef[which(bcr_jd_coef$Species == s & bcr_jd_coef$model == m),
                  leftover] = 0
      
      bcr_indices <- setdiff(which(grepl("BCR", names(coef))),
                             union(bcr_jd_indices, bcr_jd2_indices))
      bcr_coef_list <- coef[bcr_indices]
      bcr_used <- as.integer(sub(".*BCR", "", names(bcr_coef_list)))
      bcr_coef[which(bcr_coef$Species == s & bcr_coef$model == m),
               bcr_used] = coef[bcr_indices]
      leftover <- setdiff(all_bcr_used, bcr_used)
      bcr_coef[which(bcr_coef$Species == s & bcr_coef$model == m),
               leftover] = 0
    }
    rem_vcv_list[[m]][[s]] <- removal_list[[m]]$vcov
  }
}

write.table(rem_coef, 
            file = "../results/coefficients/removal.csv",
            sep = ",",
            row.names = FALSE)
write.table(bcr_coef, 
            file = "../results/coefficients/removal_bcr.csv",
            sep = ",",
            row.names = FALSE)
write.table(bcr_jd_coef, 
            file = "../results/coefficients/removal_bcr_jd.csv",
            sep = ",",
            row.names = FALSE)
write.table(bcr_jd2_coef, 
            file = "../results/coefficients/removal_bcr_jd2.csv",
            sep = ",",
            row.names = FALSE)

save(rem_vcv_list, file = "../results/var-covar/rem_vcv_list.rda")
save(rem_bic, file = "../results/bic/rem_bic.rda")

####### Distance Model Selection ##################

species <- substr(list.files(path = "data/distance"), 
                  start = 1, 
                  stop = 4)

dist_coef <- data.frame(Species = rep(species, each = n_dis_models),
                        n = NA, 
                        model = rep(seq(1,n_dis_models), times = length(species)), 
                        bic = NA,
                        intercept = NA,
                        road = NA,
                        forest = NA,
                        roadforest = NA)

dis_vcv_list <- vector(mode = "list", length = n_dis_models)
sp_list <- vector(mode = "list", length = length(species))
names(sp_list) <- species
for (m in 1:n_dis_models)
{
  dis_vcv_list[[m]] <- sp_list
}

dis_bic <- vector(mode = "list", length = length(species))
names(dis_bic) <- species

for (s in species)
{
  load(file = paste0("data/distance/", s, ".rda"))
  dist_coef[which(dist_coef$Species == s), "n"] <- nrow(distance_list[[1]]$Y)
  
  bic <- bic_multi(x = distance_list)
  dist_coef[which(dist_coef$Species == s), "bic"] <- bic
  
  bic_df <- data.frame(Model = seq(1, n_dis_models),
                       BIC = bic)
  bic_df <- bic_df[order(bic_df$BIC), ]
  bic_df$Delta_BIC <- bic_df$BIC - bic_df$BIC[1]
  dis_bic[[s]] <- bic_df
  
  for (m in 1:n_dis_models)
  {
    coef <- coef(distance_list[[m]])
    if (m == 1) {
      dist_coef[which(dist_coef$Species == s & dist_coef$model == m),"intercept"] = coef[1]
    }
    if (m == 2) {
      dist_coef[which(dist_coef$Species == s & dist_coef$model == m),"intercept"] = coef[1]
      dist_coef[which(dist_coef$Species == s & dist_coef$model == m),"road"] = coef[2]
    }
    if (m == 3) {
      dist_coef[which(dist_coef$Species == s & dist_coef$model == m),"intercept"] = coef[1]
      dist_coef[which(dist_coef$Species == s & dist_coef$model == m),"forest"] = coef[2]
    }
    if (m == 4) {
      dist_coef[which(dist_coef$Species == s & dist_coef$model == m),"intercept"] = coef[1]
      dist_coef[which(dist_coef$Species == s & dist_coef$model == m),"road"] = coef[2]
      dist_coef[which(dist_coef$Species == s & dist_coef$model == m),"forest"] = coef[3]
    }
    if (m == 5) {
      dist_coef[which(dist_coef$Species == s & dist_coef$model == m),"intercept"] = coef[1]
      dist_coef[which(dist_coef$Species == s & dist_coef$model == m),"road"] = coef[2]
      dist_coef[which(dist_coef$Species == s & dist_coef$model == m),"forest"] = coef[3]
      dist_coef[which(dist_coef$Species == s & dist_coef$model == m),"roadforest"] = coef[4]
    }    
    
    dis_vcv_list[[m]][[s]] <- distance_list[[m]]$vcov
  }
}

write.table(dist_coef, 
            file = "../results/coefficients/distance.csv",
            sep = ",",
            row.names = FALSE)

save(dis_vcv_list, file = "../results/var-covar/dis_vcv_list.rda")
save(dis_bic, file = "../results/bic/dis_bic.rda")