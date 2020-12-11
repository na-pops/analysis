####### Script Information ########################
# Brandon P.M. Edwards
# NA-POPS: analysis
# 4-model-selection.R
# Created November 2020
# Last Updated November 2020

####### Import Libraries and External Files #######

####### Set Constants #############################

# create function to sapply BIC to a list
bic_multi <- function(x)
{
  sapply(x, FUN = BIC)
}

####### Removal Model Selection ###################

species <- substr(list.files(path = "data/removal"), 
                  start = 1, 
                  stop = 4)

rem_coef <- data.frame(Species = species,
                       n = NA, 
                       best_model = NA, 
                       best_bic = NA,
                       int_best = NA,
                       tssr_best = NA,
                       tssr2_best = NA,
                       jd_best = NA,
                       jd2_best = NA,
                       full_bic = NA,
                       int_full = NA,
                       tssr_full = NA,
                       tssr2_full = NA,
                       jd_full = NA,
                       jd2_full = NA)

for (s in species)
{
  load(file = paste0("data/removal/", s, ".rda"))
  
  rem_coef[which(rem_coef$Species == s), "n"] <- nrow(removal_list[[1]]$Y)
  
  bic <- bic_multi(x = removal_list)
  min_bic <- min(bic)
  rem_coef[which(rem_coef$Species == s), "best_bic"] <- min_bic
  
  model <- which(bic == min_bic)
  rem_coef[which(rem_coef$Species == s), "best_model"] <- model
  
  coef <- coef(removal_list[[model]])
  if (model == 1) {
    rem_coef[which(rem_coef$Species == s),"int_best"] = coef[1]
  }
  if (model == 2) {
    rem_coef[which(rem_coef$Species == s),"int_best"] = coef[1]
    rem_coef[which(rem_coef$Species == s),"tssr_best"] = coef[2]
  }
  if (model == 3) {
    rem_coef[which(rem_coef$Species == s),"int_best"] = coef[1]
    rem_coef[which(rem_coef$Species == s),"jd_best"] = coef[2]
  }
  if (model == 4) {
    rem_coef[which(rem_coef$Species == s),"int_best"] = coef[1]
    rem_coef[which(rem_coef$Species == s),"tssr_best"] = coef[2]
    rem_coef[which(rem_coef$Species == s),"tssr2_best"] = coef[3]
  }
  if (model == 5) {
    rem_coef[which(rem_coef$Species == s),"int_best"] = coef[1]
    rem_coef[which(rem_coef$Species == s),"jd_best"] = coef[2]
    rem_coef[which(rem_coef$Species == s),"jd2_best"] = coef[3]
  }
  if (model == 6) {
    rem_coef[which(rem_coef$Species == s),"int_best"] = coef[1]
    rem_coef[which(rem_coef$Species == s),"tssr_best"] = coef[2]
    rem_coef[which(rem_coef$Species == s),"jd_best"] = coef[3]
  }
  if (model == 7) {
    rem_coef[which(rem_coef$Species == s),"int_best"] = coef[1]
    rem_coef[which(rem_coef$Species == s),"tssr_best"] = coef[2]
    rem_coef[which(rem_coef$Species == s),"tssr2_best"] = coef[3]
    rem_coef[which(rem_coef$Species == s),"jd_best"] = coef[4]
  }
  if (model == 8) {
    rem_coef[which(rem_coef$Species == s),"int_best"] = coef[1]
    rem_coef[which(rem_coef$Species == s),"tssr_best"] = coef[2]
    rem_coef[which(rem_coef$Species == s),"jd_best"] = coef[3]
    rem_coef[which(rem_coef$Species == s),"jd2_best"] = coef[4]
  }
  if (model == 9) {
    rem_coef[which(rem_coef$Species == s),"int_best"] = coef[1]
    rem_coef[which(rem_coef$Species == s),"tssr_best"] = coef[2]
    rem_coef[which(rem_coef$Species == s),"tssr2_best"] = coef[3]
    rem_coef[which(rem_coef$Species == s),"jd_best"] = coef[4]
    rem_coef[which(rem_coef$Species == s),"jd2_best"] = coef[5]
  }
  
  rem_coef[which(rem_coef$Species == s), "full_bic"] <- bic[9]
  coef_full <- coef(removal_list[[9]])
  rem_coef[which(rem_coef$Species == s),"int_full"] = coef_full[1]
  rem_coef[which(rem_coef$Species == s),"tssr_full"] = coef_full[2]
  rem_coef[which(rem_coef$Species == s),"tssr2_full"] = coef_full[3]
  rem_coef[which(rem_coef$Species == s),"jd_full"] = coef_full[4]
  rem_coef[which(rem_coef$Species == s),"jd2_full"] = coef_full[5]
}

write.table(rem_coef, 
            file = "../results/coefficients/removal.csv",
            sep = ",",
            row.names = FALSE)

####### Distance Model Selection ##################

species <- substr(list.files(path = "data/distance"), 
                  start = 1, 
                  stop = 4)

dist_coef <- data.frame(Species = species,
                        n = NA, 
                        best_model = NA, 
                        best_bic = NA,
                        int_best = NA,
                        road_best = NA,
                        forest_best = NA,
                        roadforest_best = NA,
                        full_bic = NA,
                        int_full = NA,
                        road_full = NA,
                        forest_full = NA,
                        roadforest_full = NA)

for (s in species)
{
  load(file = paste0("data/distance/", s, ".rda"))
  
  dist_coef[which(dist_coef$Species == s), "n"] <- nrow(distance_list[[1]]$Y)
  
  bic <- bic_multi(x = distance_list)
  min_bic <- min(bic)
  dist_coef[which(dist_coef$Species == s), "best_bic"] <- min_bic
  
  model <- which(bic == min_bic)
  dist_coef[which(dist_coef$Species == s), "best_model"] <- model
  
  coef <- coef(distance_list[[model]])
  if (model == 1) {
    dist_coef[which(dist_coef$Species == s),"int_best"] = coef[1]
  }
  if (model == 2) {
    dist_coef[which(dist_coef$Species == s),"int_best"] = coef[1]
    dist_coef[which(dist_coef$Species == s),"road_best"] = coef[2]
  }
  if (model == 3) {
    dist_coef[which(dist_coef$Species == s),"int_best"] = coef[1]
    dist_coef[which(dist_coef$Species == s),"forest_best"] = coef[2]
  }
  if (model == 4) {
    dist_coef[which(dist_coef$Species == s),"int_best"] = coef[1]
    dist_coef[which(dist_coef$Species == s),"road_best"] = coef[2]
    dist_coef[which(dist_coef$Species == s),"forest_best"] = coef[3]
  }
  if (model == 5) {
    dist_coef[which(dist_coef$Species == s),"int_best"] = coef[1]
    dist_coef[which(dist_coef$Species == s),"road_best"] = coef[2]
    dist_coef[which(dist_coef$Species == s),"forest_best"] = coef[3]
    dist_coef[which(dist_coef$Species == s),"roadforest_best"] = coef[4]
  }

  
  dist_coef[which(dist_coef$Species == s), "full_bic"] <- bic[5]
  coef_full <- coef(distance_list[[5]])
  dist_coef[which(dist_coef$Species == s),"int_full"] = coef_full[1]
  dist_coef[which(dist_coef$Species == s),"road_full"] = coef_full[2]
  dist_coef[which(dist_coef$Species == s),"forest_full"] = coef_full[3]
  dist_coef[which(dist_coef$Species == s),"roadforest_full"] = coef_full[4]
}

write.table(dist_coef, 
            file = "../results/coefficients/distance.csv",
            sep = ",",
            row.names = FALSE)
