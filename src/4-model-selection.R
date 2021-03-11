####### Script Information ########################
# Brandon P.M. Edwards
# NA-POPS: analysis
# 4-model-selection.R
# Created November 2020
# Last Updated March 2021

####### Import Libraries and External Files #######

library(detect)

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

rem_vcv_best <- vector(mode = "list", length = length(species))
rem_vcv_full <- vector(mode = "list", length = length(species))
names(rem_vcv_best) <- species
names(rem_vcv_full) <- species

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
  
  vcv_best <- rem_vcv_best[[s]] <- removal_list[[model]]$vcov
  write.table(vcv_best,
              file = paste0("../results/var-covar/removal/",
                            s,
                            "_best.csv"),
              sep = ",", row.names = FALSE)
  
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
  
  vcv_full <- rem_vcv_full[[s]] <- removal_list[[9]]$vcov
  write.table(vcv_full,
              file = paste0("../results/var-covar/removal/",
                            s,
                            "_full.csv"),
              sep = ",", row.names = FALSE)
}

write.table(rem_coef, 
            file = "../results/coefficients/removal.csv",
            sep = ",",
            row.names = FALSE)

save(rem_vcv_best, file = "../results/var-covar/rem_vcv_best.rda")
save(rem_vcv_full, file = "../results/var-covar/rem_vcv_full.rda")

####### Distance Model Selection ##################

species <- substr(list.files(path = "data/distance"), 
                  start = 1, 
                  stop = 4)

dis_vcv_best <- vector(mode = "list", length = length(species))
dis_vcv_full <- vector(mode = "list", length = length(species))
names(dis_vcv_best) <- species
names(dis_vcv_full) <- species

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
  
  vcv_best <- dis_vcv_best[[s]] <- distance_list[[model]]$vcov
  write.table(vcv_best,
              file = paste0("../results/var-covar/distance/",
                            s,
                            "_best.csv"),
              sep = ",", row.names = FALSE)
  
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
  
  vcv_full <- dis_vcv_full[[s]] <- distance_list[[5]]$vcov
  write.table(vcv_full,
              file = paste0("../results/var-covar/distance/",
                            s,
                            "_full.csv"),
              sep = ",", row.names = FALSE)
}

write.table(dist_coef, 
            file = "../results/coefficients/distance.csv",
            sep = ",",
            row.names = FALSE)

save(dis_vcv_best, file = "../results/var-covar/dis_vcv_best.rda")
save(dis_vcv_full, file = "../results/var-covar/dis_vcv_full.rda")