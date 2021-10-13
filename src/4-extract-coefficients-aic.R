####### Script Information ########################
# Brandon P.M. Edwards
# NA-POPS: analysis
# 4-extract-coefficients-aic.R
# Created November 2020
# Last Updated October 2021

####### Import Libraries and External Files #######

library(detect)

####### Set Constants #############################

# create function to sapply AIC to a list
aic_multi <- function(x)
{
  sapply(x, FUN = AIC)
}

n_rem_models <- 9
n_dis_models <- 5

####### Removal Model Coefficients ################

species <- substr(list.files(path = "data/removal"), 
                  start = 1, 
                  stop = 4)

rem_coef <- data.frame(Species = rep(species, each = n_rem_models),
                       n = NA, 
                       model = rep(seq(1,n_rem_models), times = length(species)), 
                       aic = NA,
                       intercept = NA,
                       tssr = NA,
                       tssr2 = NA,
                       jd = NA,
                       jd2 = NA)

rem_vcv_list <- vector(mode = "list", length = n_rem_models)
sp_list <- vector(mode = "list", length = length(species))
names(sp_list) <- species
for (m in 1:n_rem_models)
{
  rem_vcv_list[[m]] <- sp_list
}

rem_aic <- vector(mode = "list", length = length(species))
names(rem_aic) <- species

for (s in species)
{
  load(file = paste0("data/removal/", s, ".rda"))
  n_mod_sp <- length(removal_list)
  
  rem_coef[which(rem_coef$Species == s & 
                   rem_coef$model <= n_mod_sp), "n"] <- nrow(removal_list[[1]]$Y)
  
  aic <- aic_multi(x = removal_list)
  rem_coef[which(rem_coef$Species == s & 
                   rem_coef$model <= n_mod_sp), "aic"] <- aic
  
  aic_df <- data.frame(Model = seq(1, 9),
                       AIC = aic[1:9])
  aic_df <- aic_df[order(aic_df$AIC), ]
  aic_df$Delta_AIC <- aic_df$AIC - aic_df$AIC[1]
  rem_aic[[s]] <- aic_df
  
  for (m in 1:n_mod_sp)
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
    
    rem_vcv_list[[m]][[s]] <- removal_list[[m]]$vcov
  }
}

write.table(rem_coef, 
            file = "../results/coefficients/removal.csv",
            sep = ",",
            row.names = FALSE)

save(rem_vcv_list, file = "../results/var-covar/rem_vcv_list.rda")
save(rem_aic, file = "../results/aic/rem_aic.rda")

####### Distance Model Selection ##################

species <- substr(list.files(path = "data/distance"), 
                  start = 1, 
                  stop = 4)

dist_coef <- data.frame(Species = rep(species, each = n_dis_models),
                        n = NA, 
                        model = rep(seq(1,n_dis_models), times = length(species)), 
                        aic = NA,
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

dis_aic <- vector(mode = "list", length = length(species))
names(dis_aic) <- species

for (s in species)
{
  load(file = paste0("data/distance/", s, ".rda"))
  dist_coef[which(dist_coef$Species == s), "n"] <- nrow(distance_list[[1]]$Y)
  
  aic <- aic_multi(x = distance_list)
  dist_coef[which(dist_coef$Species == s), "aic"] <- aic
  
  aic_df <- data.frame(Model = seq(1, n_dis_models),
                       AIC = aic)
  aic_df <- aic_df[order(aic_df$AIC), ]
  aic_df$Delta_AIC <- aic_df$AIC - aic_df$AIC[1]
  dis_aic[[s]] <- aic_df
  
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
save(dis_aic, file = "../results/aic/dis_aic.rda")