####### Script Information ########################
# Brandon P.M. Edwards
# NA-POPS: analysis
# 8-simulate-phi.R
# Created January 2021
# Last Updated March 2021

####### Import Libraries and External Files #######

library(MASS)

source("../utilities/order-taxo.R")
source("../utilities/rm-non-sp.R")

####### Read Data #################################

rem <- rm_non_sp(order_taxo(read.csv("../results/coefficients/removal.csv")))
family <- read.csv("../utilities/NACC_list_species.csv")[, c("common_name",
                                                             "family")]
ibp <- read.csv("../utilities/IBP-Alpha-Codes20.csv")[, c("SPEC",
                                                          "COMMONNAME")]
load("../results/var-covar/rem_vcv_best.rda")

####### Simulate phi ##############################

jd <- seq(91, 212)
tssr_values <- seq(-2, 6)
tssr <- rep(tssr_values, each = length(jd))

sim_data <- data.frame(Species = rep(rem$Species,
                                     each = length(tssr)),
                       Model = rep(rem$best_model, each = length(tssr)),
                       Intercept = rep(1, times = length(rem$Species) * 
                                         length(tssr)),
                       TSSR = rep(tssr, length(rem$Species)),
                       JD = rep(jd, length(tssr_values) * length(rem$Species)))

phi <- NULL
phi_low <- NULL
phi_high <- NULL

design <- sim_data[which(sim_data$Species == rem$Species[1]),
                   c("Intercept", "TSSR", "JD")]
design$TSSR <- design$TSSR / 24
design$TSSR2 <- design$TSSR ^ 2
design$JD <- design$JD / 365
design$JD2 <- design$JD ^ 2

design <- design[, c("Intercept", "TSSR", "TSSR2", "JD", "JD2")]

for (sp in rem$Species)
{
  coefficients <- as.numeric(rem[which(rem$Species == sp), 
                                 c("int_best", "tssr_best", "tssr2_best", 
                                   "jd_best", "jd2_best")])
  
  zeros_indices <- which(is.na(coefficients)) - 1
  if (length(zeros_indices) > 0)
  {
    coefficients <- coefficients[-which(is.na(coefficients))]    
  }
  
  vcv <- rem_vcv_best[[sp]]
  
  # Simulate a bunch of possible coefficients
  # Try-catch surrounding for the odd non positive definite var-covar
  sim_coef <-tryCatch(
    {
      rbind(coefficients, MASS::mvrnorm(10^4, coefficients, vcv))
      #sim_coef <- rbind(coefficients, MASS::mvrnorm(10^4, coefficients, vcv))      
    },
    error = function(e)
    {
      return(NA)
    }
  )
  
  if (is.na(sim_coef))
  {
    # In the case of an error, just output the calculated phi and NA
    # for the upper and lower
    if (length(zeros_indices) > 0)
    {
      coef_zeros <- c(coefficients, rep(0, length(zeros_indices)))
      id <- c(seq_along(coefficients), zeros_indices + 0.5)
      coef_zeros <- coef_zeros[order(id)]    
    }else
    {
      coef_zeros <- coefficients
    }
    
    phi_pred <- exp(as.matrix(design) %*% (coef_zeros))
    phi <- c(phi, as.numeric(phi_pred[,1]))
    phi_low <- c(phi_low, rep(NA, nrow(design)))
    phi_high <- c(phi_high, rep(NA, nrow(design)))    
  }else
  {
    # Add columns of zeros back in to where NA coefficients were previously
    # See https://stackoverflow.com/a/1495204/5665609 for explanation
    if (length(zeros_indices) > 0)
    {
      coef_zeros <- cbind(sim_coef, matrix(0,
                                           ncol = length(zeros_indices),
                                           nrow = nrow(sim_coef)))
      id <- c(seq_along(sim_coef[1,]), zeros_indices + 0.5)
      coef_zeros <- coef_zeros[,order(id)]    
    }else
    {
      coef_zeros <- sim_coef
    }
    
    phi_pred <- exp(as.matrix(design) %*% t(coef_zeros))
    phi <- c(phi, as.numeric(phi_pred[,1]))
    
    # Calculate quantiles
    phi_pred <- phi_pred[,-1]
    phi_low <- c(phi_low,
                 as.numeric(apply(phi_pred,
                                  1,
                                  quantile,
                                  probs = c(0.025),
                                  na.rm = TRUE)))
    phi_high <- c(phi_high,
                  as.numeric(apply(phi_pred,
                                   1,
                                   quantile,
                                   probs = c(0.975),
                                   na.rm = TRUE)))       
  }
  
}

sim_data$phi <- phi
sim_data$phi_2.5 <- phi_low
sim_data$phi_97.5 <- phi_high

####### Simulate p ################################

time_values <- c(1, 3, 5, 10)
time <- rep(time_values, times = nrow(sim_data))

sim_data <- sim_data[rep(seq_len(nrow(sim_data)),
                         each = length(time_values)), ]
sim_data$Time <- time
sim_data$p <- 1 - exp(-(sim_data$Time * sim_data$phi))
sim_data$p_2.5 <- 1 - exp(-(sim_data$Time * sim_data$phi_2.5))
sim_data$p_97.5 <- 1 - exp(-(sim_data$Time * sim_data$phi_97.5))

# Save in results, and a local copy in the na-pops dashboard repo
phi_df <- sim_data
save(x = phi_df, file = "../results/simulations/phi.rda")
save(x = phi_df, file = "../napops-dashboard/data/phi.rda")
