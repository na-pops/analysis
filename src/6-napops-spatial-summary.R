####### Script Information ########################
# Brandon P.M. Edwards
# NA-POPS: analysis
# 6-napops-spatial-summary.R
# Created February 2021
# Last Updated March 2021

####### Import Libraries and External Files #######

library(bbsBayes)
library(sf)
library(tidyverse)

####### Read Data #################################

load(file = here::here("data/samples.rda"))

# Required from 5-napops-quantitative-summary.R
load("../results/quant-summary/dis_species_summary.rda")
load("../results/quant-summary/rem_species_summary.rda")

bcr <- read_sf(dsn = system.file("maps",
                                   package="bbsBayes"),
                 layer = "BBS_BCR_strata")
strat <- read_sf(dsn = system.file("maps",
                                      package="bbsBayes"),
                    layer = "BBS_usgs_strata")
state <- read_sf(dsn = system.file("maps",
                                      package="bbsBayes"),
                    layer = "BBS_ProvState_strata")

####### Analysis by BCR x State/Prov ##############

coords <- project_samples[, c("Latitude", "Longitude")]

coords <- coords[!is.na(coords$Latitude), ]
coords <- coords[!is.na(coords$Longitude), ]

coords <- st_as_sf(coords,coords = c("Longitude","Latitude"), crs = 4326)
strat <- strat %>% st_transform(st_crs(coords))

counts_strat <- st_intersects(strat, coords)
bcr_state_coverage <- data.frame(ST_12 = strat$ST_12,
                             ncounts = lengths(counts_strat),
                             sqrt_ncounts = sqrt(lengths(counts_strat)),
                             stringsAsFactors = FALSE)
bcr_state_coverage$nc_cat <- cut(bcr_state_coverage$ncounts, breaks = c(-1,0,(c(100,200,500,1000,15000))))
bcr_state_coverage <- left_join(strat, bcr_state_coverage)

####### Analysis by State Only ####################

coords <- project_samples[, c("Latitude", "Longitude")]

coords <- coords[!is.na(coords$Latitude), ]
coords <- coords[!is.na(coords$Longitude), ]

coords <- st_as_sf(coords,coords = c("Longitude","Latitude"), crs = 4326)
state <- state %>% st_transform(st_crs(coords))

counts_state <- st_intersects(state, coords)
state_coverage <- data.frame(ST_12 = state$ST_12,
                             ncounts = lengths(counts_state),
                             sqrt_ncounts = sqrt(lengths(counts_state)),
                             stringsAsFactors = FALSE)
state_coverage$nc_cat <- cut(state_coverage$ncounts, breaks = c(-1,0,(c(100,200,500,1000,15000))))
state_coverage <- left_join(state, state_coverage)

####### Analysis by BCR ####################

coords <- project_samples[, c("Latitude", "Longitude")]

coords <- coords[!is.na(coords$Latitude), ]
coords <- coords[!is.na(coords$Longitude), ]

coords <- st_as_sf(coords,coords = c("Longitude","Latitude"), crs = 4326)
bcr <- bcr %>% st_transform(st_crs(coords))

counts_bcr <- st_intersects(bcr, coords)
bcr_coverage <- data.frame(ST_12 = bcr$ST_12,
                           ncounts = lengths(counts_bcr),
                           sqrt_ncounts = sqrt(lengths(counts_bcr)),
                           stringsAsFactors = FALSE)
bcr_coverage$nc_cat <- cut(bcr_coverage$ncounts, breaks = c(-1,0,(c(100,200,500,1000,15000))))
bcr_coverage <- left_join(bcr, bcr_coverage)

####### Analysis by Species (BCR, Distance) #######

bcr_dis_coverage <- vector(mode = "list", length = length(names(dis_species_summary)))
names(bcr_dis_coverage) <- names(dis_species_summary)

for (sp in names(dis_species_summary))
{
  df <- dis_species_summary[[sp]]
  coords <- df[, c("Latitude", "Longitude")]
  
  coords <- coords[!is.na(coords$Latitude), ]
  coords <- coords[!is.na(coords$Longitude), ]
  
  coords <- st_as_sf(coords,coords = c("Longitude","Latitude"), crs = 4326)
  bcr <- bcr %>% st_transform(st_crs(coords))
  
  counts_bcr_sp <- st_intersects(bcr, coords)
  bcr_coverage_sp <- data.frame(ST_12 = bcr$ST_12,
                             ncounts = lengths(counts_bcr_sp),
                             sqrt_ncounts = sqrt(lengths(counts_bcr_sp)),
                             stringsAsFactors = FALSE)
  bcr_coverage_sp$nc_cat <- cut(bcr_coverage_sp$ncounts, breaks = c(-1,0,(c(100,200,500,1000,15000))))
  bcr_coverage_sp <- left_join(bcr, bcr_coverage_sp)
  bcr_dis_coverage[[sp]] <- bcr_coverage_sp
}

####### Analysis by Species (BCR, Removal) ########

bcr_rem_coverage <- vector(mode = "list", length = length(names(rem_species_summary)))
names(bcr_rem_coverage) <- names(rem_species_summary)

for (sp in names(rem_species_summary))
{
  df <- rem_species_summary[[sp]]
  coords <- df[, c("Latitude", "Longitude")]
  
  coords <- coords[!is.na(coords$Latitude), ]
  coords <- coords[!is.na(coords$Longitude), ]
  
  coords <- st_as_sf(coords,coords = c("Longitude","Latitude"), crs = 4326)
  bcr <- bcr %>% st_transform(st_crs(coords))
  
  counts_bcr_sp <- st_intersects(bcr, coords)
  bcr_coverage_sp <- data.frame(ST_12 = bcr$ST_12,
                                ncounts = lengths(counts_bcr_sp),
                                sqrt_ncounts = sqrt(lengths(counts_bcr_sp)),
                                stringsAsFactors = FALSE)
  bcr_coverage_sp$nc_cat <- cut(bcr_coverage_sp$ncounts, breaks = c(-1,0,(c(100,200,500,1000,15000))))
  bcr_coverage_sp <- left_join(bcr, bcr_coverage_sp)
  bcr_rem_coverage[[sp]] <- bcr_coverage_sp
}

####### Output Summary Statistics and Tables ######

save(bcr_state_coverage, file = "../results/spatial-summary/project_coverage_bcr_state.rda")
save(state_coverage, file = "../results/spatial-summary/project_coverage_state.rda")
save(bcr_coverage, file = "../results/spatial-summary/project_coverage_bcr.rda")
save(bcr_dis_coverage, file = "../results/spatial-summary/dis_coverage_bcr.rda")
save(bcr_rem_coverage, file = "../results/spatial-summary/rem_coverage_bcr.rda")
