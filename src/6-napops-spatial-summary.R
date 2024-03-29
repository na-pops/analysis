####### Script Information ########################
# Brandon P.M. Edwards
# NA-POPS: analysis
# 6-napops-spatial-summary.R
# Created February 2021
# Last Updated November 2022

####### Import Libraries and External Files #######

library(sf)
library(magrittr)

####### Set Constants #############################

laea <- 4326# sf::st_crs("+proj=laea +lat_0=40 +lon_0=-95") # Lambert equal area coord reference system
sf::sf_use_s2(FALSE)

####### Read Data #################################

# load(file = "data/combined/samples.rda")
# 
# # Required from 5-napops-quantitative-summary.R
load("../results/quant-summary/dis_species_summary.rda")
load("../results/quant-summary/rem_species_summary.rda")

load(file = "data/combined/counts.rda")
load(file = "data/combined/samples.rda")
load(file = "data/combined/dis_covars_used.rda")
load(file = "data/combined/rem_covars_used.rda")

bcr <- read_sf("../utilities/shp/bcr",
               layer = "BBS_BCR_strata")

####### Analysis by BCR ####################

coords <- project_samples[, c("Latitude", "Longitude")]

coords <- coords[!is.na(coords$Latitude), ]
coords <- coords[!is.na(coords$Longitude), ]

coords <- st_as_sf(coords,coords = c("Longitude","Latitude"), crs = laea)
bcr <- bcr %>% st_transform(st_crs(coords))

counts_bcr <- st_intersects(bcr, coords)
bcr_coverage <- data.frame(BCR = bcr$ST_12,
                           ncounts = lengths(counts_bcr),
                           sqrt_ncounts = sqrt(lengths(counts_bcr)),
                           stringsAsFactors = FALSE)
bcr_coverage$nc_cat <- cut(bcr_coverage$ncounts, breaks = c(-1,0,(c(100,200,500,1000,15000))))

####### Projects by BCR ###########################

coords <- project_samples[, c("Latitude", "Longitude")]

coords <- coords[!is.na(coords$Latitude), ]
coords <- coords[!is.na(coords$Longitude), ]

coords <- st_as_sf(coords,coords = c("Longitude","Latitude"), crs = laea)
bcr <- bcr %>% st_transform(st_crs(coords))

bcr_list <- apply(st_intersects(bcr, coords, sparse = FALSE), 2, 
                     function(col) { 
                       bcr[which(col), ]$ST_12
                     })
bcr_list[which(lengths(bcr_list) == 0)] <- NA

ps_red <- project_samples[!is.na(project_samples$Latitude), ]
ps_red$BCR <- unlist(bcr_list)
ps_red <- ps_red[-which(is.na(ps_red$BCR)), ]

# Create unique indicator for project x BCR, then get rid of duplicates
ps_red$proj_bcr <- paste0(ps_red$Project, "-", ps_red$BCR)
ps_red_unique <- ps_red[!duplicated(ps_red$proj_bcr), ]

coords <- ps_red_unique[, c("Latitude", "Longitude")]

coords <- st_as_sf(coords,coords = c("Longitude","Latitude"), crs = laea)
bcr <- bcr %>% st_transform(st_crs(coords))

counts_bcr_proj <- st_intersects(bcr, coords)
bcr_coverage_project <- data.frame(BCR = bcr$ST_12,
                           ncounts = lengths(counts_bcr_proj),
                           sqrt_ncounts = sqrt(lengths(counts_bcr_proj)),
                           stringsAsFactors = FALSE)
bcr_coverage_project$nc_cat <- cut(bcr_coverage_project$ncounts, breaks = c(-1,0,(c(10,20,30, 40,50,60,70,80,90,100))))

####### Analysis by Species (BCR, Distance) #######

bcr_dis_coverage <- vector(mode = "list", length = length(names(dis_species_summary)))
names(bcr_dis_coverage) <- names(dis_species_summary)

for (sp in names(dis_species_summary))
{
  temp <- project_counts[which(project_counts$Species == sp), ]
  temp <- temp[!duplicated(temp$Sample_ID), ]
  temp <- merge(x = temp, y = dis_covars_used, by = "Sample_ID")
  temp <- merge(x = temp, y = project_samples[, c("Sample_ID", "Latitude", "Longitude")],
                by = "Sample_ID")
  coords <- temp[, c("Latitude", "Longitude")]
  
  coords <- coords[!is.na(coords$Latitude), ]
  coords <- coords[!is.na(coords$Longitude), ]
  
  coords <- st_as_sf(coords,coords = c("Longitude","Latitude"), crs = laea)
  bcr <- bcr %>% st_transform(st_crs(coords))
  
  counts_bcr_sp <- st_intersects(bcr, coords)
  bcr_coverage_sp <- data.frame(BCR = bcr$ST_12,
                                ncounts = lengths(counts_bcr_sp),
                                sqrt_ncounts = sqrt(lengths(counts_bcr_sp)),
                                stringsAsFactors = FALSE)
  bcr_coverage_sp$nc_cat <- cut(bcr_coverage_sp$ncounts, breaks = c(-1,0,(c(100,200,500,1000,15000))))
  bcr_dis_coverage[[sp]] <- bcr_coverage_sp
}

####### Analysis by Species (BCR, Removal) ########

bcr_rem_coverage <- vector(mode = "list", length = length(names(rem_species_summary)))
names(bcr_rem_coverage) <- names(rem_species_summary)

for (sp in names(rem_species_summary))
{
  temp <- project_counts[which(project_counts$Species == sp), ]
  temp <- temp[!duplicated(temp$Sample_ID), ]
  temp <- merge(x = temp, y = rem_covars_used, by = "Sample_ID")
  temp <- merge(x = temp, y = project_samples[, c("Sample_ID", "Latitude", "Longitude")],
                by = "Sample_ID")
  coords <- temp[, c("Latitude", "Longitude")]
  
  coords <- coords[!is.na(coords$Latitude), ]
  coords <- coords[!is.na(coords$Longitude), ]
  
  coords <- st_as_sf(coords,coords = c("Longitude","Latitude"), crs = laea)
  bcr <- bcr %>% st_transform(st_crs(coords))
  
  counts_bcr_sp <- st_intersects(bcr, coords)
  bcr_coverage_sp <- data.frame(BCR = bcr$ST_12,
                                ncounts = lengths(counts_bcr_sp),
                                sqrt_ncounts = sqrt(lengths(counts_bcr_sp)),
                                stringsAsFactors = FALSE)
  bcr_coverage_sp$nc_cat <- cut(bcr_coverage_sp$ncounts, breaks = c(-1,0,(c(100,200,500,1000,15000))))
  bcr_rem_coverage[[sp]] <- bcr_coverage_sp
}

####### Output Summary Statistics and Tables ######

save(bcr_coverage, file = "../results/spatial-summary/project_coverage_bcr.rda")
save(bcr_coverage_project, file = "../results/spatial-summary/nproj_coverage_bcr.rda")
save(bcr_dis_coverage, file = "../results/spatial-summary/dis_coverage_bcr.rda")
save(bcr_rem_coverage, file = "../results/spatial-summary/rem_coverage_bcr.rda")
