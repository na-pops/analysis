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

load(file = here::here("data/counts.rda"))
load(file = here::here("data/samples.rda"))
# bcrf <- read_sf(dsn = "../utilities/shp/bcr",
#                 layer = 'bcrfinalg0_Project2_Dissolve')
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

# Create df with no BAM data, we have to deal with that separately
samples_no_bam <- project_samples[-which(grepl("bam-can_m",
                                               project_samples$Project) |
                                           grepl("bam-mn_m",
                                                 project_samples$Project)), ]

coords <- samples_no_bam[, c("Longitude", "Latitude")]
coords <- coords[which(!is.na(coords$Latitude)), ]
coords <- coords[which(!is.na(coords$Longitude)), ]
pc = st_as_sf(coords,coords = c("Longitude","Latitude"), crs = 4326)

strat1 <- strat %>% st_transform(st_crs(pc))

tt <- st_intersects(strat1,pc)
tt_df <- data.frame(ST_12 = strat1$ST_12,
                    ncounts = lengths(tt),
                    sqrt_ncounts = sqrt(lengths(tt)),
                    stringsAsFactors = FALSE)


# Wrangle Abstracted BAM Canada Data

#' Deal with abstracted BAM counts by using BCR x Prov combination
#' per sampling point

project <- "bam-can_m"
load("../project-bam-can_m/rawdata/bam-can_m_raw.rda")
load("../project-bam-can_m/output/bam-can_m_samples.rda")

# Recreate Sample ID in original data set
# Format time and date
time_formatted <- sprintf("%s:%s:%s", 
                          DT$HOUR, 
                          DT$MIN, 
                          "00")

date_formatted <- sprintf("%s-%s-%s",
                          DT$YEAR,
                          DT$MONTH,
                          DT$DAY)

DT$utc <- paste0(date_formatted,
                 "T",
                 time_formatted,
                 "Z")

DT$Sample_ID <- paste0(project, ":",
                       DT$SS, ":",
                       DT$utc)

#' Match up Sample IDs in data set with sampling DF, and extract
#' the BCR and Jurisdiction fields
bam_sample <- merge(x = sampling_df,
                    y = DT[, c("Sample_ID", "JURS", "BCR")],
                    by = "Sample_ID")
bam_sample <- bam_sample[!duplicated(bam_sample$Sample_ID),]
bam_sample$bcr_prov <- paste0(bam_sample$JURS, "-", bam_sample$BCR)
bam_canada <- bam_sample

# Wrangle Abstracted BAM MN Data

#' Deal with abstracted BAM counts by using BCR x Prov combination
#' per sampling point

project <- "bam-mn_m"
load("../project-bam-mn_m/rawdata/bam-mn_m_raw.rda")
load("../project-bam-mn_m/output/bam-mn_m_samples.rda")

# Recreate Sample ID in original data set
# Format time and date
time_formatted <- sprintf("%s:%s:%s", 
                          MN$HOUR, 
                          MN$MIN, 
                          "00")

date_formatted <- sprintf("%s-%s-%s",
                          MN$YEAR,
                          MN$MONTH,
                          MN$DAY)

MN$utc <- paste0(date_formatted,
                 "T",
                 time_formatted,
                 "Z")

MN$Sample_ID <- paste0(project, ":",
                       MN$SS, ":",
                       MN$utc)

#' Match up Sample IDs in data set with sampling DF, and extract
#' the BCR and Jurisdiction fields
bam_sample <- merge(x = sampling_df,
                    y = MN[, c("Sample_ID", "JURS", "BCR")],
                    by = "Sample_ID")
bam_sample <- bam_sample[!duplicated(bam_sample$Sample_ID),]
bam_sample$bcr_prov <- paste0(bam_sample$JURS, "-", bam_sample$BCR)
bam_mn <- bam_sample

# Combine Wrangled Data

bam_counts <- rbind(bam_canada, bam_mn)
bam <- as.data.frame(table(bam_counts$bcr_prov))
names(bam) <- c("bcr_prov", "count")

bam$bcr_prov <- as.character(bam$bcr_prov)
bam$ST_12 = paste0(c(rep("CA-",13),rep("US-",4),
                     rep("CA-",20),rep("US-",1)),bam$bcr_prov)
bam = bam[-which(bam$bcr_prov %in% c("MN-0","ON-0")),]
for(sst in as.character(bam$ST_12))
{
  ww = which(tt_df$ST_12 == sst)
  wwb = which(bam$ST_12 == sst)
  if(length(ww) > 0)
  {
    tt_df[ww,"ncounts"] <- tt_df[ww,"ncounts"]+bam[wwb,"count"]
  }
}

tt_df$nc_cat <- cut(tt_df$ncounts,breaks = c(-1,0,(c(100,200,500,1000,15000))))
tt_df$sqrt_ncounts <- sqrt(tt_df$ncounts)
bcr_state_coverage <- left_join(strat1,tt_df)

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

####### Output Summary Statistics and Tables ######

save(bcr_state_coverage, file = "../results/spatial-summary/project_coverage_bcr_state.rda")
save(state_coverage, file = "../results/spatial-summary/project_coverage_state.rda")
save(bcr_coverage, file = "../results/spatial-summary/project_coverage_bcr.rda")
