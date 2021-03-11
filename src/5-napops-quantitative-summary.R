####### Script Information ########################
# Brandon P.M. Edwards
# NA-POPS: analysis
# 5-napops-quantitative-summary.R
# Created January 2021
# Last Updated February 2021

####### Import Libraries and External Files #######

library(sf)
library(bbsBayes)

source("../utilities/get-data.R")
source("../utilities/order-taxo.R")
source("../utilities/rm-non-sp.R")

####### Create Empty Summary Data List ############

summary_stats <- list()

####### Read Data #################################

load(file = here::here("data/counts.rda"))
load(file = here::here("data/samples.rda"))

project_list <- read.table(here::here("../utilities/proj-list"))[,1]
ibp <- read.csv("../utilities/IBP-Alpha-Codes20.csv")
dis <- rm_non_sp(order_taxo(read.csv("../results/coefficients/distance.csv")))
rem <- rm_non_sp(order_taxo(read.csv("../results/coefficients/removal.csv")))

####### Summary of Samples and Counts DFs #########

summary_stats[["n_observations"]] <- nrow(project_counts)
summary_stats[["n_samples"]] <- nrow(project_samples)

utc <- project_samples$UTC
year <- unique(strptime(utc, 
                        format = "%Y-%m-%dT%H:%M:%SZ", 
                        tz = "UTC")$year + 1900)
summary_stats[["years"]] <- sort(year[!is.na(year)])

####### Number of Projects ########################

all_projects <- unique(project_samples$Project)
p_split <- strsplit(all_projects, split = ":")

# First deal with the metaprojects
p_split_meta <- p_split[which(lengths(p_split) == 2)]
project_df <- data.frame(Metaproject = sapply(p_split_meta, "[[", 1),
                      Project = sapply(p_split_meta, "[[", 2))
project_df <- project_df[order(project_df$Metaproject, project_df$Project), ]


# Now the standalone projects
p_split_alone <- sort(unlist(p_split[which(lengths(p_split) == 1)]))
project_df <- rbind(project_df,
                    data.frame(Metaproject = rep("Standalone projects not part of a larger metaproject",
                                                 length(p_split_alone)),
                               Project = p_split_alone))

summary_stats[["total_projects"]] <- nrow(project_df)

####### Number of Species & Species Table #########

codes <- union(rem$Species, dis$Species)
codes <- codes[match(ibp$SPEC, codes)]
codes <- codes[!is.na(codes)]

common_name <- ibp[which(ibp$SPEC %in% codes), "COMMONNAME"]

sci_name <- ibp[which(ibp$SPEC %in% codes), "SCINAME"]

in_removal <- ifelse(codes %in% rem$Species, 1, 0)

in_distance <- ifelse(codes %in% dis$Species, 1, 0)

species_table <- data.frame(Code = codes,
                            Common_Name = common_name,
                            Scientific_Name = sci_name,
                            removal = in_removal,
                            distance = in_distance)

summary_stats[["n_species"]] <- nrow(species_table)

####### Number of States/Provinces ################

# coords <- do.call(rbind, project_samples)[, c("Latitude", "Longitude")]
# coords <- coords[!is.na(coords$Latitude), ]
# coords <- coords[!is.na(coords$Longitude), ]
# 
# coords <- st_as_sf(coords,coords = c("Longitude","Latitude"), crs = 4326)
# state = sf::read_sf(dsn = system.file("maps",
#                                       package="bbsBayes"),
#                     layer = "BBS_ProvState_strata")
# state <- state %>% st_transform(st_crs(coords))
# 
# counts_state <- st_intersects(state, coords)
# cs_df <- data.frame(State = state$ST_12,
#                     Counts = lengths(counts_state))
# 
# summary_stats[["n_states"]] <- nrow(cs_df[which(cs_df$Counts > 0), ])

####### Output Summary Statistics and Tables ######

write.table(x = project_df, file = "../results/quant-summary/project_list.csv",
            sep = ",", row.names = FALSE)

write.table(x = species_table, file = "../results/quant-summary/species_table.csv",
            sep = ",", row.names = FALSE)

save(summary_stats, file = "../results/quant-summary/summary_statistics.rda")