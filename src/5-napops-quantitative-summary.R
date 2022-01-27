####### Script Information ########################
# Brandon P.M. Edwards
# NA-POPS: analysis
# 5-napops-quantitative-summary.R
# Created January 2021
# Last Updated January 2022

####### Import Libraries and External Files #######

source("../utilities/get-data.R")
source("../utilities/order-taxo.R")
source("../utilities/rm-non-sp.R")

####### Create Empty Summary Data List ############

summary_stats <- list()

####### Read Data #################################

load(file = "data/combined/counts.rda")
load(file = "data/combined/samples.rda")
load(file = "data/combined/dis_covars_used.rda")
load(file = "data/combined/rem_covars_used.rda")

project_list <- read.table("../utilities/proj-list")[,1]
ibp <- read.csv("../utilities/IBP-Alpha-Codes20.csv")
nacc <- read.csv("../utilities/NACC_list_species.csv")
dis <- order_taxo(read.csv("../results/coefficients/distance.csv"))
rem <- order_taxo(read.csv("../results/coefficients/removal.csv"))

####### Summary of Samples and Counts DFs #########

landbirds <- ibp[-which(ibp$LANDBIRD == "FALSE" |
                                 ibp$SP == "+"), "SPEC"]
summary_stats[["n_observations"]] <- nrow(project_counts)
summary_stats[["n_samples"]] <- nrow(project_samples)
summary_stats[["n_rem_samples"]] <- nrow(rem_covars_used)
summary_stats[["n_dis_samples"]] <- nrow(dis_covars_used)

rem_covars <- rem_covars_used[, c("OD", "TSSR", "Time_Method")]
dis_covars <- dis_covars_used[, c("ForestOnly_5x5", "roadside", "Distance_Method")]

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

dis <- dis[!duplicated(dis$Species), ]
rem <- rem[!duplicated(rem$Species), ]

codes <- union(rem$Species, dis$Species)
codes <- codes[match(ibp$SPEC, codes)]
codes <- codes[!is.na(codes)]

common_name <- ibp[which(ibp$SPEC %in% codes), "COMMONNAME"]

sci_name <- ibp[which(ibp$SPEC %in% codes), "SCINAME"]

families <- nacc[which(nacc$common_name %in% common_name), "family"]

in_removal <- ifelse(codes %in% rem$Species, 1, 0)

in_distance <- ifelse(codes %in% dis$Species, 1, 0)

species_table <- data.frame(Code = codes,
                            Common_Name = common_name,
                            Scientific_Name = sci_name,
                            Family = families,
                            removal = in_removal,
                            distance = in_distance)

summary_stats[["n_species"]] <- nrow(species_table)

####### Output Summary Statistics and Tables ######

write.table(x = project_df, file = "../results/quant-summary/project_list.csv",
            sep = ",", row.names = FALSE)

write.table(x = species_table, file = "../results/quant-summary/species_table.csv",
            sep = ",", row.names = FALSE)

save(summary_stats, file = "../results/quant-summary/summary_statistics.rda")

save(dis_covars, file = "../results/quant-summary/dis_covars.rda")
save(rem_covars, file = "../results/quant-summary/rem_covars.rda")
