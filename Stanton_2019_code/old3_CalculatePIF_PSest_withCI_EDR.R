##############################
# Calculated PIF Population Size Estimate
# Substituting EDRs estimated through na-pops
# for the MDD values used in previous estimates
#
# Data assemble Data from
# North American Breeding Bird Survey
# (run 1_AssembleBBSdatabySpecies.R first)
#
# Script uses mean count by route for last 10 years of
# data supplied
#
# Currently only outputs North American landbird species
# (specified in Species_PSEst_Adjustments.csv file)
#
# by J.C. Stanton
#
######  USER INPUTS ########
#############################
n.samp <- 100
#library(tidyverse)
# Provide Path to Input file directory - Adjustment parameter files
# Inputfiles.dir <- readline("Provide directory path to location of INPUT files: ")
# if(!file.exists(Inputfiles.dir)) stop("Provide valid directory path to location of input files. ")
Inputfiles.dir <- "Stanton_2019_code/input_files/" #paste(Inputfiles.dir, '/', sep="")

# Spp.dir <- readline("Provide directory path to location of Species files generated in step 1: ")
# if(!file.exists(Spp.dir)) stop("Provide valid directory path to location of species files. ")
#Spp.dir <- "Stanton_2019_code/input_files/species_data/"#paste(Spp.dir, '/', sep="")

# Provide Path to where model output should be written
# Output.dir <- readline("Provide directory path to store output files: ")
# if(!file.exists(Output.dir)) stop("Provide valid directory path to output files. ")
###########################################################

##### dependent datasets  #####
## Parameters used in Population Size estimation
## See associated manuscript for details and source references

# Species specific adjustment factors (Time of Day, Detection Distance, Pair)
ExpAdjs_orig <- read.csv(paste(Inputfiles.dir, 'Species_PSEst_Adjustments.csv', sep=''), stringsAsFactors=FALSE)
#
#
# Add na-pops EDRs --------------------------------------------------------
#
# library(napops)
# napops_species <- list_species() |>
#   dplyr::select(Species,Common_Name)
# ExpAdjs <- ExpAdjs |>
#   dplyr::left_join(napops_species,
#             by = c("cn" = "Common_Name"))
#
# for(i in 1:nrow(ExpAdjs)){
#   if(is.na(ExpAdjs[i,"Species"])){next}
#
#   sp <- ExpAdjs[i,"Species"]
#
#   # choose better supported model between models 1 (intercept) and 2 (roadside)
# cefs <- coef_distance(species = sp) |>
#   dplyr::filter(Model %in% c(1,2)) |>
#   dplyr::arrange(AIC)
#
# if(nrow(cefs) == 0){next} # skip if no napops output
#
# w_mod <- cefs[1,"Model"]
#
# edrt <- try(edr(species = sp,
#               model = w_mod,
#               forest = 0.5,
#               road = TRUE,
#               quantiles = c(0.1625),
#               samples = 1000),
#               silent = TRUE)
#
# if(edrt$EDR_est > 500){ # two species where the model 2 estimates of edr are silly (> 5 km!)
#   edrt <- try(edr(species = sp,
#                   model = 1,
#                   forest = 0.5,
#                   road = TRUE,
#                   quantiles = c(0.1625),
#                   samples = 1000),
#               silent = TRUE)
#   w_mod <- 1
# }
#
# if(class(edrt) == "try-error"){
#     break(paste("edr extraction did not work for",sp))
#
#     }
 %>% %>% %>%
