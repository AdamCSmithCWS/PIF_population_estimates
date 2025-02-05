#### Download and prep NA BBS data ####
# Required libraries
library(bbsBayes2) # devtools::install_github("bbsBayes/bbsBayes2")

library(tidyverse)
# Download NA BBS data up to 2023 using
# bbsBayes2 package; this only needs to be done once as the data will be stored
# and used for any future calls from bbsBayes functions
if(!bbsBayes2::have_bbs_data(quiet = TRUE)){
  bbsBayes2::fetch_bbs_data() # only downloads if doesn't already exist locally
}

# Load the full bbs dataset
all_dat <- bbsBayes2::load_bbs_data()

counts <- all_dat$birds %>%  # positive counts of each species during every BBS survey
  select(route_data_id,aou,species_total) #dropping all but the critical columns

sampling_events <- all_dat$routes %>%  # date, time, starting location, for every BBS survey since 1966
  filter(year > 2012) %>% # 10 years from 2013 - 2023, missing 2020.
  select(country,state,st_abrev,route_name,bcr,
         country_num,state_num,route,
         latitude,longitude,
         route_data_id, #this is the critical unique id for a sampling event
         year,month,day,obs_n) %>%
  mutate(strata_name = paste(st_abrev,bcr,sep = "-"),
         strata_name_numeric = paste(state_num,bcr,sep = "-"))

species_list <- all_dat$species %>% # species list
  dplyr::filter(unid_combined == TRUE,
                !(grepl("^(unid.)",english) | grepl("^(hybrid)",english) )) %>%  #this is specific to the BBS,
  select(aou,english,french,genus,species) %>%  # the TRUE option combines 13 taxonomic units that have been split or
  mutate(latin = paste(genus,species)) # lumped over the history of the BBS, this approach (== FALSE) retains the taxonomic
# units, exactly as they are stored in the BBS database


luni <- function(x){
  y <- length(unique(x))
  return(y)
}


# strata loop to limit species lists by strata ----------------------------
#setting limits on the strata that are included
more_than_n_routes_stratum <- 1 # More than this number of BBS routes in the stratum
more_than_n_surveys_stratum <- 1 # More than this number of surveys (routes * years) in the stratum

strata_list <- sampling_events %>%
  group_by(strata_name, strata_name_numeric) %>%
  summarise(.,
            n_surveys = n(),
            n_routes = luni(route_name)) %>%
  filter(n_routes > more_than_n_routes_stratum, #filter on the total number of routes (sampling locations)
         n_surveys > more_than_n_surveys_stratum) #filter on the total number of sampling events




data_all <- NULL # empty object to facilitate bind_rows() at the end of the next loop

# this strata-specific approach removes the species by strata
# combinations in the data that are exclusively 0-values.
# Species are only included in a given strata, if they meet the
# following inclusion criteria:
more_than_n_counts_species <- 2
# More than this number of surveys (routes * years) on which the species has been observed
# this is probably too low... only ~2 counts/year

for(reg in strata_list$strata_name){
  samp_ev_select <- sampling_events %>%
    filter(strata_name == reg)

  counts_select <- counts %>%
    filter(route_data_id %in% samp_ev_select$route_data_id)

  # counting number of observations of each species
  species_inc <- counts_select %>%
    group_by(aou) %>%
    summarise(.,
              n_obs = n()) #n_obs represents the number of non-zero counts across routes and year

  #identifying which species to keep (n_obs > more_than_n_counts_species)
  species_keep <- species_inc %>%
    filter(n_obs > more_than_n_counts_species) %>%
    select(aou)

  # complete set of species by surveys
  samp_ev_select <- samp_ev_select %>%
    expand_grid(.,species_keep)
  # above creates a complete set of species by survey combinations

  #zero-filling
  data_sel <- samp_ev_select %>%
    full_join(.,counts_select,
              by = c("route_data_id",
                     "aou")) %>%
    mutate(species_total = ifelse(is.na(species_total),
                                  0,
                                  species_total))
  # above creates the full set of observations with appropriate zero values for
  # the species in species_keep at all surveys conducted in the region


  data_all <- bind_rows(data_all,
                        data_sel)
  rm(data_sel)


}

data_all <- data_all %>%
  left_join(.,species_list,
            by = "aou")
n_species = length(unique(data_all$english))


saveRDS(data_all, "Stanton_2019_code/input_files/All_BBS_data_by_strata_name.rds")


#
# df <- readRDS("Stanton_2019_code/input_files/All_BBS_data_by_strata_name.rds")
# all.sp <- unique(df$aou)
# st.spp <- c()
# for(i in 1:length(all.sp)){
#   if(all.sp[i] %in% st.spp) next # Check if spp group already processed
#   # if(all.sp[i] %in% sp.group.aou){
#   #   group.i <- sp.group[[grep(all.sp[i], sp.group)]]$aou
#   #   df.i <- df[df$aou %in% group.i, ]
#   #   df.i$aou <- all.sp[i]
#   #   st.spp <- c(st.spp, group.i) # Add group to list of spp already processed
#   # } else {
#     df.i <- df[df$aou == all.sp[i], ]
#     st.spp <- c(st.spp, all.sp[i])
#   #}
#   write.csv(df.i, paste("Stanton_2019_code/input_files/species_data/", all.sp[i], '_BBS_10yr', Sys.Date(), '.csv', sep=''), row.names = FALSE)
# }



