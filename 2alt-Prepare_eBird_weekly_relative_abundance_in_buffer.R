library(sf)
library(tidyverse)
library(ebirdst)
library(terra)
library(bbsBayes2)
library(exactextractr)

#source("eBirdst.R") #script to set access key for local user (not included in archive)
# download and extract eBird seasonal relative abundance within
# BBS route path buffers
yr_ebird <- 2023 # prediction year for eBird relative abundance

if(Sys.info()[["nodename"]] == "WNCRLABN72960"){
  setwd("c:/Users/SmithAC/Documents/GitHub/PIF_population_estimates")
}

routes_buf <- readRDS("data/all_routes_buffered_400m.rds")

# run once to generate list of species with simplified english names
# and then saved as rds
# sps_list <- readRDS("data/all_counts_by_route.rds") %>%
#   filter(count > 0) %>%
#   select(aou,english,french,genus,species,route_name) %>%
#   distinct() %>%
#   group_by(aou,english,french,genus,species) %>%
#   summarise(n_routes_w_obs = n()) %>%
#   filter(n_routes_w_obs > 19,
#          !is.na(english),
#          !grepl(pattern = "^unid.",english),
#          !grepl(pattern = "^\\(",english),
#          english != c("Western Grebe (Clark's/Western)")) %>%
#   mutate(english = ifelse(grepl(pattern = " \\(",
#                                 english),
#                           str_extract(string = english,pattern = ".*(?= \\()"),
#                           english))
#
# saveRDS(sps_list,"data/all_bbs_data_species_list.rds")

sps_list <- readRDS("data/all_bbs_data_species_list.rds")
# eBird relative abundance ------------------------------------------------

#check ebird diagnostics

qual_ebird <- ebirdst_runs

re_calc <- TRUE # TRUE if want to re-calcluate the eBird-BBS overlap
      # e.g., when new route spatial data

use_weekly <- TRUE # only use TRUE if goal is to model all species
      # as weekly abundance, including breeding distributions


use_uncertainty <- FALSE # TRUE if teh uncertainty values of the relative abundance are desired

bbs_start <- as_date("2023-06-01") # Latest date for start of BBS-relevant breeding season
bbs_end <- as_date("2023-07-07") # Latest data for end of BBS-relevant breeding season



metric_used <- "mean" #options are "max", mean" or "median" (although "median" won't work for seasonal)





adjustments <- read_csv("Population_estimates_workflow_support.csv")



raw_counts_all <- readRDS("data/all_counts_by_route.rds")


ii <- which(sps_list$english %in% adjustments$common_name)

for(i in ii[c(154:404)]){

  sp_sel <- unname(unlist(sps_list[i,"english"]))
species_ebird <- ebirdst::get_species(sp_sel)

raw_counts <- raw_counts_all %>%
  filter(grepl(paste0("^",sp_sel),english))

if((file.exists(paste0("data/species_relative_abundance_2023/",
                      species_ebird,"_weekly_relative_abundance.rds")) &
   file.exists(paste0("data/species_relative_abundance_2023/",
                      species_ebird,"_route_weekly_relative_abundance.rds"))) &
   !re_calc){
  next
}

if(is.na(species_ebird)){
  sps_list[i,"available_ebird_data"] <- "Not in ebird species list"


  next

}

qual_sel <- qual_ebird[which(qual_ebird$species_code == species_ebird),]
breed_qual <- unname(unlist(qual_sel[,"breeding_quality"]))
resident_qual <- unname(unlist(qual_sel[,"resident_quality"]))
resident <- unname(unlist(qual_sel[,"is_resident"]))

breeding_start <- unname(unlist(qual_sel[,"breeding_start"]))
breeding_end <- unname(unlist(qual_sel[,"breeding_end"]))


sps_list[i,"is_resident"] <- resident
sps_list[i,"breeding_quality"] <- breed_qual
sps_list[i,"resident_quality"] <- resident_qual
sps_list[i,"species_ebird"] <- species_ebird
sps_list[i,"breeding_start"] <- breeding_start
sps_list[i,"breeding_end"] <- breeding_end

if(!resident & breed_qual == 0){
  sps_list[i,"available_ebird_data"] <- "Only nonbreeding abundance"


  next
}


down <- try(ebirdst::ebirdst_download_status(species_ebird,
                                             download_ranges = FALSE,
                                             download_abundance = TRUE,
                                             download_occurrence = FALSE,
                                             force = FALSE,
                                             pattern = "abundance_median_3km_"),
            silent = TRUE)


if(use_uncertainty){
down_lower <- try(ebirdst::ebirdst_download_status(species_ebird,
                                             download_ranges = FALSE,
                                             download_abundance = TRUE,
                                             download_occurrence = FALSE,
                                             force = FALSE,
                                             pattern = "abundance_lower_3km_"),
            silent = TRUE)


down_upper <- try(ebirdst::ebirdst_download_status(species_ebird,
                                                   download_ranges = FALSE,
                                                   download_abundance = TRUE,
                                                   download_occurrence = FALSE,
                                                   force = FALSE,
                                                   pattern = "abundance_upper_3km_"),
                  silent = TRUE)


}
if(class(down) == "try-error"){

  sps_list[i,"available_ebird_data"] <- "Download fail"


  next}

sps_list[i,"available_ebird_data"] <- "Yes"





  # Extract the weekly relative abundance for each BBS route
  # so that the model can match the week of the observation with
  # the week of the prediction



  abd_weekly_abundance <- ebirdst::load_raster(species = species_ebird,
                                               resolution = "3km",
                                               period = "weekly",
                                               product = "abundance")

  sel_weeks <- try(as_date(names(abd_weekly_abundance)),silent = TRUE)

  # match survey dates to weeks of the relative abundance surface
  species_doys <- raw_counts |>
    mutate(date = (make_date(yr_ebird,1,1)+days(doy)))

  first_day <- min(species_doys$date)
  last_day <- max(species_doys$date)

  sel_weeks <- sel_weeks[which(sel_weeks >= first_day-6 &
                                 sel_weeks <= last_day+6)]

  breed_abundance <- abd_weekly_abundance[[as.character(sel_weeks)]]






saveRDS(breed_abundance,paste0("data/species_relative_abundance_2023/",
                               species_ebird,"_weekly_relative_abundance.rds"))




# # clipping ebird data to BBS strata region --------------------------------
# # this clipping is done only to speed up the overlay between BBS and ebird
#  the rest of the ebird distribution map is unaffected for the eventual
#  predictions of density across the full species' range (i.e., see
#  saveRDS() call above)
bbs_strata <- bbsBayes2::load_map("bbs") %>%
  filter(country_code %in% c("CA","US")) |>
  st_buffer(dist = 10000) #add 10km buffer for clipping eBird data

# project boundary to match raster data
region_boundary_proj <- st_transform(bbs_strata, st_crs(breed_abundance))

bbs_strata2 <- bbsBayes2::load_map("bbs")%>%
  filter(country_code %in% c("CA","US"))
strata_proj <- st_transform(bbs_strata2, st_crs(breed_abundance))

routes_nm <- routes_buf %>%
  select(route_name,route) %>%
  st_transform(crs = st_crs(bbs_strata2)) %>%
  sf::st_join(bbs_strata2, largest = TRUE) %>%
  group_by(strata_name) %>%
  summarise() %>%
  sf::st_transform(crs = st_crs(breed_abundance))
# crop and mask to boundary of BBS
breed_abundance <- crop(breed_abundance, region_boundary_proj) |>
  mask(region_boundary_proj)


# project buffers to match raster data
routes_buf_proj <- st_transform(routes_buf, st_crs(breed_abundance))



abundance_df <- exact_extract(breed_abundance,
                                       routes_buf_proj,
                                       fun = "mean",
                                       #force_df = TRUE,
                                       append_cols = "route_name",
                              default_value = 0)


# abundance_df <- abundance_df |>
#   filter(
#   ebird_abund > 0,
#   !is.na(ebird_abund)
# ) |>
#   distinct()

saveRDS(abundance_df,
        paste0("data/species_relative_abundance_2023/",species_ebird,"_route_weekly_relative_abundance.rds"))







abundance_strat <- exact_extract(breed_abundance,
                                     strata_proj,
                                     fun = "mean",
                                     force_df = TRUE,
                                     append_cols = "strata_name")|>
  rename(mean_ebird_abundance = mean)


total_abundance_in_strata <- exact_extract(breed_abundance,
                                           strata_proj,
                                           fun = "sum",
                                           force_df = TRUE,
                                           append_cols = "strata_name")|>
  rename(total_ebird_abundance = sum)

abundance_strat <- abundance_strat |>
  left_join(total_abundance_in_strata,
            by = "strata_name")

abundance_sampled <- exact_extract(breed_abundance,
                                   routes_nm,
              fun = "mean",
              force_df = TRUE,
              append_cols = "strata_name") |>
  rename(sampled_abundance = mean)

abundance_strat <- abundance_strat |>
  left_join(abundance_sampled,
            by = "strata_name") %>%
  mutate(sampled_avail_ratio = sampled_abundance/mean_ebird_abundance,
         log_ratio = log(sampled_avail_ratio))



#
# abundance_strat[nrow(abundance_strat)+1,"strata_name"] <- "USA_CAN"
# abundance_strat[nrow(abundance_strat),"mean_ebird_abundance"] <- terra::global(breed_abundance, "mean", na.rm = TRUE)
# abundance_strat[nrow(abundance_strat),"total_ebird_abundance"] <- terra::global(breed_abundance, "sum", na.rm = TRUE)



saveRDS(abundance_strat,
        paste0("data/species_relative_abundance_2023/",species_ebird,"_bbs_strata_relative_abundance.rds"))

}
