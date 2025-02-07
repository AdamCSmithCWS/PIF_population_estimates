library(sf)
library(tidyverse)
library(ebirdst)
#ebirdst::set_ebirdst_access_key("t9el4omae1c3",overwrite = TRUE) #expire April 28 2025
#library(patchwork)
library(terra)
library(bbsBayes2)

# download and extract eBird seasonal relative abundance within
# BBS route path buffers

routes_buf <- readRDS("data/all_routes_buffered.rds")
#routes_buf <- readRDS("all_routes_buffered_1km.rds")

#
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

bbs_start <- as_date("2022-06-01") # Latest date for start of BBS-relevant breeding season
bbs_end <- as_date("2022-07-01") # Latest data for end of BBS-relevant breeding season



metric_used <- "mean" #options are "max", mean" or "median" (although "median" won't work for seasonal)

ii <- which(sps_list$english %in% c("Brown Creeper","Canyon Wren",
                                    "Black-capped Chickadee","Canada Warbler",
                                    "Connecticut Warbler"))

for(i in ii){#rev(1:nrow(sps_list))){

  sp_sel <- unname(unlist(sps_list[i,"english"]))
species_ebird <- ebirdst::get_species(sp_sel)

# paste0("data/species_relative_abundance/",
#        sp_ebird,
#        "_derived_breeding_relative_abundance.rds")


if((file.exists(paste0("data/species_relative_abundance/",
                      species_ebird,"_relative_abundance.rds")) &
   file.exists(paste0("data/species_relative_abundance/",
                      species_ebird,"_derived_breeding_relative_abundance.rds"))) &
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

if(!resident){
sel_start <- ifelse(as_date(breeding_start)>bbs_start, # if stable breeding season is later than the 1st of June (mostly only northern species)
                    bbs_start,# then use the first of June
                    as_date("2022-05-20")) # otherwise use May 20th as a compromise between stable breeding season and BBS observation dates


sel_end <- ifelse(as_date(breeding_end)>bbs_end, # if stable breeding season ending is later than the 1st week of July
                    bbs_end,# then use the last week of bbs season
                    as_date("2022-07-07")) # otherwise use July 7 as a compromise between stable breeding season end and BBS observation dates


}else{
  sel_start <- as_date("2022-05-20") # use May 20, as start of BBS observations season for resident species
  sel_end <- as_date("2022-07-07")# use July 7, as end of BBS observations season for resident species
  }

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


if(!resident & breed_qual > 0 & !use_weekly) {
  down <- try(ebirdst::ebirdst_download_status(species_ebird,
                                             download_ranges = FALSE,
                                             download_abundance = TRUE,
                                             download_occurrence = FALSE,
                                             force = FALSE,
                                             pattern = paste0("abundance_seasonal_",metric_used,"_3km_")),
            silent = TRUE)
}

if(resident | use_weekly) {

down <- try(ebirdst::ebirdst_download_status(species_ebird,
                                             download_ranges = FALSE,
                                             download_abundance = TRUE,
                                             download_occurrence = FALSE,
                                             force = FALSE,
                                             pattern = "abundance_median_3km_"),
            silent = TRUE)
}

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



if(class(down) == "try-error"){

  sps_list[i,"available_ebird_data"] <- "Download fail"


  next}

sps_list[i,"available_ebird_data"] <- "Yes"





if(resident | use_weekly){

  # The resident relative abundance is not linearly related
  # to the BBS mean observations (log-linear) (e.g., BCCH)
  # Therefore: if resident, then extract the weekly abundance data for the
  # BBS survey season and generate a breeding season abundance



  abd_weekly_abundance <- ebirdst::load_raster(species = species_ebird,
                                               resolution = "3km",
                                               period = "weekly",
                                               product = "abundance")

  sel_weeks <- try(as_date(names(abd_weekly_abundance)),silent = TRUE)
  sel_weeks <- sel_weeks[which(sel_weeks > sel_start &
                                 sel_weeks < sel_end)]

  abd_weekly_abundance_bbs_season <- abd_weekly_abundance[[as.character(sel_weeks)]]

  if(metric_used == "max"){
  breed_abundance <- terra::quantile(abd_weekly_abundance_bbs_season,
                                     probs = 1, na.rm = TRUE)
  }
  if(metric_used == "median"){
    breed_abundance <- terra::median(abd_weekly_abundance_bbs_season, na.rm = TRUE)
  }
  if(metric_used == "mean"){
    breed_abundance <- terra::mean(abd_weekly_abundance_bbs_season, na.rm = TRUE)
  }




# SD of abundance ---------------------------------------------------------

  abd_weekly_abundance_lower <- ebirdst::load_raster(species = species_ebird,
                                               resolution = "3km",
                                               period = "weekly",
                                               product = "abundance",
                                               metric = "lower")

  abd_weekly_abundance_upper <- ebirdst::load_raster(species = species_ebird,
                                                     resolution = "3km",
                                                     period = "weekly",
                                                     product = "abundance",
                                                     metric = "upper")


  abd_weekly_abundance_lower_bbs_season <- abd_weekly_abundance_lower[[as.character(sel_weeks)]]
  abd_weekly_abundance_upper_bbs_season <- abd_weekly_abundance_upper[[as.character(sel_weeks)]]

  abd_weekly_abundance_ci <- abd_weekly_abundance_upper_bbs_season - abd_weekly_abundance_lower_bbs_season

  abd_weekly_abundance_lower_bbs_season <- log(abd_weekly_abundance_lower[[as.character(sel_weeks)]]+0.0001)
  abd_weekly_abundance_upper_bbs_season <- log(abd_weekly_abundance_upper[[as.character(sel_weeks)]]+0.0001)

  abd_weekly_abundance_logci <- abd_weekly_abundance_upper_bbs_season - abd_weekly_abundance_lower_bbs_season



  ci_to_sd <- (2 * qnorm(p = 0.9))
  abd_weekly_abundance_logsd <- abd_weekly_abundance_logci/ci_to_sd
  # if(metric_used == "max"){
  #   breed_abundance <- terra::max(abd_weekly_abundance_bbs_season)
  # }
  # if(metric_used == "median"){
  #   breed_abundance <- terra::median(abd_weekly_abundance_bbs_season)
  # }
  # if(metric_used == "mean"){
  breed_abundance_sd <- terra::mean(abd_weekly_abundance_logsd, na.rm = TRUE)

  breed_abundance_ci <- terra::mean(abd_weekly_abundance_ci, na.rm = TRUE)


}else{ # if breeding not resident


  abd_seasonal_abundance <- ebirdst::load_raster(species = species_ebird,
                                                 resolution = "3km",
                                                 period = "seasonal",
                                                 product = "abundance",
                                                 metric = metric_used)  #3km high resolution

  breed_abundance <- abd_seasonal_abundance[["breeding"]]





  # SD of abundance ---------------------------------------------------------

  abd_weekly_abundance_lower <- ebirdst::load_raster(species = species_ebird,
                                                     resolution = "3km",
                                                     period = "weekly",
                                                     product = "abundance",
                                                     metric = "lower")

  abd_weekly_abundance_upper <- ebirdst::load_raster(species = species_ebird,
                                                     resolution = "3km",
                                                     period = "weekly",
                                                     product = "abundance",
                                                     metric = "upper")

  sel_weeks <- try(as_date(names(abd_weekly_abundance_upper)),silent = TRUE)
  sel_weeks <- sel_weeks[which(sel_weeks > sel_start &
                                 sel_weeks < sel_end)]

  abd_weekly_abundance_lower_bbs_season <- abd_weekly_abundance_lower[[as.character(sel_weeks)]]
  abd_weekly_abundance_upper_bbs_season <- abd_weekly_abundance_upper[[as.character(sel_weeks)]]

    abd_weekly_abundance_ci <- abd_weekly_abundance_upper_bbs_season - abd_weekly_abundance_lower_bbs_season

    abd_weekly_abundance_lower_bbs_season <- log(abd_weekly_abundance_lower[[as.character(sel_weeks)]]+0.0001)
    abd_weekly_abundance_upper_bbs_season <- log(abd_weekly_abundance_upper[[as.character(sel_weeks)]]+0.0001)

    abd_weekly_abundance_logci <- abd_weekly_abundance_upper_bbs_season - abd_weekly_abundance_lower_bbs_season



  ci_to_sd <- (2 * qnorm(p = 0.9))
  abd_weekly_abundance_logsd <- abd_weekly_abundance_logci/ci_to_sd
  # if(metric_used == "max"){
  #   breed_abundance <- terra::max(abd_weekly_abundance_bbs_season)
  # }
  # if(metric_used == "median"){
  #   breed_abundance <- terra::median(abd_weekly_abundance_bbs_season)
  # }
  # if(metric_used == "mean"){
  breed_abundance_sd <- terra::mean(abd_weekly_abundance_logsd, na.rm = TRUE)

  breed_abundance_ci <- terra::mean(abd_weekly_abundance_ci, na.rm = TRUE)

 #}


}

saveRDS(breed_abundance,paste0("data/species_relative_abundance/",species_ebird,"_derived_breeding_relative_abundance.rds"))
saveRDS(breed_abundance_sd,paste0("data/species_relative_abundance/",species_ebird,"_derived_breeding_relative_abundance_logsd.rds"))






# # clipping ebird data to BBS strata region --------------------------------
# # this clipping is done only to speed up the overlay between BBS and ebird
#  the rest of the ebird distribution map is unaffected for the eventual
#  predictions of density across the full species' range (i.e., see
#  saveRDS() call above)
bbs_strata <- bbsBayes2::load_map("bbs_usgs") %>%
  st_buffer(., dist = 10000) #add 10km buffer for clipping eBird data

# project boundary to match raster data
region_boundary_proj <- st_transform(bbs_strata, st_crs(breed_abundance))
# crop and mask to boundary of BBS
breed_abundance <- crop(breed_abundance, region_boundary_proj) |>
  mask(region_boundary_proj)

breed_abundance_sd <- crop(breed_abundance_sd, region_boundary_proj) |>
  mask(region_boundary_proj)

breed_abundance_ci <- crop(breed_abundance_ci, region_boundary_proj) |>
  mask(region_boundary_proj)

# project buffers to match raster data
routes_buf_proj <- st_transform(routes_buf, st_crs(breed_abundance))



# extract mean abundance in BBS route buffers

abundance_in_buffers <- terra::extract(breed_abundance,
                                       routes_buf_proj,
                                       fun = mean,
                                       na.rm = TRUE,
                                       ID = FALSE,
                                       exact = TRUE)

abundance_sd_in_buffers <- terra::extract(breed_abundance_sd,
                                       routes_buf_proj,
                                       fun = mean,
                                       na.rm = TRUE,
                                       ID = FALSE,
                                       exact = TRUE)

abundance_ci_in_buffers <- terra::extract(breed_abundance_ci,
                                          routes_buf_proj,
                                          fun = mean,
                                          na.rm = TRUE,
                                          ID = FALSE,
                                          exact = TRUE)

# this argument creates an area-weighted mean of the cells that are overlapped
# by the buffer
# plot(breed_abundance)
# plot(routes_buf_proj,add = TRUE)
# plot(breed_abundance,add = TRUE)

# drop NA routes that are outside of the eBird prediction surface -----------------------

# these NA routes could be dropped because they don't represent
# much when it comes to scaling to the surface...

# may as well drop routes with no estimated abundance as they are outside
# of the range of the modeled seasonal distribution

abundance_df <- data.frame(route_name = routes_buf_proj$route_name,
                             ebird_abund = abundance_in_buffers[[1]],
                           ebird_abund_logsd = abundance_sd_in_buffers[[1]],
                           ebird_abund_ci = abundance_ci_in_buffers[[1]]) %>%
  filter(
    ebird_abund > 0,
    !is.na(ebird_abund)
    ) %>%
  distinct()



# ## displaying routes with no abundance data
# tstrt <- routes_buf %>%
#   inner_join(.,abundance_df) %>%
#   filter(is.na(ebird_abund))
#
# tst <- ggplot(data = tstrt)+
#   geom_sf(data = bbs_strata)+
#   geom_sf(aes(colour = ebird_abund))
# tst

saveRDS(abundance_df,
        paste0("data/species_relative_abundance/",species_ebird,"_relative_abundance.rds"))


}
