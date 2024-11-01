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



for(i in rev(1:nrow(sps_list))){

  sp_sel <- unname(unlist(sps_list[i,"english"]))
species_ebird <- ebirdst::get_species(sp_sel)

if(is.na(species_ebird)){
  sps_list[i,"available_ebird_data"] <- "Not in ebird species list"


  next

}

qual_sel <- qual_ebird[which(qual_ebird$species_code == species_ebird),]
breed_qual <- unname(unlist(qual_sel[,"breeding_quality"]))
resident_qual <- unname(unlist(qual_sel[,"resident_quality"]))
resident <- unname(unlist(qual_sel[,"is_resident"]))

sps_list[i,"is_resident"] <- resident
sps_list[i,"breeding_quality"] <- breed_qual
sps_list[i,"resident_quality"] <- resident_qual
sps_list[i,"species_ebird"] <- species_ebird

if(!resident & breed_qual == 0){
  sps_list[i,"available_ebird_data"] <- "Only nonbreeding abundance"


  next
}
if(!resident & breed_qual > 0) {
  down <- try(ebirdst::ebirdst_download_status(species_ebird,
                                             download_ranges = FALSE,
                                             download_abundance = TRUE,
                                             download_occurrence = FALSE,
                                             force = FALSE,
                                             pattern = "abundance_seasonal_mean_3km_"),
            silent = TRUE)
}

if(resident) {

down <- try(ebirdst::ebirdst_download_status(species_ebird,
                                             download_ranges = FALSE,
                                             download_abundance = TRUE,
                                             download_occurrence = FALSE,
                                             force = FALSE,
                                             pattern = "abundance_median_3km_"),
            silent = TRUE)
}

if(class(down) == "try-error"){

  sps_list[i,"available_ebird_data"] <- "Download fail"


  next}

sps_list[i,"available_ebird_data"] <- "Yes"





if(resident){

  # The resident relative abundance is not linearly related
  # to the BBS mean observations (log-linear) (e.g., BCCH)
  # Therefore: if resident, then extract the weekly abundance data for the
  # BBS survey season and generate a breeding season abundance

  abd_weekly_abundance <- ebirdst::load_raster(species = species_ebird,
                                               resolution = "3km",
                                               period = "weekly",
                                               product = "abundance")

  bbs_weeks <- try(as_date(names(abd_weekly_abundance)),silent = TRUE)
  bbs_weeks <- bbs_weeks[which(bbs_weeks > as_date("2022-05-24") &
                                 bbs_weeks < as_date("2022-07-07"))]

  abd_weekly_abundance_bbs_season <- abd_weekly_abundance[[as.character(bbs_weeks)]]

  breed_abundance <- terra::mean(abd_weekly_abundance_bbs_season)

}else{


  abd_seasonal_abundance <- ebirdst::load_raster(species = species_ebird,
                                                 resolution = "3km",
                                                 period = "seasonal",
                                                 product = "abundance")  #3km high resolution

  breed_abundance <- abd_seasonal_abundance[["breeding"]]
}




# # clipping ebird data to BBS strata region --------------------------------
bbs_strata <- bbsBayes2::load_map("bbs_usgs") %>%
  st_buffer(., dist = 10000) #add 10km buffer for clipping eBird data

# project boundary to match raster data
region_boundary_proj <- st_transform(bbs_strata, st_crs(breed_abundance))
# crop and mask to boundary of BBS
breed_abundance <- crop(breed_abundance, region_boundary_proj) |>
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
                             ebird_abund = abundance_in_buffers[[1]]) %>%
  filter(!is.na(ebird_abund),
         ebird_abund > 0)

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
