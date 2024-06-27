library(sf)
library(tidyverse)
library(ebirdst)
#library(patchwork)
library(terra)
library(bbsBayes2)

# download and extract eBird seasonal relative abundance within
# BBS route path buffers

routes_buf <- readRDS("data/all_routes_buffered.rds")

# eBird relative abundance ------------------------------------------------

sp_sel <- "Canada Warbler"


species_ebird <- ebirdst::get_species(sp_sel)

down <- try(ebirdst::ebirdst_download_status(species_ebird,
                                             download_ranges = FALSE,
                                             download_abundance = TRUE,
                                             download_occurrence = FALSE),
            silent = TRUE)


abd_seasonal_abundance <- ebirdst::load_raster(species = species_ebird,
                                               resolution = "3km",
                                               period = "seasonal",
                                               product = "abundance")  #27km low resolution
season_sel <- ifelse(any(grepl("resident",names(abd_seasonal_abundance))),"resident","breeding")

# Need to explore the resident vs breeding season abundance patterns
# At least for BCCH the resident relative abundance is not linearly related
# to the BBS mean observations (log-linear)
# Option: if resident, then extract the weekly abundance data for the
# BBS survey season and generate a breeding season abundance

breed_abundance <- abd_seasonal_abundance[[season_sel]]


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
                             ebird_abund = abundance_in_buffers[[season_sel]]) %>%
  filter(!is.na(ebird_abund),
         ebird_abund > 0)

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


