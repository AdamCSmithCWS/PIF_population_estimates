library(sf)
library(tidyverse)
library(ebirdst)
library(patchwork)
library(terra)
library(bbsBayes2)






# BBS Routes --------------------------------------------------------------

## US routes - old unofficial shapefile created by John Sauer
comb <- readRDS("data/United_States_selected_BBS_route_paths.rds")



## Canadian routes
route_n <- function(prov,rt){
y <- vector("character",length(prov))
for(i in 1:length(prov)){
  rr <- as.integer(substr(rt[i],(nchar(rt[i])-2),nchar(rt[i])))
  y[i] <- paste(prov[i],rr,sep = "-")
}
return(y)
}

routes_can <- st_read(dsn = "data/Canada/ALLROUTES_2022.shp") %>%
  st_make_valid() #%>%

routes_can <- routes_can %>%
  mutate(route_name = route_n(Province,ProvRoute))




vie <- ggplot()+
  geom_sf(data = routes_can,
          aes(colour = Province))

vie

# Route buffering ---------------------------------------------------------



routes_can_buf <- routes_can %>%
  #filter(Province == 4) %>%
  group_by(route_name) %>%
  summarise(., do_union = FALSE)  %>%
  st_buffer(., dist = 1500)

buf_a <- as.numeric(st_area(routes_can_buf)/1e6)

routes_can_buf$area_km <- buf_a


vie <- ggplot()+
  geom_sf(data = routes_can_buf,
          aes(colour = area_km))

vie


# load route start locations ----------------------------------------------


routes_run <- bbsBayes2::load_bbs_data()[["routes"]] %>%
  filter(year > 2011) %>%
  mutate(route_name = paste(state_num,route,sep = "-"))


route_starts <- routes_run %>%
  select(route_name,latitude,longitude) %>%
  distinct() %>%
  filter(route_name %in% unique(routes_can_buf$route_name)) %>%
  st_as_sf(., coords = c("longitude",
                         "latitude"),crs = st_crs(4326)) %>%
  st_transform(.,crs = st_crs(routes_can_buf))

## routes_can_buf is now a 1km buffer of all route paths with surveys
routes_can_buf <- routes_can_buf %>%
  filter(route_name %in% route_starts$route_name)

vie <- ggplot()+
  geom_sf(data = routes_can_buf,
          aes(colour = area_km))+
  geom_sf(data = route_starts)

vie

## the area of each route is different and so the summed
## ebird abundance on each each route will vary just based on
## the length of the route path
## Therefore
## if summarising this information at the route-level
## it must be adjuste for the area of each buffered region
# BBS counts --------------------------------------------------------------

surveys <- routes_run %>%
  filter(route_name %in% unique(routes_can_buf$route_name))




bird_obs <- bbsBayes2::load_bbs_data()[["birds"]]%>%
  filter(year > 2011) %>%
  mutate(route_name = paste(state_num,route,sep = "-")) %>%
  filter(route_name %in% routes_can_buf$route_name,
         route_data_id %in% surveys$route_data_id) %>%
  select(route_data_id,aou,species_total) %>%
  distinct()

bird_incl <- bird_obs %>%
  select(aou) %>%
  distinct()

surveys_all <- surveys %>%
  expand_grid(.,bird_incl)

full_bird <- surveys_all %>%
  left_join(.,bird_obs,
            by = c("route_data_id",
                   "aou")) %>%
  mutate(count = ifelse(is.na(species_total),0,
                              species_total)) %>%
  select(route_data_id,
         route_name,
         aou,
         count,
         year) %>%
  distinct()

species_list <- bbsBayes2::load_bbs_data()[["species"]] %>%
  filter(unid_combined == FALSE)

mean_counts_route <- full_bird %>%
  group_by(route_name,aou) %>%
  summarise(n_years = n(),
            mean_count = mean(count,na.rm = TRUE),
            median_count = median(count,na.rm = TRUE),
            sd_count = sd(count,na.rm = TRUE)) %>%
  left_join(.,species_list,
            by = "aou")



raw_counts_route <- full_bird %>%
  left_join(.,species_list,
            by = "aou")

# eBird relative abundance ------------------------------------------------

sp_sel <- "American Robin"


species_ebird <- ebirdst::get_species(sp_sel)

sp_mean_counts <- mean_counts_route %>%
  filter(english == sp_sel)

sp_raw_counts <- raw_counts_route %>%
  filter(english == sp_sel)

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

breed_abundance <- abd_seasonal_abundance[[season_sel]]


# # clipping ebird data to BBS strata region --------------------------------
# bbs_strata <- bbsBayes2::load_map("bbs_usgs") %>%
#   st_buffer(., dist = 10000) #add 10km buffer for clipping eBird data
#
# # project boundary to match raster data
# region_boundary_proj <- st_transform(bbs_strata, st_crs(breed_abundance))
# # crop and mask to boundary of BBS
# # breed_abundance <- crop(breed_abundance, region_boundary_proj) |>
# #   mask(region_boundary_proj)

# project buffers to match raster data
routes_can_buf_proj <- st_transform(routes_can_buf, st_crs(breed_abundance))



# extract mean abundance in BBS route buffers

abundance_in_buffers <- terra::extract(breed_abundance,
                                routes_can_buf_proj,
                                fun = mean,
                                na.rm = TRUE,
                                ID = FALSE,
                                exact = TRUE)

abundance_join <- data.frame(route_name = routes_can_buf_proj$route_name,
                             ebird_abund = abundance_in_buffers[[season_sel]]) %>%
  inner_join(.,sp_raw_counts,
            by = "route_name")



abund_df <- abundance_join %>%
  group_by(route_name) %>%
  summarise(ebird_abund = mean(ebird_abund),
            mean_bbs_count = mean(count,na.rm = TRUE),
            n_years = n())

routes_w_ebird_bbs_abund <- routes_can_buf %>%
  full_join(abund_df,
             by = "route_name") %>%
  mutate(mean_bbs_count = ifelse(is.na(mean_bbs_count),
                                 1,
                                 mean_bbs_count+1)) %>%
  filter(ebird_abund > 0)


comp_plot  <- ggplot(data = abund_df,
                     aes(x = ebird_abund,
                         y = mean_bbs_count+1))+
  geom_point(aes(colour = n_years))+
  geom_smooth(method = "gam")+
  #scale_y_continuous(transform = "log10")+
  scale_colour_viridis_c(direction = -1)


comp_plot


tst1 <- ggplot()+
  geom_sf(data = routes_w_ebird_bbs_abund,
          aes(fill = ebird_abund,
              colour = ebird_abund))+
  scale_colour_viridis_c(aesthetics = c("colour","fill"),
                         trans = "log10")+
  labs(title = "log(Mean) eBird relative abund")
tst2 <- ggplot()+
  geom_sf(data = routes_w_ebird_bbs_abund,
          aes(fill = mean_bbs_count,
              colour = mean_bbs_count))+
  scale_colour_viridis_c(aesthetics = c("colour","fill"),
                         trans = "log10")+
  labs(title = "log(Mean) BBS counts")
tst3 <- ggplot()+
  geom_sf(data = routes_w_ebird_bbs_abund,
          aes(fill = n_years,
              colour = n_years))+
  labs(title = "Number BBS surveys")
print(tst1/tst2/tst3)# + plot_layout(guides = "collect"))



## reprojecting the abundance surface
breed_sel1 <- project(breed_abundance, crs(routes_w_ebird_bbs_abund), method = "near") |>
  # remove areas of the raster containing no data
  trim()


reg_sel <- bbsBayes2::load_map("prov_state") %>%
  filter(strata_name %in% c("YT")) %>%
  st_transform(.,crs = st_crs(routes_w_ebird_bbs_abund))

rts_sel <- routes_w_ebird_bbs_abund %>%
  st_join(.,reg_sel,
          left = FALSE)


#rts_sel <- st_transform(rts_sel, st_crs(breed_abundance))

breed_sel <- terra::crop(breed_sel1,rts_sel)



plot(breed_sel)
plot(routes_w_ebird_bbs_abund,add = TRUE, col = grey(0.7))

png("figures/AMRO_YT_example.png",
    width = 8,
    height = 8,
    units = "in",
    res = 300)
plot(breed_sel)
plot(routes_w_ebird_bbs_abund,add = TRUE, col = NA)
plot(reg_sel,add = TRUE,col = NA)
dev.off()





reg_sel <- bbsBayes2::load_map("prov_state") %>%
  filter(strata_name %in% c("NS")) %>%
  st_transform(.,crs = st_crs(routes_w_ebird_bbs_abund))

rts_sel <- routes_w_ebird_bbs_abund %>%
  st_join(.,reg_sel,
          left = FALSE)

#rts_sel <- st_transform(rts_sel, st_crs(breed_abundance))

breed_sel <- terra::crop(breed_sel1,rts_sel)



plot(breed_sel)
plot(routes_w_ebird_bbs_abund,add = TRUE, col = NA)

png("figures/AMRO_NS_example.png",
    width = 8,
    height = 8,
    units = "in",
    res = 300)
plot(breed_sel)
plot(routes_w_ebird_bbs_abund,add = TRUE, col = NA)
dev.off()




reg_sel <- bbsBayes2::load_map("prov_state") %>%
  filter(strata_name %in% c("AB")) %>%
  st_transform(.,crs = st_crs(routes_w_ebird_bbs_abund))

rts_sel <- routes_w_ebird_bbs_abund %>%
  st_join(.,reg_sel,
          left = FALSE)

#rts_sel <- st_transform(rts_sel, st_crs(breed_abundance))

breed_sel <- terra::crop(breed_sel1,rts_sel)



plot(breed_sel)
plot(routes_w_ebird_bbs_abund,add = TRUE, col = NA)

png("figures/AMRO_AB_example.png",
    width = 8,
    height = 8,
    units = "in",
    res = 300)
plot(breed_sel)
plot(routes_w_ebird_bbs_abund,add = TRUE, col = NA)
plot(reg_sel,add = TRUE,col = NA)
dev.off()



# Sum of eBird relative abundance -----------------------------------------

length(which(is.na(abundance_in_buffers$breeding)))

eBird_in_BBS_routes <- sum(abundance_join$ebird_abund,na.rm = TRUE)

# Sum of population estimates






















#
#
# # define projection
# crs_laea <- paste0("+proj=laea +lat_0=", 42.1,
#                    " +lon_0=", -120)
#
# tmpb <- tmpb %>%
#   st_transform(.,crs = st_crs(crs_laea))
#
# tmp <- tmp %>%
#   st_transform(.,crs = st_crs(crs_laea))
#
# tmpbb <- tmp %>%
#   st_buffer(., dist = 10000) %>%
#   st_bbox()
#
# # transform to the custom projection using nearest neighbor resampling
# breed_abundance_pr <- project(breed_abundance,
#                         crs_laea,
#                              method = "near")
#
#
# breed_abundance_df <- as.data.frame(breed_abundance_pr, xy = TRUE )
#
# vie <- ggplot()+
#   geom_tile(data = breed_abundance_df,
#               aes(x = x, y = y,
#                   fill = breeding))+
#   scale_fill_viridis_c()+
#   theme_void()
#   # geom_sf(data = tmpb)+
#   # geom_sf(data = tmp)+
#   # coord_sf(xlim = tmpbb[c("xmin","xmax")],
#   #          ylim = tmpbb[c("ymin","ymax")])
#
# vie
#
#
# pdf("rufu_hab_cartoon.pdf",
#     height = 22,width = 17)
# print(vie)
# dev.off()
#



