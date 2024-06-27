library(sf)
library(tidyverse)
#library(ebirdst)
library(patchwork)
library(terra)
library(bbsBayes2)

# Script to prepare the BBS count data from the last 10-years and
# combine it with the spatial information on route-paths.
# Two key objects are saved at the end of the file:
# 1 - routes_buf - a polygon simple features object representing a
#       1500 m buffer around each BBS route that has some observational
#       data during the ~10-year period
# 2 - mean_counts_route - a dataframe summarizing the mean observed counts
#       for each species on each route for the 10-year period. Also includes
#       columns for the variance of the observed counts and the number of years
#       surveys were conducted
# These two files are saved as rds files in data/ and are loaded for
# subsequent stages of the analysis.




# BBS Routes --------------------------------------------------------------


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



# US routes ---------------------------------------------------------------


## US routes - old unofficial shapefile created by John Sauer
routes_us <- readRDS("data/United_States_selected_BBS_route_paths.rds")





# Combine spatial data ----------------------------------------------------



# routes_all <- routes_can %>%
#   bind_rows(.,routes_us)

routes_all <- routes_can

vie <- ggplot()+
  geom_sf(data = routes_all,
          aes(colour = Province))

vie

# Route buffering ---------------------------------------------------------



routes_buf <- routes_all %>%
  group_by(route_name) %>%
  summarise(., do_union = FALSE)  %>%
  st_buffer(., dist = 1500)

buf_a <- as.numeric(st_area(routes_buf)/1e6)

routes_buf$area_km <- buf_a


vie <- ggplot()+
  geom_sf(data = routes_buf,
          aes(colour = area_km))

vie


# load route start locations ----------------------------------------------


routes_run <- bbsBayes2::load_bbs_data()[["routes"]] %>%
  filter(year > 2011) %>%
  mutate(route_name = paste(state_num,route,sep = "-"))


route_starts <- routes_run %>%
  select(route_name,latitude,longitude) %>%
  distinct() %>%
  filter(route_name %in% unique(routes_buf$route_name)) %>%
  st_as_sf(., coords = c("longitude",
                         "latitude"),crs = st_crs(4326)) %>%
  st_transform(.,crs = st_crs(routes_buf))

## routes_buf is now a 1km buffer of all route paths with surveys
routes_buf <- routes_buf %>%
  filter(route_name %in% route_starts$route_name)

vie <- ggplot()+
  geom_sf(data = routes_buf,
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
  filter(route_name %in% unique(routes_buf$route_name))




bird_obs <- bbsBayes2::load_bbs_data()[["birds"]]%>%
  filter(year > 2011) %>%
  mutate(route_name = paste(state_num,route,sep = "-")) %>%
  filter(route_name %in% routes_buf$route_name,
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

saveRDS(mean_counts_route,
        "data/mean_counts_by_route.rds")

raw_counts_route <- full_bird %>%
  left_join(.,species_list,
            by = "aou")

saveRDS(raw_counts_route,
        "data/all_counts_by_route.rds")

saveRDS(routes_buf,
        "data/all_routes_buffered.rds")






















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



