library(sf)
library(tidyverse)
library(patchwork)
library(terra)
library(bbsBayes2)


#############################################
#############################################
# NOTE: This script will not run
#
# The objects saved by this script are included so that subsequent scripts will run.
#
# The original GIS data required for this script are not official
# Therefore the output from this script is included in the data and code archive
# The script is included for transparency
# When the data for the route-paths is officially released, this script and the data
# archive will be updated.








# Script to prepare the BBS count data from the last 10-years and
# combine it with the spatial information on route-paths.
# Two key objects are saved at the end of the file:
# 1 - routes_buf - a polygon simple features object representing a
#       400 m buffer around each BBS route that has some observational
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
## a shapefile currated by the Canadian national BBS office bbs@ec.gc.ca
routes_can <- st_read(dsn = "data/Canada/ALLROUTES_2022.shp") %>%
  st_make_valid() #%>%

routes_can <- routes_can %>%
  mutate(route_name = route_n(Province,ProvRoute),
         prov_state = Province) %>%
select(route_name,prov_state,route)



# US routes ---------------------------------------------------------------

## sample of route-paths from unofficial gis data
## US routes - old unofficial shapefile archived online
## at https://earthworks.stanford.edu/catalog/stanford-vy474dv5024
routes_us <- readRDS("data/United_States_selected_BBS_route_paths.rds") %>%
  sf::st_transform(crs = sf::st_crs(routes_can)) %>%
  distinct() %>%
  st_make_valid()

## translate route text names to state*route numerical route_names
rts <- load_bbs_data(release = 2024)$routes %>%
  select(state_num,route,route_name,st_abrev,year) %>%
  group_by(state_num,route,route_name,st_abrev) %>%
  slice_min(order_by = year, n = 1)

## join with rts to add common route-names
routes_us <- routes_us %>%
  inner_join(rts, by = c("RTENAME" = "route_name",
                        "prov_state" = "st_abrev",
                        "first_year_surveyed" = "year")) %>%
  distinct() %>%
  mutate(route_name = paste(state_num,route,sep = "-")) %>%
  select(route_name,state_num,route) %>%
  rename(prov_state = state_num)






# Combine spatial data ----------------------------------------------------



routes_all <- routes_can %>%
  bind_rows(.,routes_us) %>%
  group_by(route_name,prov_state,route) %>%
  summarise()

lgth <- st_length(routes_all)

routes_all$length_km <- as.numeric(lgth/1000)
# routes_all <- routes_can

vie <- ggplot()+
  geom_sf(data = routes_all,
          aes(colour = prov_state))

vie

# Route buffering ---------------------------------------------------------

## 400m buffer around all BBS route paths

routes_buf <- routes_all %>%
  group_by(route_name) %>%
  summarise(., do_union = FALSE)  %>%
  st_buffer(., dist = 400) # 400 m

buf_a <- as.numeric(st_area(routes_buf)/1e6)

routes_buf$area_km <- buf_a

route_dat <- routes_all %>%
  st_drop_geometry() %>%
  distinct()

routes_buf <- routes_buf %>%
  inner_join(route_dat,by = "route_name")

vie <- ggplot()+
  geom_sf(data = routes_buf,
          aes(colour = area_km))

vie

#### some of the route paths are far too long
#### These must represent alternate route paths that were not
#### corrected.
#### These are ok: The analysis of the eBird data relies only on a
#### mean eBird relative abundance within the route buffer. So the
#### total length or area of the buffer is not critical information
####
rts_toolarge <- routes_buf %>%
  filter(area_km > 50)

rts_toolong <- routes_all %>%
  filter(route_name %in% rts_toolarge$route_name)

### Given these toolong and toolarge routes will not bias the results
### I've done nothing to correct them

# load route start locations ----------------------------------------------

routes_w_data <- load_bbs_data(release = 2024)$routes %>%
  select(-route_name) %>%
  filter(year > 2012) %>%  # 10 years from 2013 - 2023, missing 2020.
  mutate(route_name = paste(state_num,route,sep = "-")) %>%
  select(state_num,route,route_name,st_abrev,year,latitude,longitude) %>%
  group_by(state_num,route,route_name,st_abrev,latitude,longitude) %>%
  slice_min(order_by = year, n = 1) %>%
  ungroup()



route_starts <- routes_w_data %>%
  select(route_name,latitude,longitude) %>%
  distinct() %>%
  filter(route_name %in% unique(routes_buf$route_name)) %>%
  st_as_sf(., coords = c("longitude",
                         "latitude"),crs = st_crs(4326)) %>%
  st_transform(.,crs = st_crs(routes_buf))

## routes_buf is now an 400mm buffer of all route paths with surveys
routes_buf <- routes_buf %>%
  filter(route_name %in% route_starts$route_name)

# vie <- ggplot()+
#   geom_sf(data = routes_buf,
#           aes(colour = first_year))+
#   geom_sf(data = route_starts)
#
# vie

## the area of each route is different and so the summed
## ebird abundance on each each route will vary just based on
## the length of the route path
## Therefore
## if summarising this information at the route-level
## it must be adjuste for the area of each buffered region
# BBS counts --------------------------------------------------------------

surveys <- bbsBayes2::load_bbs_data(release = 2024)$routes %>%
  select(-route_name) %>%
  filter(year > 2012) %>%  # 10 years from 2013 - 2023, missing 2020.
  mutate(route_name = paste(state_num,route,sep = "-"),
         doy = lubridate::yday(paste(year,month,day,sep = "-"))) %>%
  filter(route_name %in% unique(routes_buf$route_name))


n_obs_summary <- surveys %>%
  select(route_name,obs_n,year) %>%
  distinct() %>%
  group_by(route_name,obs_n) %>%
  summarise(n_surveys = n())

n_obs_route <- n_obs_summary %>%
  group_by(route_name) %>%
  summarise(n_obs = n(),
            min_yrs_by_obs = min(n_surveys),
            max_yrs_by_obs = max(n_surveys))

n_routes_ob <- n_obs_summary %>%
  group_by(obs_n) %>%
  summarise(n_routes = n(),
            min_yrs_by_route = min(n_surveys),
            max_yrs_by_route = max(n_surveys))





bird_obs <- bbsBayes2::load_bbs_data()[["birds"]]%>%
  filter(year > 2012) %>%
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
                              species_total),
         strata_name = paste(country,st_abrev,bcr,sep ="-")) %>%
  select(route_data_id,
         route_name,
         aou,
         count,
         year,
         obs_n,
         doy,
         strata_name) %>%
  distinct()
unique(full_bird$strata_name)

species_list <- bbsBayes2::load_bbs_data()[["species"]] %>%
  filter(unid_combined == TRUE)

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
        "data/all_routes_buffered_400m.rds")
























