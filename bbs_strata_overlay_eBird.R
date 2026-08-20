library(tidyverse)
library(sf)


strata <- st_read("data/bcr2026_statprov.shp") |>
  filter(is_lake == 0) |>
  st_make_valid()

tst <- ggplot()+
  geom_sf(data = strata,
          aes(fill = bcr))

tst


centroids <- readxl::read_excel("data/LatLongs 27km pixel centroids CERW PIPL WOTH.xlsx") |>
  mutate(long = longitude,
         lat = latitude) |>
  st_as_sf( coords = c("long","lat"),crs = 4326) |>
  st_transform(crs = st_crs(strata)) |>
  st_join(strata,
          left = TRUE)


tst <- ggplot()+
  geom_sf(data = centroids,
          aes(colour = bcr))

tst

# identify non-overlap points ---------------------------------------------


miss <- centroids |>
  filter(is.na(bcr)) |>
  select(longitude, latitude, country)

nrst <- miss |>
  st_nearest_feature(strata)

strata_near <- strata[nrst,]

dist <- as.numeric(st_distance(miss,strata_near,by_element = TRUE))/1000

strata_near <- strata_near |>
  mutate(distance_to_strata_km = dist) |>
  st_drop_geometry()

miss <- miss |>
  bind_cols(strata_near)


final_set <- centroids |>
  filter(!is.na(bcr)) |>
  bind_rows(miss) |>
  mutate(distance_to_strata_km = ifelse(is.na(distance_to_strata_km),
                                        0,
                                        distance_to_strata_km))



tst <- ggplot()+
  geom_sf(data = final_set,
          aes(colour = bcr))

tst


tst2 <- ggplot()+
  geom_sf(data = final_set,
          aes(colour = distance_to_strata_km))+
  scale_color_binned(breaks = c(0,13.5,50,100,1000),
                     transform = "log1p")

tst2

final_set_out <- final_set |>
  st_drop_geometry()


write_excel_csv(final_set_out,
          "LatLongs 27km pixel centroids CERW PIPL WOTH with BCR prov state etc.csv")
