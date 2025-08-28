library(naturecounts)
library(tidyverse)


bbs <- naturecounts::nc_data_dl(collections = "BBS50-CAN",
                                username = "adam.smith",
                                info = "Analysis of stop-level BBS data")

bbs <- bbs %>%
  mutate(ObservationCount = as.integer(ObservationCount))
nrow(bbs)

locations <- bbs %>%
  select(SiteCode,latitude,longitude) %>%
  distinct()
nrow(locations)

locations_only <- bbs %>%
  select(latitude,longitude) %>%
  distinct()
nrow(locations_only)
