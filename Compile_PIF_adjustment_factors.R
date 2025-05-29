
library(tidyverse)
#remotes::install_github("bbsBayes/bbsBayes2") # allows access to BBS raw data
library(bbsBayes2)


# combine PIF adjustment factors with na-pops estimates --------------------------------------------------------


Inputfiles.dir <- "Stanton_2019_code/input_files/" #paste(Inputfiles.dir, '/', sep="")
napops_species <- napops::list_species() |>
  dplyr::select(Species,Common_Name)

# Species specific adjustment factors (Time of Day, Detection Distance, Pair)
# csv file in next line is from the supplement to Stanton et al. 2019
# https://doi.org/10.5751/ACE-01331-140104 - appendix-2
ExpAdjs <- read.csv(paste(Inputfiles.dir, 'Species_PSEst_Adjustments.csv', sep=''), stringsAsFactors=FALSE) %>%
  mutate(cn = ifelse(cn == "McCown's Longspur", #updating some species names that have changed since 2019
                     "Thick-billed Longspur",
                     cn),
         cn = ifelse(cn == "Gray Jay",
                     "Canada Jay",
                     cn),
         cn = ifelse(grepl("Le Conte's",cn), # alternate spelling of LeConte
                     gsub("Le Conte's","LeConte's",cn),
                     cn))



# appendix with species names from BBS database --------



bbs_sci_names <- bbsBayes2::load_bbs_data()[["species"]] %>%
  filter(unid_combined) %>%
  select(aou, english,family,genus, order) %>%
  mutate(english = gsub(pattern = " (all forms)",
                               x = english,
                               replacement = "", fixed = TRUE))

# append with additional species PIF adjustments from Rosenberg et --------

Rosen_pif <- read_csv("data/Rosenberg_appendix_2_aaw1313_data_s2.csv") %>%
  rename(Dist2 = `Distance Adj.`,
         Pair2 = `Pair Adj.`,
         cn = Species) %>%
  mutate(cn = ifelse(cn == "McCown's Longspur", #updating some species names that have changed since 2019
                                            "Thick-billed Longspur",
                                            cn),
                                cn = ifelse(cn == "Gray Jay",
                                            "Canada Jay",
                                            cn),
                                cn = ifelse(grepl("Le Conte's",cn), # alternate spelling of LeConte
                                            gsub("Le Conte's","LeConte's",cn),
                                            cn)) %>%
  #filter(!cn %in% ExpAdjs$cn) %>%
  select(-sci_name)

ExpAdjs <- ExpAdjs %>%
  filter(!cn %in% Rosen_pif$cn)

ExpAdjs <- ExpAdjs %>%
  bind_rows(Rosen_pif)




# Additional TOD adjustments calculated for any species > 100 BBS  --------
# from "tod_calculations.R"


tod_adjs <- read_csv("Additional_time_of_day_adjustments.csv") %>%
  filter(!cn %in% ExpAdjs$cn)

ExpAdjs <- ExpAdjs %>%
  bind_rows(tod_adjs) %>%
  left_join(bbs_sci_names,
            by = c("cn" = "english"))




# generate combined adjustment scores for each species --------------------


new_adjs <- ExpAdjs %>%
  select(cn,aou, TimeAdj.meanlog,
         Dist2,
         Pair2,
         order,
         family,
         genus)



# calculate means by genus, family, and order -----------------------------


family_scores <- new_adjs %>%
  group_by(family) %>%
  summarise(family_TimeAdj.meanlog = mean(TimeAdj.meanlog, na.rm = TRUE),
            family_Dist2 = mean(Dist2, na.rm = TRUE),
            family_Pair2 = mean(Pair2, na.rm = TRUE))
genus_scores <- new_adjs %>%
  group_by(genus) %>%
  summarise(genus_TimeAdj.meanlog = mean(TimeAdj.meanlog, na.rm = TRUE),
            genus_Dist2 = mean(Dist2, na.rm = TRUE),
            genus_Pair2 = mean(Pair2, na.rm = TRUE))
order_scores <- new_adjs %>%
  group_by(order) %>%
  summarise(order_TimeAdj.meanlog = mean(TimeAdj.meanlog, na.rm = TRUE),
            order_Dist2 = mean(Dist2, na.rm = TRUE),
            order_Pair2 = mean(Pair2, na.rm = TRUE))


# append the group means --------------------------------------------------


new_adjs <- new_adjs %>%
  left_join(family_scores,
            by = "family") %>%
  left_join(genus_scores,
            by = "genus")%>%
  left_join(order_scores,
            by = "order")

missing_dist <- new_adjs %>%
  filter(is.na(order_Dist2))


missing_pair <- new_adjs %>%
  filter(is.na(order_Pair2))

missing_tod <- new_adjs %>%
  filter(is.na(order_TimeAdj.meanlog))

## only Wood Stork is missing an order-based mean



# replace missing values with genus then family then order means ----------


new_adjs_fixed <- new_adjs %>%
  mutate(Pair2_final = ifelse(is.na(Pair2),
                        genus_Pair2,
                        Pair2),
         Dist2_final = ifelse(is.na(Dist2),
                        genus_Dist2,
                        Dist2))

new_adjs_fixed <- new_adjs_fixed %>%
  mutate(Pair2_final = ifelse(is.na(Pair2_final),
                        family_Pair2,
                        Pair2_final),
         Dist2_final = ifelse(is.na(Dist2_final),
                        family_Dist2,
                        Dist2_final))

new_adjs_fixed <- new_adjs_fixed %>%
  mutate(Pair2_final = ifelse(is.na(Pair2_final),
                        order_Pair2,
                        Pair2_final),
         Dist2_final = ifelse(is.na(Dist2_final),
                        order_Dist2,
                        Dist2_final))


# Add values for Wood Stork -----------------------------------------------
new_adjs_fixed <- new_adjs_fixed %>%
  mutate(Pair2_final = ifelse(cn == "Wood Stork",
                              1, # likely no difference in detectability by sex
                              Pair2_final),
         Dist2_final = ifelse(cn == "Wood Stork",
                              400, # likely highly detectable given large size and soaring behaviour
                              Dist2_final),
         TimeAdj.meanlog_final = TimeAdj.meanlog)


# Calculate the species-level pif adjustments -----------------------------

PIF_adjs <- new_adjs_fixed %>%
  mutate(C_t = exp(TimeAdj.meanlog_final),
         C_d = (1/(50*(Dist2_final/1000)^2*3.14159)), # area within detection distance of 50 BBS stops (square km)
         C_p = Pair2_final,
         PIF_adj_birds_km2 = C_t * C_p * C_d) # species-level adjustment using the PIF factors from Rosenberg et al. 2019 (birds/km^2))



write_csv(PIF_adjs,"PIF_adjustments_and_napops_edr_avail.csv")


