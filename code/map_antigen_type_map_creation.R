#setup page and load metadata
rm(list = ls())
library(Racmacs)
library(tidyverse)
library(meantiter)
set.seed(100)

source("./functions/forest_plot_functions.R")

# ---------------------------------- read in meta data
ag_colors_info <- read.csv("./data/map/ag_colors.csv", stringsAsFactors = FALSE)
sr_colors_info <- read.csv("./data/map/sr_colors.csv", stringsAsFactors = FALSE, sep = ";")
#alignment_map <- read.acmap("./data/map/wilks_et_al-map_ndsubset_no_outliers.ace")


# ----------------------------------------------- read in data -----------------------------------------------
readxl::read_excel("./data/omicron_folddrops.xlsx") %>% 
 filter(manuscript_data =="y") %>%
  filter(!is.na(TitersHAg)) %>%
  factorise() %>%
  pretty_plot_names() -> forest_data

forest_data$serum_name <- paste0(forest_data$standardise_encounters,"-", forest_data$vaccine_manufacturer, "-",forest_data$Standardised_sera_names, "-",  forest_data$time,"-",forest_data$Sera_details_no_time, "-", forest_data$Study)


pseudotype_sr <- forest_data %>% filter(standardised_assay == "PV") %>% pull(serum_name)
live_sr <- forest_data %>% filter(standardised_assay == "LV") %>% pull(serum_name)

# find study that's not distinguishable based on sr name and remove it
both_antigen_types <- pseudotype_sr[pseudotype_sr %in% live_sr]
pseudotype_sr <- unique(pseudotype_sr[!(pseudotype_sr %in% both_antigen_types)])
live_sr <- unique(live_sr[!(live_sr %in% both_antigen_types)])

# Map with all sera
full_map <- read.acmap("./data/map/omicron_neut_full_map.ace") 
pseudo_map <- subsetMap(full_map, sera = pseudotype_sr)
lv_map <- subsetMap(full_map, sera = live_sr)

save.acmap(pseudo_map, filename = "./data/map/omicron_neut_PV_map.ace")
save.acmap(lv_map, filename = "./data/map/omicron_neut_LV_map.ace")

                     
# Subset to only convalescent seraM
pseudo_map_conv <- subsetMap(pseudo_map, sera = srNames(pseudo_map)[grepl("conv", srNames(pseudo_map))])
pseudo_map_conv <- optimizeMap(pseudo_map_conv, number_of_dimensions = 2, number_of_optimizations = 500)
pseudo_map_conv <- realignMap(pseudo_map_conv, full_map)

lv_map_conv <- subsetMap(lv_map, sera = srNames(lv_map)[grepl("conv", srNames(lv_map))])
lv_map_conv <- optimizeMap(lv_map_conv, number_of_dimensions = 2, number_of_optimizations = 500)
lv_map_conv <- realignMap(lv_map_conv, full_map)

save.acmap(pseudo_map_conv, filename = "./data/map/omicron_neut_PV_conv_only_map.ace")
save.acmap(lv_map_conv, filename = "./data/map/omicron_neut_LV_conv_only_map.ace")

