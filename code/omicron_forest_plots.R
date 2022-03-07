##' Omicron
###' Forest plots for omicron neutralization data

# Setup workspace
rm(list = ls())
library(labbook)
library(patchwork)
library(tidyverse)
library(ggplot2)
library(meantiter)
library(grid)
library(gtable)


#----------------------------------------------- set path to save -----------------------------------------------
fileext <- "png"
path_to_save <- paste0("./figures/", fileext,"/")
dir.create(path_to_save, recursive = TRUE)



# --------------------------------------------- load required functions -----------------------------------------------

source("./functions/forest_plot_functions.R")


# ----------------------------------------------- read in data -----------------------------------------------
readxl::read_excel("./data/omicron_folddrops.xlsx") %>% mutate(
  log_fold_change = conv_to_logfold(`Titre drop`)
) -> forest_data

# ----------------------------------------------- prepare data -----------------------------------------------
forest_data <- forest_data %>% filter(shared_data =="y")

forest_data <- forest_data[!(is.na(forest_data$`Titre drop`)),]

# factorise variables and check for NAs
forest_data <- factorise(forest_data)

# get proper rowlabel with monospacing
time_forest_data <- standardise_time_to_mean_days(forest_data)
time_forest_data <- character_rowlabel_standard_time(time_forest_data)

# add time rowlabel
forest_data <- character_rowlabel(forest_data)
forest_data$rowlabel <- time_forest_data$rowlabel


# log fold change where log fold change is "large"
forest_data$log_fold_change[forest_data$`Titre drop` == "large"] <- log2(forest_data$TitersHAg[forest_data$`Titre drop` == "large"]) - log2(forest_data$TitersOmicron[forest_data$`Titre drop` == "large"])

forest_data$arrow_length <- unlist(lapply(forest_data$`Titre drop`, function(x) {
  if(length(grep(">>|large", x))>0) {
    1
  } else if(length(grep(">", x))>0) {
    0.5
  } else {
    NA
  }
}))

# add log2 titers and log2 omicron
forest_data$Log2HAg[is.na(forest_data$Log2HAg)] <- log2(forest_data$TitersHAg[is.na(forest_data$Log2HAg)]/10)

# get omicron titer by subtracting log fold change from log2 Titer of comparator antigen
forest_data$Log2Omi <- forest_data$Log2HAg - forest_data$log_fold_change

# add Omicron titers
forest_data$TitersOmicron[is.na(forest_data$TitersOmicron)] <- 2^(forest_data$Log2Omi[is.na(forest_data$TitersOmicron)])*10

# melt variant titer data into comparator antigen
forest_data <- melt_variant_comparator(forest_data, "Alpha")
forest_data <- melt_variant_comparator(forest_data, "Beta")
forest_data <- melt_variant_comparator(forest_data, "Gamma")
forest_data <- melt_variant_comparator(forest_data, "Delta")

# set titers below 10 to 5 (log2 of below 0 to -1)
forest_data$Log2HAg[forest_data$Log2HAg < 0] <- -1
forest_data$Log2Omi[forest_data$Log2Omi < 0] <- -1

# add row long for plotting purposes
forest_data["Row_long"] <- 1:nrow(forest_data)

# change to pretty plot names
forest_data <- pretty_plot_names(forest_data)


# ------------------------------------------------ do the plots -----------------------------------------------

# order by antigen
forest_data <- reorder_data(forest_data[order(forest_data$`Comparator antigen`),], rev = TRUE)

#  order by overall drop
forest_data_order<- reorder_data(forest_data[order(forest_data$log_fold_change),], rev = TRUE)

# order by fold drop within serum groups
forest_data_encounter <- reorder_data(forest_data_order[order(forest_data_order$standardise_encounters),], rev = FALSE)
save_in_scaled_format(forest_data_encounter, which_plot = "titer_drop", to_save = "/omicron_encounter", hline_by = "standardise_encounters")

# save subfigures
save_in_scaled_format_encounter(forest_data_encounter, which_plot = "titer_drop", to_save = "/omicron_encounter", hline_by = "standardise_encounters")


#order by vaccine manufacturer and 2 vs 3x vax
forest_data_vacc_enc_manuf <- reorder_data(forest_data_order[order(forest_data_order$vaccine_manufacturer, forest_data_order$standardise_encounters),], rev = FALSE)
save_in_scaled_format(forest_data_vacc_enc_manuf, which_plot = "titer_drop", to_save = "/omicron_manuf_encounter", hline_by = c("vaccine_manufacturer", "standardise_encounters"))

# split by encounter, order by vaccine type
forest_data_encounter_vacc <- reorder_data(forest_data_encounter[order(forest_data_encounter$vacc_type, forest_data_encounter$standardise_encounters),], rev = FALSE)
save_in_scaled_format(forest_data_encounter_vacc, which_plot = "titer_drop", to_save = "/omicron_vacc_encounter", hline_by = c("vacc_type", "standardise_encounters"))

# split by encounter, order by assay
forest_data_assay <- reorder_data(forest_data_encounter[order(forest_data_encounter$standardise_encounters, forest_data_encounter$standardised_assay),], rev = FALSE)
save_in_scaled_format(forest_data_assay, which_plot = "titer_drop", to_save = "/omicron_assay_encounter", hline_by = c("standardise_encounters", "standardised_assay"))

# split by encounter, order by cell

forest_data_cell <- reorder_data(forest_data_encounter[order(forest_data_encounter$standardise_encounters, forest_data_encounter$standardised_cell),], rev = FALSE)
save_in_scaled_format(forest_data_cell, which_plot = "titer_drop", to_save = "/omicron_cell_encounter", hline_by = c("standardise_encounters", "standardised_cell"), height = 14)


# -------------------- Titer plots

forest_data_titers <- forest_data %>% filter(!is.na(TitersHAg))

#  order by titers against omicron
forest_data_titer_ordered <- reorder_data(forest_data_titers[order(forest_data_titers$Log2Omi),], rev = FALSE)

# order by titers and serum group
forest_data_titer_grouped <- reorder_data(forest_data_titer_ordered[order(forest_data_titer_ordered$standardise_encounters),], rev = FALSE)
save_in_scaled_format(forest_data_titer_grouped, which_plot = "not_fold_drop", to_save = "/omicron_encounters_omi", hline_by = "standardise_encounters", order_by_omicron = TRUE)

# ordered by titer against reference antigen, not Omicron
forest_data_titer_ordered_hAG <- forest_data_titers
forest_data_titer_ordered_hAG <- reorder_data(forest_data_titer_ordered_hAG[order(forest_data_titer_ordered_hAG$Log2HAg),], rev = FALSE)

# split by encounter
forest_data_titer_hAG_grouped <- reorder_data(forest_data_titer_ordered_hAG[order(forest_data_titer_ordered_hAG$standardise_encounters),], rev = FALSE)
save_in_scaled_format(forest_data_titer_hAG_grouped, which_plot = "not_fold_drop", to_save = "omicron_encounters_hAG", hline_by = "standardise_encounters")
save_in_scaled_format_encounter(forest_data_titer_hAG_grouped, which_plot = "not_fold_drop", to_save = "omicron_encounters_hAG", hline_by = "standardise_encounters")

# split by encounter, order by assay
forest_data_titer_hAG_encounter <- character_rowlabel_assay(forest_data_titer_hAG_grouped)
forest_data_titer_hAG_assay <- reorder_data(forest_data_titer_hAG_encounter[order(forest_data_titer_hAG_encounter$standardise_encounters, forest_data_titer_hAG_encounter$standardised_assay),], rev = FALSE)
save_in_scaled_format(forest_data_titer_hAG_assay, which_plot = "not_titer_drop", to_save = "/omicron_assay_encounter_hAG", hline_by = c("standardise_encounters", "standardised_assay"))



