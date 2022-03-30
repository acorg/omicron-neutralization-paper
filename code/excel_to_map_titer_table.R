#setup page and load metadata
rm(list = ls())
library(tidyverse)
library(meantiter)
# --------------------------------------------- load required functions -----------------------------------------------

source("./functions/forest_plot_functions.R")


# ----------------------------------------------- read in data -----------------------------------------------
readxl::read_excel("./data/omicron_folddrops.xlsx") %>% mutate(
  log_fold_change = conv_to_logfold(`Titre drop`)
) -> forest_data

# ----------------------------------------------- prepare data -----------------------------------------------
forest_data <- forest_data %>% filter(manuscript_data =="y")

forest_data$`Comparator antigen`[forest_data$`Comparator antigen` %in% c("Wu-1", "WT", "WA1", "Wu-1?")] <- "WT"
forest_data$`Comparator antigen`[forest_data$`Comparator antigen` %in% c("D614G","B.1")] <- "D614G"

# factorise variables and check for NAs
forest_data <- factorise(forest_data)

# log fold change where log fold change is "large"
forest_data$log_fold_change[forest_data$`Titre drop` == "large"] <- log2(forest_data$TitersHAg[forest_data$`Titre drop` == "large"]) - log2(forest_data$TitersOmicron[forest_data$`Titre drop` == "large"])

# add log2 titers and log2 omicron
forest_data$Log2HAg[is.na(forest_data$Log2HAg)] <- log2(forest_data$TitersHAg[is.na(forest_data$Log2HAg)]/10)

# get omicron titer by subtracting log fold change from log2 Titer of comparator antigen
forest_data$Log2Omi <- forest_data$Log2HAg - forest_data$log_fold_change
forest_data$TitersOmicron[is.na(forest_data$TitersOmicron)] <- 2^(forest_data$Log2Omi[is.na(forest_data$TitersOmicron)])*10

# add variant data
forest_data <- melt_variant_comparator(forest_data, "Alpha")
forest_data <- melt_variant_comparator(forest_data, "Beta")
forest_data <- melt_variant_comparator(forest_data, "Gamma")
forest_data <- melt_variant_comparator(forest_data, "Delta")

# get omicron titre
omicron_data <- forest_data %>% filter(!(is.na(TitersOmicron)))
omicron_data$TitersHAg <- omicron_data$TitersOmicron
omicron_data$`Comparator antigen` <- "Omicron"

forest_data <- rbind(forest_data, omicron_data)

# change the sr group names
forest_data <- pretty_plot_names(forest_data)

# only keep data that has reported titers, not only fold drops
forest_data <- forest_data %>% filter(!is.na(TitersHAg))

# create serum name
forest_data$serum_name <- paste0(forest_data$standardise_encounters,"-", forest_data$vaccine_manufacturer, "-",forest_data$Standardised_sera_names, "-",  forest_data$time,"-",forest_data$`Sera details long`, "-", forest_data$Study)

# ----------------------------------------------- Create titer tables -----------------------------------------------

# If more than one omicron titer present, eg against both WT and D614G, then take the mean of these two measurements
mean_titer <- function(titers) {
  if(length(titers) == 0) {
    "*"
  } else if(length(titers) >1) {
    mean_t <- mean(log2(titers/10))
    2^mean_t*10
  } else {
    titers
  }
}

ag_names <- c("WT" ="WT", "D614G" = "D614G", "Delta" = "B.1.617.2", "Alpha" = "B.1.1.7", "Beta" = "B.1.351", "Gamma" = "P.1", "Omicron" = "B.1.1.529")

# ---------------- Full titer table -------------
forest_data %>% select(`Comparator antigen`, serum_name, TitersHAg) %>%
  pivot_wider(names_from = serum_name, values_from = TitersHAg, values_fn = mean_titer) %>%
  column_to_rownames("Comparator antigen")-> omicron_table

omicron_table <- as.data.frame(omicron_table)
omicron_table[is.na(omicron_table)] <- "*"

rownames(omicron_table) <- ag_names[match(rownames(omicron_table), names(ag_names))]
write.csv(x = omicron_table, file = "./data/titer_tables/omicron_neut_titer_table.csv")
