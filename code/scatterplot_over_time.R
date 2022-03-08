# Setup workspace
rm(list = ls())
library(patchwork)
library(tidyverse)
library(ggplot2)
library(meantiter)
library(grid)
library(ggpubr)
# --------------------------------------------- load required functions -----------------------------------------------

source("./functions/forest_plot_functions.R")
path_to_save <- "./figures/"
file_ext <- "png"

# --------------------------------------------- plot functions -----------------------------------------------

scatter_over_time <- function(data, vars = c("log_fold_change", "Log2HAg", "Log2Omi")) {
  
  names <- c("log_fold_change" = "Fold change", "Log2HAg" = "Titer", "Log2Omi" = "Titer")
  plots <- lapply(vars, function(x){
    
    data %>% 
      ggplot(aes_string(x="mean_days", y = x)) +
      geom_point() + 
      geom_smooth(method='lm') +
      stat_regline_equation(aes(label = paste(..eq.label.., ..rr.label.., sep ="~~~")),
                            label.x.npc = "right",
                            label.y.npc = "top",
                            hjust =1) +
      facet_wrap(vars(standardise_encounters), nrow = 1) +
      scale_y_continuous(name = names[x], labels = function(y) 2^y, breaks = c(0:10), limits = c(0,10)) +
      scale_x_continuous(name = "Time (months)", labels = function(d) d/30, breaks = seq(from = 0, to =360, by = 30), limits = c(0, 360)) + 
      theme(strip.background =element_rect(fill="white"))
    
  })
  
  return(plots)
  
}

scatter_over_time_to_WT <- function(data, vars = c("Alpha", "Beta", "Gamma", "Delta")) {
  
  plots <- lapply(vars, function(x){
    
    data %>% filter(`Comparator antigen` == x) %>%
      ggplot(aes_string(x="mean_days", y = "fold_drop_to_WT")) +
      geom_point() + 
      geom_smooth(method='lm') +
      stat_regline_equation(aes(label = paste(..eq.label.., ..rr.label.., sep ="~~~")),
                            label.x.npc = "right",
                            label.y.npc = "top",
                            hjust =1) +
      scale_y_continuous(name = "Fold drop", labels = function(y) round(2^y), breaks = c(0:10), limits = c(0,10)) +
      facet_wrap(vars(standardise_encounters), nrow = 1) +
      scale_x_continuous(name = "Time (months)", labels = function(d) d/30, breaks = seq(from = 0, to =360, by = 30), limits = c(0, 360)) + 
      theme(strip.background =element_rect(fill="white"))
    
  })
  names(plots) <- vars
  return(plots)
  
}

scatter_over_time_WT_titers <- function(data, vars =  c("WT", "Alpha", "Beta", "Gamma", "Delta")) {
  
  plots <- lapply(vars, function(x){
    
    data %>% filter(`Comparator antigen` == x) %>%
      ggplot(aes_string(x="mean_days", y = "Log2HAg")) +
      geom_point() + 
      geom_smooth(method='lm') +
      stat_regline_equation(aes(label = paste(..eq.label.., ..rr.label.., sep ="~~~")),
                            label.x.npc = "right",
                            label.y.npc = "top",
                            hjust =1) +
      scale_y_continuous(name = "Titer", labels = function(y) round(2^y*10) , breaks = c(0:12), limits = c(0,12)) +
      facet_wrap(vars(standardise_encounters), nrow = 1) +
      scale_y_continuous(name = "Titer", labels = function(y) round(2^y*10) , breaks = c(0:12), limits = c(0,12)) +
      scale_x_continuous(name = "Time (months)", labels = function(d) d/30, breaks = seq(from = 0, to =360, by = 30), limits = c(0, 360)) + 
      theme(strip.background =element_rect(fill="white"))
    
  })
  names(plots) <- vars
  return(plots)
  
}


# ----------------------------------------------- read in data -----------------------------------------------
readxl::read_excel("./data/omicron_folddrops.xlsx") %>% mutate(
  log_fold_change = conv_to_logfold(`Titre drop`)
) -> forest_data

# ----------------------------------------------- prepare data -----------------------------------------------
forest_data <- forest_data %>% filter(shared_data =="y")
forest_data <- forest_data[!(is.na(forest_data$`Titre drop`)),]

# factorise variables and check for NAs
forest_data <- factorise(forest_data)

# log fold change where log fold change is "large"
forest_data$log_fold_change[forest_data$`Titre drop` == "large"] <- log2(forest_data$TitersHAg[forest_data$`Titre drop` == "large"]) - log2(forest_data$TitersOmicron[forest_data$`Titre drop` == "large"])

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

# rename all WT-like antigens 
forest_data$`Comparator antigen`[forest_data$`Comparator antigen` %in% c("Wu-1", "D614G","B.1", "Wu-1?", "WT", "WA1")] <- "WT"

forest_data <- pretty_plot_names(forest_data)

# add standardised time column
forest_data <- standardise_time_to_mean_days(forest_data)

# Add fold drops to WT

mean_drop_df <- forest_data

# Add fold drop to WT for other antigens than Omicron
mean_drop_df %>%
  select(Study, Sourcelink, `Sera details long`,`Comparator antigen`, time, TitersHAg) -> fold_drop_wt

na_studies <- fold_drop_wt[is.na(fold_drop_wt$TitersHAg), "Study"]

fold_drop_wt %>% filter(!(Study %in% na_studies)) %>%
  pivot_wider(names_from = `Comparator antigen`, values_from = TitersHAg) -> wide_fold_drop

wide_fold_drop %>%
  pivot_longer(cols = c("Alpha", "Beta","Gamma", "Delta")) %>%
  mutate(drop_to_wt = as.numeric(WT)/as.numeric(value)) %>%
  filter(!is.na(drop_to_wt))-> wide_fold_drop

# match fold drop to WT for other variant antigens
mean_drop_df$fold_drop_to_WT <- log2(wide_fold_drop$drop_to_wt[match(interaction(mean_drop_df[,c("Study", "Sera details long","Comparator antigen")]), 
                                                                     interaction(wide_fold_drop[,c("Study", "Sera details long","name")]))])

# match the rows to add Omicron, WT comparison in wide format
mean_data_wt <- mean_drop_df %>% filter(`Comparator antigen` == "WT")
mean_data_wt[,c("fold_drop_to_WT", "TitersHAg", "Log2HAg")] <- mean_data_wt[,c("log_fold_change", "TitersOmicron", "Log2Omi")]
mean_data_wt$`Comparator antigen` <- "Omicron"

# bind the wide variant to WT fold drop 
mean_drop_df <- rbind(mean_drop_df, mean_data_wt)

# --------------------- Titer and fold drop vs time since infection

# removed Omicron breakthrough infection Chen(HKU)
forest_data %>% 
  filter(!(Standardised_sera_names %in% forest_data$Standardised_sera_names[grepl("+ Omicron", forest_data$Standardised_sera_names)])) %>%
  filter(!`Comparator antigen` %in% c("Alpha", "Beta", "Gamma", "Delta")) %>%
  filter(standardise_encounters %in% c("WT conv", "2x Vax", "3x Vax", "2x Vax + Inf", "Inf + 2x Vax")) %>%
  filter(standardised_assay == "LV") %>% 
  scatter_over_time() -> plots_LV

forest_data %>% 
  filter(!(Standardised_sera_names %in% forest_data$Standardised_sera_names[grepl("+ Omicron", forest_data$Standardised_sera_names)])) %>%
  filter(!`Comparator antigen` %in% c("Alpha", "Beta", "Gamma", "Delta")) %>%
  filter(standardise_encounters %in% c("WT conv", "2x Vax", "3x Vax", "2x Vax + Inf", "Inf + 2x Vax")) %>%
  filter(standardised_assay == "PV") %>% 
  scatter_over_time() -> plots_PV


plots_LV[[1]] + plots_PV[[1]] + plot_layout(nrow = 2) + plot_annotation(tag_levels = 'A') -> fold_change_lm
ggsave(paste0(path_to_save, "scatter_fold_change_time.", file_ext), fold_change_lm, width = 12, height =6, dpi = 300)

# --------- do above with WT and Omicron titers

titers_WT_PV <- mean_drop_df %>%
  filter(standardised_assay == "PV") %>%
  filter(standardise_encounters %in% c("WT conv", "2x Vax", "3x Vax", "2x Vax + Inf", "Inf + 2x Vax")) %>%
  scatter_over_time_WT_titers(., vars = c("WT", "Omicron"))

titers_WT_LV <- mean_drop_df %>%
  filter(standardised_assay == "LV") %>%
  filter(standardise_encounters %in% c("WT conv", "2x Vax", "3x Vax", "2x Vax + Inf", "Inf + 2x Vax")) %>%
  scatter_over_time_WT_titers(., vars = c("WT", "Omicron"))

ggarrange(titers_WT_LV[[1]] + rremove("xlab"), titers_WT_LV[[2]] + rremove("xlab"),
          titers_WT_PV[[1]] + rremove("xlab"), titers_WT_PV[[2]], nrow=4, labels = c("A", "B", "C", "D")) -> WT_titers

ggsave(paste0(path_to_save, "scatter_time_titers_WT_assay.", file_ext), WT_titers, width = 12, height =12, dpi = 300)

