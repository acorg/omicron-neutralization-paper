##' Data investigation
###' Impact of digitized data on mean fold drop and GMT by bootstrap shared Data
#'
#'

#setup page and load metadata
rm(list = ls())
library(tidyverse)
library(ggplot2)
library(meantiter)
library(rstatix)
library(patchwork)

# --------------------------------------------- load required functions -----------------------------------------------

source("./functions/forest_plot_functions.R")


# define required functions
do_stat_test <- function(data, encounter, variable, non_normal, 
                         formulas = c("fold_change" = "log_fold_change ~ WD_fold_change", "omi_titers" = "Log2Omi ~ WD_titers","hAG_titers"= "Log2HAg ~ WD_titers")) {
  
  
  temp <- data %>%
    filter(standardise_encounters == encounter)
  if(variable == "fold_change" & encounter == "Gamma conv") {
    return(NULL)
  }
  if(encounter %in% as.character(non_normal)) {
    
    temp %>%
      wilcox_effsize(as.formula(formulas[variable])) %>%
      select(effsize, magnitude) -> magna
    
    temp %>%
      wilcox_test(as.formula(formulas[variable]))%>%
      mutate("test" = "Wilcoxon",
             "effect_size" = magna$effsize,
             "magnitude" = as.character(magna$magnitude),
             "variable" = variable,
             "serum_group" = encounter) -> res
    return(res)
  } else {
    temp %>% 
      cohens_d(as.formula(formulas[variable]))%>%
      select(effsize, magnitude) -> magna
    
    temp %>%
      t_test(as.formula(formulas[variable])) %>%
      mutate("test" = "T-Test",
             "effect_size" = magna$effsize,
             "magnitude" = as.character(magna$magnitude),
             "variable" = variable,
             "serum_group" = encounter) -> res
    
    return(res)
  }
  
}

combine_stat_tests <- function(stat, encounters = unique(as.character(forest_data$standardise_encounters))) {
  stat <- stat[unlist(lapply(stat, function(x) !is.null(x)))]
  stat <- bind_rows(stat)
  stat$label <- paste0(stat$test, " p = ", stat$p, "\n Effect size = ", round(as.numeric(stat$effect_size),2))
  
  return(stat)
}


# ----------------------------------------------- read in data -----------------------------------------------
readxl::read_excel("./data/omicron_folddrops.xlsx") %>% mutate(
  log_fold_change = conv_to_logfold(`Titre drop`)
) -> forest_data

# ----------------------------------------------- prepare data -----------------------------------------------
forest_data <- forest_data %>% filter(manuscript_data =="y")
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

# change to pretty plot names
forest_data <- pretty_plot_names(forest_data)

# rename all WT-like antigens 
forest_data$`Comparator antigen`[forest_data$`Comparator antigen` %in% c("Wu-1", "D614G","B.1", "Wu-1?", "WT", "WA1")] <- "WT"

# Webplotdigitizer column
WD_names <- c("n" = "No", "y" = "Fold drop & titers", "y (Titers)" = "Titers", "y (titers)" = "Titers", "y (Delta titers)" = "No")
forest_data$Webplotdigitizer <- as.character(WD_names[as.character(forest_data$Webplotdigitizer)])
forest_data$Webplotdigitizer <- factor(forest_data$Webplotdigitizer, level = unique(WD_names))

forest_data <- forest_data %>%
  mutate(WD_fold_change = unlist(lapply(as.character(forest_data$Webplotdigitizer), function(x) {
    if(as.character(x) == "Titers") {
      "No"
    } else if(as.character(x) == "No") {
      "No"
    } else {
      "Yes"
    }
  })),
  WD_titers = unlist(lapply(as.character(forest_data$Webplotdigitizer), function(x) {
    if(as.character(x) == "No") {
      "No"
    } else {
      "Yes"
    }
  })))

# remove Omicron breakthrough infections
forest_data <- forest_data[!(grepl("+ Omicron", forest_data$Standardised_sera_names)),]

  
# --------------------------------- Do statistical test fold changes ------------------------


# check for normal distribution of fold changes
forest_data %>%
  filter(`Comparator antigen` == "WT") %>%
  group_by(standardise_encounters) %>%
  summarise(shapiro_p_fold_change = shapiro_test(log_fold_change),
            shapiro_p_omi_titers = shapiro_test(Log2Omi),
            shapiro_p_hAg_titers = shapiro_test(Log2HAg)
            ) -> shapiro_stat

non_normal <- lapply(shapiro_stat[2:length(shapiro_stat)], function(x) {
  shapiro_stat$standardise_encounters[x$p.value < 0.05]
})


forms <- c("fold_change" = "log_fold_change ~ WD_fold_change", "omi_titers" = "Log2Omi ~ WD_titers","hAG_titers"= "Log2HAg ~ WD_titers")

# do statistical tests
# log fold change
fc_stat <-  lapply(unique(as.character(forest_data$standardise_encounters)), function(x) {
  do_stat_test(data = forest_data, encounter = x, variable = "fold_change", non_normal = non_normal$shapiro_p_fold_change, formulas = forms)
}) %>% combine_stat_tests()
        

hAG_titers_stat <-  lapply(unique(as.character(forest_data$standardise_encounters)), function(x) {
  do_stat_test(data = forest_data, encounter = x, variable = "hAG_titers", non_normal = non_normal$shapiro_p_hAg_titers, formulas = forms)
}) %>% combine_stat_tests()

omi_titers_stat <-  lapply(unique(as.character(forest_data$standardise_encounters)), function(x) {
  do_stat_test(data = forest_data, encounter = x, variable = "omi_titers", non_normal = non_normal$shapiro_p_omi_titers, formulas = forms)
}) %>% combine_stat_tests()


stat_tests <- rbind(fc_stat, hAG_titers_stat, omi_titers_stat)

# add label
forest_data$labels_fc <-fc_stat$label[match(forest_data$standardise_encounters, fc_stat$serum_group)]
forest_data$labels_hAG <-hAG_titers_stat$label[match(forest_data$standardise_encounters, hAG_titers_stat$serum_group)]
forest_data$labels_omi<- omi_titers_stat$label[match(forest_data$standardise_encounters, omi_titers_stat$serum_group)]


# do histogram plot of fold change
forest_data %>%
  filter(`Comparator antigen` == "WT") %>%
  ggplot() +
  geom_histogram(aes(x=log_fold_change, fill = WD_fold_change), color = "grey50") +
  geom_text(aes(label = labels_fc, x = 0, y = 9), hjust = 0, size = 3) +
  scale_fill_manual(name = "Digitized data", values = c("lightgreen", "orange")) + 
  facet_wrap(~standardise_encounters, nrow =1) +
  xlab("Log fold change") + 
  ylab("Count") + 
  theme_bw() + 
  theme(strip.background =element_rect(fill="white"),
        strip.text.x = element_text(size = 10),
        legend.position = "top") -> fc_hist


# do histogram plot of Log2 WT
forest_data %>%
  filter(`Comparator antigen` == "WT") %>%
  ggplot() +
  geom_histogram(aes(x=Log2HAg, fill = WD_titers), color = "grey50") +
  geom_text(aes(label = labels_hAG, x = 0, y = 10), hjust = 0, size = 3) +
  scale_fill_manual(name = "Digitized data", values = c("lightgreen", "orange")) + 
  facet_wrap(~standardise_encounters, nrow = 1) +
  xlab("Log2 WT Titers") + 
  ylab("Count") + 
  theme_bw() + 
  theme(strip.background =element_rect(fill="white"),
        strip.text.x = element_text(size = 10),
        legend.position = "top") -> fc_hAg

# do histogram plot of Log2 Omicron
forest_data %>%
  filter(`Comparator antigen` == "WT") %>%
  ggplot() +
  geom_histogram(aes(x=Log2Omi, fill = WD_titers), color = "grey50") +
  geom_text(aes(label = labels_omi, x = 0, y = 20), hjust = 0, size = 3) +
  scale_fill_manual(name = "Digitized data", values = c("lightgreen", "orange")) + 
  facet_wrap(~standardise_encounters, nrow = 1) +
  xlab("Log2 Omicron Titers") + 
  ylab("Count") + 
  theme_bw() + 
  theme(strip.background =element_rect(fill="white"),
        strip.text.x = element_text(size = 10),
        legend.position = "top") -> fc_omi

design <- "
1
2
2
3
3
4
4"
guide_area() + fc_hist  + fc_hAg +  fc_omi + plot_layout(design=design, guides = "collect") + plot_annotation(tag_levels = 'A') -> combined_plot

ggsave("./figures/digitized_data_distributions.png", plot = combined_plot, width = 12, height = 7, dpi = 300)
