# ------------------------------------------------------------------------------
#-------------------------------- IMPORTANT ------------------------------------
# ------------------------------------------------------------------------------
# Restart R before running this script as other packages and Rmisc might interfere
# and result in inaccurate grouping 
# ------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
# ------------------------------------------------------------------------------


# Setup workspace
rm(list = ls())
library(Rmisc)
library(tidyverse)
library(ggplot2)
library(meantiter)
library(huxtable)

# function to remove NA values from table
clear_NA_from_table <- function(table) {
  levels_encounter <- levels(table$standardise_encounters)
  table <- as.data.frame(table)
  for(col in colnames(table)) {
    table[,col] <- gsub("NaN", "NA", table[,col])
    table[,col] <- gsub("NA\\\n\\(NA\\; NA\\)", "", table[,col])
    
  }
  
  table$standardise_encounters <- factor(table$standardise_encounters, levels = levels_encounter)
  return(table)
}

add_ba.1.1_stats <- function(table1, table2) {
  
  table1_rows <- match(interaction(table2$standardise_encounters, table2$OmicronVariant),
                       interaction(table1$standardise_encounters, table1$OmicronVariant))
  
  table2_rows <- 1:length(table1_rows)
  
  for(col in c("mean_fold_drop", "gmt_hAG", "gmt_Omic"))
    for(x in 1:length(table1_rows)) {
      table1[[table1_rows[x], col]] <- paste0(table1[[table1_rows[x], col]], "\n\n", table2[[table2_rows[x], col]])
    }
  
  return(table1)
}

#----------------------------------------------- set path to save -----------------------------------------------

path_to_save <- paste0("./tables/")
dir.create(path_to_save, recursive = TRUE)



# --------------------------------------------- load required functions -----------------------------------------------

source("./functions/forest_plot_functions.R")


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

# get the count of WT studies by serum group
forest_data %>% filter(`Comparator antigen` == "WT") %>%
  group_by(standardise_encounters) %>%
  count()

#-------------------- Prepare data frame for mean fold drops

mean_drop_df <- forest_data

# Add fold drop to WT for other antigens than Omicron
mean_drop_df %>%
  select(Study, Sourcelink, `Sera details long`,`Comparator antigen`, time, TitersHAg) -> fold_drop_wt

na_studies <- fold_drop_wt[is.na(fold_drop_wt$TitersHAg), "Study"]

fold_drop_wt %>% filter(!(Study %in% na_studies$Study)) %>%
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


#--------------- Calculate ratio of fold drops 2/3 and GMTs 3/2

mean_drop_df %>% 
  group_by(`Comparator antigen`,standardise_encounters) %>%
  mutate(mean_fold_drop =round(2^CI(log_fold_change)[["mean"]],1), 
         gmt = round(2^mean_titers(TitersHAg, method= "truncated_normal", dilution_stepsize = 0)$mean*10),
         gmt_o = round(2^mean_titers(TitersOmicron, method= "truncated_normal", dilution_stepsize = 0)$mean*10)) %>%
  select(standardise_encounters, `Comparator antigen`, mean_fold_drop,
         gmt, gmt_o) %>% 
  unique()-> mean_data_ratio

mean_data_ratio %>%
  pivot_wider(names_from = `Comparator antigen`, values_from = c("gmt", "gmt_o", "mean_fold_drop")) %>%
  column_to_rownames(var = "standardise_encounters")-> mean_data_ratio_wide


mean_data_ratio_wide["3x Vax", ]/mean_data_ratio_wide["2x Vax",] -> mean_3_2_ratio
mean_data_ratio_wide["2x Vax", ]/mean_data_ratio_wide["3x Vax",] -> mean_2_3_ratio


#------------------------------ Calculate mean tables ------------------------------

# -------------- Means by serum group --------------------
# not the same for WT when doing mean of fold drops and GMT WT/ GMT Omicron because some studies
# have only fold drops and no titers (e.g. Suzuki)

mean_drop_df %>% 
  group_by(`Comparator antigen`,standardise_encounters) %>%
  mutate(mean_fold_drop_m = round(2^CI(log_fold_change)[["mean"]],1), 
         mean_fold_drop_lower = round(2^CI(log_fold_change)[["lower"]],1), 
         mean_fold_drop_upper = round(2^CI(log_fold_change)[["upper"]],1), 
         mean_fold_drop_to_WT_m = round(2^CI(na.omit(fold_drop_to_WT))[["mean"]],1), 
         mean_fold_drop_to_WT_lower = round(2^CI(na.omit(fold_drop_to_WT))[["lower"]],1),
         mean_fold_drop_to_WT_upper = round(2^CI(na.omit(fold_drop_to_WT))[["upper"]],1),
         gmt = round(2^mean_titers(TitersHAg, method= "truncated_normal", dilution_stepsize = 0)$mean*10),
         gmt_lower = round(2^mean_titers(TitersHAg, method= "truncated_normal", dilution_stepsize = 0)$mean_lower*10),
         gmt_upper = round(2^mean_titers(TitersHAg, method= "truncated_normal", dilution_stepsize = 0)$mean_upper*10),
         gmt_o = round(2^mean_titers(TitersOmicron, method= "truncated_normal", dilution_stepsize = 0)$mean*10),
         gmt_lower_o = round(2^mean_titers(TitersOmicron, method= "truncated_normal", dilution_stepsize = 0)$mean_lower*10),
         gmt_upper_o = round(2^mean_titers(TitersOmicron, method= "truncated_normal", dilution_stepsize = 0)$mean_upper*10),
         mean_fold_drop = paste0(mean_fold_drop_m, "\n(", mean_fold_drop_lower, "; ",mean_fold_drop_upper,")"),
         mean_fold_drop_to_WT = paste0(mean_fold_drop_to_WT_m, "\n(", mean_fold_drop_to_WT_lower, "; ",mean_fold_drop_to_WT_upper,")"),
         gmt_hAG = paste0(gmt, "\n(", gmt_lower, "; ", gmt_upper, ")"),
         gmt_Omic = paste0(gmt_o, "\n(", gmt_lower_o, "; ", gmt_upper_o, ")")) %>%
  select(standardise_encounters, `Comparator antigen`, mean_fold_drop, mean_fold_drop_to_WT,
         gmt_hAG, gmt_Omic) %>% 
  unique()-> mean_data_encounter


mean_data_encounter <- mean_data_encounter[order(mean_data_encounter$`Comparator antigen`),]
mean_data_encounter %>%  select(!gmt_Omic) %>% 
  pivot_wider(names_from = c("Comparator antigen"), values_from = c("gmt_hAG", "mean_fold_drop", "mean_fold_drop_to_WT")) %>%
  select(!c("mean_fold_drop_to_WT_WT","mean_fold_drop_Omicron", "mean_fold_drop_to_WT_Omicron"))-> mean_encounters_wide

mean_encounters_wide <- mean_encounters_wide[rev(order(mean_encounters_wide$standardise_encounters)),]
mean_encounters_wide <- clear_NA_from_table(mean_encounters_wide)

# create pretty table
ht <- hux(mean_encounters_wide) %>%
  insert_row(colnames(mean_encounters_wide))

ht <- format_huxtable_cols_standard(ht)

ft <- as_flextable(ht)
if (TRUE) {
  my_doc <- officer::read_docx()
  my_doc <- flextable::body_add_flextable(
    my_doc, ft)
  print(my_doc, target =
          paste0(path_to_save,"/mean_encounters.docx"))
}



#-------------------- Mean by vaccine manufacturer

mean_drop_df %>% filter(standardise_encounters %in% c("2x Vax", "3x Vax")) %>%
  group_by(`Comparator antigen`, standardise_encounters, vaccine_manufacturer) %>%
  mutate(mean_fold_drop_m = round(2^CI(log_fold_change)[["mean"]],1), 
         mean_fold_drop_lower = round(2^CI(log_fold_change)[["lower"]],1), 
         mean_fold_drop_upper = round(2^CI(log_fold_change)[["upper"]],1), 
         mean_fold_drop_to_WT_m = round(2^CI(na.omit(fold_drop_to_WT))[["mean"]],1), 
         mean_fold_drop_to_WT_lower = round(2^CI(na.omit(fold_drop_to_WT))[["lower"]],1),
         mean_fold_drop_to_WT_upper = round(2^CI(na.omit(fold_drop_to_WT))[["upper"]],1),
         gmt = round(2^mean_titers(TitersHAg, method= "truncated_normal", dilution_stepsize = 0)$mean*10),
         gmt_lower = round(2^mean_titers(TitersHAg, method= "truncated_normal", dilution_stepsize = 0)$mean_lower*10),
         gmt_upper = round(2^mean_titers(TitersHAg, method= "truncated_normal", dilution_stepsize = 0)$mean_upper*10),
         gmt_o = round(2^mean_titers(TitersOmicron, method= "truncated_normal", dilution_stepsize = 0)$mean*10),
         gmt_lower_o = round(2^mean_titers(TitersOmicron, method= "truncated_normal", dilution_stepsize = 0)$mean_lower*10),
         gmt_upper_o = round(2^mean_titers(TitersOmicron, method= "truncated_normal", dilution_stepsize = 0)$mean_upper*10),
         mean_fold_drop = paste0(mean_fold_drop_m, "\n(", mean_fold_drop_lower, "; ",mean_fold_drop_upper,")"),
         mean_fold_drop_to_WT = paste0(mean_fold_drop_to_WT_m, "\n(", mean_fold_drop_to_WT_lower, "; ",mean_fold_drop_to_WT_upper,")"),
         gmt_hAG = paste0(gmt, "\n(", gmt_lower, "; ", gmt_upper, ")"),
         gmt_Omic = paste0(gmt_o, "\n(", gmt_lower_o, "; ", gmt_upper_o, ")")) %>%
  select(vaccine_manufacturer, standardise_encounters,`Comparator antigen`,
         mean_fold_drop, mean_fold_drop_to_WT,
         gmt_hAG, gmt_Omic) %>%
  unique()-> mean_data_vacc

mean_data_vacc <- mean_data_vacc[order(mean_data_vacc$`Comparator antigen`),]
mean_data_vacc %>%  select(!gmt_Omic) %>% 
  pivot_wider(names_from = c("Comparator antigen"), values_from = c("gmt_hAG", "mean_fold_drop", "mean_fold_drop_to_WT")) %>%
  select(!c("mean_fold_drop_to_WT_WT","mean_fold_drop_Omicron", "mean_fold_drop_to_WT_Omicron"))-> mean_vacc_wide

mean_vacc_wide <- mean_vacc_wide[rev(order(mean_vacc_wide$vaccine_manufacturer, mean_vacc_wide$standardise_encounters)),]
mean_vacc_wide <- mean_vacc_wide[rev(order(mean_vacc_wide$vaccine_manufacturer)),]

mean_vacc_wide <- clear_NA_from_table(mean_vacc_wide)


ht_vax <- hux(mean_vacc_wide) %>%
  insert_row(colnames(mean_vacc_wide))

ht_vax <- format_huxtable_cols_standard(ht_vax, single_grouping = FALSE, outer_group = "Vaccine")
ht_vax <- format_huxtable_row(ht_vax, data = mean_vacc_wide, column_name = "vaccine_manufacturer", rev = TRUE)

ft_vax <- as_flextable(ht_vax)
if (TRUE) {
  my_doc <- officer::read_docx()
  my_doc <- flextable::body_add_flextable(
    my_doc, ft_vax)
  print(my_doc, target =
          paste0(path_to_save,"/mean_encounters_vax.docx"))
}

#-------------------- Mean by assay
mean_drop_df %>% 
  group_by(`Comparator antigen`, standardise_encounters, standardised_assay) %>%
  mutate(mean_fold_drop_m = round(2^CI(log_fold_change)[["mean"]],1), 
         mean_fold_drop_lower = round(2^CI(log_fold_change)[["lower"]],1), 
         mean_fold_drop_upper = round(2^CI(log_fold_change)[["upper"]],1), 
         mean_fold_drop_to_WT_m = round(2^CI(na.omit(fold_drop_to_WT))[["mean"]],1), 
         mean_fold_drop_to_WT_lower = round(2^CI(na.omit(fold_drop_to_WT))[["lower"]],1),
         mean_fold_drop_to_WT_upper = round(2^CI(na.omit(fold_drop_to_WT))[["upper"]],1),
         gmt = round(2^mean_titers(TitersHAg, method= "truncated_normal", dilution_stepsize = 0)$mean*10),
         gmt_lower = round(2^mean_titers(TitersHAg, method= "truncated_normal", dilution_stepsize = 0)$mean_lower*10),
         gmt_upper = round(2^mean_titers(TitersHAg, method= "truncated_normal", dilution_stepsize = 0)$mean_upper*10),
         gmt_o = round(2^mean_titers(TitersOmicron, method= "truncated_normal", dilution_stepsize = 0)$mean*10),
         gmt_lower_o = round(2^mean_titers(TitersOmicron, method= "truncated_normal", dilution_stepsize = 0)$mean_lower*10),
         gmt_upper_o = round(2^mean_titers(TitersOmicron, method= "truncated_normal", dilution_stepsize = 0)$mean_upper*10),
         mean_fold_drop = paste0(mean_fold_drop_m, "\n(", mean_fold_drop_lower, "; ",mean_fold_drop_upper,")"),
         mean_fold_drop_to_WT = paste0(mean_fold_drop_to_WT_m, "\n(", mean_fold_drop_to_WT_lower, "; ",mean_fold_drop_to_WT_upper,")"),
         gmt_hAG = paste0(gmt, "\n(", gmt_lower, "; ", gmt_upper, ")"),
         gmt_Omic = paste0(gmt_o, "\n(", gmt_lower_o, "; ", gmt_upper_o, ")")) %>%
  select(standardise_encounters, standardised_assay, `Comparator antigen`,
         mean_fold_drop, mean_fold_drop_to_WT,
         gmt_hAG, gmt_Omic) %>%
  unique()-> mean_data_assay


mean_data_assay <- mean_data_assay[order(mean_data_assay$`Comparator antigen`),]
mean_data_assay %>%  select(!gmt_Omic) %>% 
  pivot_wider(names_from = c("Comparator antigen"), values_from = c("gmt_hAG", "mean_fold_drop", "mean_fold_drop_to_WT")) %>%
  select(!c("mean_fold_drop_to_WT_WT","mean_fold_drop_Omicron", "mean_fold_drop_to_WT_Omicron"))-> mean_assay_wide

mean_assay_wide <- mean_assay_wide[rev(order(mean_assay_wide$standardise_encounters, mean_assay_wide$standardised_assay)),]
mean_assay_wide <- clear_NA_from_table(mean_assay_wide)

ht_assay <- hux(mean_assay_wide) %>%
  insert_row(colnames(mean_assay_wide))

ht_assay <- format_huxtable_cols_standard(ht_assay, single_grouping = FALSE, outer_group = "Serum Group")
ht_assay[1,2] <- "Antigen type"
ht_assay <- format_huxtable_row(ht_assay, data = mean_assay_wide, column_name = "standardise_encounters", rev = TRUE)

ft_assay <- as_flextable(ht_assay)
if (TRUE) {
  my_doc <- officer::read_docx()
  my_doc <- flextable::body_add_flextable(
    my_doc, ft_assay)
  print(my_doc, target =
          paste0(path_to_save,"/mean_encounters_assay.docx"))
}

#-------------------- Mean by pseudovirus type
mean_drop_df %>% 
  filter(`Comparator antigen` %in% c("WT", "Delta", "Omicron")) %>%
  group_by(`Comparator antigen`, standardise_encounters, standardised_pseudo) %>%
  mutate(mean_fold_drop_m = round(2^CI(log_fold_change)[["mean"]],1), 
         mean_fold_drop_lower = round(2^CI(log_fold_change)[["lower"]],1), 
         mean_fold_drop_upper = round(2^CI(log_fold_change)[["upper"]],1), 
         mean_fold_drop_to_WT_m = round(2^CI(na.omit(fold_drop_to_WT))[["mean"]],1), 
         mean_fold_drop_to_WT_lower = round(2^CI(na.omit(fold_drop_to_WT))[["lower"]],1),
         mean_fold_drop_to_WT_upper = round(2^CI(na.omit(fold_drop_to_WT))[["upper"]],1),
         gmt = round(2^mean_titers(TitersHAg, method= "truncated_normal", dilution_stepsize = 0)$mean*10),
         gmt_lower = round(2^mean_titers(TitersHAg, method= "truncated_normal", dilution_stepsize = 0)$mean_lower*10),
         gmt_upper = round(2^mean_titers(TitersHAg, method= "truncated_normal", dilution_stepsize = 0)$mean_upper*10),
         gmt_o = round(2^mean_titers(TitersOmicron, method= "truncated_normal", dilution_stepsize = 0)$mean*10),
         gmt_lower_o = round(2^mean_titers(TitersOmicron, method= "truncated_normal", dilution_stepsize = 0)$mean_lower*10),
         gmt_upper_o = round(2^mean_titers(TitersOmicron, method= "truncated_normal", dilution_stepsize = 0)$mean_upper*10),
         mean_fold_drop = paste0(mean_fold_drop_m, "\n(", mean_fold_drop_lower, "; ",mean_fold_drop_upper,")"),
         mean_fold_drop_to_WT = paste0(mean_fold_drop_to_WT_m, "\n(", mean_fold_drop_to_WT_lower, "; ",mean_fold_drop_to_WT_upper,")"),
         gmt_hAG = paste0(gmt, "\n(", gmt_lower, "; ", gmt_upper, ")"),
         gmt_Omic = paste0(gmt_o, "\n(", gmt_lower_o, "; ", gmt_upper_o, ")")) %>%
  select(standardise_encounters, standardised_pseudo, `Comparator antigen`,
         mean_fold_drop, mean_fold_drop_to_WT,
         gmt_hAG, gmt_Omic) %>%
  unique()-> mean_data_pseudo

mean_data_pseudo <- mean_data_pseudo[order(mean_data_pseudo$`Comparator antigen`),]
mean_data_pseudo %>%  select(!gmt_Omic) %>% 
  pivot_wider(names_from = c("Comparator antigen"), values_from = c("gmt_hAG", "mean_fold_drop", "mean_fold_drop_to_WT")) %>%
  select(!c("mean_fold_drop_to_WT_WT","mean_fold_drop_Omicron", "mean_fold_drop_to_WT_Omicron"))-> mean_pseudo_wide

mean_pseudo_wide <- mean_pseudo_wide[rev(order(mean_pseudo_wide$standardise_encounters, mean_pseudo_wide$standardised_pseudo)),]
mean_pseudo_wide <- clear_NA_from_table(mean_pseudo_wide)

mean_pseudo_wide[is.na(mean_pseudo_wide)] <- ""
ht_pseudo <- hux(mean_pseudo_wide) %>%
  insert_row(colnames(mean_pseudo_wide))

ht_pseudo <- format_huxtable_cols_custom(ht_pseudo, single_grouping = FALSE, outer_group = "Serum Group", antigens = c("WT", "Delta", "Omicron"))
ht_pseudo[1,2] <- "Antigen type"
ht_pseudo <- format_huxtable_row(ht_pseudo, data = mean_pseudo_wide, column_name = "standardise_encounters", rev = TRUE)

ft_pseudo <- as_flextable(ht_pseudo)
if (TRUE) {
  my_doc <- officer::read_docx()
  my_doc <- flextable::body_add_flextable(
    my_doc, ft_pseudo)
  print(my_doc, target =
          paste0(path_to_save,"/mean_encounters_pseudo.docx"))
}


#---------------------- BA.1 BA.1.1 comparison
# look only at WT, once all and then in studies that titrated both
ba1_ba11_studies <- Reduce(intersect, list(mean_drop_df %>% filter(OmicronVariant == "BA.1") %>% pull(Study),
                                     mean_drop_df %>% filter(OmicronVariant == "BA.1+R346K") %>% pull(Study)))
mean_drop_df %>%
  filter(`Comparator antigen`=="WT") %>%
  group_by(standardise_encounters, OmicronVariant) %>%
  mutate(mean_fold_drop_m = round(2^CI(log_fold_change)[["mean"]],1), 
         mean_fold_drop_lower = round(2^CI(log_fold_change)[["lower"]],1), 
         mean_fold_drop_upper = round(2^CI(log_fold_change)[["upper"]],1), 
         mean_fold_drop_to_WT_m = round(2^CI(na.omit(fold_drop_to_WT))[["mean"]],1), 
         mean_fold_drop_to_WT_lower = round(2^CI(na.omit(fold_drop_to_WT))[["lower"]],1),
         mean_fold_drop_to_WT_upper = round(2^CI(na.omit(fold_drop_to_WT))[["upper"]],1),
         n = length(Study),
         gmt = round(2^mean_titers(TitersHAg, method= "truncated_normal", dilution_stepsize = 0)$mean*10),
         gmt_lower = round(2^mean_titers(TitersHAg, method= "truncated_normal", dilution_stepsize = 0)$mean_lower*10),
         gmt_upper = round(2^mean_titers(TitersHAg, method= "truncated_normal", dilution_stepsize = 0)$mean_upper*10),
         gmt_o = round(2^mean_titers(TitersOmicron, method= "truncated_normal", dilution_stepsize = 0)$mean*10),
         gmt_lower_o = round(2^mean_titers(TitersOmicron, method= "truncated_normal", dilution_stepsize = 0)$mean_lower*10),
         gmt_upper_o = round(2^mean_titers(TitersOmicron, method= "truncated_normal", dilution_stepsize = 0)$mean_upper*10),
         mean_fold_drop = paste0(mean_fold_drop_m," (n=", n, ")\n(", mean_fold_drop_lower, "; ",mean_fold_drop_upper,")"),
         mean_fold_drop_to_WT = paste0(mean_fold_drop_to_WT_m," (n=", n, ")\n(", mean_fold_drop_to_WT_lower, "; ",mean_fold_drop_to_WT_upper,")"),
         gmt_hAG = paste0(gmt," (n=", n, ")\n(", gmt_lower, "; ", gmt_upper, ")"),
         gmt_Omic = paste0(gmt_o, " (n=", n, ")\n(",gmt_lower_o, "; ", gmt_upper_o, ")")) %>%
  select(standardise_encounters,`Comparator antigen`, OmicronVariant,
         mean_fold_drop, gmt_hAG, gmt_Omic, n) %>%
  unique()-> mean_data_omicron_variant

mean_drop_df %>%
  filter(`Comparator antigen`=="WT") %>%
  filter(Study %in% ba1_ba11_studies) %>%
  group_by(standardise_encounters, OmicronVariant) %>%
  mutate(mean_fold_drop_m = round(2^CI(log_fold_change)[["mean"]],1), 
         mean_fold_drop_lower = round(2^CI(log_fold_change)[["lower"]],1), 
         mean_fold_drop_upper = round(2^CI(log_fold_change)[["upper"]],1), 
         mean_fold_drop_to_WT_m = round(2^CI(na.omit(fold_drop_to_WT))[["mean"]],1), 
         mean_fold_drop_to_WT_lower = round(2^CI(na.omit(fold_drop_to_WT))[["lower"]],1),
         mean_fold_drop_to_WT_upper = round(2^CI(na.omit(fold_drop_to_WT))[["upper"]],1),
         n = length(Study),
         gmt = round(2^mean_titers(TitersHAg, method= "truncated_normal", dilution_stepsize = 0)$mean*10),
         gmt_lower = round(2^mean_titers(TitersHAg, method= "truncated_normal", dilution_stepsize = 0)$mean_lower*10),
         gmt_upper = round(2^mean_titers(TitersHAg, method= "truncated_normal", dilution_stepsize = 0)$mean_upper*10),
         gmt_o = round(2^mean_titers(TitersOmicron, method= "truncated_normal", dilution_stepsize = 0)$mean*10),
         gmt_lower_o = round(2^mean_titers(TitersOmicron, method= "truncated_normal", dilution_stepsize = 0)$mean_lower*10),
         gmt_upper_o = round(2^mean_titers(TitersOmicron, method= "truncated_normal", dilution_stepsize = 0)$mean_upper*10),
         mean_fold_drop = paste0(mean_fold_drop_m," (n=", n, ")\n(", mean_fold_drop_lower, "; ",mean_fold_drop_upper,")"),
         mean_fold_drop_to_WT = paste0(mean_fold_drop_to_WT_m," (n=", n, ")\n(", mean_fold_drop_to_WT_lower, "; ",mean_fold_drop_to_WT_upper,")"),
         gmt_hAG = paste0(gmt," (n=", n, ")\n(", gmt_lower, "; ", gmt_upper, ")"),
         gmt_Omic = paste0(gmt_o, " (n=", n, ")\n(",gmt_lower_o, "; ", gmt_upper_o, ")")) %>%
  select(standardise_encounters,`Comparator antigen`, OmicronVariant,
         mean_fold_drop, gmt_hAG, gmt_Omic, n) %>%
  unique()-> mean_data_omicron_studies


mean_data_omicron <- add_ba.1.1_stats(mean_data_omicron_variant, mean_data_omicron_studies)

mean_data_omicron %>% 
  select(!n) %>%
  select(!`Comparator antigen`) %>%
  pivot_wider(names_from = c("OmicronVariant"), values_from = c("gmt_hAG","gmt_Omic", "mean_fold_drop")) -> mean_omicron_wide

mean_omicron_wide <- mean_omicron_wide[rev(order(mean_omicron_wide$standardise_encounters)),] 
mean_omicron_wide <- clear_NA_from_table(mean_omicron_wide)
mean_omicron_wide[is.na(mean_omicron_wide)] <- ""

ht <- hux(mean_omicron_wide) %>%
  insert_row(colnames(mean_omicron_wide))

ht %>%
  merge_cells(.,1, 2:5) %>%
  merge_cells(.,1, 6:7) %>%
  merge_cells(.,1:2, 1) -> ht

ht[1,2] <- "GMT"
ht[1,6] <- "Mean Omicron fold drop from WT"
ht[2, 2:7] <- c("WT BA.1", "WT BA.1+R346K", "BA.1", "BA.1+R346K", "BA.1", "BA.1+R346K")


ht <- format_huxtable_row(ht, data = mean_omicron_wide, column_name = "standardise_encounters", rev = TRUE)
ht[1,1] <- c("Serum Group")

ht %>% set_align(1, everywhere, "center") %>%
  set_align(-1, 1, ".") %>%
  set_all_padding(2) %>% 
  set_outer_padding(0.1) %>% 
  set_bold(row = 1, col = everywhere) %>% 
  set_bold(row = everywhere, col = 1) %>% 
  set_bottom_border(row = 1, col = everywhere) %>% 
  set_width(0.9) %>% 
  set_font_size(7) %>% 
  theme_article() %>%
  set_right_border(everywhere, 1, 1) %>% 
  set_right_border(everywhere, 5, 1) %>% 
  set_right_border_color(everywhere, 1, "grey")  %>%
  set_right_border_color(everywhere, 5, "grey")  -> ht

ft_omi <- as_flextable(ht)
if (TRUE) {
  my_doc <- officer::read_docx()
  my_doc <- flextable::body_add_flextable(
    my_doc, ft_omi)
  print(my_doc, target =
          paste0(path_to_save,"/mean_encounters_omicron_sublineage.docx"))
}


  

