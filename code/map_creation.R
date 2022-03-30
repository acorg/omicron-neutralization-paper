#setup page and load metadata
rm(list = ls())
library(Racmacs)
library(tidyverse)
library(meantiter)
set.seed(100)


# ---------------------------------- read in meta data
ag_colors_info <- read.csv("./data/map/ag_colors.csv", stringsAsFactors = FALSE)
sr_colors_info <- read.csv("./data/map/sr_colors.csv", stringsAsFactors = FALSE, sep = ";")
alignment_map <- read.acmap("./data/map/wilks_et_al-map_ndsubset_no_outliers.ace")


# ---------------------------------- map function
make_map <- function(table, ag_colors_info, sr_colors_info, alignment_map) {
  
  # make map
  map <- acmap(
    ag_names = rownames(table),
    sr_names = colnames(table),
    titer_table = table
  )
  dilutionStepsize(map) <- 0
  
  sr_groups_map <- unlist(lapply(colnames(titerTable(map)), function(x) {
    strsplit(x, "-")[[1]][1]
  }))
  
  
  #set serum groups
  srGroups(map) <- factor(
    sr_groups_map,
    levels = c("2x Vax","3x Vax","Inf + 2x Vax", "2x Vax + Inf", "WT conv", "Alpha conv", "Beta conv", "Gamma conv", "Delta conv")
  )
  
  #set antigen colours
  ag_colors <- ag_colors_info$Color[match(agNames(map), ag_colors_info$Antigen)]
  agFill(map) <- ag_colors
  
  #set Serum colours
  sr_colors <- sr_colors_info$Color[match(srGroups(map), sr_colors_info$Serum.group)]
  srOutline(map) <- sr_colors
  
  # Set styles
  srOutlineWidth(map) <- 1
  srSize(map) <- 10
  agSize(map) <- 18
  
  ptDrawingOrder(map) <- rev(ptDrawingOrder(map))
  
  
  map_optim <- optimizeMap(
    map,
    number_of_dimensions = 2,
    number_of_optimizations = 1000,
    minimum_column_basis = "none"
  )
  
  # realign map
  map_optim <- realignMap(map_optim, target_map = alignment_map)
  
  return(map_optim)
}

# ---------------------------------- map creation

# Map with all sera
titer_table <- read.titerTable("./data/titer_tables/omicron_neut_titer_table.csv")
map_optim <- make_map(titer_table, ag_colors_info = ag_colors_info, sr_colors_info = sr_colors_info, alignment_map = alignment_map)
save.acmap(map_optim, filename = "./data/map/omicron_neut_full_map.ace")


# Map with only convalescent sera
sr_info <- unlist(lapply(colnames(titer_table), function(x) str_split(x, "-")[[1]][1]))
titer_table_only_conv <- titer_table[,grepl("conv", sr_info)]
map_optim <- make_map(titer_table_only_conv, ag_colors_info = ag_colors_info, sr_colors_info = sr_colors_info, alignment_map = alignment_map)
save.acmap(map_optim, filename = "./data/map/omicron_neut_conv_only_map.ace")

