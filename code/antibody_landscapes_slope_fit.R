#setup page and load metadata
rm(list = ls())
library(Racmacs)
library(tidyverse)
library(meantiter)
set.seed(100)

#------ Function to fit slope per serum group, and map coordinates per serum
fit_cone_all_sera <- function(
  pars,
  ag_coords,
  log_titers,
  colbase
) {
  
  pred_titers_all <- unlist(lapply(1:length(colbase), function(ser) {
    coords <- c(pars[paste0("x",ser)], pars[paste0("y",ser)])
    ag_distances <- as.matrix(dist(rbind(coords, ag_coords)))[1, -1]
    predicted_logtiters <- colbase[ser] - ag_distances*pars["slope"]
    return(predicted_logtiters)
  }))
  
  
  sum((log_titers - pred_titers_all)^2, na.rm = T)
  
}


#--------------------- Read the base map and vaccine map
map <- read.acmap("data/map/omicron_neut_conv_only_map.ace")
map_orig <-  read.acmap("data/map/omicron_neut_full_map.ace")

sr_info <- unlist(lapply(srNames(map_orig), function(x) {
  strsplit(x, "-")[[1]][1]
}))

# extract vaccine serum group titers
map_vax <- subsetMap(map_orig, sera = srNames(map_orig)[sr_info == "3x Vax" | sr_info == "2x Vax"])

log_titers <- logtiterTable(map_vax)
colbase <- colBases(map_vax)

ag_coords <- unname(agCoords(map))



# Fit 2x Vax
log_2xvax <- log_titers[,grepl("2x Vax", srNames(map_vax))]
colb_2xvax <- colbase[grepl("2x Vax", srNames(map_vax))]

result_single_slope_2xvax <-  optim(
  par = c(
    x = ag_coords[apply(log_2xvax, 2, which.max),1],
    y = ag_coords[apply(log_2xvax, 2, which.max),2],
    slope = 1
  ),
  fn = fit_cone_all_sera,
  method = "L-BFGS-B",
  ag_coords = ag_coords,
  log_titers = log_2xvax,
  colbase = colb_2xvax
)


# Fit 3x vax
log_3xvax <- log_titers[,grepl("3x Vax", srNames(map_vax))]
colb_3xvax <- colbase[grepl("3x Vax", srNames(map_vax))]

result_single_slope_3xvax <-  optim(
  par = c(
    x = ag_coords[apply(log_3xvax, 2, which.max),1],
    y = ag_coords[apply(log_3xvax, 2, which.max),2],
    slope = 1
  ),
  fn = fit_cone_all_sera,
  method = "L-BFGS-B",
  ag_coords = ag_coords,
  log_titers = log_3xvax,
  colbase = colb_3xvax
)

saveRDS(result_single_slope_2xvax, file = "./data/landscape_fit/2xVax_sr_coord_group_slope_fit.rds")
saveRDS(result_single_slope_3xvax, file = "./data/landscape_fit/3xVax_sr_coord_group_slope_fit.rds")
