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
  
  if(length(colbase)>1) {
    pred_titers_all <- unlist(lapply(1:length(colbase), function(ser) {
      coords <- c(pars[paste0("x",ser)], pars[paste0("y",ser)])
      ag_distances <- as.matrix(dist(rbind(coords, ag_coords)))[1, -1]
      predicted_logtiters <- colbase[ser] - ag_distances*pars["slope"]
      return(predicted_logtiters)
    }))
  } else {
    coords <- c(pars["x"], pars["y"])
    ag_distances <- as.matrix(dist(rbind(coords, ag_coords)))[1, -1]
    pred_titers_all <- colbase - ag_distances*pars["slope"]
  }
  sum((log_titers - pred_titers_all)^2, na.rm = T)
  
}



#--------------------- Read the base map and vaccine map
map <- read.acmap("data/map/omicron_neut_conv_only_map.ace")
map_orig <-  read.acmap("data/map/omicron_neut_full_map.ace")

sr_info <- unlist(lapply(srNames(map_orig), function(x) {
  strsplit(x, "-")[[1]][1]
}))

# extract vaccine serum group titers
map_vax <- subsetMap(map_orig, sera = srNames(map_orig)[!(grepl("conv", srGroups(map_orig)))])

# add Omicron breakthrough sr group
sr_levels <- c(levels(srGroups(map_vax)), "Omicron breakthrough")
sr_groups <- as.character(srGroups(map_vax))
sr_groups[grepl("\\+ Omicron", srNames(map_vax))] <- "Omicron breakthrough"

srGroups(map_vax) <- factor(sr_groups, levels = sr_levels)

# get titers to fit
log_titers <- logtiterTable(map_vax)
colbase <- colBases(map_vax)

ag_coords <- unname(agCoords(map))


# fit for each serum group level
fit_results <- lapply(as.character(unique(srGroups(map_vax))), function(x) {
  log <-as.matrix(log_titers[,as.character(srGroups(map_vax)) == x])
  colb <- colbase[as.character(srGroups(map_vax)) == x]
  
  optim(
    par = c(
      x = ag_coords[apply(log, 2, which.max),1],
      y = ag_coords[apply(log, 2, which.max),2],
      slope = 1
    ),
    fn = fit_cone_all_sera,
    method = "L-BFGS-B",
    ag_coords = ag_coords,
    log_titers = log,
    colbase = colb
  )
})

names(fit_results) <- as.character(unique(srGroups(map_vax)))
saveRDS(fit_results, "./data/landscape_fit/multi_exposure_landscape_fits.rds")
