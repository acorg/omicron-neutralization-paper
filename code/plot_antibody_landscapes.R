rm(list = ls())
library(Racmacs)
library(ablandscapes)
library(tidyverse)
library(meantiter)
library(r3js)
set.seed(100)

#---------------- Load required functions ---------------------
source("./functions/landscape_remove_buttons.R")

#---------------- Set up base plane map ---------------------
map <- read.acmap("data/map/omicron_neut_conv_only_map.ace")
sr_group_colors <- read.csv("./data/map/sr_colors.csv", stringsAsFactors = FALSE, sep = ";")

plot_xlim <- c(-3, 4) + 1.5
plot_ylim <- c(-3, 4) - 0.4
plot_zlim <- c(-1, 10)

# Work out which sera are out of bounds
sr_x_coords <- srCoords(map)[,1]
sr_y_coords <- srCoords(map)[,2]
margin <- 0.2
sr_out_of_bound <- sr_x_coords < plot_xlim[1] + margin | 
  sr_x_coords > plot_xlim[2] - margin | 
  sr_y_coords < plot_ylim[1] + margin |
  sr_y_coords > plot_ylim[2] - margin

lndscp_xlim <- range(agCoords(map)[,1])
lndscp_ylim <- range(agCoords(map)[,2])

# Set viewing angles
angle <- list(
  rotation = c(-1.4370, 0.0062, -0.5350),
  translation = c(0.0344, 0.0459, 0.1175),
  zoom = 1.5646
)

#----------------------------------  Functions to plot Ab landscapes -------------------------
# serum group landscape
plot_sr_group_lndscp <- function(map, sr_group_logtiters, sr_group_cone_coords, sr_group_colbases, sr_group_cone_slopes, sr_group) {
  
  
  sr_group_mean_logtiters <- rowMeans(sr_group_logtiters, na.rm = T)
  sr_group_mean_logtiters[sr_group_mean_logtiters < plot_zlim[1]] <- plot_zlim[1]
  # Set map subset
  map_subset <- map
  
  # Plot the base plot
  data3js <- lndscp_r3jsplot(
    fit = list(acmap = map_subset),
    aspect.z = 0.5,
    show.surface = FALSE,
    show.titers = FALSE,
    output = "data3js",
    xlim = plot_xlim,
    ylim = plot_ylim,
    zlim = plot_zlim,
    show.sidegrid = TRUE,
    show.axis = FALSE,
    options = list(
      opacity.basemap.ags = 1,
      cex.basemap.ags = 3,
      cex.basemap.sr = 1.5,
      opacity.basemap.sr = 1
    )
  )
  
  # Get fitted surface
  grid_x_coords <- seq(from = lndscp_xlim[1], to = lndscp_xlim[2], by = 0.25)
  grid_y_coords <- seq(from = lndscp_ylim[1], to = lndscp_ylim[2], by = 0.25)
  grid_x_matrix <- matrix(grid_x_coords, length(grid_y_coords), length(grid_x_coords), byrow = T)
  grid_y_matrix <- matrix(grid_y_coords, length(grid_y_coords), length(grid_x_coords), byrow = F)
  grid_z_matrix <- matrix(NA, length(grid_y_coords), length(grid_x_coords))
  
  fit_lndscp_val <- function(x, y, sr_cone_coords, sr_colbases, sr_cone_slopes) {
    
    sr_distances <- as.matrix(dist(rbind(c(x, y), sr_cone_coords)))[1, -1]
    mean(sr_colbases - sr_distances*sr_cone_slopes)
    
  }
  
  grid_z_matrix[] <- vapply(
    seq_len(length(grid_z_matrix)),
    \(n) {
      fit_lndscp_val(
        grid_x_matrix[n], grid_y_matrix[n],
        sr_group_cone_coords,
        sr_group_colbases,
        sr_group_cone_slopes
      )
    }, numeric(1)
  )
  
  # Add the surface
  data3js <- r3js::surface3js(
    data3js,
    x = grid_x_matrix,
    y = grid_y_matrix,
    z = grid_z_matrix,
    col = sr_group_colors$Color[match(sr_group, sr_group_colors$Serum.group)],
    opacity = 0.8,
    toggle = "Mean landscape",
    wireframe = FALSE,
    doubleSide = TRUE
  )
  
  data3js <- r3js::surface3js(
    data3js,
    x = grid_x_matrix,
    y = grid_y_matrix,
    z = grid_z_matrix,
    col = adjustcolor(
      sr_group_colors$Color[match(sr_group, sr_group_colors$Serum.group)],
      red.f = 0.25,
      green.f = 0.25,
      blue.f = 0.25
    ),
    opacity = 0.8,
    toggle = "Mean landscape",
    wireframe = TRUE
  )
  
  # Add individual landscapes
  for (i in seq_along(sr_group_colbases)) {
    
    ## Calculate the individual surface
    sr_cone_coord <- sr_group_cone_coords[i,]
    sr_colbase <- sr_group_colbases[i]
    sr_cone_slope <- sr_group_cone_slopes
    
    grid_dists <- as.matrix(dist(rbind(
      sr_cone_coord,
      cbind(
        as.vector(grid_x_matrix),
        as.vector(grid_y_matrix)
      )
    )))[1, -1]
    
    grid_z_matrix <- matrix(NA, length(grid_y_coords), length(grid_x_coords))
    grid_z_matrix[] <- sr_colbase - grid_dists*sr_cone_slope
    
    # Add the individual surface
    data3js <- r3js::surface3js(
      data3js,
      x = grid_x_matrix,
      y = grid_y_matrix,
      z = grid_z_matrix,
      col = "grey70",
      opacity = 0.2,
      toggle = "Individual landscapes",
      wireframe = FALSE,
      doubleSide = TRUE
    )
    
    data3js <- r3js::surface3js(
      data3js,
      x = grid_x_matrix,
      y = grid_y_matrix,
      z = grid_z_matrix,
      col = "grey70",
      opacity = 0.4,
      toggle = "Individual landscapes",
      wireframe = TRUE
    )
    
  }
  # 
  # # Add the titers
  data3js <- lndscp3d_titers(
    data3js = data3js,
    object = list(
      coords = agCoords(map)[!is.na(sr_group_mean_logtiters),],
      logtiters = sr_group_mean_logtiters[!is.na(sr_group_mean_logtiters)],
      indices = which(!is.na(sr_group_mean_logtiters)),
      acmap = map
    ),
    zlim = plot_zlim,
    options = list(
      cex.titer = 1.5,
      col.impulse = "grey60"
    )
  )
  
  # Add the titer 50 plane
  x_grid <- seq(from = plot_xlim[1], to = plot_xlim[2], by = 0.5)
  y_grid <- seq(from = plot_ylim[1], to = plot_ylim[2], by = 0.5)
  z_grid <- matrix(log2(5), length(x_grid), length(y_grid))
  data3js <- r3js::surface3js(
    data3js,
    x = x_grid,
    y = y_grid,
    z = z_grid,
    col = "grey80",
    opacity = 0.2,
    toggle = "Titer 50"
  )
  
  data3js <- r3js::surface3js(
    data3js,
    x = x_grid,
    y = y_grid,
    z = z_grid,
    col = "grey40",
    opacity = 0.4,
    toggle = "Titer 50",
    wireframe = TRUE
  )
  
  # Draw border
  data3js <- r3js::lines3js(
    data3js,
    x = c(plot_xlim[1], plot_xlim[1], plot_xlim[2], plot_xlim[2], plot_xlim[1]),
    y = c(plot_ylim[1], plot_ylim[2], plot_ylim[2], plot_ylim[1], plot_ylim[1]),
    z = rep(plot_zlim[1], 5),
    lwd = 2,
    col = "grey70"
  )
  
  data3js <- remove_buttons(data3js)
  
  # Create html widget
  widget <- r3js(
    data3js = data3js,
    rotation = angle$rotation,
    translation = angle$translation,
    zoom = angle$zoom
  )
  
  htmlwidgets::onRender(
    widget,
    jsCode = paste0("function(el, x, data){
    el.style.outline = 'solid 2px #eeeeee';
    }")
  )
  
}

# GMT landscapes of multiple serum groups
plot_sr_group_lndscp_gmt <- function(map, landscape_pars, sr_groups) {
  
  sr_group <- sr_groups[1]
  
  sr_group_mean_logtiters <- rowMeans(landscape_pars[[sr_group]]$log_titers, na.rm = T)
 sr_group_mean_logtiters[sr_group_mean_logtiters < plot_zlim[1]] <- plot_zlim[1]
  # Set map subset
  map_subset <- map
  
  # Plot the base plot
  data3js <- lndscp_r3jsplot(
    fit = list(acmap = map_subset),
    aspect.z = 0.5,
    show.surface = FALSE,
    show.titers = FALSE,
    output = "data3js",
    xlim = plot_xlim,
    ylim = plot_ylim,
    zlim = plot_zlim,
    show.sidegrid = TRUE,
    show.axis = FALSE,
    options = list(
      opacity.basemap.ags = 1,
      cex.basemap.ags = 3,
      cex.basemap.sr = 1.5,
      opacity.basemap.sr = 1
    )
  )
  
  # Get fitted surface
  grid_x_coords <- seq(from = lndscp_xlim[1], to = lndscp_xlim[2], by = 0.25)
  grid_y_coords <- seq(from = lndscp_ylim[1], to = lndscp_ylim[2], by = 0.25)
  grid_x_matrix <- matrix(grid_x_coords, length(grid_y_coords), length(grid_x_coords), byrow = T)
  grid_y_matrix <- matrix(grid_y_coords, length(grid_y_coords), length(grid_x_coords), byrow = F)
  grid_z_matrix <- matrix(NA, length(grid_y_coords), length(grid_x_coords))
  
  fit_lndscp_val <- function(x, y, sr_cone_coords, sr_colbases, sr_cone_slopes) {
    
    sr_distances <- as.matrix(dist(rbind(c(x, y), sr_cone_coords)))[1, -1]
    mean(sr_colbases - sr_distances*sr_cone_slopes)
    
  }
  
  grid_z_matrix[] <- vapply(
    seq_len(length(grid_z_matrix)),
    \(n) {
      fit_lndscp_val(
        grid_x_matrix[n], grid_y_matrix[n],
        landscape_pars[[sr_group]]$sr_cone_coords,
        landscape_pars[[sr_group]]$colbases,
        landscape_pars[[sr_group]]$slope
      )
    }, numeric(1)
  )
  
  # Add the surface
  data3js <- r3js::surface3js(
    data3js,
    x = grid_x_matrix,
    y = grid_y_matrix,
    z = grid_z_matrix,
    col = sr_group_colors$Color[match(sr_group, sr_group_colors$Serum.group)],
    opacity = 0.8,
    toggle = "Mean landscape",
    wireframe = FALSE,
    doubleSide = TRUE
  )
  
  data3js <- r3js::surface3js(
    data3js,
    x = grid_x_matrix,
    y = grid_y_matrix,
    z = grid_z_matrix,
    col = adjustcolor(
      sr_group_colors$Color[match(sr_group, sr_group_colors$Serum.group)],
      red.f = 0.25,
      green.f = 0.25,
      blue.f = 0.25
    ),
    opacity = 0.8,
    toggle = "Mean landscape",
    wireframe = TRUE
  )
  
  # Add other gmt_surfaces
  for (i in 2:length(sr_groups)) {
    
    sr_group <- sr_groups[[i]]
    sr_group_mean_logtiters <- rowMeans(landscape_pars[[sr_group]]$log_titers, na.rm = T)
    sr_group_mean_logtiters[sr_group_mean_logtiters < plot_zlim[1]] <- plot_zlim[1]
    
    grid_z_matrix[] <- vapply(
      seq_len(length(grid_z_matrix)),
      \(n) {
        fit_lndscp_val(
          grid_x_matrix[n], grid_y_matrix[n],
          landscape_pars[[sr_group]]$sr_cone_coords,
          landscape_pars[[sr_group]]$colbases,
          landscape_pars[[sr_group]]$slope
        )
      }, numeric(1)
    )
    
    # Add the surface
    data3js <- r3js::surface3js(
      data3js,
      x = grid_x_matrix,
      y = grid_y_matrix,
      z = grid_z_matrix,
      col = sr_group_colors$Color[match(sr_group, sr_group_colors$Serum.group)],
      opacity = 0.8,
      toggle = "Mean landscape",
      wireframe = FALSE,
      doubleSide = TRUE
    )
    
    data3js <- r3js::surface3js(
      data3js,
      x = grid_x_matrix,
      y = grid_y_matrix,
      z = grid_z_matrix,
      col = adjustcolor(
        sr_group_colors$Color[match(sr_group, sr_group_colors$Serum.group)],
        red.f = 0.25,
        green.f = 0.25,
        blue.f = 0.25
      ),
      opacity = 0.8,
      toggle = "Mean landscape",
      wireframe = TRUE
    )
    
  }
  
  # Add the titer 50 plane
  x_grid <- seq(from = plot_xlim[1], to = plot_xlim[2], by = 0.5)
  y_grid <- seq(from = plot_ylim[1], to = plot_ylim[2], by = 0.5)
  z_grid <- matrix(log2(5), length(x_grid), length(y_grid))
  data3js <- r3js::surface3js(
    data3js,
    x = x_grid,
    y = y_grid,
    z = z_grid,
    col = "grey80",
    opacity = 0.2,
    toggle = "Titer 50"
  )
  
  data3js <- r3js::surface3js(
    data3js,
    x = x_grid,
    y = y_grid,
    z = z_grid,
    col = "grey40",
    opacity = 0.4,
    toggle = "Titer 50",
    wireframe = TRUE
  )
  
  # Draw border
  data3js <- r3js::lines3js(
    data3js,
    x = c(plot_xlim[1], plot_xlim[1], plot_xlim[2], plot_xlim[2], plot_xlim[1]),
    y = c(plot_ylim[1], plot_ylim[2], plot_ylim[2], plot_ylim[1], plot_ylim[1]),
    z = rep(plot_zlim[1], 5),
    lwd = 2,
    col = "grey70"
  )
  
  data3js <- remove_buttons(data3js)
  
  # Create html widget
  widget <- r3js(
    data3js = data3js,
    rotation = angle$rotation,
    translation = angle$translation,
    zoom = angle$zoom
  )
  
  htmlwidgets::onRender(
    widget,
    jsCode = paste0("function(el, x, data){
    el.style.outline = 'solid 2px #eeeeee';
    }")
  )
  
}

#-------------------------------------------- Do the plots----------------

# Load map with 2x Vax and 3x Vax tritations
map_orig <-  read.acmap("data/map/omicron_neut_full_map.ace")

sr_info <- unlist(lapply(srNames(map_orig), function(x) {
  strsplit(x, "-")[[1]][1]
}))

# subset to only vaccine sera
map_vax <- subsetMap(map_orig, sera = srNames(map_orig)[sr_info == "3x Vax" | sr_info == "2x Vax"])

log_titers <- logtiterTable(map_vax)
colbase <- colBases(map_vax)

sr_groups <- c("2x Vax" = "2xVax", "3x Vax" = "3xVax")

# load parameters
landscape_paramters <- lapply(names(sr_groups), function(x){
  fit <- readRDS(paste0("./data/landscape_fit/",sr_groups[x],"_sr_coord_group_slope_fit.rds"))
  colb <- colbase[grepl(x, srNames(map_vax))]
  log <- log_titers[,grepl(x, srNames(map_vax))]
  
  pars <- fit$par
  slope <- pars["slope"]
  
  sr_cone_coords <- as.matrix(cbind(pars[grepl("x", names(pars))], pars[grepl("y", names(pars))]))
  
  list("sr_cone_coords" = sr_cone_coords, "colbases" = colb, "log_titers" = log, "slope" = slope)
})
names(landscape_paramters) <- names(sr_groups)

#- plot GMT landscapes
plot_sr_group_lndscp_gmt(map, landscape_paramters, names(sr_groups))

# Plot single serum group landscapes
# 2x vax
plot_sr_group_lndscp(map, landscape_paramters$`2x Vax`$log_titers, landscape_paramters$`2x Vax`$sr_cone_coords, landscape_paramters$`2x Vax`$colbases, landscape_paramters$`2x Vax`$slope, "2x Vax")

# 3x vax
plot_sr_group_lndscp(map, landscape_paramters$`3x Vax`$log_titers, landscape_paramters$`3x Vax`$sr_cone_coords, landscape_paramters$`3x Vax`$colbases, landscape_paramters$`3x Vax`$slope, "3x Vax")

