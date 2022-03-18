
# Setup workspace
rm(list = ls())
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(Racmacs)

# Read the map
map_full <- read.acmap("./data/map/omicron_neut_full_map.ace")
map_conv <- read.acmap("./data/map/omicron_neut_conv_only_map.ace")

# Get map plotting limits
lims <- Racmacs:::mapPlotLims(map_conv, sera = T, padding = 0.5)

# Setup plotting function
doplot <- function(map, label) {
  
  # Setup the plot
  par(mar = c(0.1,0.5,0.1,0.5))
  
  # Plot the regular map
  srOutlineWidth(map) <- 2
  srOutline(map) <- adjustcolor(srOutline(map), alpha.f = 0.8)
  plot(map, xlim = lims$xlim, ylim = lims$ylim, fill.alpha = 0.9)
  
  # Plot labels
  label_adjustments <- matrix(0, numAntigens(map), 2)
  rownames(label_adjustments) <- agNames(map)
  label_adjustments["B.1.351",] <- c(0.5, 0.7)
  label_adjustments["P.1",] <- c(0.7, -0.6)
  label_adjustments["B.1.1.529",] <- c(-0.4, -0.6)
  label_adjustments["B.1.1.7",] <- c(0, 1)
  label_adjustments["D614G",] <- c(0, -0.5)
  label_adjustments["WT",] <- c(-0.5, -0.2)
  label_adjustments["B.1.617.2",] <- c(0, -0.6)
 
  
  labels <- agNames(map)
  names(labels) <- agNames(map)
  labels["B.1.351"] <- "Beta\n(B.1.351)"
  labels["P.1"] <- "Gamma\n(P.1)"
  labels["B.1.617.2"] <- "Delta\n(B.1.617.2)"
  labels["B.1.1.529"] <- "Omicron\n(B.1.1.529/BA.1)"
  labels["B.1.1.7"] <- "Alpha\n(B.1.1.7)"
  
  
  label_size <- rep(1.2, numAntigens(map))
  names(label_size) <- agNames(map)
  
  text(
    agCoords(map) + label_adjustments,
    cex = label_size,
    label = labels,
    font = 1
  )
  
  text(
    x = lims$xlim[1]+0.5, 
    y = lims$ylim[2]-0.5,
    cex = 3,
    label = label,
    font = 1
  )
  
}

# Do the plots
pdf("./figures/map/omicron_neut_full_gmt.pdf", 7*1.4*1.1, 3*1.5*1.1)
layout(matrix(c(1,2), 1, 2))
doplot(map_full, label = "A")
doplot(map_conv, label = "B")
dev.off()

png("./figures/map/omicron_neut_full_gmt.png", 7*1.4*1.1, 3*1.5*1.1, units = "in", res = 300)
layout(matrix(c(1,2), 1, 2))
doplot(map_full, label = "A")
doplot(map_conv, label = "B")
dev.off()
