
# Setup workspace
rm(list = ls())
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(Racmacs)

# Read the map
map_wilks <- read.acmap("data/map/wilks_et_al-map_ndsubset_no_outliers.ace")
map <- read.acmap("./data/map/omicron_neut_full_map.ace")

# Get map plotting limits
lims <- Racmacs:::mapPlotLims(map_wilks, sera = T)#, padding = 0.3)
lims$xlim[1] <- lims$xlim[1] + 2
lims$ylim[1] <- lims$ylim[1] + 2
lims$ylim[2] <- lims$ylim[2] - 2

# Setup plotting function
doplot <- function(map, label) {
  
  # Setup the plot
  par(mar = rep(0.5, 4))
  
  # Plot the regular map
  srOutlineWidth(map) <- 2
  srOutline(map) <- adjustcolor(srOutline(map), alpha.f = 0.8)
  plot(map, xlim = lims$xlim, ylim = lims$ylim, fill.alpha = 0.9)
  
  # Plot labels
  label_adjustments <- matrix(0, numAntigens(map), 2)
  rownames(label_adjustments) <- agNames(map)
  label_adjustments["B.1.351",] <- c(0, 0.5)
  label_adjustments["P.1",] <- c(-0.7, 0)
  label_adjustments["B.1.1.7+E484K",] <- c(0, -0.3)
  label_adjustments["B.1.526+E484K",] <- c(0.9, 0)
  label_adjustments["B.1.621",] <- c(0.6, 0)
  label_adjustments["B.1.1.529",] <- c(0, -0.5)
  label_adjustments["C.37",] <- c(0, 0.5)
  label_adjustments["B.1.617.1",] <- c(0, 0.5)
  label_adjustments["B.1.1.7",] <- c(0, 0.5)
  label_adjustments["D614G",] <- c(-0.6, 0)
  label_adjustments["B.1.429",] <- c(0, 0.5)
  label_adjustments["B.1.617.2",] <- c(-0.6, -0.3)
  label_adjustments["B.1.617.2 (AY.3)+E484Q",] <- c(0.2, -0.3)
  label_adjustments["B.1.617.2 (AY.1)+K417N",] <- c(-0.95, 0)
  label_adjustments["B.1.617.2 (AY.2)+K417N",] <- c(0.6, -0.3)
  label_adjustments["B.1.617.2+K417N",] <- c(-0.8, 0)
  
  labels <- agNames(map)
  names(labels) <- agNames(map)
  labels["B.1.351"] <- "B.1.351\n(Beta)"
  labels["P.1"] <- "P.1\n(Gamma)"
  labels["B.1.621"] <- "B.1.621\n(Mu)"
  labels["B.1.1.529"] <- "B.1.1.529\n(Omicron)"
  labels["B.1.617.1"] <- "B.1.617.1\n(Kappa)"
  labels["C.37"] <- "C.37\n(Lambda)"
  labels["B.1.1.7"] <- "B.1.1.7\n(Alpha)"
  labels["B.1.429"] <- "B.1.429\n(Epsilon)"
  labels["B.1.617.2"] <- "B.1.617.2\n(Delta)"
  labels["B.1.526+E484K"] <- "B.1.526+E484K\n(Iota)"
  labels["B.1.617.2 (AY.3)+E484Q"] <- "B.1.617.2 (AY.3)+E484Q"
  labels["B.1.617.2 (AY.1)+K417N"] <- "B.1.617.2 (AY.1)+K417N"
  labels["B.1.617.2 (AY.2)+K417N"] <- "B.1.617.2 (AY.2)+K417N"
  labels["B.1.617.2+K417N"] <- "B.1.617.2+K417N"
  labels["B.1.1.7+E484K"] <- "B.1.1.7+E484K"
  
  label_size <- rep(0.7, numAntigens(map))
  names(label_size) <- agNames(map)
  label_size["B.1.1.7+E484K"] <- 0.6
  label_size["B.1.617.2 (AY.1)+K417N"] <- 0.6
  label_size["B.1.617.2 (AY.2)+K417N"] <- 0.6
  label_size["B.1.617.2 (AY.3)+E484Q"] <- 0.6
  label_size["B.1.617.2+K417N"] <- 0.6
  
  text(
    agCoords(map) + label_adjustments,
    cex = label_size,
    label = labels,
    font = 1
  )
  
  text(
    x = lims$xlim[1]+0.7, 
    y = lims$ylim[2]-0.7,
    cex = 3,
    label = label,
    font = 1
  )
}

# Setup plotting function
doplot_proc <- function(map, label) {
  
  # Setup the plot
  par(mar = rep(0.5, 4))
  
  # Plot the regular map
  srOutlineWidth(map) <- 2
  srOutline(map) <- adjustcolor(srOutline(map), alpha.f = 0.8)
  plot(map, xlim = lims$xlim, ylim = lims$ylim, fill.alpha = 0.9)
  
  text(
    x = lims$xlim[1]+0.7, 
    y = lims$ylim[2]-0.7,
    cex = 3,
    label = label,
    font = 1
  )
}


# Do the plots
pdf("figures/map/map_comparison.pdf", 7*1.4*1.1, 3*1.5*1.1)
layout(matrix(c(1,2), 1, 2))
doplot(map_wilks, "A")
doplot_proc(procrustesMap(map_wilks, map), "B")
dev.off()

png("figures/map/map_comparison.png", 7*1.4*1.1, 3*1.5*1.1, units = "in", res = 300)
layout(matrix(c(1,2), 1, 2))
doplot(map_wilks, "A")
doplot_proc(procrustesMap(map_wilks, map), "B")
dev.off()
