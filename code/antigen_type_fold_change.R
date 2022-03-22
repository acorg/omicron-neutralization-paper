# Setup workspace
rm(list = ls())
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(Racmacs)

mapValues <- function(val_fn, name_fn) {
  function(map) {
    values <- val_fn(map)
    names(values) <- name_fn(map)
    values
  }
}

calc_fold_changes <- function(map, titer_table) {
  # Setup to store results
  all_group_results <- tibble(
    sr_group = character(0),
    ag = character(0),
    diff = numeric(0),
    diff_upper = numeric(0),
    diff_lower = numeric(0),
    homologous = logical(0)
  )
  
  # Append results for each group
  for (sr_group in names(homologous_ags)) {
    
    sr <- which(srGroups(map) == sr_group)
    homo_ag_name <- homologous_ags[sr_group]
    homo_ag <- match(homo_ag_name, agNames(map))
    sr_group_results <- tibble(
      sr_group = rep(sr_group, numAntigens(map)),
      ag = character(numAntigens(map)),
      diff = numeric(numAntigens(map)),
      diff_upper = numeric(numAntigens(map)),
      diff_lower = numeric(numAntigens(map)),
      homologous = logical(numAntigens(map))
    )
    
    for (ag in seq_len(numAntigens(map))) {
      homologous_titers <- titer_table[homo_ag, sr]
      ag_titers <- titer_table[ag, sr]
      titer_diff_est <- meantiter::mean_titer_diffs(
        titers1 = homologous_titers,
        titers2 = ag_titers,
        method = "truncated_normal",
        dilution_stepsize = 0
      )
      
      # Populate results
      sr_group_results$ag[ag] <- agNames(map)[ag]
      if (ag == homo_ag) {
        sr_group_results$diff[ag] <- 0
        sr_group_results$diff_upper[ag] <- 0
        sr_group_results$diff_lower[ag] <- 0
        sr_group_results$homologous[ag] <- TRUE
      } else {
        sr_group_results$diff[ag] <- titer_diff_est$mean_diff
        sr_group_results$diff_upper[ag] <- titer_diff_est$mean_diff_upper
        sr_group_results$diff_lower[ag] <- titer_diff_est$mean_diff_lower
        sr_group_results$homologous[ag] <- FALSE
      }
      
    }
    
    # Append the results
    all_group_results <- bind_rows(all_group_results, sr_group_results)
    
  }
  
  # Remove NA diffs
  all_group_results %>% 
    filter(
      !is.na(diff)
    ) -> all_group_results
  
  return(all_group_results)
  
}


map_lv<- read.acmap("./data/map/omicron_neut_LV_map.ace")
map_pv<- read.acmap("./data/map/omicron_neut_PV_map.ace")

agFillValues <- mapValues(agFill, agNames)

homologous_ags <- c(
  "3x Vax" = "WT",
  "2x Vax" = "WT",
  "WT conv" = "WT",
  "Alpha conv" = "B.1.1.7",
  "Beta conv" = "B.1.351",
  "Gamma conv" = "P.1",
  "Delta conv" = "B.1.617.2"
)


all_groups_lv <- calc_fold_changes(map_lv, titerTable(map_lv)) %>%
  mutate("Antigen type" = "Live-virus")
all_groups_pv <- calc_fold_changes(map_pv, titerTable(map_pv)) %>%
  mutate("Antigen type" = "Pseudovirus")

combo <- rbind(all_groups_lv, all_groups_pv)

# Cycle through serum groups
foldchange <- \(x) {
  
  xabs <- abs(x)
  foldchange <- 2^xabs
  foldchange[x < 0] <- -foldchange[x < 0]
  as.character(round(foldchange, 1))
  
}

plots <- list()
for (sr_group_name in unique(combo$sr_group)) {
  
  
  combo %>%
    filter(
      sr_group == sr_group_name
    ) -> sr_group_results
  
  temp_ag_levels <- unique(sr_group_results$ag[order(-sr_group_results$diff)])
  
  sr_group_results$x <- unlist(lapply(sr_group_results$ag, function(x) grep(x, temp_ag_levels)))
  sr_group_results %>%
    ggplot(
      aes(
        x = ifelse(`Antigen type` == "Live-virus", x - 0.15, x + 0.15),
        y = diff,
        ymin = diff_lower,
        ymax = diff_upper,
        color = ag,
        shape = `Antigen type`
      )
    ) + 
    geom_pointrange(
      show.legend = FALSE
    ) +
    scale_x_continuous(
      breaks = 1:length(temp_ag_levels),
      labels = temp_ag_levels
    ) + 
    scale_color_manual(
      values = agFillValues(map_lv)
    ) + 
    scale_y_continuous(
      breaks = 2:min(floor(combo$diff_lower), na.rm=T),
      labels = foldchange
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey40"
    ) + 
    coord_cartesian(
      ylim = c(min(floor(combo$diff_lower), na.rm=T), 2.5)
    ) +
    labs(
      x = "",
      y = "Fold drop compared to homologous",
      title = sr_group_name
    ) + 
    theme_bw() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      panel.grid.minor = element_blank()
    ) -> gp
  
  gp <- gp + 
    geom_text(data = sr_group_results %>% filter(`Antigen type` == "Live-virus"),
              aes(
                x = x - 0.2,
                y = 2,
                label = paste0("LV: ",foldchange(diff))
              ),
              size = 2.5,
              color = "grey20"
    ) +
    geom_text(data = sr_group_results %>% filter(`Antigen type` == "Pseudovirus"),
              aes(
                x = x + 0.2,
                y = 2,
                label = paste0("PV: ",foldchange(diff))
              ),
              size = 2.5,
              color = "grey20"
    )
  
  # Output the plot
  num_ags <- length(unique(sr_group_results$ag))
  plots <- c(plots, list(gp))
  
  
}

plots

