# Setup workspace
rm(list = ls())
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(Racmacs)
library(patchwork)

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
srGroupValues <- mapValues(srOutline, srGroups)

#View(titerTable(map_pv)[,srGroups(map_pv)=="3x Vax"])
# sr colors
#sr_colors <- unique(srOutline(map_lv))
#names(sr_colors) <- unique(as.character(srGroups(map_lv)))

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

do_fold_change_plot <- function(combo, sr_group_name) {
  combo %>%
    filter(
      sr_group == sr_group_name
    ) -> sr_group_results
  
  temp_ag_levels <- unique(sr_group_results$ag[order(-sr_group_results$diff)])
  
  sr_group_results$x <- unlist(lapply(sr_group_results$ag, function(x) grep(x, temp_ag_levels)))
  
  sr_group_results %>%
    ggplot(
      aes(
        x = x,#ifelse(`Antigen type` == "Live-virus", x - 0.15, x + 0.15),
        y = diff,
        ymin = diff_lower,
        ymax = diff_upper,
        shape = `Antigen type`,
        group = `Antigen type`
      )
    ) + 
    geom_line(aes(color = sr_group,
                  linetype = `Antigen type`)) +
    geom_pointrange(
      aes(color = ag),
      show.legend = FALSE,
      size = 0.5
    ) +
    scale_x_continuous(
      breaks = 1:length(temp_ag_levels),
      labels = temp_ag_levels
    ) + 
    scale_color_manual(
      values = c(agFillValues(map_lv), srGroupValues(map_lv))
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
      axis.title.y = element_text(
        size = 8
      ),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    ) -> gp
  
  gp <- gp + 
    geom_text(data = sr_group_results %>% filter(`Antigen type` == "Live-virus"),
              aes(
                x = x,
                y = 2,
                label = paste0("LV: ",foldchange(diff))
              ),
              size = 2,
              color = "grey20"
    ) +
    geom_text(data = sr_group_results %>% filter(`Antigen type` == "Pseudovirus"),
              aes(
                x = x,
                y = 2-0.4,
                label = paste0("PV: ",foldchange(diff))
              ),
              size = 2,
              color = "grey20"
    )
  
  return(gp)
}

do_ratio_plot <- function(combo, sr_group_name) {
  combo %>%
    filter(
      sr_group == sr_group_name
    ) %>%
    mutate(fc = as.numeric(foldchange(diff)))-> sr_group_results
  
  temp_ag_levels <- unique(sr_group_results$ag[order(-sr_group_results$diff)])
  
  sr_group_results$x <- unlist(lapply(sr_group_results$ag, function(x) grep(x, temp_ag_levels)))
  
  sr_group_results %>%
    select(sr_group, ag, `Antigen type`, fc, x) %>%
    pivot_wider(names_from = `Antigen type`, values_from = fc) %>%
    mutate(ratio_full = `Live-virus`/Pseudovirus,
           ratio = log2(abs(ratio_full))) -> ratio_df
          # ratio = ifelse(ratio_full < 0, -ratio, ratio)) -> ratio_df
  
  ratio_df %>%
    ggplot(
      aes(
        x = x,#ifelse(`Antigen type` == "Live-virus", x - 0.15, x + 0.15),
        y = ratio
      )
    ) + 
    geom_line(aes(color = sr_group)) +
    geom_point(
      aes(color = ag),
      show.legend = FALSE,
      size = 2
    ) +
    scale_x_continuous(
      breaks = 1:length(temp_ag_levels),
      labels = temp_ag_levels
    ) + 
    scale_color_manual(
      values = c(agFillValues(map_lv), srGroupValues(map_lv))
    ) + 
    scale_y_continuous(
      labels = function(x) round(2^x,1),
      breaks = seq(-2.5,2.5,0.5),
      limits = c(-2.7,2.7)
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey40"
    ) + 
    labs(
      x = "",
      y = "Fold drop LV/PV"
    ) + 
    theme_bw() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      axis.title.y = element_text(
        size = 8
      ),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    ) -> gp_ratio
  
  return(gp_ratio)
  
}

plots <- list()
for (sr_group_name in unique(combo$sr_group)) {
  
  gp <- do_fold_change_plot(combo, sr_group_name)
  gp_diff <- do_ratio_plot(combo, sr_group_name)
  
  combo_plot <- (gp + theme(axis.title.x=element_blank(),
                                       axis.text.x=element_blank(),
                                       axis.ticks.x=element_blank())) /gp_diff +
    plot_layout(heights = c(1.5, 1)) 

  plots <- c(plots, list(combo_plot))
  
}

ggpubr::ggarrange(plotlist = plots, labels = c("A", "B", "C", "D", "E", "F", "G")) -> all_plots
ggsave("./figures/fold_change_from_homologous.png", plot = all_plots, width = 12, height = 14)

