# Setup workspace
rm(list = ls())
library(tidyverse)
library(Racmacs)

map <- read.acmap("./data/map/220221_omicron_neut_full_map.ace")

# Run the dimensionality testing
mapDims <- dimensionTestMap(
  map,
  dimensions_to_test = 1:5,
  test_proportion = 0.1,
  minimum_column_basis = "none",
  fixed_column_bases = colBases(map),
  number_of_optimizations = 500,
  replicates_per_dimension = 100,
  options = list()
)

# Plot the figure
df <- data.frame(dimensions=c(1, 2, 3, 4, 5),
                 rmse=c(mapDims$mean_rmse_detectable))


ggplot(data=df, aes(x=dimensions, y=rmse)) +
  geom_line()+
  geom_point() +
  theme_bw() +
  xlab('Dimension') +
  ylab('Mean RMSE of detectable titers') +
  theme(strip.background = element_blank()) ->dp
dp

saveRDS(df, "./data/map/220221_omicron_neut_full_map-dim_test.rds")

