library(tidyverse)
library(ComplexHeatmap)

go_annotation <- read_rds('autoreg/data/go_annotation.rds')
factor_targets <- read_rds('autoreg/data/factor_targets.rds') 
factor_occupancy <- read_rds('autoreg/data/factor_occupancy.rds')

target_lists <- factor_targets %>%
  mutate(annotation = str_split(annotation, ' \\(', simplify = TRUE)[, 1]) %>%
  filter(annotation == 'Promoter',
         geneId %in% go_annotation$SYMBOL,
         factor %in% c('PPARG', 'MED1'),
         group != 'non') %>%
  group_by(factor, group) %>%
  summarise(targets = list(unique(geneId)))

comb_mat <- target_lists$targets %>%
  set_names(paste(target_lists$factor, target_lists$group, sep = '_')) %>%
  list_to_matrix() %>%
  make_comb_mat(mode = 'intersect')
UpSet(comb_mat[comb_degree(comb_mat) >= 2])
