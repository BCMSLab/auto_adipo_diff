# loading required libraries
library(tidyverse)
library(SummarizedExperiment)

# loading data
transformed_counts <- read_rds('autoreg/data/transformed_counts.rds')
go_annotation <- read_rds('autoreg/data/go_annotation.rds')
tf_annotation <- read_rds('autoreg/data/tf_annotation.rds')

# define variables
tf <- c('Cebpb', 'Pparg', 'Rxrg', 'Ep300', 'Med1')

# generating figure
(list('All Genes' = rownames(transformed_counts),
     'Autophagy Genes' = rownames(transformed_counts) %in% go_annotation$SYMBOL,
     'Autophagy TF' = rownames(transformed_counts) %in% intersect(go_annotation$SYMBOL, tf_annotation$SYMBOL),
     "Adipogenic TF" = rownames(transformed_counts) %in% tf) %>%
  map(function(x) {
    cmdscale(dist(t(assay(transformed_counts)[x,]))) %>%
      as.data.frame() %>%
      mutate(group = transformed_counts$group,
             group = case_when(group == 'non' ~ 'Non',
                               group == 'early' ~ 'Early',
                               group == 'late' ~ 'Late'),
             group = factor(group, levels = c('Non', 'Early', 'Late')))
  }) %>%
  bind_rows(.id = 'type') %>%
  mutate(type = factor(type,
                       levels = c('All Genes', 'Autophagy Genes', 'Autophagy TF', 'Adipogenic TF'))) %>%
  ggplot(aes(x = V1, y = V2, color = group)) +
  geom_point(size = 1.5, alpha = .7) +
  facet_wrap(~type, nrow = 1, scales = 'free') +
  theme_bw() +
  theme(legend.position = 'top',
        panel.grid = element_blank(),
        strip.background = element_blank()) +
  labs(x = 'Dimension 1', y = 'Dimension 2', color = 'Stage of Differentiation')) %>%
  ggsave(plot = .,
         filename = 'manuscript/figures/mds_gene_group.png',
         height = 7, width = 20, units = 'cm')
