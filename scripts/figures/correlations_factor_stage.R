# loading required libraries
library(tidyverse)
library(reshape2)

# loading and cleaning data
deg_res <- read_rds('autoreg/data/deg_res.rds') %>%
  filter(contrast != 'late_vs_early') %>%
  dplyr::select(group = contrast, Gene1 = row, fc = log2FoldChange) %>%
  mutate(group = ifelse(group == 'early_vs_non', 'early', 'late'),
         dir = case_when(fc > 1 ~ 'UP',
                         fc < -1 ~ 'Down',
                         TRUE ~ 'None'))

ddcor <- read_rds('autoreg/data/ddcor.rds') %>%
  map(function(x) {
    list(early = filter(x$early) %>% dplyr::select(Gene1, Gene2, cor = early_cor),
         late  = filter(x$late) %>% dplyr::select(Gene1, Gene2, cor = late_cor)) %>%
      bind_rows(.id = 'group')
  }) %>%
  bind_rows()

# generating figure
(left_join(ddcor, deg_res) %>%
  na.omit() %>%
  filter(Gene2 != 'Ctcf') %>%
  mutate(group = case_when(group == 'early' ~ 'Early',
                           group == 'late' ~ 'Late')) %>%
  ggplot(aes(x = group, y = abs(cor), fill = dir)) +
  geom_boxplot() +
  lims(y = c(0, 1)) +
  facet_wrap(~Gene2, nrow = 1) +
  theme_bw() +
  theme(legend.position = 'top',
        panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(0,"null"))+
  labs(y = "Pearson's Correlation",
       x = 'Stage of Differentiation',
       fill = 'Regulation')) %>%
  ggsave(plot = .,
         'manuscript/figures/correlations_factor_stage.png',
         width = 20, height = 8, units = 'cm')
  