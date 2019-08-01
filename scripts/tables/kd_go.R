# loading required libraries
library(tidyverse)
library(xtable)

# loading data
kd_go <- read_rds('autoreg/data/kd_go_res.rds')

kd_go %>%
  filter(grepl('lipid|autophag', term),
         pvalue < .05) %>%
  mutate(category = ifelse(grepl('lipid', term), 'Lipid Metabolism', 'Autophagy'),
         term = paste0(term, '(', go_id, ')'),
         time = case_when(
           factor == 'cebpb' ~ paste0(time, 'h'),
           factor == 'pparg' ~ paste0(time, 'd')
         )) %>%
  group_by(category, term, factor) %>%
  summarise(times = paste0(time, collapse = '/')) %>%
  spread(factor, times) %>%
  ungroup() %>%
  mutate(category = ifelse(duplicated(category), '', category)) %>%
  setNames(c('Category', 'GO Term (ID)', 'Cebpb KD', 'Pparg KD')) %>%
  xtable(align = 'cllll') %>%
  print(floating = FALSE,
        include.rownames = FALSE,
        booktabs = TRUE,
        sanitize.text.function = identity,
        comment = FALSE,
        add.to.row = list(pos = list(7),
                          command = '\\midrule '),
        file = 'manuscript/tables/kd_go.tex')
