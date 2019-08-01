# loading required libraries
library(tidyverse)
library(xtable)

# loading data
deg_go_res <- read_rds('autoreg/data/deg_go_res.rds')

# table header
header <- paste0("\\multirow{2}{*}{Category} & \\multirow{2}{*}{GO Term (ID)} &",
                 "\\multicolumn{2}{c}{Early vs Non} & \\multicolumn{2}{c}{Late vs Non} & \\multicolumn{2}{c}{Late vs Early}\\\\",
                 "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}",
                 "&& Ratio & P-value & Ratio & P-value & Ratio & P-value\\\\")

# generating tabls
deg_go_res %>%
  filter(grepl('lipid|autophag', term),
         over_represented_pvalue < .05,
         ontology != 'MF') %>%
  mutate(type = ifelse(grepl('lipid', term), 'Lipid Metabolism', 'Autophagy'),
         term = paste0(term, ' (', category, ')'),
         ratio  = round(numDEInCat/numInCat, 2),
         pvalue = ifelse(over_represented_pvalue < .01, '$<0.01$', round(over_represented_pvalue, 2))) %>%
  unite(value, ratio, pvalue) %>%
  dplyr::select(type, term, contrast, value) %>%
  spread(contrast, value) %>%
  separate(early_vs_non, c('en_ratio', 'en_pvalue'), sep = '_') %>%
  separate(late_vs_non, c('ln_ratio', 'ln_pvalue'), sep = '_') %>%
  separate(late_vs_early, c('le_ratio', 'le_pvalue'), sep = '_') %>%
  mutate(type = ifelse(duplicated(type), '', type)) %>%
  xtable(align = 'clp{.4\\textwidth}cccccc') %>%
  print(floating = FALSE,
        include.rownames = FALSE,
        booktabs = TRUE,
        sanitize.text.function = identity,
        comment = FALSE,
        include.colnames=FALSE,
        add.to.row = list(pos = list(0, 3),
                          command = c(header, '\\midrule ')),
        file = 'manuscript/tables/deg_go.tex')
