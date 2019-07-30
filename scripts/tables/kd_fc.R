# loading required libraries
library(tidyverse)
library(xtable)

# loading data
cebpb_kd_res <- read_rds('autoreg/data/cebpb_kd_res.rds')
pparg_kd_res <- read_rds('autoreg/data/pparg_kd_res.rds')

# defining variables
targets <- list('Adipogenic TF' = c('Cebpa', 'Cebpb', 'Pparg'),
                'Lipogenesis' = c('Lpl', 'Acly', 'Fasn'),
                'Autophagy TF' = c('Foxo1', 'Tfeb', 'Xbp1'),
                'Autophagy Gene' = c('Map1lc3b', 'Becn1', 'Sqstm1')) %>%
  reshape2::melt() %>%
  setNames(c('row', 'category'))

cebpb <- cebpb_kd_res %>%
  filter(row %in% targets$row) %>%
  filter(padj < .2) %>%
  dplyr::select(time, row, log2FoldChange, lfcSE) %>%
  mutate_at(vars(log2FoldChange, lfcSE), ~round(.x, 2)) %>%
  unite(change, log2FoldChange, lfcSE) %>%
  spread(time, change) %>%
  separate(`0`, c('cebpb_fc_0', 'cebpb_se_0'), sep = '_') %>%
  separate(`4`, c('cebpb_fc_4', 'cebpb_se_4'), sep = '_') %>%
  left_join(targets) %>%
  dplyr::select(category, everything())

pparg <- pparg_kd_res %>%
  filter(gene_id %in% as.character(targets$row)) %>%
  filter(adj.P.Val < .2,
         time %in% c(0, 2, 5)) %>%
  dplyr::select(time, row = gene_id, logFC, se) %>%
  mutate_at(vars(logFC, se), ~round(.x, 2)) %>%
  unite(change, logFC, se) %>%
  spread(time, change) %>%
  separate(`0`, c('pparg_fc_0', 'pparg_se_0'), sep = '_') %>%
  separate(`2`, c('pparg_fc_2', 'pparg_se_2'), sep = '_') %>%
  separate(`5`, c('pparg_fc_5', 'pparg_se_5'), sep = '_') %>%
  left_join(targets) %>%
  dplyr::select(category, everything())

header <- paste0("\\multirow{3}{*}{Category} & \\multirow{3}{*}{Gene} &",
                 "\\multicolumn{4}{c}{Cebpb KD vs Control} & \\multicolumn{6}{c}{Pparg KD vs Control}\\\\",
                 "\\cmidrule(lr){3-6}\\cmidrule(lr){7-12}",
                 "& & \\multicolumn{2}{c}{0 h} & \\multicolumn{2}{c}{4 h} & \\multicolumn{2}{c}{0 d} & \\multicolumn{2}{c}{2 d} & \\multicolumn{2}{c}{5 d}\\\\",
                 "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}\\cmidrule(lr){9-10}\\cmidrule(lr){11-12}",
                 "& & FC & SE & FC & SE & FC & SE & FC & SE & FC & SE\\\\")

full_join(cebpb, pparg) %>%
  arrange(category) %>%
  mutate(category = ifelse(duplicated(category), "", category)) %>%
  xtable(align = 'cllcccccccccc') %>%
  print(floating = FALSE,
        include.rownames = FALSE,
        booktabs = TRUE,
        sanitize.text.function = identity,
        comment = FALSE,
        include.colnames=FALSE,
        add.to.row = list(pos = list(0, 3, 5, 7),
                          command = c(header, rep('\\midrule ', 3))),
        file = 'manuscript/tables/kd_fc.tex')
