library(tidyverse)
library(reshape2)
library(SummarizedExperiment)
library(xtable)


transformed_counts <- read_rds('autoreg/data/transformed_counts.rds')
go_annotation <- read_rds('autoreg/data/go_annotation.rds')
tf_annotation <- read_rds('autoreg/data/tf_annotation.rds')
tf <- c('Cebpb', 'Pparg', 'Rxrg', 'Ep300', 'Med1')

list('All Genes' = rownames(transformed_counts),
      'Autophagy Genes' = rownames(transformed_counts) %in% go_annotation$SYMBOL,
      'Autophagy TF' = rownames(transformed_counts) %in% intersect(go_annotation$SYMBOL, tf_annotation$SYMBOL),
      "Adipogenic TF" = rownames(transformed_counts) %in% tf) %>%
    map(function(x) {
      mds <- cmdscale(dist(t(assay(transformed_counts)[x,])),
                      eig = TRUE,
                      add = TRUE,
                      k = 2) 
      (mds$eig*100/sum(mds$eig))[1:2]
    }) %>%
  melt() %>%
  setNames(c('value', 'Category')) %>%
  mutate(value = as.integer(round(value)),
         dim = rep(paste('Dim', 1:2), 4),
         Category = factor(Category, levels = c('All Genes', 'Autophagy Genes', 'Adipogenic TF', 'Autophagy TF'))) %>%
  arrange(Category) %>%
  spread(dim, value) %>%
  xtable(align = 'clcc') %>%
  print(floating = FALSE,
        include.rownames = FALSE,
        booktabs = TRUE,
        sanitize.text.function = identity,
        comment = FALSE,
        file = 'manuscript/tables/variance_explained.tex')

#caption = 'Percent of gene expression variance explained by the stage of differentiation.'
#label = 'tab:variance_explained'

