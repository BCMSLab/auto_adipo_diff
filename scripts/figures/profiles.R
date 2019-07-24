library(tidyverse)
library(reshape2)
library(SummarizedExperiment)

# expression and transcripiton
gene_counts <- read_rds('autoreg/data/gene_counts.rds')
occupancy <- read_rds('autoreg/data/factor_occupancy.rds')
binding_data <- read_rds('autoreg/data/binding_data.rds')

# defining variables
pd <- tibble(time = gene_counts$time,
             group = gene_counts$group,
             id = gene_counts$id)

targets <- list(adipogenic_tf = c('Pparg', 'Cebpb'),
                autophagy_tf = c('Foxo1', 'Tfeb', 'Xbp1'),
                autophagy_genes = c('Becn1', 'Map1lc3b', 'Sqstm1'))

# generating figures
imap(targets, function(x, .y) {
  filename <- paste0('manuscript/figures/profile_expression_', .y, '.png')
  
  (list('Gene Expression' = assay(gene_counts)[rownames(gene_counts) %in% x,] %>%
         melt() %>%
         setNames(c('geneId', 'id', 'count')) %>%
         left_join(pd),
       'Transcription Activity' = occupancy %>%
         filter(geneId %in% x) %>%
         dplyr::select(-factor)) %>%
    bind_rows(.id = 'count_type') %>%
    mutate(group = factor(group, levels = c('non', 'early', 'late'))) %>%
    na.omit() %>%
    mutate(count = log2(count + 1)) %>%
    group_by(count_type, geneId, group) %>%
    mutate(ave = mean(count), sd = sd(count)) %>%
    ggplot(aes(x = group, y = count)) +
    geom_jitter(width = .2, color = 'darkgray', alpha = .5) +
    geom_point(aes(y = ave), color = 'red') +
    geom_linerange(aes(ymin = ave - sd, ymax = ave + sd), color = 'red') +
    facet_grid(count_type~geneId) +
    expand_limits(y = 0) +
    scale_x_discrete(labels = c('Non', 'Early', 'Late')) +
    theme_bw() + 
    labs(x = 'Stage of Differentiation',
         y = 'Reads Count (Log 2)') +
    theme(strip.background = element_blank(),
          panel.grid = element_blank(),
          panel.spacing = unit(0,"null"))) %>%
    ggsave(plot = .,
           filename = filename,
           width = 12, height = 8, units = 'cm')
})

imap(targets, function(x, .y) {
  filename <- paste0('manuscript/figures/profile_peaks_', .y, '.png')
  
  (binding_data %>%
    map(function(x) {
      pd <- tibble(Var2 = x$id,
                   group = x$group,
                   time = x$time)
      fd <- tibble(Var1 = rownames(x),
                   geneId = mcols(x)$geneId)
      assay(x) %>%
        melt() %>%
        left_join(pd) %>%
        left_join(fd)
    }) %>%
    bind_rows(.id = 'factor') %>%
    filter(geneId %in% x,
           factor %in% c('PPARG', 'CEBPB'),
           !is.na(group)) %>%
    mutate(group = factor(group, levels = c('non', 'early', 'late'))) %>%
    group_by(group, factor, geneId, Var1) %>%
    summarise(value = log2(mean(value) + 1)) %>%
    group_by(group, factor, geneId) %>%
    mutate(ave = mean(value), sd = sd(value)) %>%
    ggplot(aes(x = group, y = value)) +
    geom_jitter(width = .2, color = 'darkgray', alpha = .5) +
    geom_point(aes(y = ave), color = 'red') +
    geom_linerange(aes(ymin = ave - sd, ymax = ave + sd), color = 'red') +
    facet_grid(factor~geneId) +
    expand_limits(y = 0) +
    scale_x_discrete(labels = c('Non', 'Early', 'Late')) + 
    theme_bw() + 
    labs(x = 'Stage of Differentiation',
         y = 'Reads in Peaks Count (Log 2)') +
    theme(strip.background = element_blank(),
          panel.grid = element_blank(),
          panel.spacing = unit(0,"null"))) %>%
    ggsave(plot = .,
           filename = filename,
           width = 12, height = 8, units = 'cm')
})
