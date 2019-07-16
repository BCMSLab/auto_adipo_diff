library(tidyverse)

go_annotation <- read_rds('autoreg/data/go_annotation.rds')
occupancy_res <- read_rds('autoreg/data/occupancy_res.rds')
factor_targets <- read_rds('autoreg/data/factor_targets.rds') %>%
  filter(n > 1,
         geneId %in% go_annotation$SYMBOL) %>%
  mutate(factor = factor(factor)) %>%
  with(split(.$geneId, .$factor)) %>%
  map(unique)

map(c('PPARG', 'CEBPB'), function(x) {
  (occupancy_res$tf %>%
    as_tibble() %>%
    mutate(contrast = ifelse(contrast == 'early', 'Early vs Non', 'Late vs Non')) %>%
    filter(row %in% factor_targets[[x]]) %>%
    dplyr::select(factor, contrast, row, log2FoldChange) %>%
    spread(factor, log2FoldChange) %>%
    gather(cofac, value, MED1, RXRG, EP300) %>%
    ggplot(aes_string(x = x, y = 'value')) +
    geom_vline(xintercept = 0, lty = 2, color = 'red') +
    geom_hline(yintercept = 0, lty = 2, color = 'red') +
    geom_point(color = 'darkgray', alpha = .7) +
    facet_grid(contrast ~ cofac) +
    labs(y = 'Co-factor Occupancy Change', 
         x = paste(x, 'Occupancy Change')) +
    theme_bw() +
    theme(strip.background = element_blank(),
          panel.grid = element_blank(),
          panel.spacing = unit(0,"null"))) %>%
    ggsave(plot = .,
           filename = paste0('manuscript/figures/', x, '_occupancy_change.png'),
           width = 12, height = 8, units = 'cm')
})

map(c('PPARG', 'CEBPB'), function(x) {
  tf <- occupancy_res$tf %>%
    as_tibble() %>%
    filter(factor == x) %>%
    mutate(contrast = ifelse(contrast == 'early', 'Early vs Non', 'Late vs Non'))
  
  hm <- occupancy_res$hm %>%
    as_tibble() %>%
    mutate(contrast = ifelse(contrast == 'early', 'Early vs Non', 'Late vs Non'))
  
  (rbind(tf, hm) %>%
      filter(row %in% factor_targets[[x]]) %>%
      dplyr::select(factor, contrast, row, log2FoldChange) %>%
      spread(factor, log2FoldChange) %>%
      gather(cofac, value, H3K27ac, H3K4me1, H3K4me3) %>%
      ggplot(aes_string(x = x, y = 'value')) +
      geom_vline(xintercept = 0, lty = 2, color = 'red') +
      geom_hline(yintercept = 0, lty = 2, color = 'red') +
      geom_point(color = 'darkgray', alpha = .7) +
      facet_grid(contrast ~ cofac) +
      labs(y = 'Histone Tags Change', 
           x = paste(x, 'Occupancy Change')) +
      theme_bw() +
      theme(strip.background = element_blank(),
            panel.grid = element_blank(),
            panel.spacing = unit(0,"null"))) %>%
    ggsave(plot = .,
           filename = paste0('manuscript/figures/', x, '_hm_occupancy_change.png'),
           width = 12, height = 8, units = 'cm')
})



