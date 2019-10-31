# loading required librars
library(tidyverse)

# loadind data
data_tracks_tidy <- read_rds('autoreg/data/data_tracks_tidy.rds')
data_tracks_tissue_tidy <- read_rds('autoreg/data/data_tracks_tissue_tidy.rds')

# defining variables
targets <- list(adipogenic_tf = c('Pparg', 'Cebpb'),
                autophagy_tf = c('Foxo1', 'Tfeb', 'Xbp1'),
                autophagy_genes = c('Becn1', 'Map1lc3b', 'Sqstm1'))
factors <- list(PPARG = 'PPARG', CEBPB = 'CEBPB')
controls <- list('Cidec', 'Klf5')

# generating figures
map2(factors, controls, function(f, c) {
  
  imap(targets, function(x, .y) {
    filename <- paste0('manuscript/figures/signal_', f, '_', .y, '.png')
    df <- data_tracks_tidy %>%
      filter(gene_id == c, factor == f) %>%
      group_by(time) %>%
      summarise(max = max(score))
    
    (data_tracks_tidy %>%
      filter(factor == f,
             gene_id %in% c(x, c)) %>%
        left_join(df) %>%
        mutate(width = end - start,
               width = ifelse(width < 10, 10, width),
               score = score/max,
               gene_id = factor(gene_id, levels = c(x, c))) %>%
        ggplot(aes(x = pos, y = score, width = width)) +
      geom_col() +
      geom_vline(xintercept = 0, color = 'red', lty = 2) +
      facet_grid(time ~ gene_id, scales = 'free') +
      labs(x = '', y = 'Score (Relative Peak Enrichment)') +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +
      theme_classic() +
      theme(panel.grid = element_blank(),
            strip.background = element_blank(),
            axis.text = element_text(size = 8)) +
      scale_x_continuous(breaks = c(-3000, 0, 3000),
                         labels = c('-3kb', 'TSS', '3kb'))) %>%
      ggsave(plot = .,
             filename = filename,
             width = 18, height = 9, units = 'cm')
  })
})

map2(factors, controls, function(f, c) {
  imap(targets, function(x, .y) {
    filename <- paste0('manuscript/figures/signal_tissue_', f, '_', .y, '.png')
    df <- data_tracks_tissue_tidy %>%
      filter(gene_id == c, factor == f) %>%
      group_by(tissue, rep) %>%
      summarise(max = max(score))
    
    (data_tracks_tissue_tidy %>%
        filter(factor == f,
               gene_id %in% c(x, c)) %>%
        left_join(df) %>%
        mutate(width = end - start,
               width = ifelse(width < 10, 10, width),
               score = score/max,
               gene_id = factor(gene_id, levels = c(x, c))) %>%
        ggplot(aes(x = pos, y = score, width = width)) +
        geom_col() +
      geom_vline(xintercept = 0, color = 'red', lty = 2) +
      facet_grid(tissue + rep ~ gene_id, scales = 'free') +
      labs(x = '', y = 'Score (Relative Peak Enrichment)') +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +
      theme_classic() +
      theme(panel.grid = element_blank(),
            strip.background = element_blank(),
            axis.text = element_text(size = 8)) +
      scale_x_continuous(breaks = c(-3000, 0, 3000),
                         labels = c('-3kb', 'TSS', '3kb'))) %>%
      ggsave(plot = .,
             filename = filename,
             width = 18, height = 9, units = 'cm')
  })
})
