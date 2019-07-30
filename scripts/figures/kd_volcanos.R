# loading required libraries
library(tidyverse)
library(reshape2)

# loading data
cebpb_kd_res <- read_rds('autoreg/data/cebpb_kd_res.rds')
pparg_kd_res <- read_rds('autoreg/data/pparg_kd_res.rds')

# defining variables
targets <- list('Adipogenic TF' = c('Cebpa', 'Cebpb', 'Pparg'),
                'Lipogenesis' = c('Lpl', 'Acly', 'Fasn'),
                'Autophagy TF' = c('Foxo1', 'Tfeb', 'Xbp1'),
                'Autophagy Gene' = c('Map1lc3b', 'Becn1', 'Sqstm1')) %>%
  melt() %>%
  setNames(c('row', 'category'))

(cebpb_kd_res %>%
  ggplot(aes(x = log2FoldChange, y = -log10(pvalue))) +
  facet_wrap(~time) +
  geom_point(color = 'darkgray', alpha = .5) +
  geom_vline(xintercept = c(-1,1), lty = 2) +
  geom_hline(yintercept = 5, lty = 2) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(0,"null"))+
  labs(y = "P-value (-Log 10)",
       x = 'Fold-Change (Log 2)')) %>%
  ggsave(plot = .,
         filename = 'manuscript/figures/volcanos_kd_cebpb.png',
         width = 10, height = 7, units = 'cm')

(pparg_kd_res %>%
    filter(time  %in% c(0, 2, 5)) %>%
    ggplot(aes(x = logFC, y = -log10(P.Value))) +
    facet_wrap(~time) +
    geom_point(color = 'darkgray', alpha = .5) +
    geom_vline(xintercept = c(-1,1), lty = 2) +
    geom_hline(yintercept = 5, lty = 2) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          panel.spacing = unit(0,"null"))+
    labs(y = "P-value (-Log 10)",
         x = 'Fold-Change (Log 2)')) %>%
  ggsave(plot = .,
         filename = 'manuscript/figures/volcanos_kd_pparg.png',
         width = 12, height = 7, units = 'cm')
