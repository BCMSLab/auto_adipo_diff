library(tidyverse)
library(cowplot)

deg_res <- read_rds('autoreg/data/deg_res.rds')
go_annotation <- read_rds('autoreg/data/go_annotation.rds')
tf_annotation <- read_rds('autoreg/data/tf_annotation.rds')
targets <- intersect(go_annotation$SYMBOL, tf_annotation$SYMBOL)

(deg_res %>%
  mutate(contrast = case_when(
    contrast == 'early_vs_non' ~ 'Early vs Non',
    contrast == 'late_vs_non' ~ 'Late vs Non',
    contrast == 'late_vs_early' ~ 'Late vs Early')) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(color = 'darkgray', alpha = .5) +
  geom_vline(xintercept = c(-1,1), lty = 2) +
  geom_hline(yintercept = 5, lty = 2) +
  facet_wrap(~contrast, nrow = 1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(0,"null")) +
  labs(y = "P-value (-Log 10)",
       x = 'Fold-Change (Log 2)')) %>%
  ggsave(plot = .,
         filename = 'manuscript/figures/volcanos_all.png',
         width = 12, height = 7, units = 'cm')

(deg_res %>%
  mutate(contrast = case_when(
    contrast == 'early_vs_non' ~ 'Early vs Non',
    contrast == 'late_vs_non' ~ 'Late vs Non',
    contrast == 'late_vs_early' ~ 'Late vs Early')) %>%
  filter(row %in% go_annotation$SYMBOL) %>%
  mutate(name = ifelse(row %in% targets, row, '')) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(color = 'darkgray', alpha = .5) +
  geom_vline(xintercept = c(-1,1), lty = 2) +
  geom_hline(yintercept = 5, lty = 2) +
  geom_text(aes(label = name), color = 'magenta', size = 3) +
  facet_wrap(~contrast, nrow = 1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(0,"null"))+
  labs(y = "P-value (-Log 10)",
       x = 'Fold-Change (Log 2)')) %>%
  ggsave(plot = .,
         filename = 'manuscript/figures/volcanos_autophagy.png',
         width = 12, height = 7, units = 'cm')
