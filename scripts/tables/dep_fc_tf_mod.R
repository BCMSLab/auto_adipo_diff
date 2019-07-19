library(tidyverse)
library(xtable)
library(SummarizedExperiment)

binding_data <- read_rds('autoreg/data/binding_data.rds')
dep_res <- read_rds('autoreg/data/dep_res.rds')

ind <- c('Ctcf', 'Cebpb', 'Pparg', 'Rxrg', 'Ep300', 'Med1')


peak_symbol <- map(binding_data, function(x){
  tibble(row = mcols(x)$name,
         gene_id = mcols(x)$geneId)
}) %>%
  bind_rows() %>%
  unique() %>%
  filter(gene_id %in% ind)

header <- paste0("\\multirow{2}[3]{*}{Factor} & \\multirow{2}[3]{*}{Gene} &",
                 " \\multicolumn{4}{c}{Early vs Non} & \\multicolumn{4}{c}{Late vs Non} &  \\multicolumn{4}{c}{Late vs Early} \\\\",
                 " \\cmidrule(lr){3-6} \\cmidrule(lr){7-10} \\cmidrule(lr){11-14}",
                 "&& (N) & Range & Ave & SD & (N) & Range & Ave & SD & (N) & Range & Ave & SD \\\\")

cat_fac <- list(factor = c('CTCF', 'CEBPB', 'PPARG'),
                co_factor = c('EP300', 'MED1', 'RXRG'),
                hm = c('H3K27ac', 'H3K4me3'))
fac <- tibble(cat = factor(rep(c('Factor', 'Cofactor', 'Histone Marker'), times = c(3,3,2)),
                           levels = c('Factor', 'Cofactor', 'Histone Marker')),
              factor = factor(unlist(cat_fac, use.names = FALSE),
                              levels = unlist(cat_fac, use.names = FALSE)))

peak_symbol %>%
  inner_join(dep_res) %>%
  filter(factor %in% c('PPARG', 'CEBPB')) %>%
  filter(gene_id %in% c('Pparg', 'Cebpb')) %>%
  group_by(factor, contrast, gene_id) %>%
  summarise(n = n(),
            range = ifelse(n() == 1, as.character(round(log2FoldChange, 2)),
                           paste(round(min(log2FoldChange), 2), round(max(log2FoldChange), 2), sep = '/')),
            ave = round(mean(log2FoldChange), 2),
            sd = ifelse(is.na(sd(log2FoldChange)), '', as.character(round(sd(log2FoldChange), 2)))) %>%
  unite(values, n, range, ave, sd, sep = '_') %>%
  spread(contrast, values) %>%
  ungroup() %>%
  inner_join(fac) %>%
  arrange(cat) %>%
  dplyr::select(factor, gene = gene_id, early_vs_non, late_vs_non, late_vs_early) %>%
  separate(early_vs_non, sep = '_', into = c('n1', 'range1', 'ave1', 'sd1')) %>%
  separate(late_vs_non, sep = '_', into = c('n2', 'range2', 'ave2', 'sd2')) %>%
  separate(late_vs_early, sep = '_', into = c('n3', 'range3', 'ave3', 'sd3')) %>%
  mutate(factor = ifelse(duplicated(factor), '', factor)) %>%
  xtable(align = 'cllcccccccccccc') %>%
  print(floating = FALSE,
        include.rownames = FALSE,
        booktabs = TRUE,
        sanitize.text.function = identity,
        comment = FALSE,
        include.colnames=FALSE,
        add.to.row = list(pos = list(0, 2),
                          command = c(header, '\\midrule ')),
        file = 'manuscript/tables/dep_fc_tf_mod.tex')

#caption = 'Significant peaks of adipogenic factors on adipogenic transcription factor genes.',
#label = 'tab:dep_fc_tf',
