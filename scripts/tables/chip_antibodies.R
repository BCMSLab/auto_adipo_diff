library(tidyverse)
library(SummarizedExperiment)
library(xtable)

id_antibody <- read_csv('autoreg/data/id_antibody.csv')
peak_counts <- read_rds('autoreg/data/peak_counts.rds')

tfs <- c('CEBPB', 'PPARG', 'POLR2A', 'RXRG', 'EP300', 'MED1')
hms <- c('H3K27ac', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K9me3')

pd <- colData(peak_counts)

pd[, !colnames(pd) %in% 'qc'] %>%
  as_tibble() %>%
  left_join(id_antibody) %>%
  filter(!is.na(factor),
         factor %in% c(tfs, hms)) %>%
  select(factor, antibody, study, bibtexkey) %>%
  group_by(factor) %>%
  summarise(antibody = paste(unique(antibody), collapse = '; '),
            study = paste(unique(study), collapse = ', '),
            bibtexkey = paste(unique(bibtexkey), collapse = ',')) %>%
  mutate(bibtexkey = paste0('\\cite{', bibtexkey, '}')) %>%
  ungroup() %>%
  setNames(c('Factor', 'Antibody', 'Study', 'Ref.')) %>%
  xtable(align = 'clp{.3\\textwidth}p{.3\\textwidth}l') %>%
  print(floating = FALSE,
        include.rownames = FALSE,
        booktabs = TRUE,
        sanitize.text.function = identity,
        comment = FALSE,
        file = 'manuscript/tables/chip_antibodies.tex')

#caption = 'ChIP antibodies for transcription factors, co-factors and histone markers.',
#label = 'tab:chip_antibodies',
