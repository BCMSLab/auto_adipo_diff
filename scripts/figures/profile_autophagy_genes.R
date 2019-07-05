library(tidyverse)
library(reshape2)
library(cowplot)
library(SummarizedExperiment)
library(Gviz)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

# expression and transcripiton
gene_counts <- read_rds('autoreg/data/gene_counts.rds')
occupancy <- read_rds('autoreg/data/factor_occupancy.rds')

targets <- c('Becn1', 'Map1lc3b', 'Sqstm1')

pd <- tibble(time = gene_counts$time,
             group = gene_counts$group,
             id = gene_counts$id)

profiles <- list('Gene Expression' = assay(gene_counts)[rownames(gene_counts) %in% targets,] %>%
                   melt() %>%
                   setNames(c('geneId', 'id', 'count')) %>%
                   left_join(pd),
                 'Transcription Activity' = occupancy %>%
                   filter(geneId %in% targets) %>%
                   dplyr::select(-factor)) %>%
  bind_rows(.id = 'count_type') %>%
  mutate(group = factor(group, levels = c('non', 'early', 'late')))

(profiles %>%
    na.omit() %>%
    mutate(count = log2(count + 1)) %>%
    group_by(count_type, geneId, group) %>%
    mutate(ave = mean(count), sd = sd(count)) %>%
    ggplot(aes(x = group, y = count)) +
    geom_jitter(width = .2, color = 'darkgray', alpha = .5) +
    geom_point(aes(y = ave), color = 'red') +
    geom_linerange(aes(ymin = ave - sd, ymax = ave + sd), color = 'red') +
    facet_grid(count_type~geneId) +
    theme_bw() + 
    labs(x = 'Stage of Differentiation',
         y = 'Log2 Count') +
    theme(strip.background = element_blank(),
          panel.grid = element_blank(),
          panel.spacing = unit(0,"null"))) %>%
  ggsave(plot = .,
         filename = 'manuscript/figures/profile_autophagy_genes_a.png',
         width = 10, height = 6, units = 'cm')

# peaks
binding_data <- read_rds('autoreg/data/binding_data.rds')

df <- binding_data %>%
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
  bind_rows(.id = 'factor')

(df %>%
    filter(geneId %in% targets,
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
    theme_bw() + 
    labs(x = 'Stage of Differentiation',
         y = 'Log2 Count in Peaks') +
    theme(strip.background = element_blank(),
          panel.grid = element_blank(),
          panel.spacing = unit(0,"null"))) %>%
  ggsave(plot = .,
         filename = 'manuscript/figures/profile_autophagy_genes_b.png',
         width = 10, height = 6, units = 'cm')

# signal tracks
track_signal <- read_rds('autoreg/data/data_tracks.rds')

gene_id <- list(Becn1 = 56208,
                Map1lc3b = 67443,
                Sqstm1 = 18412,
                Cidec = 14311)

# generate figures
## PPARG
png(filename = 'manuscript/figures/profile_autophagy_genes_c.png',
    width = 12, height = 6, units = 'cm', res = 300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, length(gene_id))))

imap(gene_id, function(x, .y) {
  # get promoter region
  prom_gr <- promoters(TxDb.Mmusculus.UCSC.mm10.knownGene,
                       filter = list(gene_id = x),
                       upstream = 3000,
                       downstream = 3000) %>%
    GenomicRanges::reduce()
  pushViewport(viewport(layout.pos.col = unname(which(gene_id == x)),
                        layout.pos.row = 1))
  # extract gr info
  gr_genome <- unique(as.character(genome(prom_gr)))
  gr_chrom <- unique(as.character(seqnames(prom_gr)))
  gr_strand <- unique(as.character(strand(prom_gr)))
  gr_start <- min(start(prom_gr))
  gr_end <- max(end(prom_gr))
  
  title_add = ifelse(.y == 'Becn1', TRUE, FALSE)
  
  # PPARG
  plotTracks(list(track_signal$PPARG$`0`,
                  track_signal$PPARG$`24`,
                  track_signal$PPARG$`48`,
                  track_signal$PPARG$`72`,
                  track_signal$PPARG$`144`,
                  track_signal$PPARG$`168`),
             type = 'h',
             chromosome = gr_chrom,
             from = gr_start,
             to = gr_end,
             add = TRUE,
             main = .y,
             showTitle = title_add,
             cex.main = .7,
             cex.axis = .5,
             cex.title = .5)
  popViewport(1)
  return(NULL)
})

dev.off()

dt <- read_rds('autoreg/data/data_tracks_tissue.rds')

png(filename = 'manuscript/figures/profile_autophagy_genes_e.png',
    width = 16, height = 8, units = 'cm', res = 300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, length(gene_id))))

imap(gene_id, function(x, .y) {
  # get promoter region
  prom_gr <- promoters(TxDb.Mmusculus.UCSC.mm10.knownGene,
                       filter = list(gene_id = x),
                       upstream = 3000,
                       downstream = 3000) %>%
    GenomicRanges::reduce()
  pushViewport(viewport(layout.pos.col = unname(which(gene_id == x)),
                        layout.pos.row = 1))
  # extract gr info
  gr_genome <- unique(as.character(genome(prom_gr)))
  gr_chrom <- unique(as.character(seqnames(prom_gr)))
  gr_strand <- unique(as.character(strand(prom_gr)))
  gr_start <- min(start(prom_gr))
  gr_end <- max(end(prom_gr))
  
  title_add = ifelse(.y == 'Becn1', TRUE, FALSE)
  
  # PPARG
  plotTracks(list(
    dt$Adipocyte_PPARG_1,
    dt$Adipocyte_PPARG_2,
    dt$Fibroblast_PPARG_1,
    dt$Fibroblast_PPARG_2,
    dt$Macrophage_PPARG_1,
    dt$Macrophage_PPARG_2),
    type = 'h',
    chromosome = gr_chrom,
    from = gr_start,
    to = gr_end,
    add = TRUE,
    main = .y,
    showTitle = title_add,
    cex.main = .7,
    cex.axis = .5,
    cex.title = .4)
  popViewport(1)
  return(NULL)
})
dev.off()

# CEBPB
gene_id <- list(Becn1 = 56208,
                Map1lc3b = 67443,
                Sqstm1 = 18412,
                Klf5 = 12224)

png(filename = 'manuscript/figures/profile_autophagy_genes_d.png',
    width = 12, height = 6, units = 'cm', res = 300)

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, length(gene_id))))

imap(gene_id, function(x, .y) {
  # get promoter region
  prom_gr <- promoters(TxDb.Mmusculus.UCSC.mm10.knownGene,
                       filter = list(gene_id = x),
                       upstream = 3000,
                       downstream = 3000) %>%
    GenomicRanges::reduce()
  pushViewport(viewport(layout.pos.col = unname(which(gene_id == x)),
                        layout.pos.row = 1))
  # extract gr info
  gr_genome <- unique(as.character(genome(prom_gr)))
  gr_chrom <- unique(as.character(seqnames(prom_gr)))
  gr_strand <- unique(as.character(strand(prom_gr)))
  gr_start <- min(start(prom_gr))
  gr_end <- max(end(prom_gr))
  
  title_add = ifelse(.y == 'Becn1', TRUE, FALSE)
  
  # CEBPB
  plotTracks(list(track_signal$CEBPB$`0`,
                  track_signal$CEBPB$`2`,
                  track_signal$CEBPB$`4`,
                  track_signal$CEBPB$`6`,
                  track_signal$CEBPB$`48`),
             type = 'h',
             chromosome = gr_chrom,
             from = gr_start,
             to = gr_end,
             add = TRUE,
             main = .y,
             showTitle = title_add,
             cex.main = .7,
             cex.axis = .5,
             cex.title = .5)
  popViewport(1)
  return(NULL)
})
dev.off()

png(filename = 'manuscript/figures/profile_autophagy_genes_f.png',
    width = 14, height = 7, units = 'cm', res = 300)

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, length(gene_id))))

imap(gene_id, function(x, .y) {
  # get promoter region
  prom_gr <- promoters(TxDb.Mmusculus.UCSC.mm10.knownGene,
                       filter = list(gene_id = x),
                       upstream = 3000,
                       downstream = 3000) %>%
    GenomicRanges::reduce()
  pushViewport(viewport(layout.pos.col = unname(which(gene_id == x)),
                        layout.pos.row = 1))
  # extract gr info
  gr_genome <- unique(as.character(genome(prom_gr)))
  gr_chrom <- unique(as.character(seqnames(prom_gr)))
  gr_strand <- unique(as.character(strand(prom_gr)))
  gr_start <- min(start(prom_gr))
  gr_end <- max(end(prom_gr))
  
  title_add = ifelse(.y == 'Becn1', TRUE, FALSE)
  
  # CEBPB
  plotTracks(list(
    dt$Adipocyte_CEBPB_1,
    dt$Adipocyte_CEBPB_2,
    dt$Fibroblast_CEPBP_1,
    dt$Fibroblast_CEBPB_2,
    dt$Macrophage_CEBPB_1),
    type = 'h',
    chromosome = gr_chrom,
    from = gr_start,
    to = gr_end,
    add = TRUE,
    main = .y,
    showTitle = title_add,
    cex.main = .7,
    cex.axis = .5,
    cex.title = .5)
  popViewport(1)
  return(NULL)
})

dev.off()
