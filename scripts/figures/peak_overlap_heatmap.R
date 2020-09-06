# loading required libraries
library(tidyverse)
library(ComplexHeatmap)
library(reshape2)
library(circlize)

# loading data
peak_counts <- read_rds('autoreg/data/peak_counts.rds')

# defining variables
hms <- c('H3K27ac', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K9me3')
tfs <- c('CEBPB', 'PPARG', 'RXRG', 'EP300', 'MED1')

pd <- tibble(
  factor = peak_counts$factor,
  group = peak_counts$group,
  time = peak_counts$time,
  id = peak_counts$id
) %>%
  filter(!is.na(factor)) %>%
  filter(factor %in% c(tfs, hms)) %>%
  mutate(file = paste0('autoreg/data/peaks/', id, '_peaks.xls'))

peak_overlaps <- read_rds('autoreg/data/peak_overlaps.rds')

col_fun <- colorRamp2(c(0, 1), c('white', 'darkblue'))

# generating figure
hms <- split(pd$id, pd$group) %>%
  set_names(c('Non', 'Early', 'Late')) %>%
  imap(function(x, .y) {
    mat <- peak_overlaps %>%
      filter(qSample %in% x, tSample %in% x) %>%
      mutate(frac = N_OL / (qLen + tLen)) %>%
      acast(qSample ~ tSample, value.var = 'frac')
    
    fac <- pd$factor[which(pd$id %in% x)]
    ra <- rowAnnotation(Factor1 = anno_mark(at = which(!duplicated(fac)),
                                            labels = unique(fac)))
    ca <- columnAnnotation(Factor2 = anno_mark(at = which(!duplicated(fac)),
                                               labels = unique(fac)))
    
    Heatmap(mat,
            show_column_names = FALSE,
            show_row_names = FALSE,
            show_heatmap_legend = FALSE,
            col = col_fun,
            show_column_dend = FALSE,
            show_row_dend = FALSE,
            column_names_side = 'top',
            top_annotation = ca,
            right_annotation = ra,
            column_title = .y,
            column_title_side = 'bottom')
  })

png(filename = 'manuscript/figures/peak_overlap_heatmap.png',
    width = 24, height = 9, units = 'cm', res = 300)

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 3)))
pushViewport(viewport(layout.pos.col = 1,
                      layout.pos.row = 1))

draw(hms$Non, newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.col = 2,
                      layout.pos.row = 1))
draw(hms$Early, newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.col = 3,
                      layout.pos.row = 1))
draw(hms$Late, newpage = FALSE)
upViewport()

dev.off()
