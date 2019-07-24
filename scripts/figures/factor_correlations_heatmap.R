# loading required libraries
library(tidyverse)
library(ComplexHeatmap)
library(SummarizedExperiment)
library(circlize)

# loadind data
peak_counts <- read_rds('autoreg/data/peak_counts.rds')
go_annotation <- read_rds('autoreg/data/go_annotation.rds')

# defining variables
tfs <- c('CEBPB', 'PPARG', 'RXRG', 'EP300', 'MED1')
hms <- c('H3K27ac', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K9me3')

# subset the object
ind <- mcols(peak_counts)$geneId %in% unique(go_annotation$SYMBOL)

se <- peak_counts[ind, ]
se <- se[, !is.na(se$factor)]
se <- se[, !is.na(se$group)]
se <- se[, se$factor %in% c(tfs, hms)]

ses <- map(unique(se$group), function(x) {
  se[, se$group == x]
}) %>%
  set_names(c('Non', 'Early', 'Late'))

corr <- map(ses, function(x) {
  cor(assay(x))
})

col_fun <- colorRamp2(c(0, 1), c('white', 'darkblue'))

# generating figure
hms <- imap(ses, function(x, .y) {
  ra <- rowAnnotation(Factor1 = anno_mark(at = which(!duplicated(x$factor)),
                                          labels = unique(x$factor)))
  ca <- columnAnnotation(Factor2 = anno_mark(at = which(!duplicated(x$factor)),
                                             labels = unique(x$factor)))
  
  Heatmap(cor(assay(x)),
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

png(filename = 'manuscript/figures/factor_correlations_heatmap.png',
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

