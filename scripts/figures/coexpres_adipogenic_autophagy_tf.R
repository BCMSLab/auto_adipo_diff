library(tidyverse)
library(reshape2)
library(ComplexHeatmap)
library(circlize)

ddcor <- read_rds('autoreg/data/dgca.rds')
go_annotation <- read_rds('autoreg/data/go_annotation.rds')
tf_annotation <- read_rds('autoreg/data/tf_annotation.rds')

targets <- intersect(go_annotation$SYMBOL, tf_annotation$SYMBOL)

hms <- ddcor %>%
  filter(Gene2 != 'Ctcf') %>%
  filter(Gene1 %in% targets) %>%
  dplyr::select(Gene1, Gene2, group, cor) %>%
  unique() %>%
  with(split(dplyr::select(., -Gene2), .$Gene2)) %>%
  imap(function(x, .y) {
    
    mat <- acast(x, Gene1 ~ group, value.var = 'cor')
    col_fun <- colorRamp2(c(-1, 0, 1), c('darkgreen', 'white', 'darkblue'))
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(as.character(round(mat[i, j], 1)), x, y, gp = gpar(fontsize = 10))}
    
      Heatmap(mat,
              cluster_columns = FALSE,
              column_order = c(3,1,2),
              column_labels = c('Early', 'Late', 'Non'),
              show_heatmap_legend = FALSE,
              cell_fun = cell_fun,
              col = col_fun,
              column_names_rot = 0,
              column_names_centered = TRUE,
              column_title = .y,
              row_title_gp = gpar(fontsize = 12),
              column_title_gp = gpar(fontsize = 12),
              column_names_gp = gpar(fontsize = 12))
  })

hms_list <- hms$Pparg + hms$Cebpb + hms$Med1 + hms$Rxrg + hms$Ep300

png(filename = 'manuscript/figures/coexpres_adipogenic_autophagy_tf.png',
    width = 24, height = 7, units = 'cm', res = 300)
draw(hms_list)
dev.off()
