# load required libraries
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(reshape2)

eset1 <- read_rds('autoreg/data/primary_adipocyte.rds')

eset2 <- read_rds('autoreg/data/primary_adipocyte_insulin.rds')
exprs(eset2) <- log2(exprs(eset2) + 1)

targets <- list('Lipogenic' = c('LPL', 'FASN', 'ACLY'),
                'Adipogenic' = c('PPARG', 'CEBPB'),
                'Autophagy' = c('FOXO1', 'TFEB', 'XBP1', 'BECN1', 'MAP1LC3B', 'SQSTM1'))
length(unlist(targets))
ind1 <- featureNames(eset1) %in% unlist(targets)
ind2 <- featureNames(eset2) %in% unlist(targets)

mat <- cbind(exprs(eset1)[ind1,], exprs(eset2)[ind2,])[unlist(targets),]
row_slices <- factor(rep(names(targets), times = lengths(targets)),
                     levels = names(targets))
col_slices <- factor(toupper(c(eset1$group, eset2$group)),
                     levels = c('NDF', 'DF', 'CONTROL', 'OIS', 'OIR'))

peaks <- read_rds('autoreg/data/primary_adipocyte_chip.rds')
peaks <- map(peaks, function(x) {
  gr <- x[x$SYMBOL %in% unlist(targets)]
  df <- data.frame(
    gene = gr$SYMBOL,
    annotation = str_split(gr$annotation, ' ', simplify = TRUE)[, 1]
  ) %>%
    filter(annotation %in% c("5'",  "Distal", "Exon", "Intron", "Promoter"))
  
  mat <- left_join(tibble(gene = rownames(mat)), df) %>%
    acast(gene ~ annotation)
  
  mat <- mat[, 1:5]
  mat[mat != 0] <- '1'
  mat[unlist(targets),]
}) %>%
  setNames(c('CEBPB', 'PPARG'))

marks <- columnAnnotation(cell = anno_mark(
  labels = c('Pre-adipocyte\nNo MDI',
             'Pre-adipocyte\n+ MDI',
             'Subcutanous\nLean',
             'Subcutanous\nObese\nInsuline-sensitive',
             'Subcutanous\nObese\nInsuline-resistant'),
  at = c(1, 13, 25, 35, 55))
  )
hm <- Heatmap(scale(mat, scale = FALSE),
              column_split = col_slices,
              cluster_columns = FALSE,
              cluster_column_slices = FALSE,
              cluster_row_slices = FALSE,
              show_column_names = FALSE,
              show_column_dend = FALSE,
              show_row_dend = FALSE,
              column_title = rep('', 5),
              top_annotation = marks)

marks_hp1 <- columnAnnotation(cell = anno_mark(
  labels = c('PPARG\n hMADS\n Day 19'), at = 1
))

hp1 <- Heatmap(peaks$PPARG,
               col = c('white', 'black'),
               width = unit(3, 'cm'),
               column_names_rot = 45,
               top_annotation = marks_hp1)

marks_hp2 <- columnAnnotation(cell = anno_mark(
  labels = c('CEBPB\n hMSC\n 6 Hours'), at = 1
))

hp2 <- Heatmap(peaks$CEBPB,
               col = c('white', 'black'),
               width = unit(3, 'cm'),
               column_names_rot = 45,
               top_annotation = marks_hp2)

png(filename = 'manuscript/figures/primary_adipocyte.png',
    width = 21, height = 17, units = 'cm', res = 300)
draw(hp1 + hp2 + hm, 
     row_split = row_slices,
     cluster_rows = FALSE,
     show_heatmap_legend = FALSE)
dev.off()

# \clearpage  % To be removed
# 
# \begin{figure}
# \centering
# \includegraphics[width=\textwidth]{figures/primary_adipocyte.png}
# \caption[Gene expression and binding patterns of adipogenic transcription factors on key genes in primary human adipocytes]{Gene expression and binding patterns of adipogenic transcription factors on key genes in primary human adipocytes.
#   Probe intensity of selected genes from microarray samples (n = 48) from lean and obese human subjects were shown as color values (low, blue; high, red). Probe intensity of selected genes from microarray samples (n = 24) from MDI-induced and non-induced pre-adipocytes from human subjects were shown as color values (low, blue; high, red). Binding peaks-when present - of PPARG and CEBPB on selected genes in MDI-induced hMADS (19 days) and hMSC (6 hours) are shown as binary colors (black, peak; white, no peak) and stratified by genomic location.
# }
# \label{fig:primary_adipocyte}
# \end{figure}
# \clearpage  % To be removed

