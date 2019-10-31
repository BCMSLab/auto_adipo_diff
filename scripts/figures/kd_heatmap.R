# loading required libraries
library(tidyverse)
library(circlize)
library(ComplexHeatmap)
library(DESeq2)

# loading data
se <- read_rds('autoreg/data/counts_cebpb_kd.rds')
dds <- DESeqDataSet(se, design = ~group*time)
transformed_se <- varianceStabilizingTransformation(dds)

# defining variables
targets <- list('Adipogenic TF' = c('Cebpa', 'Cebpb', 'Pparg'),
                'Lipogenesis' = c('Lpl', 'Acly', 'Fasn'),
                'Autophagy TF' = c('Foxo1', 'Tfeb', 'Xbp1'),
                'Autophagy Gene' = c('Map1lc3b', 'Becn1', 'Sqstm1')) %>%
  reshape2::melt() %>%
  setNames(c('gene', 'category'))

mat <- assay(transformed_se)[rownames(transformed_se) %in% targets$gene,]
mat2 <- mat
mat2 <- t(apply(mat, 1, scale))

colnames(mat2) <- se$id

column_ord <- order(paste(se$time, se$group))
row_ord <- match(rownames(mat2), targets$gene)

col_fun <- colorRamp2(c(-2, 0, 2), c('darkgreen', 'white', 'darkblue'))

png(filename = 'manuscript/figures/kd_heatmap_cebpb.png',
    width = 10, height = 12, units = 'cm', res = 500)
Heatmap(mat2,
        col = col_fun,
        cluster_columns = FALSE,
        cluster_column_slices = FALSE,
        column_split = se$time,
        column_labels = se$group,
        column_order = column_ord,
        row_order = row_ord,
        cluster_rows = FALSE,
        row_split = targets$category[row_ord],
        row_title_gp = gpar(fontsize = 8),
        column_names_rot = 45,
        show_heatmap_legend = FALSE)

dev.off()

# pparg
eset <- read_rds('autoreg/data/arrays_pparg_kd.rds')
targets <- filter(targets, gene %in% rownames(eset))
eset <- eset[rownames(eset) %in% targets$gene, eset$time %in% c(0, 2, 5)] 
mat <- exprs(eset)
#mat2 <- log2(mat + 1)

mat2 <- t(apply(mat, 1, scale))
colnames(mat2) <- eset$geo_accession

column_ord <- order(paste(eset$time, eset$group))
row_ord <- match(rownames(mat2), targets$gene)

col_fun <- colorRamp2(c(-2, 0, 2), c('darkgreen', 'white', 'darkblue'))

png(filename = 'manuscript/figures/kd_heatmap_pparg.png',
    width = 18, height = 14, units = 'cm', res = 500)

Heatmap(mat2,
        col = col_fun,
        cluster_columns = FALSE,
        cluster_column_slices = FALSE,
        column_split = eset$time,
        column_labels = eset$group,
        column_order = column_ord,
        row_order = row_ord,
        cluster_rows = FALSE,
        row_split = targets$category[row_ord],
        column_names_rot = 45,
        row_title_gp = gpar(fontsize = 8),
        show_heatmap_legend = FALSE)

dev.off()
