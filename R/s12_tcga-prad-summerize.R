source('./R/routine_tasks.R')
folder_check('./results/tcga-prad')
library(readxl)
library(writexl)

dgea <- readRDS('./data/tcga-prad/dgea-edgeR.rds')
dgea$gene <- rownames(dgea)
rownames(dgea) <- 1:nrow(dgea)
df <- dgea[order(abs(dgea$logFC), decreasing = TRUE), ]
rownames(df) <- 1:nrow(df)
df <- df[, c('gene', 'logFC', 'logCPM', 'LR', 'PValue', 'FDR')]
write.csv(df, './results/tcga-prad/dgea-edgeR.csv', quote = FALSE)

pa <- readRDS('./data/tcga-prad/pathway-gsea.rds')
pa$leadingEdge <- NULL
write.csv(pa, './results/tcga-prad/pathway-gsea.csv', quote = FALSE)

dgea <- dgea[dgea$FDR <= 0.05, ]
pa <- pa[pa$padj <= 0.05,]
pa$leadingEdge <- NULL

dgea <- dgea[,c('gene',colnames(dgea)[1:ncol(dgea)-1])]

l <- list(DGEA = dgea,
          GSEA = pa)
write_xlsx(l,'./results/Supplementary-Data-3_PRAD_DGEA_GSEA.xlsx')
