library(readxl)
library(writexl)

dgea <- readRDS('./data/metabric-brca/dgea_limma.rds')
df <- dgea[order(abs(dgea$logFC), decreasing = TRUE), ]
df$FDR <- df$adj.P.Val
df$adj.P.Val <- NULL
df$Significance <- NULL
df$t <- NULL
df$B <- NULL
df$name <- NULL
df <- df[, c('gene', 'logFC', 'AveExpr', 'P.Value', 'FDR')]
write.csv(df, './results/metabric-brca/dgea_limma.csv', quote = FALSE)

pa <- readRDS('./data/metabric-brca/pathway_gsea.rds')
pa$leadingEdge <- NULL
write.csv(pa, './results/metabric-brca/pathway_gsea.csv', quote = FALSE)

dgea <- dgea[dgea$adj.P.Val <= 0.05, ]
pa <- pa[pa$padj <= 0.05,]
pa$leadingEdge <- NULL
dgea$Significance <- NULL
dgea <- dgea[,c('gene',colnames(dgea)[1:ncol(dgea)-1])]

l <- list(DGEA = dgea,
          GSEA = pa)
write_xlsx(l,'./results/Supplementary-Data-2_BRCA_diff_Pathway.xlsx')
