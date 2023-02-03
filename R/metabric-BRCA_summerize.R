library(readxl)
library(writexl)

dgea <- readRDS('./data/metabric-brca/DGEA_limma.rds')
pa <- readRDS('./data/metabric-brca/pathway_GSEA.rds')

dgea <- dgea[dgea$adj.P.Val <= 0.05, ]
pa <- pa[pa$padj <= 0.05,]
pa$leadingEdge <- NULL
dgea$Significance <- NULL
dgea <- dgea[,c('gene',colnames(dgea)[1:ncol(dgea)-1])]

l <- list(DGEA = dgea,
          GSEA = pa)
write_xlsx(l,'./results/Supplementary-Data-1_BRCA_diff_Pathway.xlsx')
