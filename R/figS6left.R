# volcano and box
library(ggplot2)
library(ggrepel)
library(fgsea)
library(writexl)
set.seed(123)
source('./R/routine_tasks.R')

volcanoPlot <- function(df, title, FCt = 2, FDRt = 0.001) {
    df$logFDR <- as.numeric(df$logFDR)
    df$logFC <- as.numeric(df$logFC)
    df$label <- rownames(df)
    df$Significance <- 'Not Significant'
    df$Significance[df$logFC > FCt & df$logFDR > -log10(FDRt)] <- 'Over Expressed'
    df$Significance[df$logFC < -FCt & df$logFDR > -log10(FDRt)] <- 'Under Expressed'
    df$logFDR <- as.numeric(df$logFDR)
    df$logFC <- as.numeric(df$logFC)
    #df2 <- df[df$label %in% top_genes,]
    limma <- readRDS('./results/metabric-brca/limma.rds')
    sharedSigGenes <- intersect(rownames(df)[df$logFDR > -log10(FDRt)], rownames(limma)[limma$adj.P.Val < FDRt])
    df2 <- df[df$label %in% sharedSigGenes,]
    df2 <- df2[order(abs(df2$logFC), decreasing = TRUE),]
    df2 <- df2[1:5,]
    color <- 'red'
    return(
        ggplot(data = df,
               aes(x = logFC, y = logFDR, color = Significance)) +
            geom_point(size = 1) +
            scale_color_manual(values = c('grey60', '#ef8a62','#67a9cf')) + 
            geom_hline(yintercept=-log10(FDRt), linetype="dashed") +
            geom_vline(xintercept=-FCt, linetype="dashed") +
            geom_vline(xintercept=FCt, linetype="dashed") +
            ylab(expression('-log'[10]*'FDR')) +
            xlab(expression('log'[2]*'FC')) +
            # add notable genes
            scale_y_continuous(breaks = c(seq(0,40, 5))) +
            scale_x_continuous(breaks = c(seq(-6,8, 1))) +
            theme_classic() + 
            geom_label_repel(
              data = df2[df2$logFC < 0, ], mapping = aes(label = label) 
              , size = 4, min.segment.length = 0
              , label.size = NA, fill = NA, color = color
              , xlim = c(-4,-4.1), nudge_y = 10, nudge_x = 0
            ) + 
            geom_label_repel(
              data = df2[df2$logFC > 0, ], mapping = aes(label = label) 
              , size = 4, nudge_x = -0.1, nudge_y = 30, min.segment.length = 0 #1.75
              ,  label.size = NA, fill = NA, color = color
              , xlim = c(4,4.1)
            ) + 
            guides(color = 'none') +
            theme(legend.title=element_text(size=12),
                  legend.text=element_text(size=12),
                  legend.position = c(0.25, 0.8),
                  axis.text.x = element_text(angle = 45, hjust = 1)
                  )
      )
}

p <- readRDS('./results/TCGA_PRAD/edgeR.rds')

#svg('./figures/figS6left.svg', width = 6, height = 8, onefile = TRUE)
pdf('./figures/figS6left.pdf', width = 6, height = 8, onefile = TRUE)
plot(rmbg(volcanoPlot(p)))
dev.off()


p_rank <- na.omit(setNames(p$logFC, rownames(p)))
calc_gsea(p_rank[order(p_rank, decreasing = TRUE)],
  , rdsout = './results/TCGA_PRAD/pathway.rds'
  , xlsxout = './results/TCGA_PRAD/pathway_enrichment_brca_metabric.xlsx')
print('done')
