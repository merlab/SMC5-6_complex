library(ggplot2)
library(ggpubr)
library(Biobase)
library(tidyr)
library(cowplot)
library(scales)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(ggrepel)
library(fgsea)
library(writexl)
set.seed(123)
genes <- c("SMC5", "SMC6", "NSMCE1", "NSMCE2", "NSMCE3", "NSMCE4A", "EID3")
source('./R/R_rainclouds.R')
source('./R/routine_tasks.R')

volcanoPlot <- function(df, title, FCt = 0.15, FDRt = 0.001) {
    df$logFDR <- as.numeric(df$logFDR)
    df$logFC <- as.numeric(df$logFC)
    df$label <- rownames(df)
    df$Significance <- 'Not Significant'
    df$Significance[df$logFC > FCt & df$logFDR > -log10(FDRt)] <- 'Over Expressed'
    df$Significance[df$logFC < -FCt & df$logFDR > -log10(FDRt)] <- 'Under Expressed'
    df$logFDR <- as.numeric(df$logFDR)
    df$logFC <- as.numeric(df$logFC)
    #print(top_genes)

    df2 <- df[df$label %in% top_genes,]
    color <- 'red'

    return(
        ggplot(data = df,
               aes(x = logFC, y = logFDR, color = Significance)) +
            geom_point(size = 1) +
            scale_color_manual(values = c('grey60', '#ef8a62','#67a9cf')) + 
            geom_hline(yintercept=-log10(0.001), linetype="dashed") +
            geom_vline(xintercept=-0.15, linetype="dashed") +
            geom_vline(xintercept=0.15, linetype="dashed") +
            ylab(expression('-log'[10]*'FDR')) +
            xlab(expression('log'[2]*'FC')) +
            # add notable genes
            geom_point(data = df2, size = 0.75, color = color, shape = 15) +
            geom_label_repel(
              data = df2[df2$logFC < 0, ], mapping = aes(label = label) 
              , size = 4, min.segment.length = 0
              , label.size = NA, fill = NA, color = color
              , xlim = c(-1,-1.1), nudge_y = 30, nudge_x = 0
            ) + 
            geom_label_repel(
              data = df2[df2$logFC > 0, ], mapping = aes(label = label) 
              , size = 4, nudge_x = -0.1, nudge_y = 30, min.segment.length = 0 #1.75
              ,  label.size = NA, fill = NA, , color = color
              #, xlim = c(-0.2,-0.7)
              , xlim = c(0.9,1)
            ) + 
            xlim(min(df$logFC), 1.25) +
            scale_x_continuous(breaks = round(c(min(df$logFC),seq(-1.5, 1.5, 0.25)), digits = 2) 
                        , limits = c(min(df$logFC), 1.25)) +
            ylim(0,135) +
            scale_y_continuous(breaks = c(seq(0,125, 25), 135)) +
            theme_classic() + 
            guides(color = guide_legend(override.aes = list(size = 2))) +
            theme(legend.title=element_text(size=12),
                  legend.text=element_text(size=12),
                  legend.position = c(0.25, 0.8),
                  axis.text.x = element_text(angle = 45, hjust = 1)
                  )
      )
}

b <- readRDS('./data/metabric-brca/dgea-limma.rds')

breastVolcano <- rmbg(volcanoPlot(b))

pdf('./figures/fig4a.pdf', width = 6, height = 6, onefile = TRUE)
plot(breastVolcano)
dev.off()

b_mRNA <- readRDS('./data/metabric-brca/microarray-metagx.rds')
b_mutation <- obtain_mut_from_mRNA(b_mRNA)
b_mRNA <- b_mRNA[,rownames(b_mutation)]
breastRainCloud <- rainCloudPlot(b_mRNA[top_genes, ], b_mutation, type = 'microarray')
p <- rmbg(ggarrange(plotlist = breastRainCloud[-5], nrow = 4, ncol = 2, align = 'hv'))
pdf('./figures/fig4b.pdf', width = 6, height = 12, onefile = TRUE)
plot(p)
dev.off()
print('done')
