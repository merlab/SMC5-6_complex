# config
trans <- 0.3
# program
source("./R/routine_tasks.R")
library(ggplot2)
library(ggpubr)
set.seed(123)

# NOTE: making edge & node data.fames
# TODO: fix names. never use PDF
df <- as.data.frame(readRDS("./data/tcga-prad/pathway-gsea.rds"))
rownames(df) <- df$pathway
df$pathway <- gsub("_", " ", df$pathway)
df <- df[df$padj < 0.15, ]
df$logFDR <- -log10(df$padj)
df$leadingEdge <- NULL
df <- df[order(df$ES, decreasing = FALSE), ]
df$pathway <- gsub("AND", "&", df$pathway)
df$pathway <- factor(df$pathway, levels = df$pathway)


pdf("./figures/figS5c.pdf", height = 8, width = 8)
p <- ggplot(df, aes(x = ES, y = pathway, size = size, color = logFDR)) +
    geom_point() +
    scale_colour_gradient(low = "blue", high = "red") +
    xlab("Enrichment score (ES)") +
    ylab("") +
    labs(colour = expression("-log"[10] * "FDR"), size = "Pathway Size") +
    scale_colour_gradient(
        breaks = seq(1, 15, 2),
        low = "blue", high = "red"
    ) +
    theme_bw() +
    theme(
        axis.text.y = element_text(size = 10),
        legend.position = "bottom",
        legend.box = "vertical",
        panel.background = element_rect(fill = "transparent", color = NA) # bg of the panel
        , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
        , legend.background = element_rect(fill = "transparent", color = NA) # get rid of legend bg
        , legend.box.background = element_rect(fill = "transparent", color = NA) # get rid of legend panel bg
    )
plot(p)
dev.off()

