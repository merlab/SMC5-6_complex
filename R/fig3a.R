# purpose: make plot of the heterogenity of tcga
source("./R/routine_tasks.R")
source("./R/R_rainclouds.R")
complexGenes <- c("SMC5", "SMC6", "NSMCE1", "NSMCE2", "NSMCE3", "NSMCE4A", "EID3")
colorsFinal <- c("#98969D", "#D76937", "#E89150", "#3B3B41", "#3B3C41", "#792825", "#BE4732", "#BE4732")
folder_check("./reviewer-addressing/tcga/")
library(ggpubr)
library(ggplot2)
library(dplyr)
library(rstatix)


plotDf <- readRDS("./results/tcga-vaf-data.rds")
plotDf$gene <- factor(plotDf$gene, levels = complexGenes)
p <- (
    ggplot(plotDf, aes(x = gene, y = variationFreq, fill = gene)) + # x = type
        geom_flat_violin(
            position = position_nudge(x = .15, y = 0),
            adjust = 2, trim = TRUE, width = .5
        ) +
        geom_point(
            position = position_jitter(width = .15),
            mapping = aes(x = as.numeric(gene) - 0.25),
            size = .25
        ) +
        geom_boxplot(
            outlier.shape = NA, alpha = 0.3, width = .1, colour = "black",
        ) +
        geom_hline(yintercept = 0.5, color = "red", linetype = "dashed") +
        scale_fill_manual(
            values = colorsFinal,
            name = "Gene"
        ) +
        xlab("Genes") +
        ylab("Allele Variation Frequency") +
        guides(fill = "none") +
        scale_y_continuous(expand = c(0, 0), limits = c(0, 1), breaks = seq(0, 1, .1)) +
        theme_classic() +
        theme(
            legend.position = "top",
            legend.direction = "horizontal"
        )
)
print(length(plotDf[plotDf$variationFreq >= 0.5, ]) / nrow(plotDf) * 100)

for (gene in complexGenes) {
    s <- plotDf[plotDf$gene == gene, ]
    pval <- (t.test(s$variationFreq, mu = 0.5, alternative = "less"))$p.value
    pval <- signif(pval, 2)
    if (pval < 1e-4) {
        pval.text <- "****"
    } else if (pval < 1e-3) {
        pval.text <- "***"
    } else if (pval < 1e-2) {
        pval.text <- "**"
    } else if (pval < .05) {
        pval.text <- "*"
    } else {
        pval.text <- "ns"
    }
    # https://stackoverflow.com/questions/41501561/how-to-write-numbers-in-scientific-notation-in-r
    p <- p + geom_text(x = gene, y = .95, label = pval.text)
}

p <- rmbg(p)

pdf("./figures/fig3a.pdf", height = 5, width = 8)
plot(p)
dev.off()
print("done")
