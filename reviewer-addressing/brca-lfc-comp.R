library(limma)
library(ggplot2)
library(DEGreport)
library(ggpubr)
source("./R/routine_tasks.R")
comp <- readRDS("./metabric-brca-complex.rds")
inst <- readRDS("./metabric-brca-instability.rds")
sharedGenes <- intersect(rownames(comp), rownames(inst))
comp <- comp[sharedGenes, ]
inst <- inst[sharedGenes, ]

plotdf <- data.frame(instabilityLFC = inst$logFC, complexLFC = comp$logFC)
pdf("./lfc-comp-inst2.pdf", width = 8, height = 8)
plot(ggplot(plotdf, aes(x = instabilityLFC, y = complexLFC)) +
    geom_point(color = "black", alpha = .5) +
    # xlab("lfc of instability score >1 or <= 1") +
    xlab("lfc of instability score >2 or <= 2") +
    # scale_y_continuous(limits = c(-1.5, 1.5)) +
    # scale_x_continuous(limits = c(-1.5, 1.5)) +
    geom_cor(method = "pearson") +
    # geom_cor(method = "spearman", ypos = 1.4) +
    # geom_cor(method = "kendall", ypos = 1.3) +
    theme_classic())
dev.off()
print("done")
cor.test(comp$logFC, inst$logFC, method = "pearson")
cor.test(comp$logFC, inst$logFC, method = "spearman")
cor.test(comp$logFC, inst$logFC, method = "kendall")
