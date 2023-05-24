# purpose: make plot of the heterogenity of tcga
source("./R/routine_tasks.R")
source("./R/R_rainclouds.R")
complexGenes <- c("SMC5", "SMC6", "NSMCE1", "NSMCE2", "NSMCE3", "NSMCE4A", "EID3")
colorsShamayita <- c("#98969D", "#D76937", "#E89150", "#3B3B41", "#3B3C41", "#792825", "#BE4732", "#BE4732")
folder_check("./reviewer-addressing/tcga/")
library(TCGAbiolinks)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(viridis)
library(rstatix)

colors1 <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f")
colors2 <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d")
colors3 <- c("#8dd3c7", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5")
colors4 <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#e5c494", "#b3b3b3")
colors5 <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#a65628", "#f781bf")

projects <- rev(getGDCprojects()$id)
projects <- grep("TCGA-", projects, value = TRUE)

plotDf <- data.frame()
for (project in projects) {
    print(project)
    f <- sprintf("./reviewer-addressing/tcga/%s-mut.rds", project)
    if (!file.exists(f)) next()
    mat <- readRDS(f)
    mat <- mat[mat$Hugo_Symbol %in% complexGenes, ]
    if (nrow(mat) == 0) next()
    t <- data.frame(type = project, variationFreq = mat$t_alt_count / (mat$t_ref_count + mat$t_alt_count), gene = mat$Hugo_Symbol)
    plotDf <- rbind(plotDf, t)
}
plotDf$gene <- factor(plotDf$gene, levels = complexGenes)
p <- (
    ggplot(plotDf, aes(x = gene, y = variationFreq, fill = gene)) + # x = type
        # geom_boxplot(alpha = .75, outline.color = NA) + # , width = .7) +
        geom_flat_violin(
            # position = position_nudge(x = .25, y = 0),
            # position = position_nudge(x = .35, y = 0),
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
        # stat_compare_means(
        #     hide.ns = FALSE,
        #     label = "p.signif",
        #     method = "t.test",
        #     label.y = .9
        # ) +
        # https://datavizpyr.com/how-to-make-grouped-boxplot-with-jittered-data-points-in-ggplot2/
        # geom_point(position = position_jitterdodge(dodge.width = .05)) +
        # scale_fill_viridis(discrete = TRUE, name = "Gene") +
        scale_fill_manual(
            values = colorsShamayita,
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

sink("./reviewer-addressing/plot/hetero-sum-pvals.txt")
for (gene in complexGenes) {
    s <- plotDf[plotDf$gene == gene, ]
    pval <- (t.test(s$variationFreq, mu = 0.5, alternative = "less"))$p.value
    pval <- signif(pval, 2)
    print(gene)
    print(pval)
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
sink(NULL)

p <- rmbg(p)

pdf("./reviewer-addressing/plot/hetero-sum.pdf", height = 5, width = 8)
plot(p)
dev.off()
print("done")
