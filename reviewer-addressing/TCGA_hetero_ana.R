# purpose: download RNASeq of TCGA PRAD from TCGA biolinker
source("./R/routine_tasks.R")
source("./R/R_rainclouds.R")
# complexGenes <- c("NSMCE2", "SMC6", "SMC5", "NSMCE1", "NSMCE3", "NSMCE4A", "EID3")
complexGenes <- c("SMC5", "SMC6", "NSMCE1", "NSMCE2", "NSMCE3", "NSMCE4A", "EID3")
folder_check("./reviewer-addressing/tcga/")
library(TCGAbiolinks)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(viridis)
library(rstatix)

allTypes <- c(
    "LAML", "ACC", "BLCA", "LGG", "BRCA", "CESC", "CHOL", "LCML", "COAD",
    "CNTL", "ESCA", "FPPP", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LIHC",
    "LUAD", "LUSC", "DLBC", "MESO", "MISC", "OV", "PAAD", "PCPG", "PRAD",
    "READ", "SARC", "SKCM", "STAD", "TGCT", "THYM", "THCA", "UCS", "UCEC",
    "UVM"
)


for (type in allTypes) {
    tryCatch(
        {
            f <- sprintf("./reviewer-addressing/tcga/%s-mut.rds", type)
            if (file.exists(f)) next()
            query <- GDCquery(
                project = paste0("TCGA-", type),
                data.category = "Simple Nucleotide Variation",
                access = "open",
                legacy = FALSE,
                data.type = "Masked Somatic Mutation"
                # workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
            )
            wd <- getwd()
            setwd("~/.cache/TCGAbiolinks/")
            GDCdownload(
                query = query
            )
            dataPrep <- GDCprepare(
                query = query,
                save = TRUE
            )
            setwd(wd)
            saveRDS(dataPrep, f)
        },
        error = function(cond) {
            message(cond)
        }
    )
}

plotDf <- data.frame()
for (type in allTypes) {
    print(type)
    f <- sprintf("./reviewer-addressing/tcga/%s-mut.rds", type)
    if (!file.exists(f)) next()
    mat <- readRDS(f)
    mat <- mat[mat$Hugo_Symbol %in% complexGenes, ]
    if (nrow(mat) == 0) next()
    t <- data.frame(type = type, variationFreq = mat$t_alt_count / (mat$t_ref_count + mat$t_alt_count), gene = mat$Hugo_Symbol)
    plotDf <- rbind(plotDf, t)
}
plotDf$gene <- factor(plotDf$gene, levels = complexGenes)
p <- (
    ggplot(plotDf, aes(x = gene, y = variationFreq, fill = gene)) + # x = type
        # geom_boxplot(alpha = .75, outline.color = NA) + # , width = .7) +
        geom_flat_violin(
            # position = position_nudge(x = .25, y = 0),
            position = position_nudge(x = .35, y = 0),
            adjust = 2, trim = TRUE, width = .5
        ) +
        geom_point(
            position = position_jitter(width = .15),
            size = .25
        ) +
        geom_boxplot(aes(x = as.numeric(gene) + 0.25),
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
        # scale_color_viridis(discrete = TRUE) +
        scale_fill_viridis(discrete = TRUE, name = "Gene") +
        # xlab("Cancer Types") +
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


pdf("./reviewer-addressing/plot/hetero-sum.pdf", height = 5, width = 8)
plot(p)
dev.off()
