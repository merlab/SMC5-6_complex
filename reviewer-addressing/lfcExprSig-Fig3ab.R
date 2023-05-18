# load libs
source("./R/R_rainclouds.R")
library(ggplot2)
library(ggpubr)
library(cowplot)
dir <- getwd()
source("./R/routine_tasks.R")
tissue <- "Breast"
colors <- c("#DE3B1C", "#707176")
# read data
complexGenes <- c("SMC5", "SMC6", "NSMCE1", "NSMCE2", "NSMCE3", "NSMCE4A", "EID3", "SLF1", "SLF2")
#
dgea <- readRDS("./data/metabric-brca/dgea-limma.rds")
dgea <- dgea[order(abs(dgea$logFC), decreasing = TRUE), ]
geneList <- dgea$gene
#
# load TCGA BRCA data
expmat <- readRDS("./reviewer-addressing/tcga/BRCA-expr.rds")
# print(head(colnames(expmat)))
colnames(expmat) <- gsub("-01[A-Z]-.*$", "", colnames(expmat))
# print(head(colnames(expmat)))
# print(range(nchar(colnames(expmat))))
expmat <- t(expmat)
signatureScore <- rowMeans(expmat[
    ,
    colnames(expmat) %in% geneList[seq_len(10)]
])
#
df <- readRDS("./data/cbioportal/formatted.rds")
df <- df[rownames(expmat), ]
df$sig <- signatureScore

df$aneuploidyScore <- as.numeric(df$aneuploidyScore)
plot_df <- df[, c("sig", "aneuploidyScore", "ploidy", "isalt")]
plot_df$sig <- as.numeric(plot_df$sig)


a_plot_df <- na.omit(data.frame(
    aneuploidyScore = plot_df$aneuploidyScore,
    # altStat = as.numeric(plot_df$isalt)
    altStat = plot_df$sig
))

a_plot_df$altStat <- as.factor(ifelse(
    a_plot_df$altStat >= median(a_plot_df$altStat),
    "High expr", "Low Expr"
))
# a_plot_df$altStat <- as.numeric(a_plot_df$altStat >= median(a_plot_df$altStat))

# p_plot_df <- na.omit(data.frame(
#     ploidy = plot_df$ploidy,
#     altStat = plot_df$sig
#     # altStat = as.numeric(plot_df$isalt)
# ))
#
#
# p_plot_df$altStat <- ifelse(
#     p_plot_df$altStat >= median(p_plot_df$altStat),
#     "High expr", "Low Expr"
# )
# table(p_plot_df$altStat)

pa <- ggplot(a_plot_df, aes(
    x = altStat, y = aneuploidyScore,
    fill = altStat, colour = altStat
)) +
    geom_flat_violin(
        position = position_nudge(x = .25, y = 0),
        adjust = 2, trim = TRUE
    ) +
    geom_point(position = position_jitter(width = .15), size = .25) +
    geom_boxplot(aes(x = as.numeric(altStat) + 0.25, y = aneuploidyScore),
        outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK"
    ) +
    ylab("Aneuploidy score") +
    xlab("") +
    ylim(0, 40) +
    scale_colour_manual(values = colors) +
    scale_fill_manual(
        values = colors,
        name = "Signature\nExpression"
    ) +
    theme_cowplot() +
    guides(color = "none") +
    # wilcox t test raw output. comment to add custom p value
    stat_compare_means(method = "wilcox.test") +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "plain"))
# ploidy
# pp <- ggplot(p_plot_df, aes(
#     x = altStat, y = ploidy,
#     fill = altStat, colour = altStat
# )) +
#     geom_flat_violin(
#         position = position_nudge(x = .25, y = 0),
#         adjust = 2, trim = TRUE
#     ) +
#     geom_point(position = position_jitter(width = .15), size = .25) +
#     geom_boxplot(aes(x = as.numeric(altStat) + 0.25, y = ploidy),
#         outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK"
#     ) +
#     ylab("Ploidy score") +
#     xlab("") +
#     ylim(1, 7) +
#     # ggtitle(i) +
#     scale_colour_manual(values = colors) +
#     scale_fill_manual(values = colors) +
#     theme_cowplot() +
#     # wilcox t test raw output. comment to add custom p value
#     stat_compare_means(method = "wilcox.test") +
#     theme(plot.title = element_text(hjust = 0.5, size = 12, face = "plain"))


# pdf("./reviewer-addressing/plot/lfcExprSig-Fig3ab.pdf", width = 7, height = 5)
pdf("./reviewer-addressing/plot/lfcExprSig-Fig3ab.pdf", width = 4, height = 5)
plot(pa)
# plot(ggarrange(pa, pp, nrow = 1, ncol = 2, align = "hv"))
dev.off()
print("done")
